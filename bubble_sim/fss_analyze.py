#!/usr/bin/env python3
"""
Analyse runs/fss/results*.jsonl produced by fss_scan.py (all shards).

Per (anisotropy a, size L):
  * <h_d>, var(h_d)          -> var(h_d) ~ L^(-2/nu)  (fit per a)
  * <h_d>(L)                 -> h_d(L) = h_ref + C ln(L / L_max)   (fit L >= lmin)
        <h_d> is not monotonic in L: the Laplace term ~ 1/R0 lowers it at small L,
        but at large L it grows again (the wall has to cross the strongest pinning
        region on its way to the border, and that grows with the swept area). An
        extrapolation h_inf + B/L is therefore not meaningful; we report the local
        drift C = d<h_d>/d ln L, the fitted value h_ref at the largest L, the
        minimum of <h_d> over L, and chi^2/dof of both models.
  * threshold                -> <h_d> = h_c + B/L - A L^(-1/nu)   (fit L >= lmin-hc)
        B/L: Laplace pressure of a bubble of radius ~ L/6 (raises h_d at small L);
        -A L^(-1/nu): usual finite-size shift of the sample threshold in an L x L box
        (h_d is set by the strongest barrier the wall must cross; approaches h_c from
        below). nu is fixed to the var(h_d) fit; the systematic error on h_c is half
        the spread of h_c over nu +- its error.
  * shape anisotropy of the critical bubble:
        eps = (aspect - 1) / (aspect + 1),  theta = long-axis orientation
        Q_b = < eps cos 2(theta - theta_exp) >   theta_exp = burgers_angle + 90
              (signed: > 0 elongated perpendicular to b, as predicted)
        |Q| = | < eps exp(2 i theta) > |         (any common axis)
  * fraction of multivalued angles (overhangs / pinch-offs) and u_o / L

Writes runs/fss/summary.csv and runs/fss/fss_summary.png, then (if shapes.npz or
configs/ exist) bubbles.png and Sq.png via fss_shapes.py.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def bootstrap_var_err(x, n=500, rng=np.random.default_rng(0)):
    if len(x) < 3:
        return np.nan
    v = [np.var(rng.choice(x, len(x)), ddof=1) for _ in range(n)]
    return float(np.std(v))


def wls(A, y, sigma):
    """Weighted least squares: returns (coefficients, covariance)."""
    W = 1 / np.asarray(sigma) ** 2
    M = A.T @ (A * W[:, None])
    cov = np.linalg.inv(M)
    return cov @ (A.T @ (W * y)), cov


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=os.path.join(HERE, "runs", "fss"))
    ap.add_argument("--lmin", type=int, default=512,
                    help="fit <h_d> and var(h_d) only for L >= lmin (small bubbles still remember "
                         "the initial circle and are dominated by the Laplace term)")
    ap.add_argument("--lmin-hc", type=int, default=256,
                    help="fit <h_d> = h_c + B/L - A L^(-1/nu) for L >= lmin-hc")
    args = ap.parse_args()
    import glob
    rows = [json.loads(l) for fn in sorted(glob.glob(os.path.join(args.dir, "results*.jsonl")))
            for l in open(fn)]
    ok = [r for r in rows if r["status"] == "pinned" and "aspect" in r]
    print(f"{len(rows)} runs, {len(ok)} usable")

    As = sorted({r["aniso"] for r in ok})
    Ls = sorted({r["L"] for r in ok})
    summ = []
    for a in As:
        for L in Ls:
            g = [r for r in ok if r["aniso"] == a and r["L"] == L]
            if not g:
                continue
            hd = np.array([r["h_d"] for r in g])
            asp = np.array([r["aspect"] for r in g])
            th = np.deg2rad([r["orientation_deg"] for r in g])
            eps = (asp - 1) / (asp + 1)
            th_exp = np.deg2rad(g[0]["burgers_angle"] + 90)
            qb = eps * np.cos(2 * (th - th_exp))
            qz = eps * np.exp(2j * th)
            mv = np.array([r["frac_multivalued"] for r in g])
            uo = np.array([r["u_o"] for r in g]) / L
            n = len(g)
            summ.append(dict(aniso=a, L=L, n=n,
                             hd_mean=hd.mean(), hd_err=hd.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan,
                             hd_var=hd.var(ddof=1) if n > 1 else np.nan, hd_var_err=bootstrap_var_err(hd),
                             aspect_mean=asp.mean(), Qb=qb.mean(), Qb_err=qb.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan,
                             Qabs=abs(qz.mean()), multival=mv.mean(), multival_err=mv.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan,
                             frac_with_overhang=(mv > 0).mean(), uo_over_L=uo.mean()))

    keys = list(summ[0].keys())
    with open(os.path.join(args.dir, "summary.csv"), "w") as f:
        f.write(",".join(keys) + "\n")
        for s in summ:
            f.write(",".join(f"{s[k]:.6g}" if isinstance(s[k], float) else str(s[k]) for k in keys) + "\n")

    # fits
    fits = {}
    print(f"\n{'a':>5} {'L':>5} {'n':>4} {'<h_d>':>10} {'var(h_d)':>10} {'Q_b':>7} {'|Q|':>6} {'<aspect>':>8} {'multival':>8} {'P(ovh)':>6}")
    for s in summ:
        print(f"{s['aniso']:5.2f} {s['L']:5d} {s['n']:4d} {s['hd_mean']:10.5f} {s['hd_var']:10.3e} "
              f"{s['Qb']:7.3f} {s['Qabs']:6.3f} {s['aspect_mean']:8.3f} {s['multival']:8.3f} {s['frac_with_overhang']:6.2f}")
    print()
    for a in As:
        S = [s for s in summ if s["aniso"] == a and s["n"] > 2 and s["L"] >= args.lmin]
        if len(S) < 2:
            continue
        L = np.array([s["L"] for s in S], float)
        # var(h_d) ~ L^slope: weighted least squares in log space (sigma = relative error)
        v = np.array([s["hd_var"] for s in S])
        sv = np.maximum(np.array([s["hd_var_err"] for s in S]) / v, 1e-3)
        (ic, sl), cv = wls(np.c_[np.ones_like(L), np.log(L)], np.log(v), sv)
        nu = -2 / sl
        nu_err = 2 / sl ** 2 * np.sqrt(cv[1, 1])
        # <h_d> = h_ref + C ln(L / L_max), weighted by the standard errors
        m = np.array([s["hd_mean"] for s in S])
        sm = np.array([s["hd_err"] for s in S])
        Lref = L.max()
        A_log = np.c_[np.ones_like(L), np.log(L / Lref)]
        (href, C), ch = wls(A_log, m, sm)
        dof = max(1, len(L) - 2)
        chi2_log = float(np.sum(((m - A_log @ [href, C]) / sm) ** 2)) / dof
        A_inv = np.c_[np.ones_like(L), 1 / L]
        c_inv, _ = wls(A_inv, m, sm)
        chi2_inv = float(np.sum(((m - A_inv @ c_inv) / sm) ** 2)) / dof
        Sall = [s for s in summ if s["aniso"] == a]
        # <h_d> = h_c + B/L - A L^(-1/nu), nu fixed from var(h_d)
        H = [s for s in Sall if s["n"] > 2 and s["L"] >= args.lmin_hc]
        LH = np.array([s["L"] for s in H], float)
        mH = np.array([s["hd_mean"] for s in H]); sH = np.array([s["hd_err"] for s in H])

        def hc_fit(nu_):
            A_ = np.c_[np.ones_like(LH), 1 / LH, -LH ** (-1 / nu_)]
            c_, C_ = wls(A_, mH, sH)
            return c_, C_, float(np.sum(((mH - A_ @ c_) / sH) ** 2)) / max(1, len(LH) - 3)

        hc = dict(hc=np.nan, hc_err=np.nan, hc_syst=np.nan, hcB=np.nan, hcA=np.nan, chi2_hc=np.nan, Lhc=LH)
        if len(LH) >= 4 and nu > 0:
            c3, C3, chi3 = hc_fit(nu)
            alt = [hc_fit(max(0.5, nu + d))[0][0] for d in (-nu_err, nu_err)]
            hc.update(hc=c3[0], hc_err=np.sqrt(C3[0, 0]), hc_syst=0.5 * abs(alt[1] - alt[0]),
                      hcB=c3[1], hcA=c3[2], chi2_hc=chi3)
        imin = int(np.argmin([s["hd_mean"] for s in Sall]))
        fits[a] = dict(slope=sl, ic=ic, nu=nu, nu_err=nu_err, href=href, href_err=np.sqrt(ch[0, 0]),
                       C=C, C_err=np.sqrt(ch[1, 1]), Lref=Lref, chi2_log=chi2_log, chi2_inv=chi2_inv,
                       hmin=Sall[imin]["hd_mean"], hmin_err=Sall[imin]["hd_err"], Lmin_h=Sall[imin]["L"], Lfit=L,
                       **hc)
        print(f"a={a:4.2f} (L>={args.lmin}): var(h_d) ~ L^{sl:.2f} -> nu = {nu:.2f} +- {nu_err:.2f};  "
              f"<h_d> = {href:.5f}({np.sqrt(ch[0, 0]) * 1e5:.0f}) + {C * 1e4:+.2f}({np.sqrt(ch[1, 1]) * 1e4:.2f})e-4 "
              f"ln(L/{Lref:.0f})  chi2/dof {chi2_log:.1f} (h_inf+B/L: {chi2_inv:.1f});  "
              f"min <h_d> {fits[a]['hmin']:.5f} at L={fits[a]['Lmin_h']}")
        if np.isfinite(hc["hc"]):
            print(f"        (L>={args.lmin_hc}, nu={nu:.2f}): h_c = {hc['hc']:.5f} +- {hc['hc_err']:.5f} (stat) "
                  f"+- {hc['hc_syst']:.5f} (nu)   B={hc['hcB']:.2f}  A={hc['hcA']:.3f}  chi2/dof {hc['chi2_hc']:.1f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm
    cols = {a: cm.viridis(i / max(1, len(As) - 1)) for i, a in enumerate(As)}
    fig, axs = plt.subplots(2, 3, figsize=(15, 9))

    ax = axs[0, 0]
    for a in As:
        S = [s for s in summ if s["aniso"] == a]
        L = np.array([s["L"] for s in S])
        f = fits.get(a)
        ax.errorbar(L, [s["hd_mean"] for s in S], [s["hd_err"] for s in S], fmt="o", color=cols[a],
                    label=f"a={a}" + (f", C={f['C'] * 1e4:+.1f}±{f['C_err'] * 1e4:.1f}e-4" if f else ""))
        if f and np.isfinite(f["hc"]):
            x = np.geomspace(f["Lhc"].min() / 1.1, f["Lhc"].max() * 1.1, 80)
            ax.plot(x, f["hc"] + f["hcB"] / x - f["hcA"] * x ** (-1 / f["nu"]), "-", color=cols[a], lw=0.8)
            ax.axhline(f["hc"], color=cols[a], ls=":", lw=0.7)
    ax.axvspan(args.lmin_hc, max(Ls) * 1.1, color="0.93", zorder=-1)
    ax.set_xscale("log"); ax.set_xlabel("L"); ax.set_ylabel(r"$\langle h_d\rangle$"); ax.legend(fontsize=7)
    ax.set_title(rf"$\langle h_d\rangle = h_c + B/L - A L^{{-1/\nu}}$, $L\geq{args.lmin_hc}$ (dotted: $h_c$)", fontsize=10)

    ax = axs[0, 1]
    for a in As:
        S = [s for s in summ if s["aniso"] == a and s["n"] > 2]
        L = np.array([s["L"] for s in S])
        ax.errorbar(L, [s["hd_var"] for s in S], [s["hd_var_err"] for s in S], fmt="o", color=cols[a],
                    label=f"a={a}" + (f", ν={fits[a]['nu']:.2f}" if a in fits else ""))
        if a in fits:
            Lf = fits[a]["Lfit"]
            ax.plot(Lf, np.exp(fits[a]["ic"]) * Lf ** fits[a]["slope"], "-", color=cols[a], lw=0.8)
    if Ls:
        L = np.array(Ls, float)
        for nu, ls in ((4 / 3, ":"), (1.76, "--")):
            ref = [s for s in summ if s["L"] == Ls[0]]
            y0 = np.mean([s["hd_var"] for s in ref])
            ax.plot(L, y0 * (L / L[0]) ** (-2 / nu), "k" + ls, lw=0.8, label=f"ν={nu:.2f}")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("L"); ax.set_ylabel(r"var($h_d$)")
    ax.axvspan(args.lmin, max(Ls) * 1.1, color="0.93", zorder=-1)
    ax.legend(fontsize=7); ax.set_title(rf"sample-to-sample fluctuations, $\sim L^{{-2/\nu}}$, fit $L\geq{args.lmin}$")

    ax = axs[0, 2]
    if fits:
        aa = sorted(fits)
        ax.errorbar(aa, [fits[a]["nu"] for a in aa], [fits[a]["nu_err"] for a in aa], fmt="o-",
                    label=r"$\nu$ from var($h_d$)")
        ax.axhline(4 / 3, color="k", ls=":", lw=0.8, label="qEW / white (1.33)")
        ax.axhline(1.76, color="k", ls="--", lw=0.8, label="letter, correlated (1.76)")
        ax2 = ax.twinx()
        ax2.errorbar(aa, [fits[a]["hc"] for a in aa],
                     [np.hypot(fits[a]["hc_err"], fits[a]["hc_syst"]) for a in aa], fmt="s-",
                     color="C3", label=r"$h_c$ (stat $\oplus$ $\nu$ syst)")
        ax2.errorbar(aa, [fits[a]["hmin"] for a in aa], [fits[a]["hmin_err"] for a in aa], fmt="^--",
                     color="C1", label=r"min over L of $\langle h_d\rangle$")
        ax2.set_ylabel(r"$\langle h_d\rangle$", color="C3")
        h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="center left")
        ax.set_xlabel("anisotropy a"); ax.set_ylabel(r"$\nu$")
        ax.set_title("exponent and threshold vs anisotropy")

    ax = axs[1, 0]
    for L in Ls:
        S = [s for s in summ if s["L"] == L]
        ax.errorbar([s["aniso"] for s in S], [s["Qb"] for s in S], [s["Qb_err"] for s in S], fmt="o-",
                    label=f"L={L}", color=cm.plasma(Ls.index(L) / max(1, len(Ls) - 1) * 0.8))
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("anisotropy a"); ax.set_ylabel(r"$Q_b=\langle\epsilon\cos 2(\theta-\theta_b-90°)\rangle$")
    ax.set_title("bubble elongation perpendicular to b"); ax.legend(fontsize=8)

    ax = axs[1, 1]
    for a in As:
        S = [s for s in summ if s["aniso"] == a]
        ax.errorbar([s["L"] for s in S], [s["multival"] for s in S], [s["multival_err"] for s in S],
                    fmt="o-", color=cols[a], label=f"a={a}")
    ax.set_xscale("log"); ax.set_xlabel("L"); ax.set_ylabel("fraction of multivalued angles")
    ax.set_title("overhangs at depinning"); ax.legend(fontsize=8)

    ax = axs[1, 2]
    for a in As:
        g = [r for r in ok if r["aniso"] == a and r["L"] == max(Ls)]
        th = np.array([r["orientation_deg"] for r in g])
        ax.hist(th, bins=np.linspace(0, 180, 13), histtype="step", color=cols[a], lw=1.5, label=f"a={a}")
    ax.axvline((ok[0]["burgers_angle"] + 90) % 180, color="k", ls="--", lw=0.8, label="⊥ b")
    ax.set_xlabel("orientation of long axis [deg]"); ax.set_ylabel("count")
    ax.set_title(f"orientation distribution, L={max(Ls)}"); ax.legend(fontsize=7)

    fig.tight_layout()
    out = os.path.join(args.dir, "fss_summary.png")
    fig.savefig(out, dpi=130)
    print("saved", out, "and summary.csv")

    # bubble shapes and S(q): needs shapes.npz (fss_shapes.py extracts it from configs/)
    import glob as _g
    if os.path.exists(os.path.join(args.dir, "shapes.npz")) or _g.glob(os.path.join(args.dir, "configs", "*.npz")):
        import fss_shapes
        sys.argv = [sys.argv[0], "--dir", args.dir]
        fss_shapes.main()
    else:
        print("no shapes.npz / configs: run `python fss_shapes.py --dir", args.dir, "--no-plot` where the scan ran")


if __name__ == "__main__":
    main()
