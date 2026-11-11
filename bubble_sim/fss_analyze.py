#!/usr/bin/env python3
"""
Analyse runs/fss/results*.jsonl produced by fss_scan.py (all shards).

Per (anisotropy a, size L):
  * <h_d>, var(h_d)          -> var(h_d) ~ L^(-2/nu)  (fit per a)
  * <h_d>(L)                 -> h_d(L) = h_inf + B / L  (bubble: Laplace pressure ~ 1/R_c ~ 1/L)
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=os.path.join(HERE, "runs", "fss"))
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
        S = [s for s in summ if s["aniso"] == a and s["n"] > 2]
        if len(S) < 2:
            continue
        L = np.array([s["L"] for s in S], float)
        v = np.array([s["hd_var"] for s in S])
        ve = np.array([s["hd_var_err"] for s in S])
        wts = 1 / np.maximum(ve / v, 1e-3)                  # weights in log space
        sl, ic = np.polyfit(np.log(L), np.log(v), 1, w=wts)
        nu = -2 / sl
        m = np.array([s["hd_mean"] for s in S])
        A = np.c_[np.ones_like(L), 1 / L]
        (hinf, B), *_ = np.linalg.lstsq(A, m, rcond=None)
        fits[a] = dict(slope=sl, ic=ic, nu=nu, hinf=hinf, B=B)
        print(f"a={a:4.2f}: var(h_d) ~ L^{sl:.2f}  -> nu = {nu:.2f};   <h_d> = {hinf:.4f} + {B:.2f}/L")

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
        ax.errorbar(1 / L, [s["hd_mean"] for s in S], [s["hd_err"] for s in S], fmt="o", color=cols[a], label=f"a={a}")
        if a in fits:
            x = np.linspace(0, 1 / L.min() * 1.05, 50)
            ax.plot(x, fits[a]["hinf"] + fits[a]["B"] * x, "-", color=cols[a], lw=0.8)
    ax.set_xlabel("1/L"); ax.set_ylabel(r"$\langle h_d\rangle$"); ax.set_xlim(left=0); ax.legend(fontsize=8)
    ax.set_title(r"mean depinning field, fit $h_\infty + B/L$")

    ax = axs[0, 1]
    for a in As:
        S = [s for s in summ if s["aniso"] == a and s["n"] > 2]
        L = np.array([s["L"] for s in S])
        ax.errorbar(L, [s["hd_var"] for s in S], [s["hd_var_err"] for s in S], fmt="o", color=cols[a],
                    label=f"a={a}" + (f", ν={fits[a]['nu']:.2f}" if a in fits else ""))
        if a in fits:
            ax.plot(L, np.exp(fits[a]["ic"]) * L ** fits[a]["slope"], "-", color=cols[a], lw=0.8)
    if Ls:
        L = np.array(Ls, float)
        for nu, ls in ((4 / 3, ":"), (1.76, "--")):
            ref = [s for s in summ if s["L"] == Ls[0]]
            y0 = np.mean([s["hd_var"] for s in ref])
            ax.plot(L, y0 * (L / L[0]) ** (-2 / nu), "k" + ls, lw=0.8, label=f"ν={nu:.2f}")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("L"); ax.set_ylabel(r"var($h_d$)")
    ax.legend(fontsize=7); ax.set_title(r"sample-to-sample fluctuations, $\sim L^{-2/\nu}$")

    ax = axs[0, 2]
    if fits:
        aa = sorted(fits)
        ax.plot(aa, [fits[a]["nu"] for a in aa], "o-", label=r"$\nu$ from var($h_d$)")
        ax.axhline(4 / 3, color="k", ls=":", lw=0.8, label="qEW / white (1.33)")
        ax.axhline(1.76, color="k", ls="--", lw=0.8, label="letter, correlated (1.76)")
        ax2 = ax.twinx()
        ax2.plot(aa, [fits[a]["hinf"] for a in aa], "s-", color="C3", label=r"$h_\infty$")
        ax2.set_ylabel(r"$h_\infty$", color="C3")
        h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="center left")
        ax.set_xlabel("anisotropy a"); ax.set_ylabel(r"$\nu$")
        ax.set_title("exponent and extrapolated threshold vs anisotropy")

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
