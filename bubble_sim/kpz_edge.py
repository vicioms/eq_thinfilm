#!/usr/bin/env python3
"""
qKPZ tests on the edge sides (wall normal || b) of the critical bubbles.

Background (Tang, Kardar & Dhar, PRL 74, 920 (1995), cond-mat/9409002; Amaral,
Barabasi & Stanley, PRL 73, 62 (1994); Mukerjee, Bonachela, Munoz & Wiese,
arXiv:2207.08341):
  * in an anisotropic medium the threshold depends on the wall tilt s; around
    the hard direction (maximum threshold) F_c(s) - F_c(0) ~ -|s|^(1/nu(1-zeta))
    and a KPZ term lambda (grad u)^2 survives at depinning -> qKPZ / directed-
    percolation depinning: zeta = 0.633, nu = 1.733, z = 1, v(F_c, s) ~ |s|;
  * the sign of lambda v matters: positive -> DPD; negative -> faceting.
Here the edge orientation is the hard direction (the edge sides barely move,
the screw tips run), so the edge sides are the candidates for positive qKPZ.

Tests (L >= 2048, all a; a = 0 is the control):
  1. large-scale exponent of the edge sides: local zeta_eff(q) of the intrinsic
     side spectrum S_u (sides_stats.npz), down to the lowest q, with bootstrap
     errors; qKPZ predicts a drift towards 0.63 at large scales.
  2. up-down asymmetry (a KPZ term breaks u -> -u): skewness of the wall's
     excess curvature kappa = d(dpsi)/ds, smoothed over sigma = 2..32, on each
     kind of side. kappa > 0 = bulge outwards. Positive-lambda KPZ growth makes
     rounded hills and sharp valleys -> negative skewness of kappa.
  3. advance vs wall direction near the hard direction: support function of the
     wall about the nucleation centre on a fine angular grid; a facet along the
     hard direction (kinetic Wulff shape of v(s) ~ |s|) gives a linear kink,
     h(phi) - h(0) ~ |phi|, a smooth shape gives ~ phi^2.

    python kpz_edge.py --dir runs/fss_D0.04
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fss_contours import Q_EDGES, resample_closed, smooth_closed, normals_angle, tangent_angle  # noqa: E402
from sides import side_extent, SIG_SIDE  # noqa: E402

QC = np.sqrt(Q_EDGES[1:] * Q_EDGES[:-1])
KSIG = (2.0, 4.0, 8.0, 16.0, 32.0)
NPHI = 1440
B = 300
RNG = np.random.default_rng(77)


def skew(x):
    x = x - x.mean()
    s = x.std()
    return np.mean(x ** 3) / s ** 3 if s > 0 else np.nan


def one(args):
    xy, L, a, seed, R0 = args
    pts, P = resample_closed(xy)
    ds = P / len(pts)
    sm = smooth_closed(pts, max(SIG_SIDE, P / 20), P)
    phn = normals_angle(sm, P)
    psi = tangent_angle(pts)[0]
    dpsi = psi - tangent_angle(sm)[0]
    dpsi -= 2 * np.pi * np.round(np.median(dpsi) / (2 * np.pi))
    out = dict(L=L, aniso=a, seed=seed)
    n = len(dpsi)
    qf = 2 * np.pi * np.fft.rfftfreq(n, d=ds)
    F = np.fft.rfft(dpsi)
    kap = {}
    for sg in KSIG:     # excess curvature smoothed over sg (periodic, spectral derivative)
        kap[sg] = np.fft.irfft(1j * qf * F * np.exp(-0.5 * (qf * sg) ** 2), n)
    for kind, dirs in (("edge", (0.0, np.pi)), ("screw", (np.pi / 2, -np.pi / 2))):
        sk = {sg: [] for sg in KSIG}
        for dr in dirs:
            ext = side_extent(sm, phn, dr)
            if ext is None or len(ext[2]) * ds < 256:
                continue
            idx = ext[2]
            core = idx[int(0.1 * len(idx)):int(0.9 * len(idx))]     # drop the ends of the side
            for sg in KSIG:
                sk[sg].append(skew(kap[sg][core]))
        for sg in KSIG:
            out[f"{kind}_skew_{sg:g}"] = np.mean(sk[sg]) if sk[sg] else np.nan
    # whole wall, for reference
    for sg in KSIG:
        out[f"all_skew_{sg:g}"] = skew(kap[sg])
    # fine support function about the nucleation centre
    c0 = np.array([(L - 1) / 2, (L - 1) / 2])
    ang = np.linspace(0, 2 * np.pi, NPHI, endpoint=False)
    nh = np.c_[np.cos(ang), np.sin(ang)]
    sup = ((pts - c0) @ nh.T).max(axis=0) - R0
    # fold the four symmetric copies onto [0, 90] deg
    k = np.arange(NPHI)
    f = 0.25 * (sup + sup[(-k) % NPHI] + sup[(k + NPHI // 2) % NPHI] + sup[(NPHI // 2 - k) % NPHI])
    out["support"] = f[: NPHI // 4 + 1]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=os.path.join(HERE, "runs", "fss_D0.04"))
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--lmin", type=int, default=2048)
    args = ap.parse_args()
    d = args.dir
    import fss_shapes
    C = fss_shapes.load_contours(d)
    R0 = {}
    for fn in glob.glob(os.path.join(d, "results*.jsonl")):
        for line in open(fn):
            r = json.loads(line)
            R0[(r["L"], r["aniso"], r["seed"])] = r["R0"]
    jobs = [(C["contours"][i][C["main"][i]], int(C["L"][i]), float(C["aniso"][i]), int(C["seed"][i]),
             R0[(int(C["L"][i]), float(C["aniso"][i]), int(C["seed"][i]))])
            for i in range(len(C["L"])) if C["L"][i] >= args.lmin]
    with ProcessPoolExecutor(args.workers) as ex:
        rows = list(ex.map(one, jobs, chunksize=2))
    S = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    SS = dict(np.load(os.path.join(d, "sides_stats.npz")))
    report(d, S, SS)


def boot_mean(x):
    x = x[np.isfinite(x)]
    if len(x) < 3:
        return np.nan, np.nan
    return x.mean(), np.std([x[RNG.integers(0, len(x), len(x))].mean() for _ in range(B)])


def report(d, S, SS):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm
    As = sorted(set(S["aniso"].tolist()))
    Ls = sorted(set(S["L"].tolist()))
    ca = {a: cm.viridis(i / max(1, len(As) - 1)) for i, a in enumerate(As)}
    out = []
    fig, axs = plt.subplots(2, 3, figsize=(17, 10))

    # 1. local zeta of edge (and screw) sides
    for col, kind in enumerate(("edge", "screw")):
        ax = axs[0, col]
        for a in As:
            for L, ls in zip(Ls, ("--", "-")):
                m_sel = (SS["L"] == L) & (SS["aniso"] == a)
                X = SS[f"{kind}_iSu"][m_sel]
                X = X[np.any(np.isfinite(X), axis=1)]
                ln = np.nanmean(SS[f"{kind}_ilen"][m_sel])
                if len(X) < 8 or not np.isfinite(ln):
                    continue
                k = 2
                ok = (QC[k:-k] > 2 * np.pi * 4 / ln) & (QC[k:-k] < 0.12)
                if ok.sum() < 2:
                    continue

                def zl(Y):
                    m = np.nanmean(Y, 0)
                    lq, lS = np.log(QC), np.log(m)
                    return -(1 + (lS[2 * k:] - lS[:-2 * k]) / (lq[2 * k:] - lq[:-2 * k])) / 2
                z0 = zl(X)
                zb = np.array([zl(X[RNG.integers(0, len(X), len(X))]) for _ in range(B)])
                ze = np.nanstd(zb, 0)
                ax.errorbar(QC[k:-k][ok], z0[ok], ze[ok], fmt="o" + ls, color=ca[a], ms=3, lw=0.8,
                            label=f"a={a}, L={L}")
                if L == Ls[-1]:
                    out.append(f"local zeta, {kind:5s} a={a:4.2f} L={L}: " +
                               "  ".join(f"q={qq:.4f}:{zz:.2f}±{ee:.2f}" for qq, zz, ee in zip(QC[k:-k][ok], z0[ok], ze[ok])))
        ax.axhline(1.25, color="k", ls="--", lw=0.6, label="qEW 1.25")
        ax.axhline(0.633, color="C3", ls="--", lw=0.8, label="qKPZ 0.633")
        ax.set_xscale("log"); ax.set_ylim(0, 2); ax.set_xlabel("q (arc length)")
        ax.set_ylabel(r"local $\zeta_{eff}(q)$ of the side")
        ax.set_title(f"{kind} sides: local exponent (bootstrap errors)"); ax.legend(fontsize=6, ncol=2)

    # 2. skewness of excess curvature
    ax = axs[0, 2]
    out.append("")
    out.append("skewness of the excess curvature kappa (positive = bulge out), mean over sides, bootstrap error:")
    for kind, mk, clr in (("edge", "o", "C0"), ("screw", "s", "C3"), ("all", "^", "0.4")):
        for L, ls in zip(Ls, ("--", "-")):
            for sg_i, sg in enumerate((4.0, 16.0)):
                vals = [boot_mean(S[f"{kind}_skew_{sg:g}"][(S["L"] == L) & (S["aniso"] == a)]) for a in As]
                m, e = np.array(vals).T
                if sg == 16.0:
                    ax.errorbar(np.array(As) + 0.01 * (Ls.index(L)), m, e, fmt=mk + ls, color=clr, ms=4,
                                label=f"{kind}, L={L}")
                out.append(f"  {kind:5s} L={L} sigma={sg:g}: " + "  ".join(f"a={a}: {mm:+.3f}±{ee:.3f}" for a, mm, ee in zip(As, m, e)))
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("anisotropy a"); ax.set_ylabel(r"skewness of $\kappa$ smoothed over 16 sites")
    ax.set_title("up-down asymmetry of the wall (KPZ: ≠ 0)"); ax.legend(fontsize=6.5, ncol=2)
    ax = axs[1, 2]
    for kind, clr in (("edge", "C0"), ("screw", "C3")):
        for a, ls in ((As[0], ":"), (As[-1], "-")):
            m = [boot_mean(S[f"{kind}_skew_{sg:g}"][(S["L"] == Ls[-1]) & (S["aniso"] == a)]) for sg in KSIG]
            mm, ee = np.array(m).T
            ax.errorbar(KSIG, mm, ee, fmt="o" + ls, color=clr, ms=4, label=f"{kind}, a={a}")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xscale("log"); ax.set_xlabel("smoothing scale σ of the curvature"); ax.set_ylabel(r"skewness of $\kappa$")
    ax.set_title(f"curvature skewness vs scale, L={Ls[-1]}"); ax.legend(fontsize=7)

    # 3. support function near the hard direction
    ang = np.rad2deg(np.linspace(0, 2 * np.pi, NPHI, endpoint=False))[: NPHI // 4 + 1]
    ax = axs[1, 0]
    ax2 = axs[1, 1]
    out.append("")
    out.append("support h(phi) - h(0) near the hard direction (phi = 0, normal || b): fit ~ |phi|^p over 1..15 deg")
    for a in As:
        for L, ls in zip(Ls, ("--", "-")):
            sup = S["support"][(S["L"] == L) & (S["aniso"] == a)]
            m = sup.mean(0)
            dh = m - m[0]
            ok = (ang > 0) & (ang <= 30)
            ax.plot(ang[ok], dh[ok], ls, color=ca[a], lw=1, label=f"a={a}, L={L}" if L == Ls[-1] else None)
            fit = (ang >= 1) & (ang <= 15) & (dh > 0)

            def pfit(Y):
                mm = Y.mean(0)
                dd = mm - mm[0]
                f2 = fit & (dd > 0)
                return np.polyfit(np.log(ang[f2]), np.log(dd[f2]), 1)[0] if f2.sum() > 5 else np.nan
            p0 = pfit(sup)
            pb = np.nanstd([pfit(sup[RNG.integers(0, len(sup), len(sup))]) for _ in range(B)])
            out.append(f"  a={a:4.2f} L={L}: h(0) - R0 = {m[0]:7.1f}, h(5)-h(0) = {np.interp(5, ang, dh):6.1f}, "
                       f"h(15)-h(0) = {np.interp(15, ang, dh):6.1f}, exponent p = {p0:.2f} ± {pb:.2f}")
            if L == Ls[-1]:
                ax2.errorbar([a], [p0], [pb], fmt="o", color=ca[a], ms=6)
    ax.set_xlabel("direction φ to b [deg]"); ax.set_ylabel(r"$h(\phi)-h(0)$")
    ax.set_xscale("log"); ax.set_yscale("log")
    x = np.geomspace(1, 30, 10)
    ax.plot(x, 3 * x, "k:", lw=0.8, label="∝ φ (facet)")
    ax.plot(x, 0.3 * x ** 2, "k--", lw=0.8, label="∝ φ² (smooth)")
    ax.set_title("advance of the wall vs direction, near the hard direction"); ax.legend(fontsize=7)
    ax2.axhline(1, color="k", ls=":", lw=0.8, label="facet (kink)")
    ax2.axhline(2, color="k", ls="--", lw=0.8, label="smooth")
    ax2.set_xlabel("anisotropy a"); ax2.set_ylabel(r"exponent p of $h(\phi)-h(0)\sim\phi^p$, 1°–15°")
    ax2.set_title(f"kink exponent, L={Ls[-1]}"); ax2.legend(fontsize=7)
    fig.tight_layout()
    fn = os.path.join(d, "kpz_edge.png")
    fig.savefig(fn, dpi=110)
    txt = "\n".join(out)
    print(txt)
    with open(os.path.join(d, "kpz_edge.txt"), "w") as f:
        f.write(txt + "\n")
    print("saved", fn)


if __name__ == "__main__":
    main()
