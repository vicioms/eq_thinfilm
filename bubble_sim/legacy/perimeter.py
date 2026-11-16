"""Length of the critical front vs system size: P(L) of the flat walls, P(R_eff) of the bubbles.

    python perimeter.py [--flat D1 D2 ...] [--bubbles B1 B2 ...] [--out D]

P = arc length of the phi = 0 contour, resampled every 0.5 site (as in fss_flat.py /
fss_contours.py), the same sub-pixel marching-squares length as NB7_perimeter.ipynb
(P vs <R> of the experimental bubbles).
  * flat walls: P over one period vs L, per orientation and a; P / L -> const if the wall
    is self-affine with zeta < 1; P / L ~ L^(zeta - 1) once the slopes grow with scale (zeta > 1)
    or the wall facets;
  * bubbles: P vs R_eff = sqrt(area / pi) of every critical bubble, pooled over L and binned
    in log R_eff; P / (2 pi R_eff) as in NB7.
alpha = d ln P / d ln (L or R): fit over all points, and local (between successive L / R bins).
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # bubble_sim/ (this file is in legacy/)
sys.path[:0] = [os.path.dirname(os.path.abspath(__file__)), HERE]
RNG = np.random.default_rng(7)
B = 300


def flat_lengths(d):
    import fss_flat as F
    W = F.load_walls(d)
    P = np.array([F.periodic_resample(c, int(L))[1] for c, L in zip(W["curve"], W["L"])])
    name = {0.0: "screw wall (line ∥ b)", 90.0: "edge wall (line ⊥ b)"}.get(W["burgers_angle"],
                                                                          f"b at {W['burgers_angle']}°")
    return name, W["L"].astype(float), W["aniso"].astype(float), P


def bubble_lengths(d):
    """(L, a, R_eff, P, R0) per critical bubble; from contour_stats.npz if present, else contours.npz."""
    import json, glob
    import fss_contours as C
    st = os.path.join(d, "contour_stats.npz")
    if os.path.exists(st):
        z = np.load(st)
        return z["L"].astype(float), z["aniso"].astype(float), z["R_eff"], z["P"], z["R0"]
    z = np.load(os.path.join(d, "contours.npz"))
    st_, off, start, sc = z["steps"].astype(np.int64), z["offsets"], z["start"].astype(np.int64), float(z["scale"])
    n = len(z["L"])
    P, R = np.empty(n), np.empty(n)
    for i in range(n):
        xy = np.vstack([start[i], start[i] + np.cumsum(st_[off[i]:off[i + 1]], axis=0)]) / sc
        pts, P[i] = C.resample_closed(xy)
        x, y = pts[:, 0], pts[:, 1]
        R[i] = np.sqrt(abs(0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)) / np.pi)
    R0 = {}
    for fn in glob.glob(os.path.join(d, "results*.jsonl")):
        for line in open(fn):
            r = json.loads(line)
            R0[(r["L"], r["aniso"], r["seed"])] = r.get("R0") or np.nan
    r0 = np.array([R0.get((int(L), float(a), int(s)), np.nan) for L, a, s in zip(z["L"], z["aniso"], z["seed"])])
    return z["L"].astype(float), z["aniso"].astype(float), R, P, r0


def slope(x, y, w=None):
    return np.polyfit(np.log(x), np.log(y), 1, w=w)[0]


def flat_alpha(L, P, Ls):
    """alpha over all L (fit of the mean P per L) and between successive L, bootstrap over walls."""
    groups = [P[L == l] for l in Ls]
    m = np.array([g.mean() for g in groups])

    def all_(m):
        return slope(Ls, m)
    loc = np.log(m[1:] / m[:-1]) / np.log(Ls[1:] / Ls[:-1])
    bs = []
    for _ in range(B):
        mb = np.array([g[RNG.integers(0, len(g), len(g))].mean() for g in groups])
        bs.append(np.r_[all_(mb), np.log(mb[1:] / mb[:-1]) / np.log(Ls[1:] / Ls[:-1])])
    e = np.std(bs, 0)
    return (all_(m), e[0]), list(zip(loc, e[1:])), m, np.array([g.std(ddof=1) / np.sqrt(len(g)) for g in groups])


def bubble_bins(R, P, nmin=8, nb=14):
    edges = np.geomspace(R.min() * 0.999, R.max() * 1.001, nb + 1)
    i = np.digitize(R, edges) - 1
    out = []
    for k in range(nb):
        s = i == k
        if s.sum() >= nmin:
            out.append((np.exp(np.log(R[s]).mean()), P[s].mean(), P[s].std(ddof=1) / np.sqrt(s.sum()),
                        (P[s] / (2 * np.pi * R[s])).mean(), (P[s] / (2 * np.pi * R[s])).std(ddof=1) / np.sqrt(s.sum()),
                        s.sum()))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flat", nargs="*", default=["runs/flat_D0.04_b0_A2", "runs/flat_D0.04_b90_A2"])
    ap.add_argument("--bubbles", nargs="*", default=["runs/fss_D0.04", "runs/fss_D0.04_R32"])
    ap.add_argument("--out", default="runs/perimeter")
    args = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm
    os.makedirs(args.out, exist_ok=True)
    txt = []
    ncol = max(len(args.flat), len(args.bubbles))
    fig, axs = plt.subplots(2, 2 * ncol, figsize=(5.2 * 2 * ncol, 9), squeeze=False)

    for j, d in enumerate(args.flat):
        name, L, a, P = flat_lengths(os.path.join(HERE, d) if not os.path.isabs(d) else d)
        txt.append(f"== flat, {name}  ({d})")
        Ls = np.unique(L)
        for ak in np.unique(a):
            s = a == ak
            (al, ale), loc, m, e = flat_alpha(L[s], P[s], Ls[np.isin(Ls, np.unique(L[s]))])
            Lk = Ls[np.isin(Ls, np.unique(L[s]))]
            c = cm.viridis(ak)
            axs[0, 2 * j].errorbar(Lk, m, e, fmt="o-", color=c, label=f"a={ak:g}: α={al:.3f}±{ale:.3f}")
            axs[0, 2 * j + 1].errorbar(Lk, m / Lk, e / Lk, fmt="o-", color=c, label=f"a={ak:g}")
            txt.append(f"  a={ak:4.2f}: " + "  ".join(f"L={l:.0f}: P/L={mm / l:.3f}" for l, mm in zip(Lk, m))
                       + f"   alpha = {al:.3f} ± {ale:.3f};  local " + " ".join(f"{x:.2f}±{y:.2f}" for x, y in loc))
        for k, (lab, yl) in enumerate((("P", "P (one period)"), ("P/L", "P / L"))):
            ax = axs[0, 2 * j + k]
            ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("L"); ax.set_ylabel(yl)
            ax.set_title(f"flat: {name}", fontsize=10); ax.legend(fontsize=8)
        Lr = np.array([Ls.min(), Ls.max()])
        axs[0, 2 * j].plot(Lr, Lr * P[L == Ls.min()].mean() / Ls.min(), "k:", lw=0.8)

    for j, d in enumerate(args.bubbles):
        L, a, R, P, R0 = bubble_lengths(os.path.join(HERE, d) if not os.path.isabs(d) else d)
        fixed = np.allclose(np.nanstd(R0), 0) if np.isfinite(R0).any() else False
        lab = f"R0 = {np.nanmean(R0):.0f}" if fixed else "R0 = L/6"
        txt.append(f"== bubbles, {lab}  ({d}, {len(P)} bubbles)")
        As = np.unique(a)
        for ak in As:
            s = (a == ak) & np.isfinite(P) & (R > 0)
            bb = bubble_bins(R[s], P[s])
            c = cm.viridis(ak)
            # alpha: per-bubble fit (bootstrap) and fit of the binned means
            al = slope(R[s], P[s])
            alb = [slope(R[s][ii], P[s][ii]) for ii in (RNG.integers(0, s.sum(), s.sum()) for _ in range(B))]
            ok = bb[:, 0] > 2 * np.nanmin(R0[s]) if fixed else np.ones(len(bb), bool)
            alg = slope(bb[ok, 0], bb[ok, 1], w=bb[ok, 1] / bb[ok, 2]) if ok.sum() >= 3 else np.nan
            axs[1, 2 * j].errorbar(bb[:, 0], bb[:, 1], bb[:, 2], fmt="o-", ms=3, color=c,
                                   label=f"a={ak:g}: α={al:.3f}±{np.std(alb):.3f}")
            axs[1, 2 * j + 1].errorbar(bb[:, 0], bb[:, 3], bb[:, 4], fmt="o-", ms=3, color=c, label=f"a={ak:g}")
            txt.append(f"  a={ak:4.2f}: alpha (all bubbles) = {al:.3f} ± {np.std(alb):.3f}"
                       f"   alpha (binned{', R > 2 R0' if fixed else ''}) = {alg:.3f}   "
                       + "  ".join(f"R={r:.0f}: P/2piR={p:.3f}" for r, p in zip(bb[:, 0], bb[:, 3])))
            if not fixed:
                for l in np.unique(L[s]):
                    t = s & (L == l)
                    txt.append(f"      L={l:.0f}: <R_eff>={R[t].mean():.0f}  <P>/(2 pi <R_eff>) = "
                               f"{P[t].mean() / (2 * np.pi * R[t].mean()):.3f}  (n={t.sum()})")
        for k, yl in enumerate(("P", r"$P / (2\pi R_{eff})$")):
            ax = axs[1, 2 * j + k]
            ax.set_xscale("log"); ax.set_xlabel(r"$R_{eff} = \sqrt{A/\pi}$"); ax.set_ylabel(yl)
            if k == 0:
                ax.set_yscale("log")
            ax.set_title(f"critical bubbles, {lab}", fontsize=10); ax.legend(fontsize=8)
    for ax in axs.ravel():
        if not ax.has_data():
            ax.axis("off")
    fig.tight_layout()
    fn = os.path.join(args.out, "perimeter.png")
    fig.savefig(fn, dpi=100)
    with open(os.path.join(args.out, "perimeter.txt"), "w") as f:
        f.write("\n".join(txt) + "\n")
    print("\n".join(txt))
    print("saved", fn)


if __name__ == "__main__":
    main()
