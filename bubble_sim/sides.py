#!/usr/bin/env python3
"""
Edge and screw sides of the critical bubbles as separate, flat-like interfaces.

With b along x, a bubble at a > 0 is elongated along y. Its wall then has two
kinds of sides:
  * edge sides   (normal || b, wall line perpendicular to b): left and right,
                 described as Cartesian height profiles x(y), cutting the
                 bubble with lines parallel to b (rows y = const);
  * screw sides  (normal _|_ b, wall line parallel to b): top and bottom,
                 described as y(x), cutting with columns x = const.
At a = 0 both kinds are statistically identical (control).

Two estimators per side (both on the same stretch of wall):
  * intrinsic (main): the arc of the contour that forms the side; tangent angle
    psi(s) minus that of the large-scale shape, its spectrum S_psi ~ q^(1-2 zeta),
    and the integrated-slope height u(s) = int dpsi ds with S_u ~ q^-(1+2 zeta).
    Well defined through overhangs (as in fss_contours.py).
  * Cartesian (KFR Eq. 9 style, for comparison): x(y) or y(x) as the mean of the
    crossings of the side with the lattice lines. On these super-rough walls a
    sizeable fraction of lines cross a side more than once, and the projection
    then flattens S(q), exactly as the polar R(theta) did.

Construction, per bubble (contours.npz):
  * crossings of the wall contour with the lines y = j + 1/2 (or x = i + 1/2),
    split at the centroid into the two halves;
  * profile u(t) = mean of the crossings of that half on line t (KFR Eq. 9),
    multiplicity m(t) kept (m > 1: overhang);
  * the extent of each side comes from the large-scale shape: on the wall
    smoothed over max(64, P/20) along its length, the run of wall around
    the extreme point in that direction whose normal is within 45 deg of
    +x / -x (edge sides) or +y / -y (screw sides); the side is the set of lines
    inside that run's extent (at least MIN_SIDE lines; the rough wall itself
    is not required to have small slope);
  * S_u(q) of each side: quadratic detrend (removes the bubble's curvature),
    Hann window, S = |sum w u e^{-iqt}|^2 / sum w^2  -> q^-(1+2 zeta);
    also the slope spectrum S_du(q) of du/dt (same window), ~ q^(1-2 zeta).
Spectra are binned in log q and averaged over sides and samples; zeta from a
weighted fit, errors by bootstrap over samples. `--selftest` checks the
estimator on synthetic self-affine profiles with a curved mean shape.

    python sides.py --dir runs/fss_D0.04
"""
from __future__ import annotations

import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fss_contours import Q_EDGES, bin_log  # noqa: E402

QC = np.sqrt(Q_EDGES[1:] * Q_EDGES[:-1])
MIN_SIDE = 64
SIG_SIDE = 64.0
Q_MAX = 0.1
K_MIN = 6            # fits start at q = 2 pi K_MIN / (side length)
B = 300
RNG = np.random.default_rng(2024)
# robustness checks: SIDES_SUB_SIGMA = scale of the shape subtracted from psi
# (default: the same P/20 smoothing that defines the sides); SIDES_TAG names outputs
SUB_SIGMA = float(os.environ.get("SIDES_SUB_SIGMA", "0"))
TAG = os.environ.get("SIDES_TAG", "")


# ---------------------------------------------------------------------------
def crossings(xy, axis):
    """Crossings of the closed polygon with lines (axis=1: y = j+1/2 -> returns
    (j, x); axis=0: x = i+1/2 -> (i, y))."""
    if not np.allclose(xy[0], xy[-1]):
        xy = np.vstack([xy, xy[:1]])
    a = xy[:, axis]            # coordinate across the lines
    o = xy[:, 1 - axis]        # coordinate along the lines (the height)
    a0, a1 = a[:-1], a[1:]
    lo, hi = np.minimum(a0, a1), np.maximum(a0, a1)
    k0 = np.ceil(lo - 0.5 - 1e-9).astype(np.int64)
    k1 = np.floor(hi - 0.5 - 1e-9).astype(np.int64)
    cnt = np.maximum(k1 - k0 + 1, 0)
    cnt[hi - lo == 0] = 0
    seg = np.repeat(np.arange(len(a0)), cnt)
    k = np.repeat(k0, cnt) + (np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt))
    c = k + 0.5 + 1e-9
    f = (c - a0[seg]) / (a1[seg] - a0[seg])
    h = o[:-1][seg] + f * (o[1:][seg] - o[:-1][seg])
    return k, h


def spectra(u):
    """Hann-windowed, quadratically detrended S_u and S of the slope du/dt."""
    n = len(u)
    t = np.arange(n)
    r = u - np.polyval(np.polyfit(t, u, 2), t)
    w = np.hanning(n)
    q = 2 * np.pi * np.fft.rfftfreq(n)
    Su = np.abs(np.fft.rfft(r * w)) ** 2 / np.sum(w ** 2)
    du = np.gradient(u)
    du = du - np.polyval(np.polyfit(t, du, 1), t)
    Sd = np.abs(np.fft.rfft(du * w)) ** 2 / np.sum(w ** 2)
    return q[1:], Su[1:], Sd[1:]


def side_extent(sm, phn, direction):
    """Coordinate range (along the side) of the run of the smoothed wall whose
    outward normal is within 45 deg of `direction` (radians), around the extreme point."""
    dvec = np.array([np.cos(direction), np.sin(direction)])
    good = np.cos(phn - direction) > np.cos(np.pi / 4)
    i0 = int(np.argmax(sm @ dvec))
    if not good[i0]:
        return None
    n = len(good)
    lo, hi = i0, i0
    while good[(lo - 1) % n] and (i0 - lo) < n:
        lo -= 1
    while good[(hi + 1) % n] and (hi - i0) < n:
        hi += 1
    idx = np.arange(lo, hi + 1) % n
    along = 1 if abs(np.cos(direction)) > 0.5 else 0      # edge sides: extent in y; screw: in x
    return sm[idx, along].min(), sm[idx, along].max(), idx


def analyse(args):
    from fss_contours import resample_closed, smooth_closed, normals_angle
    xy, L, a, seed, R0 = args
    pts, P = resample_closed(xy)
    c0 = np.array([(L - 1) / 2, (L - 1) / 2])
    x, y = pts[:, 0], pts[:, 1]
    cr = x * np.roll(y, -1) - np.roll(x, -1) * y
    A = 0.5 * cr.sum()
    xc = np.sum((x + np.roll(x, -1)) * cr) / (6 * A)
    yc = np.sum((y + np.roll(y, -1)) * cr) / (6 * A)
    from fss_contours import tangent_angle
    sm = smooth_closed(pts, max(SIG_SIDE, P / 20), P)
    phn = normals_angle(sm, P)
    ds = P / len(pts)
    psi = tangent_angle(pts)[0]
    psi_sm = tangent_angle(smooth_closed(pts, SUB_SIGMA, P) if SUB_SIGMA > 0 else sm)[0]
    dpsi_all = psi - psi_sm
    dpsi_all -= 2 * np.pi * np.round(np.median(dpsi_all) / (2 * np.pi))
    res = dict(L=L, aniso=a, seed=seed)
    # edge sides: rows y = j+1/2, heights x (right: +x, left: -x)
    # screw sides: columns x = i+1/2, heights y (top: +y, bottom: -y)
    for kind, axis, cen_h, dirs in (("edge", 1, xc, (0.0, np.pi)), ("screw", 0, yc, (np.pi / 2, -np.pi / 2))):
        k, h = crossings(xy, axis)
        Su_b, Sd_b, lens, mults = [], [], [], []
        iSu_b, iSp_b, ilens, trav = [], [], [], []
        for dr in dirs:
            ext = side_extent(sm, phn, dr)
            if ext is None:
                continue
            # intrinsic: the arc itself
            idx = ext[2]
            # advance along the side normal over its middle (|transverse offset| < R0/2), from the
            # initial circle: |x - x0| - sqrt(R0^2 - (y - y0)^2) for edge sides, x <-> y for screw sides
            hcoord = 0 if kind == "edge" else 1
            dpt = pts[idx] - c0
            mid = np.abs(dpt[:, 1 - hcoord]) < R0 / 2
            if mid.any():
                trav.append(np.mean(np.abs(dpt[mid, hcoord]) - np.sqrt(R0 ** 2 - dpt[mid, 1 - hcoord] ** 2)))
            if len(idx) * ds >= MIN_SIDE:
                dp = dpsi_all[idx]
                n = len(dp)
                tt = np.arange(n)
                w = np.hanning(n)
                qi = 2 * np.pi * np.fft.rfftfreq(n, d=ds)
                r = dp - np.polyval(np.polyfit(tt, dp, 1), tt)
                Sp = ds * np.abs(np.fft.rfft(r * w)) ** 2 / np.sum(w ** 2)
                uu = np.cumsum(dp - dp.mean()) * ds
                r = uu - np.polyval(np.polyfit(tt, uu, 2), tt)
                Sui = ds * np.abs(np.fft.rfft(r * w)) ** 2 / np.sum(w ** 2)
                iSp_b.append(bin_log(qi[1:], Sp[1:])); iSu_b.append(bin_log(qi[1:], Sui[1:])); ilens.append(n * ds)
            sgn = 1 if (np.cos(dr) + np.sin(dr)) > 0 else -1
            m = ((h - cen_h) * sgn > 0) & (k + 0.5 >= ext[0]) & (k + 0.5 <= ext[1])
            if m.sum() < MIN_SIDE:
                continue
            kk, hh = k[m], h[m]
            kmin, kmax = kk.min(), kk.max()
            nl = kmax - kmin + 1
            ssum = np.bincount(kk - kmin, hh, nl)
            cnt = np.bincount(kk - kmin, minlength=nl)
            if nl < MIN_SIDE or (cnt == 0).mean() > 0.01:
                continue
            t = np.arange(nl)
            okc = cnt > 0
            u = np.interp(t, t[okc], ssum[okc] / cnt[okc])
            q, Su, Sd = spectra(u)
            Su_b.append(bin_log(q, Su)); Sd_b.append(bin_log(q, Sd))
            lens.append(nl); mults.append(np.mean(cnt > 1))
        res[f"{kind}_Su"] = np.nanmean(Su_b, 0) if Su_b else np.full(len(QC), np.nan)
        res[f"{kind}_Sd"] = np.nanmean(Sd_b, 0) if Sd_b else np.full(len(QC), np.nan)
        res[f"{kind}_len"] = np.mean(lens) if lens else np.nan
        res[f"{kind}_nsides"] = len(lens)
        res[f"{kind}_multi"] = np.mean(mults) if mults else np.nan
        res[f"{kind}_iSu"] = np.nanmean(iSu_b, 0) if iSu_b else np.full(len(QC), np.nan)
        res[f"{kind}_iSp"] = np.nanmean(iSp_b, 0) if iSp_b else np.full(len(QC), np.nan)
        res[f"{kind}_ilen"] = np.mean(ilens) if ilens else np.nan
        res[f"{kind}_travel"] = np.mean(trav) if trav else np.nan
    return res


# ---------------------------------------------------------------------------
def fit_zeta(X, qlo, qhi, kind="u"):
    """zeta and bootstrap error from per-sample binned spectra X (n, bins)."""
    X = X[np.all(np.isfinite(X[:, (QC >= qlo) & (QC <= qhi)]), axis=1)]
    ok = (QC >= qlo) & (QC <= qhi)
    if ok.sum() < 3 or len(X) < 4:
        return np.nan, np.nan, len(X)

    def one(Y):
        m = Y[:, ok].mean(0)
        sl = np.polyfit(np.log(QC[ok]), np.log(m), 1)[0]
        return -(1 + sl) / 2 if kind == "u" else (1 - sl) / 2
    z0 = one(X)
    zb = [one(X[RNG.integers(0, len(X), len(X))]) for _ in range(B)]
    return z0, np.std(zb, ddof=1), len(X)


def synth_profile(n, zeta, rng, curv=0.0):
    """Self-affine profile (S ~ q^-(1+2 zeta), cut-off at q ~ 0.3) cut out of a 4x longer periodic one."""
    N = 4 * n
    q = 2 * np.pi * np.fft.rfftfreq(N)
    S = np.zeros_like(q)
    S[1:] = q[1:] ** (-(1 + 2 * zeta)) * np.exp(-(q[1:] / 0.3) ** 2)
    c = np.sqrt(S / 2) * (rng.normal(size=len(q)) + 1j * rng.normal(size=len(q)))
    u = np.fft.irfft(c, N)[:n]
    t = np.arange(n) - n / 2
    return u / u.std() * 20 + curv * t ** 2


def selftest():
    rng = np.random.default_rng(3)
    print("synthetic sides (n = 1500, 200 profiles), with a parabolic mean shape u = 4e-4 t^2:")
    for zeta in (1.0, 1.25, 1.43):
        Su, Sd = [], []
        for _ in range(200):
            q, a, b = spectra(synth_profile(1500, zeta, rng, curv=4e-4))
            Su.append(bin_log(q, a)); Sd.append(bin_log(q, b))
        zu = fit_zeta(np.array(Su), 2 * np.pi * K_MIN / 1500, Q_MAX, "u")
        zd = fit_zeta(np.array(Sd), 2 * np.pi * K_MIN / 1500, Q_MAX, "d")
        print(f"  input zeta = {zeta}: from S_u {zu[0]:.3f} ± {zu[1]:.3f};  from slope spectrum {zd[0]:.3f} ± {zd[1]:.3f}")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=os.path.join(HERE, "runs", "fss_D0.04"))
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--lmin", type=int, default=1024)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    import fss_shapes
    C = fss_shapes.load_contours(args.dir)
    import glob
    import json
    R0 = {}
    for fn in glob.glob(os.path.join(args.dir, "results*.jsonl")):
        for line in open(fn):
            r = json.loads(line)
            R0[(r["L"], r["aniso"], r["seed"])] = r["R0"]
    jobs = [(C["contours"][i][C["main"][i]], int(C["L"][i]), float(C["aniso"][i]), int(C["seed"][i]),
             R0[(int(C["L"][i]), float(C["aniso"][i]), int(C["seed"][i]))])
            for i in range(len(C["L"])) if C["L"][i] >= args.lmin]
    print(f"cutting {len(jobs)} bubbles into sides ...")
    with ProcessPoolExecutor(args.workers) as ex:
        rows = list(ex.map(analyse, jobs, chunksize=4))
    S = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    np.savez_compressed(os.path.join(args.dir, f"sides_stats{TAG}.npz"), **S)
    plot(args.dir, S)


def plot(d, S):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm
    As = sorted(set(S["aniso"].tolist()))
    Ls = sorted(set(S["L"].tolist()))
    ca = {a: cm.viridis(i / max(1, len(As) - 1)) for i, a in enumerate(As)}
    q_ov = 2 * np.pi / 339.0
    out = []

    def sel(L, a):
        return (S["L"] == L) & (S["aniso"] == a)

    def zfits(kind, L, a):
        ln = np.nanmean(S[f"{kind}_ilen"][sel(L, a)])
        if not np.isfinite(ln):
            return None
        qlo = 2 * np.pi * K_MIN / ln
        r = dict(len=ln, nside=int(S[f"{kind}_nsides"][sel(L, a)].sum()))
        r["iuA"] = fit_zeta(S[f"{kind}_iSu"][sel(L, a)], qlo, Q_MAX, "u")
        r["iuB"] = fit_zeta(S[f"{kind}_iSu"][sel(L, a)], max(q_ov, qlo), Q_MAX, "u")
        r["ipB"] = fit_zeta(S[f"{kind}_iSp"][sel(L, a)], max(q_ov, qlo), Q_MAX, "d")
        lc = np.nanmean(S[f"{kind}_len"][sel(L, a)])
        r["cB"] = fit_zeta(S[f"{kind}_Su"][sel(L, a)], max(q_ov, 2 * np.pi * K_MIN / lc), Q_MAX, "u") if np.isfinite(lc) else (np.nan, np.nan, 0)
        r["multi"] = np.nanmean(S[f"{kind}_multi"][sel(L, a)])
        return r

    F = {(k, L, a): zfits(k, L, a) for k in ("edge", "screw") for L in Ls for a in As}
    fig, axs = plt.subplots(2, 4, figsize=(22, 10))
    Lm = Ls[-1]
    for col, kind in enumerate(("edge", "screw")):
        ax = axs[0, col]
        for a in As:
            X = S[f"{kind}_iSu"][sel(Lm, a)]
            X = X[np.any(np.isfinite(X), axis=1)]
            if len(X) < 4:
                continue
            m = np.nanmean(X, 0); e = np.nanstd(X, 0, ddof=1) / np.sqrt(len(X))
            ln = np.nanmean(S[f"{kind}_ilen"][sel(Lm, a)])
            ok = QC > 2 * np.pi * 3 / ln
            ax.errorbar(QC[ok], (m * QC ** 3.5)[ok], (e * QC ** 3.5)[ok], fmt="o-", ms=3, lw=0.9, color=ca[a],
                        label=f"a={a} ({len(X)} bubbles, arc ≈ {ln:.0f})")
        ax.axvspan(Q_MAX, 3, color="0.93", zorder=-1)
        ax.axvline(q_ov, color="C2", ls="--", lw=0.8)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("q (arc length)")
        ax.set_ylabel(r"$q^{3.5}\,S_u(q)$, intrinsic  (flat for ζ = 1.25)")
        ttl = "edge sides: normal ∥ b, line ⊥ b" if kind == "edge" else "screw sides: normal ⊥ b, line ∥ b"
        ax.set_title(f"{ttl}, L={Lm}", fontsize=10); ax.legend(fontsize=7)
    ax = axs[0, 2]
    for kind, mk, clr in (("edge", "o", "C0"), ("screw", "s", "C3")):
        for L, ls in zip(Ls, (":", "--", "-")):
            z = [F[(kind, L, a)]["iuB"][:2] if F[(kind, L, a)] else (np.nan, np.nan) for a in As]
            z, ze = np.array(z).T
            ax.errorbar(np.array(As) + (0.012 if kind == "screw" else 0), z, ze, fmt=mk + ls, color=clr, ms=4,
                        label=f"{kind}, L={L}")
    for kind, clr in (("edge", "C0"), ("screw", "C3")):
        z = [F[(kind, Lm, a)]["cB"][:2] if F[(kind, Lm, a)] else (np.nan, np.nan) for a in As]
        z, ze = np.array(z).T
        ax.errorbar(np.array(As) - 0.012, z, ze, fmt="x", color=clr, ms=6, label=f"{kind}, Cartesian x(y)/y(x), L={Lm}")
    ax.axhline(1.25, color="k", ls="--", lw=0.6, label="qEW 1.25")
    ax.axhline(2 - 1 / 1.76, color="k", ls=":", lw=0.6, label="2−1/1.76 = 1.43")
    ax.set_xlabel("anisotropy a"); ax.set_ylabel(r"ζ ($q_{ov}$ < q < 0.1)")
    ax.set_title("roughness of the two kinds of side"); ax.legend(fontsize=6, ncol=2)
    ax = axs[1, 0]
    for kind, mk, clr in (("edge", "o", "C0"), ("screw", "s", "C3")):
        for L, ls in zip(Ls, (":", "--", "-")):
            ax.plot(As, [np.nanmean(S[f"{kind}_ilen"][sel(L, a)]) / L for a in As], mk + ls, color=clr, label=f"{kind}, L={L}")
    ax.set_xlabel("anisotropy a"); ax.set_ylabel("side arc length / L"); ax.set_title("length of each kind of side")
    ax.legend(fontsize=7)
    ax = axs[1, 1]
    for kind, mk, clr in (("edge", "o", "C0"), ("screw", "s", "C3")):
        for L, ls in zip(Ls, (":", "--", "-")):
            ax.plot(As, [F[(kind, L, a)]["multi"] if F[(kind, L, a)] else np.nan for a in As], mk + ls, color=clr, label=f"{kind}, L={L}")
    ax.set_xlabel("anisotropy a"); ax.set_ylabel("fraction of lattice lines crossing the side > once")
    ax.set_title("why the Cartesian projection fails"); ax.legend(fontsize=7)
    ax = axs[1, 2]
    band = (QC > 0.03) & (QC < 0.1)
    for L, ls in zip(Ls, (":", "--", "-")):
        r = []
        for a in As:
            Xe = S["edge_iSp"][sel(L, a)]; Xs = S["screw_iSp"][sel(L, a)]
            me = np.nanmean(Xe[:, band], 0); ms = np.nanmean(Xs[:, band], 0)
            r.append(np.exp(np.nanmean(np.log(ms / me))))
        ax.plot(As, r, "o" + ls, color="k", label=f"L={L}")
        out.append(f"S_psi(screw)/S_psi(edge), 0.03<q<0.1, L={L}: " + "  ".join(f"a={a}: {v:.2f}" for a, v in zip(As, r)))
    ax.axhline(1, color="k", lw=0.5)
    ax.set_yscale("log"); ax.set_xlabel("anisotropy a"); ax.set_ylabel(r"$S_\psi$(screw) / $S_\psi$(edge), 0.03 < q < 0.1")
    ax.set_title("slope-fluctuation amplitude, screw over edge sides"); ax.legend(fontsize=7)
    ax = axs[0, 3]
    for kind, mk, clr in (("edge", "o", "C0"), ("screw", "s", "C3")):
        for L, ls in zip(Ls, (":", "--", "-")):
            ax.plot(As, [np.nanmean(S[f"{kind}_travel"][sel(L, a)]) for a in As], mk + ls, color=clr, label=f"{kind}, L={L}")
    ax.set_yscale("log"); ax.set_xlabel("anisotropy a"); ax.set_ylabel("advance of the side along its normal (middle part)")
    ax.set_title("how far each kind of side moved from the start"); ax.legend(fontsize=7)
    ax = axs[1, 3]
    for kind, mk, clr in (("edge", "o", "C0"), ("screw", "s", "C3")):
        for L in Ls:
            for a in As:
                r = F[(kind, L, a)]
                if r is None or not np.isfinite(r["iuB"][0]):
                    continue
                tv = np.nanmean(S[f"{kind}_travel"][sel(L, a)])
                ax.errorbar([tv], [r["iuB"][0]], [r["iuB"][1]], fmt=mk, color=clr, mfc=ca[a], ms=6 if L == Lm else 4,
                            alpha=1 if L == Lm else 0.6)
                out.append(f"travel {kind:5s} L={L:5d} a={a:4.2f}: d = {tv:7.1f}  zeta(B) = {r['iuB'][0]:.3f} ± {r['iuB'][1]:.3f}")
    ax.axhline(1.25, color="k", ls="--", lw=0.6)
    ax.set_xscale("log"); ax.set_xlabel("distance travelled by the side"); ax.set_ylabel(r"ζ ($q_{ov}$ < q < 0.1)")
    ax.set_title("ζ vs distance travelled (fill colour: a; big: L=4096)", fontsize=10)
    fig.tight_layout()
    fn = os.path.join(d, f"sides{TAG}.png")
    fig.savefig(fn, dpi=110)
    out.insert(0, "")
    hdr = ["zeta per kind of side (bootstrap over bubbles).",
           "  iu: intrinsic integrated-slope height S_u;  ip: intrinsic slope spectrum S_psi;  c: Cartesian mean-crossing profile",
           "  A: 2 pi 6/len < q < 0.1;  B: q_ov = 0.0185 < q < 0.1"]
    rows = []
    for kind in ("edge", "screw"):
        for L in Ls:
            for a in As:
                r = F[(kind, L, a)]
                if r is None:
                    continue
                rows.append(f"  {kind:5s} L={L:5d} a={a:4.2f}: sides={r['nside']:4d} arc={r['len']:6.0f}  "
                            f"iu A {r['iuA'][0]:.3f}±{r['iuA'][1]:.3f}  iu B {r['iuB'][0]:.3f}±{r['iuB'][1]:.3f}  "
                            f"ip B {r['ipB'][0]:.3f}±{r['ipB'][1]:.3f}  | c B {r['cB'][0]:.3f}±{r['cB'][1]:.3f} (multi {r['multi']:.3f})")
    txt = "\n".join(hdr + rows + out)
    print(txt)
    with open(os.path.join(d, f"sides{TAG}.txt"), "w") as f:
        f.write(txt + "\n")
    print("saved", fn)


if __name__ == "__main__":
    main()
