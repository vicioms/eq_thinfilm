#!/usr/bin/env python3
"""
Critical flat walls (fss_scan.py --geometry flat): extraction and analysis.

1. extract (where the scan ran, needs <dir>/configs):  python fss_flat.py --dir D --extract
   For each critical configuration: the wall = top boundary of the phi > 0
   domain attached to the bottom (islands ignored, holes left out), as the
   phi = 0 contour over exactly one period in x (contourpy on the field tiled
   3 times in x), stored as start point + int8 steps at 1/8 site; and the KFR
   Eq. 9 profile u(x) = mean of the crossings per column with its multiplicity.
   -> <dir>/walls.npz (tracked by git).

2. analyse (anywhere):  python fss_flat.py --dirs D1 D2 ... [--out D]
   Per orientation (burgers angle), anisotropy a and L:
     * intrinsic roughness: tangent angle psi(s) of the wall along its arc length
       (periodic over one period), S_psi(q) = P |c_n|^2 ~ q^(1 - 2 zeta), valid
       through overhangs and for zeta > 1 (as fss_contours.py for bubbles);
     * Cartesian KFR profile: S_u(q) of u(x) ~ q^-(1 + 2 zeta), for comparison;
     * spectra pooled at fixed (L, a): the modes n >= 2 of all walls in one array,
       binned in log q with bins merged until each holds >= one mode per wall
       (class Pooled); bin position = mean log q of its modes; fits weighted by
       sqrt(modes per bin), errors by bootstrap over walls with the same bins;
     * local zeta_eff(q) (sliding factor-2 window), zeta fits (modes n >= 2 up to
       --qmax, default 0.1, the noise core; the periodic walls show no low-mode
       depletion, zeta is the same from n >= 1 to n >= 8);
     * memory check: the same, split at the median distance travelled by the
       wall before it got pinned for the last time;
     * h_d: mean and var(h_d) ~ L^(-2/nu) (L >= lmin).
   Reference values: qEW zeta = 1.25, nu = 4/3; qKPZ / DPD zeta = 0.633, nu = 1.733
   (Tang, Kardar & Dhar, PRL 74, 920 (1995)).
   Also writes <out>/walls.png: example critical walls at the largest L per orientation and a,
   with the arc-length distribution of the wall's tangent angle (facets) for every L.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
SCALE = 8
CFG_RE = re.compile(r"L(\d+)_a([0-9.]+)_s(\d+)\.npz$")
Q_MAX = 0.1
N_MIN = 2                      # lowest mode fitted (n = 1, the box-scale mode, left out)
DS = 0.5
B = 300
RNG = np.random.default_rng(5)


# ---------------------------------------------------------------------------
# extraction
# ---------------------------------------------------------------------------
def wall_curve(phi):
    """One period of the wall line (n, 2) array (x, y), starting where it crosses x = L."""
    import contourpy
    import bubble_depinning as bd
    L = phi.shape[1]
    m = bd.wall_mask(phi)
    f = np.where(m, phi, -np.abs(phi))
    F = np.hstack([f, f, f])
    lines = contourpy.contour_generator(z=F).lines(0.0)
    if not lines:
        return None
    c = max(lines, key=lambda a: np.ptp(a[:, 0]))
    if np.ptp(c[:, 0]) < 2 * L:
        return None
    if c[0, 0] > c[-1, 0]:
        c = c[::-1]
    i1 = int(np.argmax(c[:, 0] >= L))
    target = c[i1] + np.array([L, 0.0])
    rest = np.arange(i1 + 1, len(c))
    i2 = rest[np.argmin(np.hypot(*(c[rest] - target).T))]
    seg = c[i1:i2] - np.array([L, 0.0])
    return seg


def _extract_one(fn):
    import bubble_depinning as bd
    L, a, seed = CFG_RE.search(fn).groups()
    z = np.load(fn)
    phi = z["phi_c"].astype(np.float64)
    seg = wall_curve(phi)
    an = bd.analyse_flat(phi)
    if seg is None or an is None:
        return None
    q = np.rint(seg * SCALE).astype(np.int32)
    st = np.diff(q, axis=0)
    if np.abs(st).max(initial=0) >= 128:
        return None
    return dict(L=int(L), aniso=float(a), seed=int(seed), h_d=float(z["h_d"]),
                start=q[0].astype(np.int32), steps=st.astype(np.int8),
                u=an["u"].astype(np.float32), mult=np.minimum(an["mult"], 127).astype(np.int8))


def extract(d, workers):
    files = sorted(glob.glob(os.path.join(d, "configs", "*.npz")))
    if not files:
        sys.exit(f"no configurations in {d}/configs")
    with ProcessPoolExecutor(workers) as ex:
        rows = [r for r in ex.map(_extract_one, files, chunksize=4) if r is not None]
    rows.sort(key=lambda r: (r["aniso"], r["L"], r["seed"]))
    out = {k: np.array([r[k] for r in rows]) for k in ("L", "aniso", "seed", "h_d")}
    out["start"] = np.array([r["start"] for r in rows])
    out["offsets"] = np.r_[0, np.cumsum([len(r["steps"]) for r in rows])]
    out["steps"] = np.concatenate([r["steps"] for r in rows]) if rows else np.zeros((0, 2), np.int8)
    for L in sorted({r["L"] for r in rows}):
        out[f"u_L{L}"] = np.stack([r["u"] for r in rows if r["L"] == L])
        out[f"mult_L{L}"] = np.stack([r["mult"] for r in rows if r["L"] == L])
    out["n_configs"] = len(files)
    out["scale"] = SCALE
    fn = os.path.join(d, "walls.npz")
    np.savez_compressed(fn, **out)
    print(f"saved {fn} ({os.path.getsize(fn) / 1e6:.1f} MB, {len(rows)} of {len(files)} walls)")


def load_walls(d):
    z = np.load(os.path.join(d, "walls.npz"))
    W = {k: z[k] for k in ("L", "aniso", "seed", "h_d")}
    st, off, start = z["steps"].astype(np.int64), z["offsets"], z["start"].astype(np.int64)
    W["curve"] = [np.vstack([start[i], start[i] + np.cumsum(st[off[i]:off[i + 1]], axis=0)]) / float(z["scale"])
                  for i in range(len(W["L"]))]
    W["u"], W["mult"] = [None] * len(W["L"]), [None] * len(W["L"])
    for L in np.unique(W["L"]):
        idx = np.nonzero(W["L"] == L)[0]
        for i, a, b in zip(idx, z[f"u_L{L}"], z[f"mult_L{L}"]):
            W["u"][i], W["mult"][i] = a, b
    meta = {}
    for fn in glob.glob(os.path.join(d, "results*.jsonl")):
        for line in open(fn):
            r = json.loads(line)
            meta[(r["L"], r["aniso"], r["seed"])] = r
    W["travel"] = np.array([meta.get((int(L), float(a), int(s)), {}).get("travel", np.nan)
                            for L, a, s in zip(W["L"], W["aniso"], W["seed"])])
    W["Ly"] = np.array([meta.get((int(L), float(a), int(s)), {}).get("Ly", 2 * int(L))
                        for L, a, s in zip(W["L"], W["aniso"], W["seed"])])
    angles = {r.get("burgers_angle") for r in meta.values()}
    W["burgers_angle"] = angles.pop() if len(angles) == 1 else np.nan
    return W


# ---------------------------------------------------------------------------
# configurations
# ---------------------------------------------------------------------------
def plot_walls(sets, outdir, n_show=4):
    """walls.png: example critical walls at the largest L (phi > 0 domain grey, red arrow = b),
    one row per (orientation, a), and the arc-length distribution of the tangent angle
    theta of the wall (0 = along x) for every L: facets show up as peaks at +- theta."""
    import matplotlib.pyplot as plt
    from matplotlib import cm
    rows = []
    for name, W in sets:
        for a in sorted({float(x) for x in W["aniso"]}):
            if a == 0 and any(r[2] == 0 for r in rows):      # isotropic: shown once
                continue
            rows.append((name, W, a))
    Ls = sorted({int(L) for _, W in sets for L in W["L"]})
    cL = {L: cm.plasma(i / max(1, len(Ls) - 1) * 0.85) for i, L in enumerate(Ls)}
    bins = np.linspace(-90, 90, 37)
    fig, axs = plt.subplots(len(rows), n_show + 1, figsize=(2.3 * n_show + 4.2, 4.2 * len(rows)),
                            squeeze=False, gridspec_kw=dict(width_ratios=[1] * n_show + [1.9]))
    for r, (name, W, a) in enumerate(rows):
        sel = W["aniso"] == a
        Lm = int(W["L"][sel].max())
        idx = np.nonzero(sel & (W["L"] == Lm))[0][:n_show]
        ang = np.deg2rad(W["burgers_angle"])
        for c in range(n_show):
            ax = axs[r, c]
            ax.set_xticks([]); ax.set_yticks([])
            if c >= len(idx):
                ax.axis("off"); continue
            i = idx[c]
            cv, Ly = W["curve"][i], int(W["Ly"][i])
            # the period runs from x = 0 to x = L (wall_curve starts at the x = L crossing)
            x = np.r_[cv[:, 0], Lm, Lm, 0, 0]
            y = np.r_[cv[:, 1], cv[0, 1], 0, 0, cv[0, 1]]
            ax.fill(x, y, color="0.55", lw=0)
            ax.plot(cv[:, 0], cv[:, 1], color=cm.viridis(a) if a < 1 else "goldenrod", lw=0.6)
            ax.set_xlim(0, Lm); ax.set_ylim(0, Ly); ax.set_aspect("equal")
            ax.annotate("", xy=(0.12 * Lm + 0.1 * Lm * np.cos(ang), 0.06 * Ly + 0.1 * Lm * np.sin(ang)),
                        xytext=(0.12 * Lm, 0.06 * Ly), arrowprops=dict(arrowstyle="->", color="red", lw=1.5))
            ax.set_title(f"h_d={W['h_d'][i]:.4f}", fontsize=8)
            if c == 0:
                ax.set_ylabel(f"{name.split(' (')[0]}, a = {a:g}\nL = {Lm}", fontsize=10)
        ax = axs[r, n_show]
        for L in Ls:
            ii = np.nonzero(sel & (W["L"] == L))[0]
            if len(ii) == 0:
                continue
            th = []
            for i in ii:
                # chords over 4 sites of arc length, so the contour's lattice steps do not dominate
                c2, _ = periodic_resample(W["curve"][i], L, ds=0.5)
                d = np.roll(c2, -8, axis=0) - c2
                d[-8:, 0] += L
                th.append(np.degrees(np.arctan(d[:, 1] / np.where(d[:, 0] == 0, 1e-12, d[:, 0]))))
            h, _ = np.histogram(np.concatenate(th), bins, density=True)
            ax.plot(0.5 * (bins[1:] + bins[:-1]), h, color=cL[L], label=f"L={L} (n={len(ii)})")
        ax.set_xlabel("tangent angle θ (deg, 0 = along x)"); ax.set_ylabel("P(θ)")
        ax.set_xlim(-90, 90); ax.set_xticks([-90, -45, 0, 45, 90]); ax.legend(fontsize=7)
        ax.set_title("arc-length distribution of θ", fontsize=9)
    fig.suptitle("critical flat walls (grey: phi > 0 domain attached to the bottom; red arrow: Burgers vector b)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fn = os.path.join(outdir, "walls.png")
    fig.savefig(fn, dpi=110)
    plt.close(fig)
    print("saved", fn)


# ---------------------------------------------------------------------------
# spectra
# ---------------------------------------------------------------------------
def periodic_resample(curve, L, ds=DS):
    """Uniform arc-length resampling of one period of a wall that repeats with shift (L, 0)."""
    c = np.vstack([curve, curve[:1] + np.array([L, 0.0])])
    seg = np.hypot(*np.diff(c, axis=0).T)
    s = np.r_[0, np.cumsum(seg)]
    P = s[-1]
    n = max(32, int(round(P / ds)))
    t = np.arange(n) * (P / n)
    return np.c_[np.interp(t, s, c[:, 0]), np.interp(t, s, c[:, 1])], P


def wall_spectra(curve, u, L, m=2):
    """Raw spectra of one wall, modes n >= 1 as (q, S, n): the tangent angle psi(s) along
    the arc length, S = P |c_n|^2, and the KFR profile u(x), S = L |c_n|^2."""
    pts, P = periodic_resample(curve, L)
    n = len(pts)
    ext = np.vstack([pts[-m:] - [L, 0], pts, pts[:m] + [L, 0]])
    d = ext[2 * m:] - ext[:-2 * m]
    psi = np.unwrap(np.arctan2(d[:, 1], d[:, 0]))
    c = np.fft.rfft(psi - psi.mean()) / n
    k = np.arange(1, len(c))
    sp = (2 * np.pi * k / P, P * np.abs(c[1:]) ** 2, k)
    ok = np.isfinite(u)
    uu = np.interp(np.arange(L), np.nonzero(ok)[0], u[ok], period=L) if ok.sum() > 1 else np.zeros(L)
    cu = np.fft.rfft(uu - uu.mean()) / L
    ku = np.arange(1, len(cu))
    su = (2 * np.pi * ku / L, L * np.abs(cu[1:]) ** 2, ku)
    return sp, su, P


FINE = np.geomspace(1e-3, 3.0, 241)


class Pooled:
    """Spectrum of a set of walls at fixed L (and a): the modes n >= N_MIN of all walls in one
    array, binned in log q. Adjacent fine bins are merged until each holds at least as many
    modes as there are walls (about one mode per wall), so at low q, where the mode spacing
    2 pi / P exceeds the bin width, no bin is left to the few walls whose P happens to put a
    mode there; at high q the bins stay fine and hold many modes. Bin position = mean log q
    of its modes. Errors: bootstrap over walls with the same bins."""

    def __init__(self, spectra):
        self.parts = [(q[k >= N_MIN], S[k >= N_MIN]) for q, S, k in spectra]
        self.n = len(self.parts)
        if self.n == 0:
            self.qc = np.zeros(0); return
        q = np.concatenate([p[0] for p in self.parts])
        cnt = np.histogram(q, FINE)[0]
        edges, acc = [FINE[0]], 0
        for i, ci in enumerate(cnt):
            acc += ci
            if acc >= max(self.n, 8):
                edges.append(FINE[i + 1]); acc = 0
        if acc and len(edges) > 2:
            edges[-1] = FINE[-1]
        self.edges = np.array(edges)
        i = np.digitize(q, self.edges) - 1
        ok = (i >= 0) & (i < len(self.edges) - 1)
        with np.errstate(invalid="ignore"):
            self.qc = np.exp(np.bincount(i[ok], np.log(q[ok]), len(self.edges) - 1)
                             / np.bincount(i[ok], minlength=len(self.edges) - 1))

    def mean(self, idx=None):
        parts = self.parts if idx is None else [self.parts[j] for j in idx]
        q = np.concatenate([p[0] for p in parts]); S = np.concatenate([p[1] for p in parts])
        i = np.digitize(q, self.edges) - 1
        ok = (i >= 0) & (i < len(self.qc))
        c = np.bincount(i[ok], minlength=len(self.qc))
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.bincount(i[ok], S[ok], len(self.qc)) / c, c

    def boot(self, nb):
        return [self.mean(RNG.integers(0, self.n, self.n)) for _ in range(nb)]


def _slope(qc, S, c, ok):
    ok = ok & (c > 0) & np.isfinite(S) & (S > 0)
    if ok.sum() < 3:
        return np.nan
    return np.polyfit(np.log(qc[ok]), np.log(S[ok]), 1, w=np.sqrt(c[ok]))[0]   # log-mean error ~ 1/sqrt(count)


def _zeta(sl, kind):
    return (1 - sl) / 2 if kind == "psi" else -(1 + sl) / 2


def fit(pool, qhi, kind):
    """zeta from a weighted power-law fit of the pooled spectrum, modes n >= N_MIN up to qhi."""
    if pool.n < 4:
        return np.nan, np.nan
    win = pool.qc <= qhi
    z0 = _zeta(_slope(pool.qc, *pool.mean(), win), kind)
    zb = [_zeta(_slope(pool.qc, S, c, win), kind) for S, c in pool.boot(B)]
    return z0, np.nanstd(zb, ddof=1)


def local_zeta(pool, kind, half=np.log(2) / 2, nb=100):
    """zeta_eff(q): weighted fits over a sliding window of a factor 2 in q."""
    lq = np.log(pool.qc)

    def one(S, c):
        return np.array([_zeta(_slope(pool.qc, S, c, np.abs(lq - l) <= half), "psi") for l in lq])
    return one(*pool.mean()), np.nanstd([one(S, c) for S, c in pool.boot(nb)], 0)


# ---------------------------------------------------------------------------
def analyse(dirs, outdir, lmin):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm
    os.makedirs(outdir, exist_ok=True)
    sets = []
    for d in dirs:
        W = load_walls(d)
        ang = W["burgers_angle"]
        name = {0.0: "screw wall (line ∥ b)", 90.0: "edge wall (line ⊥ b)"}.get(ang, f"b at {ang}°")
        spec = [wall_spectra(c, u, int(L)) for c, u, L in zip(W["curve"], W["u"], W["L"])]
        W["sp"] = [s[0] for s in spec]; W["su"] = [s[1] for s in spec]
        W["P"] = np.array([s[2] for s in spec])
        W["multi"] = np.array([np.mean(m > 1) for m in W["mult"]])
        sets.append((name, W))
    plot_walls(sets, outdir)
    As = sorted({float(a) for _, W in sets for a in W["aniso"]})
    Ls = sorted({int(L) for _, W in sets for L in W["L"]})
    cL = {L: cm.plasma(i / max(1, len(Ls) - 1) * 0.85) for i, L in enumerate(Ls)}
    out = []
    fig, axs = plt.subplots(len(sets), 4, figsize=(21, 4.8 * len(sets)), squeeze=False)
    for row, (name, W) in enumerate(sets):
        sel = lambda L, a: (W["L"] == L) & (W["aniso"] == a)
        pools = {}

        def pool(L, a, kind="psi", mask=None):
            m = sel(L, a) if mask is None else sel(L, a) & mask
            key = (L, a, kind, None if mask is None else m.tobytes())
            if key not in pools:
                pools[key] = Pooled([W["sp" if kind == "psi" else "su"][i] for i in np.nonzero(m)[0]])
            return pools[key]
        out.append(f"== {name}")
        # compensated intrinsic spectra at the largest L, all a
        ax = axs[row, 0]
        Lm = max(int(L) for L in W["L"])
        for a in As:
            pl = pool(Lm, a)
            if pl.n < 2:
                continue
            m = pl.mean()[0]
            e = np.nanstd([S for S, _ in pl.boot(100)], 0)
            ax.errorbar(pl.qc, m * pl.qc ** 1.5, e * pl.qc ** 1.5, fmt="o-", ms=3, lw=0.8,
                        color=cm.viridis(As.index(a) / max(1, len(As) - 1)), label=f"a={a} (n={pl.n})")
        ax.axvspan(Q_MAX, 2, color="0.93", zorder=-1)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("q (arc length)")
        ax.set_ylabel(r"$q^{1.5}S_\psi$ (flat: ζ=1.25; slope +0.6·2: ζ=0.63)")
        ax.set_title(f"{name}, L={Lm}", fontsize=10); ax.legend(fontsize=7)
        # local zeta at the largest a, all L
        ax = axs[row, 1]
        for a, mk in ((As[0], "o"), (As[-1], "s")):
            for L in Ls:
                pl = pool(L, a)
                if pl.n < 4:
                    continue
                z0, ze = local_zeta(pl, "psi")
                ok = pl.qc < 0.3
                ax.errorbar(pl.qc[ok], z0[ok], ze[ok], fmt=mk + ("-" if a == As[-1] else ":"), color=cL[L],
                            ms=3, lw=0.8, label=f"a={a}, L={L}")
        for zz, ls, lab in ((1.25, "--", "qEW"), (0.633, "-.", "qKPZ")):
            ax.axhline(zz, color="k", ls=ls, lw=0.7, label=lab)
        ax.axvspan(Q_MAX, 0.3, color="0.93", zorder=-1)
        ax.set_xscale("log"); ax.set_ylim(0, 2); ax.set_xlabel("q"); ax.set_ylabel(r"local $\zeta_{eff}$ (intrinsic)")
        ax.set_title("local exponent", fontsize=10); ax.legend(fontsize=6, ncol=2)
        # zeta vs a: intrinsic, Cartesian, and split by travel
        ax = axs[row, 2]
        out.append(f"  zeta (window 2 pi 2/P .. {Q_MAX:g}): intrinsic S_psi | Cartesian u(x) | intrinsic, travel < median / > median")
        for L, ls in zip(Ls, (":", "-.", "--", "-")):
            zi, zc = [], []
            for a in As:
                pl = pool(L, a)
                if pl.n < 4:
                    zi.append((np.nan, np.nan)); zc.append((np.nan, np.nan)); continue
                Pm = W["P"][sel(L, a)].mean()
                z1 = fit(pl, Q_MAX, "psi")
                z2 = fit(pool(L, a, "u"), Q_MAX, "u")
                tr = W["travel"]
                med = np.nanmedian(tr[sel(L, a)])
                zlo = fit(pool(L, a, mask=tr < med), Q_MAX, "psi")
                zhi = fit(pool(L, a, mask=tr >= med), Q_MAX, "psi")
                zi.append(z1); zc.append(z2)
                out.append(f"  L={L:5d} a={a:4.2f} n={pl.n:4d}: {z1[0]:.3f}±{z1[1]:.3f} | {z2[0]:.3f}±{z2[1]:.3f} "
                           f"(multi {np.mean(W['multi'][sel(L, a)]):.3f}) | travel<{med:.0f}: {zlo[0]:.3f}±{zlo[1]:.3f}"
                           f"  >: {zhi[0]:.3f}±{zhi[1]:.3f}   P/L={Pm / L:.3f}")
            z, ze = np.array(zi).T
            ax.errorbar(As, z, ze, fmt="o" + ls, color=cL[L], label=f"intrinsic, L={L}")
            z, ze = np.array(zc).T
            ax.plot(As, z, "x" + ls, color=cL[L], alpha=0.6)
        for zz, ls in ((1.25, "--"), (0.633, "-.")):
            ax.axhline(zz, color="k", ls=ls, lw=0.7)
        ax.set_xlabel("anisotropy a"); ax.set_ylabel("ζ (o intrinsic, x Cartesian KFR Eq. 9)")
        ax.set_title("roughness exponent", fontsize=10); ax.legend(fontsize=7)
        # h_d and nu
        ax = axs[row, 3]
        out.append("  h_d: mean, var -> nu (L >= %d)" % lmin)
        for a in As:
            Lx = [L for L in Ls if sel(L, a).sum() > 3]
            v = np.array([np.var(W["h_d"][sel(L, a)], ddof=1) for L in Lx])
            ax.plot(Lx, v, "o-", color=cm.viridis(As.index(a) / max(1, len(As) - 1)), label=f"a={a}")
            Lf = [L for L in Lx if L >= lmin]
            nu = nue = np.nan
            if len(Lf) >= 2:
                vf = np.array([np.var(W["h_d"][sel(L, a)], ddof=1) for L in Lf])
                nu = -2 / np.polyfit(np.log(Lf), np.log(vf), 1)[0]
                nb = []
                for _ in range(B):
                    vb = [np.var(W["h_d"][sel(L, a)][RNG.integers(0, sel(L, a).sum(), sel(L, a).sum())], ddof=1) for L in Lf]
                    if min(vb) > 0:            # a resample of identical h_d has no variance
                        nb.append(-2 / np.polyfit(np.log(Lf), np.log(vb), 1)[0])
                nue = np.std(nb)
            out.append(f"  a={a:4.2f}: " + "  ".join(f"L={L}: <h_d>={W['h_d'][sel(L, a)].mean():.5f}" for L in Lx)
                       + f"   nu = {nu:.2f} ± {nue:.2f}")
        Lx = np.array(Ls, float)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("L"); ax.set_ylabel("var(h_d)")
        ax.set_title("threshold fluctuations (ν: qEW 4/3, qKPZ 1.73)", fontsize=10); ax.legend(fontsize=7)
    fig.tight_layout()
    fn = os.path.join(outdir, "flat_walls.png")
    fig.savefig(fn, dpi=110)
    txt = "\n".join(out)
    print(txt)
    with open(os.path.join(outdir, "flat_walls.txt"), "w") as f:
        f.write(txt + "\n")
    print("saved", fn)


def main():
    global Q_MAX
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", help="run directory to extract")
    ap.add_argument("--extract", action="store_true")
    ap.add_argument("--dirs", nargs="+", help="run directories to analyse (one per orientation)")
    ap.add_argument("--out", help="where to write the figure (default: first of --dirs)")
    ap.add_argument("--lmin", type=int, default=256)
    ap.add_argument("--qmax", type=float, default=Q_MAX,
                    help="upper end of the zeta fit window (default 0.1; e.g. 0.03 keeps only the "
                         "scales beyond the qEW -> qKPZ crossover of the edge walls)")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    args = ap.parse_args()
    Q_MAX = args.qmax
    if args.extract:
        extract(args.dir, args.workers)
    if args.dirs:
        analyse(args.dirs, args.out or args.dirs[0], args.lmin)


if __name__ == "__main__":
    main()
