#!/usr/bin/env python3
"""
Shape and roughness of the critical bubbles of an FSS scan (fss_scan.py).

Two stages:

1. extract  (needs <dir>/configs/*.npz, i.e. run it where the scan ran)
   For every critical configuration: rays from the bubble centre give R(theta)
   on n_theta(L) angles (>= 2 rays per lattice unit of wall) (mean of the phi = 0 crossings) and the number of crossings
   per ray (> 1 = overhang). A few bubble masks per (a, L) are kept for display.
   Everything goes to <dir>/shapes.npz (a few MB; tracked by git, unlike the
   configs), so the plots can be made on a laptop after `git pull`.

2. plot     (needs only shapes.npz)
   <dir>/bubbles.png    critical bubbles for every a at the largest L,
                        mean shape <R(theta)>/R_eff and overhang probability vs theta
   <dir>/Sq.png         interface structure factor S(q) of the whole bubble,
                        collapse over L, roughness exponent zeta vs a, and the
                        arc-resolved S(q) of the edge (normal || b) and screw
                        (normal _|_ b) quadrants of the wall

   Definitions: delta R(theta) = R(theta) - <R>_theta, arc length s = <R> theta,
   l = 2 pi <R>.  c_n = (1/N) sum_theta delta R e^{-i n theta},
   S(q_n) = l |c_n|^2 with q_n = n / <R>  ->  S(q) ~ q^-(1 + 2 zeta).
   Mode n = 1 (translation) is dropped; n = 2 is the elliptic deformation.

   python fss_shapes.py --dir runs/fss_cpu              # extract if needed, then plot
   python fss_shapes.py --dir runs/fss_cpu --extract    # force re-extraction
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
N_THETA_MIN = 1024


def n_theta(L):
    """Rays per bubble: >= 2 per lattice unit of wall at the typical critical
    radius R ~ L/4 (2 pi L/4 * 2 = pi L), rounded up to a power of 2."""
    return max(N_THETA_MIN, 2 ** int(np.ceil(np.log2(np.pi * L))))


N_MASKS = 6                # masks stored per (a, L)
MASK_MAX = 1024            # stored masks are subsampled to at most this size
ZETA_REF = {"qEW (1.25)": 1.25, "letter (1.43)": 2 - 1 / 1.76}
CFG_RE = re.compile(r"L(\d+)_a([0-9.]+)_s(\d+)\.npz$")


# ---------------------------------------------------------------------------
# 1. extraction
# ---------------------------------------------------------------------------
def _extract_one(fn):
    sys.path.insert(0, HERE)
    import bubble_depinning as bd
    L, a, seed = CFG_RE.search(fn).groups()
    L, a, seed = int(L), float(a), int(seed)
    z = np.load(fn)
    phi = z["phi_c"].astype(np.float64)
    an = bd.analyse(phi, n_theta=n_theta(L))
    if an is None:
        return None
    out = dict(L=L, aniso=a, seed=seed, h_d=float(z["h_d"]), R_eff=an["R_eff"],
               xc=an["xc"], yc=an["yc"], aspect=an["aspect"], orientation_deg=an["orientation_deg"],
               R_theta=an["R_theta"].astype(np.float32),
               mult=np.minimum(an["mult"], 127).astype(np.int8))
    if seed < N_MASKS:
        st = max(1, L // MASK_MAX)
        out["mask"] = np.packbits(bd.bubble_mask(phi)[::st, ::st])
    return out


def extract(d, workers):
    files = sorted(glob.glob(os.path.join(d, "configs", "*.npz")))
    if not files:
        sys.exit(f"no configurations in {d}/configs: run the extraction where the scan ran")
    print(f"extracting {len(files)} configurations ...")
    with ProcessPoolExecutor(workers) as ex:
        rows = [r for r in ex.map(_extract_one, files, chunksize=8) if r is not None]
    rows.sort(key=lambda r: (r["aniso"], r["L"], r["seed"]))
    out = {k: np.array([r[k] for r in rows]) for k in
           ("L", "aniso", "seed", "h_d", "R_eff", "xc", "yc", "aspect", "orientation_deg")}
    for L in sorted({r["L"] for r in rows}):          # ray count depends on L
        out[f"R_theta_L{L}"] = np.stack([r["R_theta"] for r in rows if r["L"] == L])
        out[f"mult_L{L}"] = np.stack([r["mult"] for r in rows if r["L"] == L])
    out["n_configs"] = len(files)
    for r in rows:
        if "mask" in r:
            out[f"mask_L{r['L']}_a{r['aniso']}_s{r['seed']}"] = r["mask"]
    fn = os.path.join(d, "shapes.npz")
    np.savez_compressed(fn, **out)
    print(f"saved {fn} ({os.path.getsize(fn) / 1e6:.1f} MB, {len(rows)} bubbles)")


def load_shapes(d):
    z = np.load(os.path.join(d, "shapes.npz"))
    S = {k: z[k] for k in z.files if not k.startswith(("mask_", "R_theta", "mult"))}
    # per-sample ray arrays as object arrays (their length depends on L);
    # old files store one (n_samples, 1024) array
    n = len(S["L"])
    for key in ("R_theta", "mult"):
        obj = np.empty(n, dtype=object)
        if key in z.files:
            for i, row in enumerate(z[key]):
                obj[i] = row
        else:
            for L in np.unique(S["L"]):
                idx = np.nonzero(S["L"] == L)[0]
                for i, row in zip(idx, z[f"{key}_L{L}"]):
                    obj[i] = row
        S[key] = obj
    masks = {}
    for k in z.files:
        if k.startswith("mask_"):
            L, a, s = re.match(r"mask_L(\d+)_a([0-9.]+)_s(\d+)", k).groups()
            n = int(L) // max(1, int(L) // MASK_MAX)
            masks[(int(L), float(a), int(s))] = np.unpackbits(z[k])[: n ** 2].reshape(n, n).astype(bool)
    return S, masks


# ---------------------------------------------------------------------------
# 2. structure factors
# ---------------------------------------------------------------------------
class Rays:
    """Per-sample rays; R[i] is one sample's array, R[mask] stacks the selected
    samples (which must share L, i.e. the ray count) into a 2-D array."""
    def __init__(self, obj):
        self.obj = obj

    def __getitem__(self, k):
        x = self.obj[k]
        return np.stack(list(x)) if isinstance(x, np.ndarray) and x.dtype == object else x


def th_of(n):
    return np.linspace(0, 2 * np.pi, n, endpoint=False)


def fill_gaps(R):
    """Rays without a crossing (R = nan, rare) -> periodic linear interpolation."""
    R = [row.astype(float) for row in R]
    for row in R:
        bad = np.isnan(row)
        if bad.any() and not bad.all():
            th = np.arange(len(row))
            row[bad] = np.interp(th[bad], th[~bad], row[~bad], period=len(row))
    obj = np.empty(len(R), dtype=object)
    obj[:] = R
    return obj


def bubble_Sq(R):
    """R: (n_samples, n_theta). Returns q (n_samples, n_modes), S (same)."""
    Rm = R.mean(1, keepdims=True)
    c = np.fft.rfft(R - Rm, axis=1) / R.shape[1]
    S = 2 * np.pi * Rm * np.abs(c) ** 2
    n = np.arange(c.shape[1])
    q = n[None, :] / Rm
    return q[:, 2:], S[:, 2:]


def low_modes(R, nmax=2):
    """Per-sample Fourier modes n <= nmax of R(theta): mean radius, centre
    offset and the elliptic deformation, i.e. the bubble's overall shape."""
    c = np.fft.rfft(R, axis=1)
    c[:, nmax + 1:] = 0
    return np.fft.irfft(c, n=R.shape[1], axis=1)


def arc_Sq(R, centres_deg, half_width_deg=30.0, nmax_shape=2):
    """Periodogram of the wall roughness on arcs centred on the given angles.
    Each bubble's own shape (Fourier modes n <= nmax_shape of R(theta)) is
    subtracted first, so that the curvature of an elongated bubble does not
    leak into the arcs near its tips; then linear detrending and a Hann window
    on each arc. q is per unit arc length (local mean radius of each arc).
    The removed modes have q <= nmax_shape / R, below the arcs' lowest mode."""
    nth = R.shape[1]
    dR = R - low_modes(R, nmax_shape)
    hw = int(round(half_width_deg / 360 * nth))
    j = np.arange(-hw, hw + 1)
    win = np.hanning(len(j))
    S_all, q_all = [], []
    for c in centres_deg:
        idx = (int(round(c / 360 * nth)) + j) % nth
        seg = dR[:, idx]
        Rloc = R[:, idx].mean(1, keepdims=True)
        ds = Rloc * 2 * np.pi / nth
        # linear detrending (per sample)
        A = np.c_[np.ones_like(j), j].astype(float)
        coef = np.linalg.lstsq(A, seg.T, rcond=None)[0]
        dev = seg - (A @ coef).T
        F = np.fft.rfft(dev * win[None, :], axis=1) * ds
        S = np.abs(F) ** 2 / (np.sum(win ** 2) * ds)
        m = np.arange(F.shape[1])
        q = 2 * np.pi * m[None, :] / (len(j) * ds)
        S_all.append(S[:, 1:]); q_all.append(q[:, 1:])
    return np.concatenate(q_all), np.concatenate(S_all)


def binned_mean(q, S, bins):
    q, S = q.ravel(), S.ravel()
    k = np.digitize(q, bins) - 1
    qm = np.full(len(bins) - 1, np.nan); Sm = qm.copy(); Se = qm.copy()
    for b in range(len(bins) - 1):
        sel = k == b
        if sel.sum() >= 3:
            qm[b] = np.exp(np.log(q[sel]).mean())
            Sm[b] = S[sel].mean()
            Se[b] = S[sel].std(ddof=1) / np.sqrt(sel.sum())
    ok = ~np.isnan(qm)
    return qm[ok], Sm[ok], Se[ok]


def fit_zeta(q, S, qmin, qmax):
    sel = (q > qmin) & (q < qmax)
    if sel.sum() < 3:
        return np.nan
    sl = np.polyfit(np.log(q[sel]), np.log(S[sel]), 1)[0]
    return (-sl - 1) / 2


# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------
def plot(d, w):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm

    S, masks = load_shapes(d)
    R = Rays(fill_gaps(S["R_theta"]))
    MULT = Rays(S["mult"])
    Ls = sorted(set(S["L"].tolist())); As = sorted(set(S["aniso"].tolist()))
    Lmax = Ls[-1]
    tb = 0.0                                          # burgers angle (fss_scan default)
    import json
    rows = [json.loads(l) for fn in glob.glob(os.path.join(d, "results*.jsonl")) for l in open(fn)]
    if rows:
        tb = rows[0].get("burgers_angle", 0.0)
    acol = {a: cm.viridis(i / max(1, len(As) - 1)) for i, a in enumerate(As)}
    lcol = {L: cm.plasma(i / max(1, len(Ls) - 1) * 0.85) for i, L in enumerate(Ls)}
    sel = lambda a=None, L=None: (((S["aniso"] == a) if a is not None else True)
                                  & ((S["L"] == L) if L is not None else True))

    # ---------------- bubbles ---------------------------------------------
    nshow = min(N_MASKS, 5)
    fig = plt.figure(figsize=(2.3 * (nshow + 2), 2.3 * len(As) + 0.6))
    gs = fig.add_gridspec(len(As), nshow + 2)
    for i, a in enumerate(As):
        pmax = 0.0
        for s in range(nshow):
            ax = fig.add_subplot(gs[i, s])
            m = masks.get((Lmax, a, s))
            if m is None:
                ax.axis("off"); continue
            ax.imshow(m, cmap="Greys", origin="lower", interpolation="nearest", vmin=0, vmax=1.6,
                      extent=(-0.5, Lmax - 0.5, -0.5, Lmax - 0.5))          # lattice units
            k = np.nonzero(sel(a, Lmax) & (S["seed"] == s))[0]
            if len(k):
                k = k[0]
                th = th_of(len(R[k]))
                ax.plot(S["xc"][k] + R[k] * np.cos(th), S["yc"][k] + R[k] * np.sin(th), "-",
                        color=acol[a], lw=0.6)
                ax.set_title(f"h_d={S['h_d'][k]:.4f}", fontsize=7)
            ax.set_xticks([]); ax.set_yticks([])
            if s == 0:
                ax.set_ylabel(f"a = {a}", fontsize=10)
            L0 = Lmax
            ax.annotate("", xy=(0.12 * L0 + 0.12 * L0 * np.cos(np.radians(tb)), 0.1 * L0 + 0.12 * L0 * np.sin(np.radians(tb))),
                        xytext=(0.12 * L0, 0.1 * L0), arrowprops=dict(arrowstyle="->", color="C3", lw=1.2))
        # mean shape (b fixed -> no alignment needed)
        ax = fig.add_subplot(gs[i, nshow], projection="polar")
        for L in Ls:
            k = sel(a, L)
            if k.sum() == 0:
                continue
            r = (R[k] / S["R_eff"][k, None]).mean(0)
            th = th_of(len(r))
            ax.plot(np.r_[th, th[:1]], np.r_[r, r[:1]], color=lcol[L], lw=1, label=f"L={L}")
        ax.plot([np.radians(tb)] * 2, [0, 1.3], color="C3", lw=1)
        ax.set_ylim(0, 1.35); ax.set_yticks([0.5, 1]); ax.tick_params(labelsize=6)
        if i == 0:
            ax.set_title(r"$\langle R(\theta)\rangle/R_{eff}$", fontsize=9)
            ax.legend(fontsize=6, loc="lower left", bbox_to_anchor=(-0.35, -0.25))
        # overhang probability vs theta
        ax = fig.add_subplot(gs[i, nshow + 1], projection="polar")
        for L in Ls:
            k = sel(a, L)
            if k.sum() == 0:
                continue
            p = (MULT[k] > 1).mean(0)
            th = th_of(len(p))
            nb = 32
            pmax = max(pmax, p.reshape(nb, -1).mean(1).max())
            p = p.reshape(nb, -1).mean(1); t = th.reshape(nb, -1).mean(1)
            ax.plot(np.r_[t, t[:1]], np.r_[p, p[:1]], color=lcol[L], lw=1)
        top = np.ceil(pmax * 20) / 20
        ax.set_ylim(0, top); ax.set_yticks([top / 2, top]); ax.set_rlabel_position(100)
        ax.plot([np.radians(tb)] * 2, [0, top], color="C3", lw=1)
        ax.tick_params(labelsize=6)
        if i == 0:
            ax.set_title("P(overhang | θ)", fontsize=9)
    fig.suptitle(f"critical bubbles, L = {Lmax}  (red arrow / line: Burgers vector b; "
                 f"contour = mean crossing radius R(θ))", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(d, "bubbles.png"), dpi=130)
    plt.close(fig)

    # ---------------- S(q) ------------------------------------------------
    qbins = np.geomspace(0.01, np.pi, 40)
    qmax_fit = 1 / w                                  # below the noise correlation length w
    qmin_of = lambda Rbar: 4 / Rbar                   # above the n <= 3 (shape) modes

    def refs(ax, x0, y0):
        x = np.geomspace(x0, x0 * 10, 10)
        for (lab, z), ls in zip(ZETA_REF.items(), (":", "--")):
            ax.plot(x, y0 * (x / x0) ** (-(1 + 2 * z)), "k" + ls, lw=1, label=f"ζ {lab}")
    fig, axs = plt.subplots(2, 3, figsize=(15, 9))
    zeta = {}
    for a in As:
        for L in Ls:
            k = sel(a, L)
            if k.sum() < 3:
                continue
            q, Sq = bubble_Sq(R[k])
            qm, Sm, Se = binned_mean(q, Sq, qbins)
            zeta[(a, L)] = fit_zeta(qm, Sm, qmin_of(R[k].mean()), qmax_fit)
            if a in (As[0], As[-1]):
                ax = axs[0, 0] if a == As[0] else axs[0, 1]
                ax.errorbar(qm, Sm, Se, fmt="o-", ms=3, lw=0.8, color=lcol[L], label=f"L={L}, ζ={zeta[(a, L)]:.2f}")
            if L == Lmax:
                axs[0, 2].plot(qm, Sm, "o-", ms=3, lw=0.8, color=acol[a], label=f"a={a}")
    kL = sel(None, Lmax)
    q0 = qmin_of(R[kL].mean())
    _, Sref = bubble_Sq(R[kL]); qref, _ = bubble_Sq(R[kL])
    y0 = 3 * np.mean(Sref[np.abs(qref - q0) < 0.2 * q0])
    for ax, ttl in ((axs[0, 0], f"a = {As[0]}"), (axs[0, 1], f"a = {As[-1]}"), (axs[0, 2], f"all a, L = {Lmax}")):
        refs(ax, q0, y0)
        ax.axvspan(q0, qmax_fit, color="0.92", zorder=-1)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("q  (n / ⟨R⟩)")
        ax.set_ylabel("S(q)"); ax.legend(fontsize=7); ax.set_title(f"bubble S(q), {ttl}")

    ax = axs[1, 0]
    for L in Ls:
        aa = [a for a in As if (a, L) in zeta]
        ax.plot(aa, [zeta[(a, L)] for a in aa], "o-", color=lcol[L], label=f"L={L}")
    for (lab, z), ls in zip(ZETA_REF.items(), (":", "--")):
        ax.axhline(z, color="k", ls=ls, lw=0.8, label=lab)
    ax.set_xlabel("anisotropy a"); ax.set_ylabel("ζ")
    ax.set_title(f"roughness exponent, fit 4/⟨R⟩ < q < 1/w = {qmax_fit:.2f}")
    ax.legend(fontsize=7)

    # arc-resolved: edge = normal || b (theta = tb, tb+180), screw = normal _|_ b
    zarc = {}
    for j, a in enumerate((As[0], As[-1])):
        ax = axs[1, 1 + j]
        k = sel(a, Lmax)
        larc = np.pi / 2 * R[k].mean()                    # quadrant arc length
        for kind, cen, ls in (("edge (normal ∥ b)", (tb, tb + 180), "-"),
                              ("screw (normal ⊥ b)", (tb + 90, tb + 270), "--")):
            q, Sq = arc_Sq(R[k], cen, half_width_deg=45)
            qm, Sm, Se = binned_mean(q, Sq, qbins)
            zarc[(a, kind)] = fit_zeta(qm, Sm, 2 * 2 * np.pi / larc, qmax_fit)
            ax.errorbar(qm, Sm, Se, fmt="o" + ls, ms=3, lw=0.9, color="C0" if ls == "-" else "C3",
                        label=f"{kind}, ζ={zarc[(a, kind)]:.2f}")
        ax.axvspan(2 * 2 * np.pi / larc, qmax_fit, color="0.92", zorder=-1)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("q (arc length)")
        ax.set_ylabel("S(q)"); ax.legend(fontsize=7)
        ax.set_title(f"quadrants (±45°), a = {a}, L = {Lmax}")
    fig.tight_layout()
    fig.savefig(os.path.join(d, "Sq.png"), dpi=130)
    plt.close(fig)

    print(f"\n{'a':>5} " + " ".join(f"{'zeta L=' + str(L):>11}" for L in Ls))
    for a in As:
        print(f"{a:5.2f} " + " ".join(f"{zeta.get((a, L), np.nan):11.3f}" for L in Ls))
    for (a, kind), z in zarc.items():
        print(f"a={a}: {kind:20s} zeta = {z:.3f}")
    print("saved", os.path.join(d, "bubbles.png"), "and", os.path.join(d, "Sq.png"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=os.path.join(HERE, "runs", "fss"))
    ap.add_argument("--extract", action="store_true", help="force re-extraction from configs/")
    ap.add_argument("--no-plot", action="store_true", help="only extract (e.g. on the cluster)")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--w", type=float, default=3.0, help="noise correlation length (sets the fit window)")
    args = ap.parse_args()
    fn = os.path.join(args.dir, "shapes.npz")
    n_cfg = len(glob.glob(os.path.join(args.dir, "configs", "*.npz")))
    stale = os.path.exists(fn) and n_cfg and int(np.load(fn)["n_configs"]) != n_cfg
    if args.extract or not os.path.exists(fn) or stale:
        extract(args.dir, args.workers)
    if not args.no_plot:
        plot(args.dir, args.w)


if __name__ == "__main__":
    main()
