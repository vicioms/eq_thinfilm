#!/usr/bin/env python3
"""Soliton fit of the bubble wall along rays (Caballero, Agoritsas, Lecomte, Giamarchi,
PRB 102, 104204 (2020), Sec. III, in polar form).

For every critical configuration phi_c and every ray theta from the area centroid of the
bubble, phi(r) is fitted in a window around the phi = 0 crossing with

    phi(r) = a tanh((rho - r) / w) + b                 (phi = +1 inside)

giving the wall position rho(theta), its width w(theta), the amplitude a and the offset b
(the field h shifts the two minima). Rays crossing the wall more than once (overhangs) are
flagged: the single-crossing fit does not apply there.

Two subcommands:

  extract RUN [--workers N]   reads RUN/configs/*.npz (the cluster side, where fss_scan.py
                              stored phi_c), writes RUN/rayfit.npz
  analyze RUN [--eps0 --c]    reads RUN/rayfit.npz, writes RUN/rayfit.png and RUN/rayfit.txt:
                              width distribution and its dependence on the wall orientation
                              relative to b, mode spectrum S(q) of u = rho - <rho>, and the
                              mean pinning force along the wall from the exact area identity
                              <F>_wall = 2 pi sigma / P - 2 h (valid at a metastable state)

Continuum values used for comparison (phi^4 kink, Sec. 1 of the derivation):
  w0 = sqrt(2 c / eps0),  sigma0 = (2 sqrt 2 / 3) sqrt(c eps0),  sigma_fit = (4/3) c a^2 / w
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
sys.path.insert(0, HERE)
from analysis.geometry import centroid, n_rays, ray_crossings, signed_area  # noqa: E402

CFG_RE = re.compile(r"L(\d+)_a([0-9.]+)_s(\d+)\.npz$")
DR = 0.25            # radial sampling step of the fit window (lattice units)
HALF = 6.0           # half width of the fit window (lattice units)
N_IT = 40            # Levenberg-Marquardt iterations


# ---------------------------------------------------------------------------
# wall and rays
# ---------------------------------------------------------------------------
def main_contour(phi):
    """Largest phi = 0 contour enclosing phi > 0 (the bubble wall), counter-clockwise."""
    import contourpy
    best, bestA = None, 0.0
    for c in contourpy.contour_generator(z=phi).lines(0.0):
        if len(c) < 8 or not np.allclose(c[0], c[-1]):
            continue
        A = signed_area(c)
        if abs(A) <= bestA:
            continue
        cx, cy = centroid(c)
        i, j = int(round(cy)), int(round(cx))
        if 0 <= i < phi.shape[0] and 0 <= j < phi.shape[1] and phi[i, j] <= 0:
            continue                 # a hole, or a horseshoe whose centroid falls outside
        best, bestA = (c if A > 0 else c[::-1]), abs(A)
    return best


def fit_rays(phi, centre, theta, r_guess, w0, half=HALF, dr=DR, n_it=N_IT, order=3):
    """Vectorised Levenberg-Marquardt fit of a tanh((rho - r)/w) + b on every ray.
    phi is sampled with a cubic spline (order=3): bilinear sampling broadens a wall of width
    ~1.5 sites by ~6%. Returns rho, w, a, b, rms residual, converged flag (arrays over rays)."""
    from scipy.ndimage import map_coordinates
    s = np.arange(-half, half + dr / 2, dr)
    r = r_guess[:, None] + s[None, :]
    x = centre[0] + r * np.cos(theta)[:, None]
    y = centre[1] + r * np.sin(theta)[:, None]
    f = map_coordinates(phi, [y.ravel(), x.ravel()], order=order, mode="nearest").reshape(r.shape)

    n = len(theta)
    p = np.c_[r_guess, np.full(n, w0), np.ones(n), np.zeros(n)]          # rho, w, a, b
    lam = np.full(n, 1e-2)

    def resid(p):
        z = (p[:, :1] - r) / p[:, 1:2]
        t = np.tanh(z)
        return p[:, 2:3] * t + p[:, 3:4] - f, z, t

    res, z, t = resid(p)
    cost = np.sum(res ** 2, 1)
    for _ in range(n_it):
        sech2 = 1 - t ** 2
        J = np.stack([p[:, 2:3] * sech2 / p[:, 1:2],                    # d/d rho
                      -p[:, 2:3] * sech2 * z / p[:, 1:2],               # d/d w
                      t, np.ones_like(t)], axis=-1)                     # d/d a, d/d b
        JTJ = np.einsum("nki,nkj->nij", J, J)
        g = np.einsum("nki,nk->ni", J, res)
        A = JTJ + lam[:, None, None] * (np.eye(4) * np.diagonal(JTJ, axis1=1, axis2=2)[:, None, :]
                                        + 1e-9 * np.eye(4))
        dp = -np.linalg.solve(A, g[..., None])[..., 0]
        pn = p + dp
        pn[:, 1] = np.clip(pn[:, 1], 0.05, 4 * half)
        rn, zn, tn = resid(pn)
        cn = np.sum(rn ** 2, 1)
        ok = cn < cost
        p[ok], res[ok], z[ok], t[ok], cost[ok] = pn[ok], rn[ok], zn[ok], tn[ok], cn[ok]
        lam = np.where(ok, lam / 3, lam * 4)
    rms = np.sqrt(cost / len(s))
    conv = (np.abs(p[:, 0] - r_guess) < half / 2) & (p[:, 1] < half) & (p[:, 2] > 0)
    return p[:, 0], p[:, 1], p[:, 2], p[:, 3], rms, conv


def fit_config(phi, w0=np.sqrt(2.0), n_theta=None):
    """All ray fits of one configuration (phi = +1 inside). None if no bubble wall."""
    phi = np.asarray(phi, np.float64)
    c = main_contour(phi)
    if c is None:
        return None
    cen = centroid(c)
    A = signed_area(c)
    P = float(np.sum(np.hypot(*np.diff(c, axis=0).T)))
    n_theta = n_theta or n_rays(phi.shape[1])
    theta = 2 * np.pi * np.arange(n_theta) / n_theta
    mult, mean, outer, _ = ray_crossings(c, cen, n_theta)
    guess = np.where(np.isfinite(outer), outer, np.sqrt(A / np.pi))
    rho, w, a, b, rms, conv = fit_rays(phi, cen, theta, guess, w0)
    return dict(xc=cen[0], yc=cen[1], A=A, P=P, rho=rho, w=w, a=a, b=b, rms=rms,
                ok=conv & (mult == 1), mult=mult, r_contour=guess)


# ---------------------------------------------------------------------------
# extract
# ---------------------------------------------------------------------------
def _one(fn):
    L, a, seed = CFG_RE.search(fn).groups()
    z = np.load(fn)
    out = fit_config(z["phi_c"])
    return int(L), float(a), int(seed), float(z["h_d"]), out


def extract(run, workers):
    files = sorted(glob.glob(os.path.join(run, "configs", "*.npz")))
    if not files:
        sys.exit(f"no configurations in {run}/configs: run the extraction where the scan ran")
    with ProcessPoolExecutor(workers) as ex:
        rows = [r for r in ex.map(_one, files, chunksize=2) if r[4] is not None]
    rows.sort(key=lambda r: (r[1], r[0], r[2]))
    out = {k: np.array([r[i] for r in rows]) for i, k in enumerate(("L", "aniso", "seed", "h_d"))}
    for k in ("xc", "yc", "A", "P"):
        out[k] = np.array([r[4][k] for r in rows])
    for L in sorted(set(out["L"])):
        sel = [r[4] for r in rows if r[0] == L]
        for k, dt in (("rho", np.float32), ("w", np.float32), ("a", np.float32), ("b", np.float32),
                      ("rms", np.float32), ("ok", bool), ("mult", np.int8), ("r_contour", np.float32)):
            out[f"{k}_L{L}"] = np.array([s[k] for s in sel]).astype(dt)
    fn = os.path.join(run, "rayfit.npz")
    np.savez_compressed(fn, **out)
    print(f"saved {fn} ({len(rows)} configurations)")


# ---------------------------------------------------------------------------
# analyze
# ---------------------------------------------------------------------------
def wall_angle_to_b(theta, rho, b_deg):
    """Angle between the local wall line and b, folded to [0, 90] deg (0: screw, 90: edge)."""
    n = len(theta)
    drho = (np.roll(rho, -1) - np.roll(rho, 1)) * n / (4 * np.pi)
    psi = theta - np.arctan2(drho, rho)              # direction of the outward normal
    a = np.mod(np.degrees(psi) + 90.0 - b_deg, 180.0)
    return np.minimum(a, 180.0 - a)


def normal_factor(rho):
    """cos of the angle between ray and wall normal, 1 / sqrt(1 + (rho' / rho)^2), per ray."""
    n = rho.shape[-1]
    d = (np.roll(rho, -1, -1) - np.roll(rho, 1, -1)) * n / (4 * np.pi)
    return 1 / np.sqrt(1 + (d / rho) ** 2)


def spectrum(rho, ok):
    """u = rho - <rho>, overhang rays filled by linear interpolation; |u_m|^2 for m >= 0."""
    th = np.arange(len(rho))
    if ok.sum() < 0.9 * len(rho):
        return None
    r = np.interp(th, th[ok], rho[ok], period=len(rho))
    u = r - r.mean()
    um = np.fft.rfft(u) / len(u)
    return np.abs(um) ** 2, r.mean()


def analyze(run, eps0, c, b_deg):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    z = np.load(os.path.join(run, "rayfit.npz"))
    L, an, h, P = z["L"], z["aniso"], z["h_d"], z["P"]
    w0 = np.sqrt(2 * c / eps0)
    sig0 = 2 * np.sqrt(2) / 3 * np.sqrt(c * eps0)
    As, Ls = sorted(set(an)), sorted(set(L))
    cmap = plt.get_cmap("viridis")
    col = {a: cmap(i / max(1, len(As) - 1)) for i, a in enumerate(As)}
    lines = [f"ray soliton fit, {run}", f"continuum: w0 = {w0:.3f}, sigma0 = {sig0:.4f} "
             f"(eps0 = {eps0}, c = {c})", ""]

    fig, ax = plt.subplots(2, 2, figsize=(11, 9))
    bins_w = np.linspace(0, 3 * w0, 61)
    bins_ang = np.linspace(0, 90, 10)
    sig_fit = np.zeros(len(L))
    idx = {Lv: np.flatnonzero(L == Lv) for Lv in Ls}
    pos = {i: (Lv, k) for Lv in Ls for k, i in enumerate(idx[Lv])}

    lines.append(f"{'a':>5} {'L':>6} {'n':>4} {'<w>':>7} {'sd w':>6} {'<a>':>6} {'<b>':>7} "
                 f"{'ok frac':>7} {'sigma_fit':>9} {'h_d':>8} {'-<F>/2 (s0)':>11} {'-<F>/2 (fit)':>12}")
    for a in As:
        Wa, Aa, pool = [], [], {}
        for Lv in Ls:
            ii = np.flatnonzero((L == Lv) & (an == a))
            if not len(ii):
                continue
            ks = [pos[i][1] for i in ii]
            ok = z[f"ok_L{Lv}"][ks]
            w = z[f"w_L{Lv}"][ks]
            am = z[f"a_L{Lv}"][ks]
            bm = z[f"b_L{Lv}"][ks]
            rho = z[f"rho_L{Lv}"][ks]
            n_t = rho.shape[1]
            theta = 2 * np.pi * np.arange(n_t) / n_t
            w = w * normal_factor(rho)               # ray width -> width along the wall normal
            for j, i in enumerate(ii):
                o = ok[j]
                sig_fit[i] = 4 / 3 * c * np.median(am[j][o] ** 2 / w[j][o]) if o.any() else np.nan
                Wa.append(w[j][o])
                Aa.append(np.c_[wall_angle_to_b(theta, rho[j], b_deg)[o], w[j][o]])
                sp = spectrum(rho[j], o)
                if sp is not None:
                    pool.setdefault(Lv, []).append(sp)
            Fs0 = 2 * np.pi * sig0 / P[ii] - 2 * h[ii]
            Ff = 2 * np.pi * sig_fit[ii] / P[ii] - 2 * h[ii]
            lines.append(f"{a:5.2f} {Lv:6d} {len(ii):4d} {np.mean(w[ok]):7.3f} {np.std(w[ok]):6.3f} "
                         f"{np.mean(am[ok]):6.3f} {np.mean(bm[ok]):7.4f} {ok.mean():7.3f} "
                         f"{np.nanmean(sig_fit[ii]):9.4f} {np.mean(h[ii]):8.5f} "
                         f"{-np.mean(Fs0) / 2:11.5f} {-np.nanmean(Ff) / 2:12.5f}")
            ax[1, 1].errorbar(Lv * (1 + 0.03 * As.index(a)), -np.nanmean(Ff) / 2,
                              np.nanstd(Ff) / 2 / np.sqrt(len(ii)), fmt="o", color=col[a],
                              label=f"a = {a:g}" if Lv == Ls[0] else None)
            ax[1, 1].plot(Lv * (1 + 0.03 * As.index(a)), np.mean(h[ii]), "x", color=col[a])
        Wa = np.concatenate(Wa)
        ax[0, 0].hist(Wa, bins_w, density=True, histtype="step", color=col[a], label=f"a = {a:g}")
        Aa = np.vstack(Aa)
        k = np.digitize(Aa[:, 0], bins_ang) - 1
        m = [np.mean(Aa[k == q, 1]) if np.any(k == q) else np.nan for q in range(len(bins_ang) - 1)]
        ax[0, 1].plot(0.5 * (bins_ang[1:] + bins_ang[:-1]), m, "o-", color=col[a], label=f"a = {a:g}")
        if a == As[-1] or a == As[0]:
            for Lv, sp in pool.items():
                S = np.mean([s[0] for s in sp], 0)
                R = np.mean([s[1] for s in sp])
                mm = np.arange(len(S))
                q = mm[2:] / R
                keep = q <= 2.0                      # wavelengths >= pi sites; rays oversample the wall
                q, S = q[keep], S[2:][keep]
                ax[1, 0].loglog(q, S * 2 * np.pi * R, ".", ms=2,
                                color=col[a], alpha=0.3 + 0.7 * Ls.index(Lv) / max(1, len(Ls) - 1))
    lines += ["", "-<F>/2: mean pinning force along the wall per unit 2h, from the area identity "
              "<F> = 2 pi sigma / P - 2 h at the last metastable state",
              "(s0: continuum sigma; fit: sigma from the fitted profile a^2/w). Compare with h_d.",
              "Widths are along the wall normal. With b along x the orientation axis is also the angle",
              "to the lattice axes: a dependence common to all a is lattice (discrete Laplacian) anisotropy.",
              "Expected offset b = h / (2 eps0) (shift of the minima by the field)."]

    ax[0, 0].axvline(w0, color="k", ls="--", lw=1, label="w0")
    ax[0, 0].set(xlabel="fitted width w", ylabel="pdf", title="wall width (single-crossing rays)")
    ax[0, 0].legend(fontsize=8)
    ax[0, 1].axhline(w0, color="k", ls="--", lw=1)
    ax[0, 1].set(xlabel="angle between wall line and b [deg] (0 screw, 90 edge)", ylabel="<w>",
                 title="width vs wall orientation")
    ax[0, 1].legend(fontsize=8)
    qq = np.logspace(-1.5, 0.3, 10)
    for zeta, ls in ((1.25, "--"), (1.0, ":")):
        ax[1, 0].loglog(qq, 1e-3 * qq ** (-(1 + 2 * zeta)), "k" + ls, lw=1, label=f"zeta = {zeta}")
    ax[1, 0].set(xlabel="q = m / R", ylabel="S(q) = 2 pi R |u_m|^2",
                 title=f"mode spectrum (a = {As[0]:g}, {As[-1]:g}; darker = larger L)")
    ax[1, 0].legend(fontsize=8)
    ax[1, 1].set_xscale("log")
    ax[1, 1].set(xlabel="L", ylabel="field", title="-<F>/2 from area identity (o) vs h_d (x)")
    ax[1, 1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(run, "rayfit.png"), dpi=130)
    open(os.path.join(run, "rayfit.txt"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("saved", os.path.join(run, "rayfit.png"), os.path.join(run, "rayfit.txt"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extract")
    e.add_argument("run")
    e.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    n = sub.add_parser("analyze")
    n.add_argument("run")
    n.add_argument("--eps0", type=float, default=1.0)
    n.add_argument("--c", type=float, default=1.0, help="mean bond coupling (cx + cy) / 2")
    n.add_argument("--burgers-angle", type=float, default=None,
                   help="default: from the run's results*.jsonl")
    args = ap.parse_args()
    if args.cmd == "extract":
        extract(args.run, args.workers)
    else:
        b = args.burgers_angle
        if b is None:
            from analysis.data import load_results
            meta = load_results(args.run)
            b = float(next(iter(meta.values())).get("burgers_angle", 0.0)) if meta else 0.0
        analyze(args.run, args.eps0, args.c, b)


if __name__ == "__main__":
    main()
