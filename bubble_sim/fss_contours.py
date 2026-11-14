#!/usr/bin/env python3
"""
Geometry of the critical bubble walls from contours.npz (fss_shapes.py --contours).

Why contours and not R(theta)
-----------------------------
R(theta) (mean of the ray crossings, the polar analogue of Kolton-Ferrero-Rosso
Eq. 9) is single valued: where a ray crosses the wall 3, 5, ... times it jumps,
and those jumps add a q^-2 tail to S(q) that grows with the number of overhangs,
i.e. with L. The contour is the wall itself, so everything below is computed on
the phi = 0 line (the bubble's own wall: the largest positive-area contour).

Per critical bubble (one row per sample):

 1. size, memory of the start
      A, P (area, length of the wall), R_eff = sqrt(A/pi), advance R_eff - R0,
      rms radial roughness w_R (modes n >= 3) -> w_R / (R_eff - R0) says whether
      the wall travelled far enough to forget the initial circle.

 2. overhangs (polar, from rays out of the centroid, KFR Sec. III.2)
      mult(theta): crossings per ray; rho = P(any ray has mult > 1)  [KFR Fig. 5]
      u_o^2 = < var_n R_n(theta) >_theta                              [KFR Eq. 11]
      n_ov  = number of separate multivalued angular runs; l_ov = P / n_ov is the
              wall length per overhang (KFR's l_o is the length scale of rho).
      R_mean (KFR Eq. 9), R_out (outermost crossing, the "envelope" of
      Zhou & Zheng, PRE 82, 031139) for the polar S(q).

 3. roughness, intrinsic (no single-valuedness needed)
      psi(s) = tangent angle vs arc length s, dpsi = psi - 2 pi s / P.
      For a small-slope wall dpsi = du/ds, so
          S_psi(q) = q^2 S_u(q) ~ q^(1 - 2 zeta)
      (zeta = 1.25 qEW -> q^-1.5; zeta = 1 -> 1/q, a self-similar wall whose
      winding angle variance grows like ln l, as for percolation hulls;
      zeta = 0.5 -> flat). Also the winding structure function
          G(l) = <[dpsi(s + l) - dpsi(s)]^2> ~ l^(2 (zeta - 1))   (zeta > 1)
      S(q) and these slope statistics stay valid for zeta > 1, where the local
      width only gives zeta_loc = 1 (Lopez et al.; Guyonnet et al.,
      arXiv:1904.11726 Sec. II.2.2).
      Divider (ruler) length P(eps) ~ eps^(1 - d_f).

 4. anisotropy (Burgers vector b along x)
      local outward normal angle phi_n of the wall smoothed over sigma along s;
      length-weighted P(phi_n) and E = <cos 2 phi_n> (+1: all wall has normal
      parallel to b = edge, -1: normal perpendicular to b = screw; 0 for a
      circle and, by symmetry, for any square-lattice artefact).
      Support function about the nucleation centre: advance of the most
      advanced wall point in each direction, h(phi) - R0.
      Every wall point is classified edge / screw / mixed by the normal of the
      wall smoothed over 32 sites; slope fluctuations below a scale sigma,
      var(psi - psi_sigma), are compared between classes for sigma = 4..32
      (no polar quadrants and no mean-shape subtraction; at sigma = 32 the
      curvature of the bubble tips, radius ~ 45 at a = 1, starts to leak in).
      Intrinsic overhangs: fraction of wall turned back by more than 90 deg
      relative to the wall smoothed over sigma.

Writes <dir>/contour_stats.npz and <dir>/contours_{overhangs,roughness,anisotropy}.png.
`--selftest` checks the estimators on synthetic walls of known zeta.

    python fss_contours.py --dir runs/fss_D0.04
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

DS = 0.5                    # arc-length resampling step (lattice units)
TAN_HALF = 2                # tangent from the chord p[k+m] - p[k-m]: 2 m DS = 2 sites
SIGMAS = (2.0, 8.0, 32.0)   # smoothing scales for the normal-angle statistics
N_PHI_BINS = 18             # P(phi_n) on [0, 180)
N_SUPPORT = 72              # support function directions on [0, 360)
EPS_RULER = np.geomspace(1, 512, 19)
ELL_G = np.geomspace(1, 2048, 23)
CLASS_HALF = 30.0           # edge: normal within 30 deg of b; screw: within 30 deg of b-perp
Q_EDGES = np.geomspace(1e-3, 2.0, 34)
W_CORE = 3.0                # noise correlation length w (spectrum cut off as exp(-k^2 w^2))
Q_MAX_FIT = 0.1             # above q ~ 0.12 zeta_eff rises steeply (core cut-off): fits stop here
N_MIN_FIT = 12              # fits start at mode n = 12 (the lowest ~5 modes are depleted)
HP_SIGMA = 16.0             # high-pass scale for the orientation-resolved slope variance
CLASS_SIGMA = 32.0          # wall points are classified edge/screw by the normal of the curve smoothed here
HP_SIGMAS = (4.0, 8.0, 16.0, 32.0)
OV_SIGMAS = (32.0, 128.0)   # intrinsic overhangs: |psi - psi_smoothed(sigma)| > pi / 2


# ---------------------------------------------------------------------------
# curve primitives
# ---------------------------------------------------------------------------
def resample_closed(xy, ds=DS):
    """Uniform arc-length resampling of a closed polygon; returns (pts, P), CCW."""
    xy = np.asarray(xy, float)
    if np.allclose(xy[0], xy[-1]):
        xy = xy[:-1]
    x, y = xy[:, 0], xy[:, 1]
    if 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y) < 0:
        xy = xy[::-1]
    seg = np.hypot(*np.diff(np.vstack([xy, xy[:1]]), axis=0).T)
    s = np.r_[0, np.cumsum(seg)]
    P = s[-1]
    n = max(16, int(round(P / ds)))
    t = np.arange(n) * (P / n)
    closed = np.vstack([xy, xy[:1]])
    pts = np.c_[np.interp(t, s, closed[:, 0]), np.interp(t, s, closed[:, 1])]
    return pts, P


def tangent_angle(pts, m=TAN_HALF):
    """Unwrapped tangent angle and total turning (2 pi for a simple CCW curve)."""
    d = np.roll(pts, -m, axis=0) - np.roll(pts, m, axis=0)
    raw = np.arctan2(d[:, 1], d[:, 0])
    step = (np.diff(np.r_[raw, raw[:1]]) + np.pi) % (2 * np.pi) - np.pi
    return np.unwrap(raw), float(step.sum())


def smooth_closed(pts, sigma, P):
    """Gaussian smoothing along arc length (periodic), via FFT."""
    n = len(pts)
    q = 2 * np.pi * np.fft.rfftfreq(n, d=P / n)
    g = np.exp(-0.5 * (q * sigma) ** 2)
    return np.c_[np.fft.irfft(np.fft.rfft(pts[:, 0]) * g, n), np.fft.irfft(np.fft.rfft(pts[:, 1]) * g, n)]


def normals_angle(pts, P):
    """Outward normal angle of a CCW closed curve (spectral derivative)."""
    n = len(pts)
    q = 2 * np.pi * np.fft.rfftfreq(n, d=P / n)
    dx = np.fft.irfft(1j * q * np.fft.rfft(pts[:, 0]), n)
    dy = np.fft.irfft(1j * q * np.fft.rfft(pts[:, 1]), n)
    return np.arctan2(dy, dx) - np.pi / 2


@njit(cache=True)
def ruler_length(pts, eps):
    """Divider length of a closed polyline (densely sampled) with ruler eps."""
    n = pts.shape[0]
    x0, y0 = pts[0, 0], pts[0, 1]
    i = 0
    steps = 0
    while True:
        j = i + 1
        while j < n and (pts[j, 0] - x0) ** 2 + (pts[j, 1] - y0) ** 2 < eps * eps:
            j += 1
        if j >= n:
            # remainder back to the start
            return steps * eps + np.sqrt((pts[0, 0] - x0) ** 2 + (pts[0, 1] - y0) ** 2)
        # intersection of the circle |p - p0| = eps with segment (j-1, j)
        ax, ay = pts[j - 1, 0] - x0, pts[j - 1, 1] - y0
        bx, by = pts[j, 0] - pts[j - 1, 0], pts[j, 1] - pts[j - 1, 1]
        A = bx * bx + by * by
        B = 2 * (ax * bx + ay * by)
        C = ax * ax + ay * ay - eps * eps
        disc = max(B * B - 4 * A * C, 0.0)
        t = (-B + np.sqrt(disc)) / (2 * A) if A > 0 else 1.0
        x0, y0 = pts[j - 1, 0] + t * bx, pts[j - 1, 1] + t * by
        i = j - 1
        steps += 1


def ray_crossings(xy, centre, n_theta):
    """All crossings of the closed polygon with n_theta rays from `centre`.
    Returns per-ray multiplicity, mean, outermost, innermost radius and variance."""
    xy = np.asarray(xy, float)
    if not np.allclose(xy[0], xy[-1]):
        xy = np.vstack([xy, xy[:1]])
    d = xy - centre
    th = np.arctan2(d[:, 1], d[:, 0])
    r = np.hypot(d[:, 0], d[:, 1])
    dth = np.diff(th)
    dth = (dth + np.pi) % (2 * np.pi) - np.pi
    T = np.r_[th[0], th[0] + np.cumsum(dth)]          # unwrapped angle along the contour
    dt = 2 * np.pi / n_theta
    a, b = T[:-1], T[1:]
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    k0 = np.ceil(lo / dt).astype(np.int64)
    k1 = np.ceil(hi / dt).astype(np.int64)            # rays k0 .. k1-1 cross the segment
    cnt = np.maximum(k1 - k0, 0)
    seg = np.repeat(np.arange(len(a)), cnt)
    k = np.repeat(k0, cnt) + (np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt))
    tk = k * dt
    f = (tk - a[seg]) / np.where(b[seg] != a[seg], b[seg] - a[seg], 1.0)
    rk = r[:-1][seg] + f * (r[1:][seg] - r[:-1][seg])
    kk = np.mod(k, n_theta)
    mult = np.bincount(kk, minlength=n_theta)
    s1 = np.bincount(kk, rk, n_theta)
    s2 = np.bincount(kk, rk * rk, n_theta)
    out = np.full(n_theta, -np.inf); np.maximum.at(out, kk, rk)
    inn = np.full(n_theta, np.inf); np.minimum.at(inn, kk, rk)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = s1 / mult
        var = np.maximum(s2 / mult - mean ** 2, 0)
    return mult, mean, out, inn, var


def runs_of(mask):
    """Number of separate cyclic runs of True."""
    if mask.all():
        return 1
    if not mask.any():
        return 0
    return int(np.sum(mask & ~np.roll(mask, 1)))


def psd_periodic(f, P):
    """S(q_n) = P |c_n|^2, c_n = (1/N) sum f e^{-i q_n s}, q_n = 2 pi n / P (n >= 1)."""
    n = len(f)
    c = np.fft.rfft(f - f.mean()) / n
    q = 2 * np.pi * np.arange(len(c)) / P
    return q[1:], P * np.abs(c[1:]) ** 2


def psd_window(f, ds):
    """Hann-windowed, linearly detrended PSD of a non-periodic segment (same normalisation)."""
    n = len(f)
    s = np.arange(n)
    f = f - np.polyval(np.polyfit(s, f, 1), s)
    w = np.hanning(n)
    c = np.fft.rfft(f * w)
    q = 2 * np.pi * np.fft.rfftfreq(n, d=ds)
    return q[1:], ds * np.abs(c[1:]) ** 2 / np.sum(w ** 2)


def bin_log_err(q, S, err, edges=Q_EDGES):
    """Bin means of S and the error of the bin mean from per-point errors."""
    idx = np.digitize(q, edges) - 1
    ok = (idx >= 0) & (idx < len(edges) - 1)
    c = np.bincount(idx[ok], minlength=len(edges) - 1)
    m = np.bincount(idx[ok], S[ok], len(edges) - 1)
    v = np.bincount(idx[ok], err[ok] ** 2, len(edges) - 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return m / c, np.sqrt(v) / c


def bin_log(q, S, edges=Q_EDGES):
    """Mean of S in log-q bins (NaN where empty)."""
    idx = np.digitize(q, edges) - 1
    ok = (idx >= 0) & (idx < len(edges) - 1)
    s = np.bincount(idx[ok], S[ok], len(edges) - 1)
    c = np.bincount(idx[ok], minlength=len(edges) - 1)
    with np.errstate(invalid="ignore"):
        return s / c


# ---------------------------------------------------------------------------
# one bubble
# ---------------------------------------------------------------------------
def analyse_wall(xy, L, R0, n_theta, b_angle=0.0):
    pts, P = resample_closed(xy)
    x, y = pts[:, 0], pts[:, 1]
    A = 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)
    cen = np.array([x.mean(), y.mean()])
    # area centroid (exact for the polygon)
    cr = x * np.roll(y, -1) - np.roll(x, -1) * y
    cen = np.array([np.sum((x + np.roll(x, -1)) * cr), np.sum((y + np.roll(y, -1)) * cr)]) / (6 * A)
    R_eff = np.sqrt(A / np.pi)
    out = dict(P=P, A=A, R_eff=R_eff, xc=cen[0], yc=cen[1], advance=R_eff - R0)

    # polar crossings (KFR Eq. 9, 11) about the centroid
    mult, Rm, Ro, Ri, var = ray_crossings(pts, cen, n_theta)
    good = mult > 0
    Rm_f = np.where(good, Rm, np.nanmean(Rm[good]))
    Ro_f = np.where(good, Ro, np.nanmean(Ro[good]))
    multi = mult > 1
    out.update(frac_multi=multi.mean(), any_ov=bool(multi.any()), n_ov=runs_of(multi),
               u_o=np.sqrt(var[good].mean()), frac_noray=1 - good.mean())
    cR = np.fft.rfft(Rm_f) / n_theta
    cO = np.fft.rfft(Ro_f) / n_theta
    out["cR"] = cR[:n_theta // 2].astype(np.complex64)     # keep complex modes: connected S(q)
    out["cO"] = cO[:n_theta // 2].astype(np.complex64)
    out["w_R"] = np.sqrt(2 * np.sum(np.abs(cR[3:]) ** 2))  # rms of R(theta) without n = 0, 1, 2

    # turning angle
    psi, tot = tangent_angle(pts)
    s = np.arange(len(pts)) * (P / len(pts))
    dpsi = psi - psi[0] - 2 * np.pi * s / P
    out["winding_total"] = tot
    q, Sp = psd_periodic(dpsi, P)
    out["Spsi_bins"] = bin_log(q, Sp)
    # winding structure function via circular autocorrelation
    F = np.fft.rfft(dpsi - dpsi.mean())
    C = np.fft.irfft(np.abs(F) ** 2, len(dpsi)) / len(dpsi)
    lag = np.rint(ELL_G / (P / len(pts))).astype(int)
    G = np.where(lag < len(pts) // 4, 2 * (C[0] - C[np.minimum(lag, len(pts) - 1)]), np.nan)
    out["G"] = G

    # divider length
    fine, _ = resample_closed(xy, ds=0.25)
    out["ruler"] = np.array([ruler_length(fine, e) if e < R_eff else np.nan for e in EPS_RULER])

    # normal-angle statistics, relative to b
    phis = {}
    for sig in SIGMAS:
        sm = smooth_closed(pts, sig, P)
        phn = normals_angle(sm, P) - np.deg2rad(b_angle)
        phis[sig] = phn
        out[f"E_{sig:g}"] = float(np.mean(np.cos(2 * phn)))
        out[f"Pphi_{sig:g}"] = np.histogram(np.mod(phn, np.pi), bins=N_PHI_BINS, range=(0, np.pi))[0] / len(phn)

    # support function about the nucleation centre
    c0 = np.array([(L - 1) / 2, (L - 1) / 2])
    ang = np.linspace(0, 2 * np.pi, N_SUPPORT, endpoint=False) + np.deg2rad(b_angle)
    nh = np.c_[np.cos(ang), np.sin(ang)]
    out["support"] = ((pts - c0) @ nh.T).max(axis=0) - R0

    # local class of every wall point from the normal of the curve smoothed at HP_SIGMA
    def psi_of(curve):
        return tangent_angle(curve)[0]

    def classify(phn):
        d = np.rad2deg(np.mod(phn, np.pi))                     # 0 / 180: normal || b (edge)
        dist_b = np.minimum(d, 180 - d)
        return np.where(dist_b < CLASS_HALF, 0, np.where(abs(d - 90) < CLASS_HALF, 1, 2))

    names = ("edge", "screw", "mixed")
    # intrinsic overhangs: the wall turns back by more than 90 deg w.r.t. its smoothed direction
    for sig in OV_SIGMAS:
        sm = smooth_closed(pts, sig, P)
        dps = psi - psi_of(sm)
        dps -= 2 * np.pi * np.round(np.median(dps) / (2 * np.pi))
        back = np.abs(dps) > np.pi / 2
        cls = classify(phis.get(sig, normals_angle(sm, P) - np.deg2rad(b_angle)))
        out[f"ovi_{sig:g}"] = back.mean()
        out[f"ovi_runs_{sig:g}"] = runs_of(back)
        for ci, c in enumerate(names):
            m = cls == ci
            out[f"ovi_{sig:g}_{c}"] = back[m].mean() if m.any() else np.nan
            out[f"len_{sig:g}_{c}"] = m.mean()
    # high-passed turning angle (removes the bubble's own shape at scales > HP_SIGMA)
    cls_pt = classify(phis[CLASS_SIGMA])
    out["psi_hp_var"] = np.full((len(HP_SIGMAS), 3), np.nan)    # [sigma, class]: var of psi - psi_sigma
    out["n_class"] = np.bincount(cls_pt, minlength=3)
    for si, sh in enumerate(HP_SIGMAS):
        hp = psi - psi_of(smooth_closed(pts, sh, P))
        hp -= 2 * np.pi * np.round(np.median(hp) / (2 * np.pi))
        for ci in range(3):
            if (cls_pt == ci).any():
                out["psi_hp_var"][si, ci] = np.var(hp[cls_pt == ci])
        if sh == HP_SIGMA:
            psi_hp = hp
    return out


def _worker(args):
    xy, L, R0, nt = args
    return analyse_wall(xy, L, R0, nt)


# ---------------------------------------------------------------------------
# whole run
# ---------------------------------------------------------------------------
def n_theta(L):
    return max(1024, 2 ** int(np.ceil(np.log2(np.pi * L))))


def compute(d, workers):
    import fss_shapes
    C = fss_shapes.load_contours(d)
    R0 = {}
    for fn in glob.glob(os.path.join(d, "results*.jsonl")):
        for line in open(fn):
            r = json.loads(line)
            R0[(r["L"], r["aniso"], r["seed"])] = r["R0"]
    n = len(C["L"])
    jobs = []
    for i in range(n):
        key = (int(C["L"][i]), float(C["aniso"][i]), int(C["seed"][i]))
        jobs.append((C["contours"][i][C["main"][i]], key[0], R0[key], n_theta(key[0])))
    print(f"analysing {n} walls ...")
    with ProcessPoolExecutor(workers) as ex:
        rows = list(ex.map(_worker, jobs, chunksize=4))
    out = {k: np.array([C[k][i] for i in range(n)]) for k in ("L", "aniso", "seed", "h_d")}
    out["R0"] = np.array([j[2] for j in jobs])
    out["n_holes"] = np.array([int((np.asarray(a) < 0).sum()) for a in C["area"]])
    out["n_islands"] = np.array([int((np.asarray(a) > 0).sum()) - 1 for a in C["area"]])
    for k in rows[0]:
        v = [r[k] for r in rows]
        if isinstance(v[0], np.ndarray) and k in ("cR", "cO"):
            for L in np.unique(out["L"]):
                out[f"{k}_L{L}"] = np.stack([v[i] for i in range(n) if out["L"][i] == L])
        else:
            out[k] = np.array(v)
    fn = os.path.join(d, "contour_stats.npz")
    np.savez_compressed(fn, **out)
    print("saved", fn)
    return out


# ---------------------------------------------------------------------------
# aggregation helpers
# ---------------------------------------------------------------------------
def groups(S):
    As = sorted(set(S["aniso"].tolist()))
    Ls = sorted(set(S["L"].tolist()))
    return As, Ls


def mean_err(X):
    X = np.asarray(X, float)
    n = np.sum(np.isfinite(X), axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        m = np.nanmean(X, axis=0)
        e = np.nanstd(X, axis=0, ddof=1) / np.sqrt(n)
    return m, e


def fit_power(q, S, Se, qlo, qhi):
    ok = np.isfinite(S) & (S > 0) & (q >= qlo) & (q <= qhi) & np.isfinite(Se) & (Se > 0)
    if ok.sum() < 3:
        return np.nan, np.nan
    x, y = np.log(q[ok]), np.log(S[ok])
    sy = Se[ok] / S[ok]
    W = 1 / sy ** 2
    A = np.c_[np.ones_like(x), x]
    cov = np.linalg.inv(A.T @ (A * W[:, None]))
    c = cov @ (A.T @ (W * y))
    return c[1], np.sqrt(cov[1, 1])


def connected_polar_S(S, key, L, a):
    """Connected S(q_n) = 2 pi <R> (<|c_n|^2> - |<c_n>|^2) of a polar profile, q_n = n / <R>."""
    sel = S["aniso"][S["L"] == L] == a
    c = S[f"{key}_L{L}"][sel].astype(complex)
    R = np.abs(c[:, 0]).mean()
    full = (np.abs(c) ** 2).mean(0)
    conn = full - np.abs(c.mean(0)) ** 2
    n = np.arange(c.shape[1])
    # error of the connected mean from sample-to-sample spread of |c_n - <c_n>|^2
    dev = np.abs(c - c.mean(0)) ** 2
    err = dev.std(0, ddof=1) / np.sqrt(len(c))
    return n[3:] / R, 2 * np.pi * R * conn[3:], 2 * np.pi * R * full[3:], 2 * np.pi * R * err[3:], R


# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------
def plot(d, S):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm
    As, Ls = groups(S)
    cL = {L: cm.plasma(i / max(1, len(Ls) - 1) * 0.85) for i, L in enumerate(Ls)}
    ca = {a: cm.viridis(i / max(1, len(As) - 1)) for i, a in enumerate(As)}
    qc = np.sqrt(Q_EDGES[1:] * Q_EDGES[:-1])
    names = ("edge", "screw", "mixed")
    report = []

    def sel(L=None, a=None):
        m = np.ones(len(S["L"]), bool)
        if L is not None:
            m &= S["L"] == L
        if a is not None:
            m &= S["aniso"] == a
        return m

    def ratio_of_means(x, y):
        (mx, ex), (my, ey) = mean_err(x), mean_err(y)
        r = mx / my
        return r, abs(r) * np.hypot(ex / mx, ey / my)

    # l_ov at the largest sizes, all a (sets the overhang wavenumber)
    big = S["L"] >= Ls[-2]
    l_ov = S["P"][big].sum() / max(1, S["n_ov"][big].sum())
    q_ov = 2 * np.pi / l_ov

    def spsi(L, a):
        m, e = mean_err(S["Spsi_bins"][sel(L, a)])
        return m, e, S["P"][sel(L, a)].mean()

    def zeta_fit(L, a, qlo=None, qhi=Q_MAX_FIT):
        m, e, Pm = spsi(L, a)
        lo = 2 * np.pi * N_MIN_FIT / Pm if qlo is None else max(qlo, 2 * np.pi * N_MIN_FIT / Pm)
        if qhi / lo < 2:
            return np.nan, np.nan
        sl, se = fit_power(qc, m, e, lo, qhi)
        return (1 - sl) / 2, se / 2

    # ================================================================== overhangs & memory
    fig, axs = plt.subplots(2, 4, figsize=(20, 9))
    ax = axs[0, 0]
    for a in As:
        ax.plot(Ls, [np.mean(S["any_ov"][sel(L, a)]) for L in Ls], "o-", color=ca[a], label=f"a={a}")
    ax.set_xscale("log"); ax.set_xlabel("L"); ax.set_ylabel(r"$\rho$ = P(some ray crosses the wall > once)")
    ax.set_title("polar overhang probability (KFR Fig. 5 analogue)"); ax.legend(fontsize=8)
    ax = axs[0, 1]
    for a in As:
        lov = [S["P"][sel(L, a)].sum() / S["n_ov"][sel(L, a)].sum() if S["n_ov"][sel(L, a)].sum() else np.nan
               for L in Ls]
        ax.plot(Ls, lov, "o-", color=ca[a], label=f"a={a}")
    ax.axhline(l_ov, color="k", ls=":", lw=0.8, label=rf"$\ell_{{ov}}$={l_ov:.0f} (L ≥ {Ls[-2]})")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("L")
    ax.set_ylabel(r"$\ell_{ov}$ = wall length / number of polar overhangs")
    ax.set_title("wall length per overhang"); ax.legend(fontsize=8)
    ax = axs[0, 2]
    for a in As:
        m, e = zip(*[mean_err(S["u_o"][sel(L, a)]) for L in Ls])
        ax.errorbar(Ls, m, e, fmt="o-", color=ca[a], label=f"a={a}")
    ax.plot(Ls, np.array(Ls) / Ls[-1] * np.nanmean(S["u_o"][sel(Ls[-1])]), "k:", lw=0.8, label="∝ L (KFR)")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("L"); ax.set_ylabel(r"$u_o$ (KFR Eq. 11, polar)")
    ax.set_title("overhang size"); ax.legend(fontsize=8)
    ax = axs[0, 3]
    for a in As:
        for sig, ls in zip(OV_SIGMAS, ("-", "--")):
            m, e = zip(*[mean_err(S[f"ovi_{sig:g}"][sel(L, a)]) for L in Ls])
            ax.errorbar(Ls, m, e, fmt="o" + ls, color=ca[a], ms=3, label=f"a={a}" if sig == OV_SIGMAS[0] else None)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("L")
    ax.set_ylabel(r"fraction of wall with $|\psi-\psi_\sigma|>\pi/2$")
    ax.set_title(f"intrinsic overhangs (σ={OV_SIGMAS[0]:g} solid, {OV_SIGMAS[1]:g} dashed)"); ax.legend(fontsize=8)
    ax = axs[1, 0]
    for a in As:
        m, e = zip(*[mean_err(S["advance"][sel(L, a)]) for L in Ls])
        ax.errorbar(Ls, m, e, fmt="o-", color=ca[a], label=f"a={a}")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("L"); ax.set_ylabel(r"$R_{eff}-R_0$")
    ax.set_title("distance travelled by the wall"); ax.legend(fontsize=8)
    ax = axs[1, 1]
    for a in As:
        m, e = zip(*[ratio_of_means(S["w_R"][sel(L, a)], S["advance"][sel(L, a)]) for L in Ls])
        ax.errorbar(Ls, m, e, fmt="o-", color=ca[a], label=f"a={a}")
    ax.set_xscale("log"); ax.set_xlabel("L")
    ax.set_ylabel(r"$\langle w_R\rangle/\langle R_{eff}-R_0\rangle$  ($w_R$: rms of R(θ), n ≥ 3)")
    ax.set_title("roughness / distance travelled"); ax.legend(fontsize=8)
    ax = axs[1, 2]
    for a in As:
        m, e = zip(*[mean_err(S["P"][sel(L, a)] / (2 * np.pi * S["R_eff"][sel(L, a)])) for L in Ls])
        ax.errorbar(Ls, m, e, fmt="o-", color=ca[a], label=f"a={a}")
    ax.set_xscale("log"); ax.set_xlabel("L"); ax.set_ylabel(r"$P/(2\pi R_{eff})$")
    ax.set_title("excess wall length (resampled at 0.5 sites)"); ax.legend(fontsize=8)
    ax = axs[1, 3]
    sig = OV_SIGMAS[0]
    for a, ls in ((As[0], "--"), (As[-1], "-")):
        for c, clr in (("edge", "C0"), ("screw", "C3")):
            m, e = zip(*[mean_err(S[f"ovi_{sig:g}_{c}"][sel(L, a)]) for L in Ls])
            ax.errorbar(Ls, m, e, fmt="o" + ls, color=clr, ms=3, label=f"{c}, a={a}")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("L")
    ax.set_ylabel(f"intrinsic overhang fraction (σ={sig:g})")
    ax.set_title("overhangs on edge vs screw wall"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(d, "contours_overhangs.png"), dpi=110); plt.close(fig)

    report.append(f"wall length per polar overhang at L >= {Ls[-2]}: l_ov = {l_ov:.0f}  (q_ov = 2 pi / l_ov = {q_ov:.4f})")
    report.append(f"{'a':>5} {'L':>5} {'n':>4} {'rho':>5} {'l_ov':>7} {'u_o':>7} {'adv':>7} {'w_R/adv':>7} {'P/2piR':>7} "
                  f"{'ovi32':>7} {'ovi128':>7} {'tip/half':>8}")
    for a in As:
        for L in Ls:
            m = sel(L, a)
            if not m.any():
                continue
            lov = S["P"][m].sum() / max(1, S["n_ov"][m].sum())
            tip = np.mean(S["R0"][m] + S["support"][m].max(axis=1)) / (L / 2)
            report.append(f"{a:5.2f} {L:5d} {m.sum():4d} {S['any_ov'][m].mean():5.2f} {lov:7.0f} {S['u_o'][m].mean():7.2f} "
                          f"{S['advance'][m].mean():7.1f} {S['w_R'][m].mean() / S['advance'][m].mean():7.2f} "
                          f"{np.mean(S['P'][m] / (2 * np.pi * S['R_eff'][m])):7.3f} {S['ovi_32'][m].mean():7.4f} "
                          f"{S['ovi_128'][m].mean():7.4f} {tip:8.2f}")
    report.append("  tip/half = <R0 + max_phi support> / (L/2): how close the farthest wall point gets to the box edge")

    # ================================================================== roughness
    fig, axs = plt.subplots(2, 3, figsize=(16, 9.5))
    for col, a in enumerate((As[0], As[-1])):
        ax = axs[0, col]
        for L in Ls:
            m, e, Pm = spsi(L, a)
            z, ze = zeta_fit(L, a)
            ax.errorbar(qc, m, e, fmt=".-", color=cL[L], lw=0.8, ms=4, label=f"L={L}: ζ={z:.2f}±{ze:.2f}")
        x = np.geomspace(2e-3, 0.1, 10)
        m4, _, _ = spsi(Ls[-1], a)
        y0 = np.nanmean(m4[(qc > 0.03) & (qc < 0.04)])
        for z, ls in ((1.25, "--"), (1.0, ":")):
            ax.plot(x, y0 * (x / 0.035) ** (1 - 2 * z), "k" + ls, lw=0.8, label=f"ζ={z}")
        ax.axvspan(Q_MAX_FIT, 2, color="0.93", zorder=-1)
        ax.axvline(q_ov, color="C2", lw=0.8, ls="--", label=r"$q_{ov}=2\pi/\ell_{ov}$")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("q (arc length)")
        ax.set_ylabel(r"$S_\psi(q)$  ~  $q^{1-2\zeta}$"); ax.legend(fontsize=6.5, ncol=2)
        ax.set_title(f"turning-angle spectrum of the wall, a = {a} (fit: n ≥ {N_MIN_FIT}, q < {Q_MAX_FIT})", fontsize=10)
    ax = axs[0, 2]
    for a, mk, ls in ((As[0], "o", "-"), (As[-1], "s", "--")):
        for L in Ls[2:]:
            m, e, Pm = spsi(L, a)
            k = 2
            lq, lS, rel = np.log(qc), np.log(m), e / m
            zl = (1 - (lS[2 * k:] - lS[:-2 * k]) / (lq[2 * k:] - lq[:-2 * k])) / 2
            ze = 0.5 * np.hypot(rel[2 * k:], rel[:-2 * k]) / (lq[2 * k:] - lq[:-2 * k])
            ok = (qc[k:-k] > 2 * np.pi * N_MIN_FIT / Pm) & (qc[k:-k] < 0.3)
            ax.errorbar(qc[k:-k][ok], zl[ok], ze[ok], fmt=mk + ls, color=cL[L], ms=3, lw=0.8, label=f"a={a}, L={L}")
    for z in (1.25, 1.0):
        ax.axhline(z, color="k", lw=0.5, ls=":")
    ax.axvspan(Q_MAX_FIT, 0.3, color="0.93", zorder=-1)
    ax.axvline(q_ov, color="C2", lw=0.8, ls="--")
    ax.set_xscale("log"); ax.set_ylim(0.4, 2.2); ax.set_xlabel("q (arc length)")
    ax.set_ylabel(r"local $\zeta_{eff}(q)$ from $S_\psi$ (slope over a factor 2.4 in q)")
    ax.set_title(r"local exponent (circles a=0, squares a=1); dashed green: $q_{ov}$", fontsize=10)
    ax.legend(fontsize=6.5, ncol=2)
    ax = axs[1, 0]
    for a, ls in ((As[0], "-"), (As[-1], "--")):
        for L in Ls:
            m, e = mean_err(S["G"][sel(L, a)])
            ax.errorbar(ELL_G, m, e, fmt="o" + ls, color=cL[L], ms=3, lw=0.8, label=f"a={a}, L={L}" if L in (Ls[0], Ls[-1]) else None)
    x = np.geomspace(10, 300, 10)
    m4, _ = mean_err(S["G"][sel(Ls[-1], As[0])])
    g0 = np.interp(10, ELL_G, m4)
    ax.plot(x, g0 * (x / 10) ** 0.5, "k--", lw=0.8, label=r"$\ell^{2(\zeta-1)}$, ζ=1.25")
    ax.axvline(l_ov, color="C2", lw=0.8, ls="--", label=r"$\ell_{ov}$")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(r"arc length $\ell$")
    ax.set_ylabel(r"$G(\ell)=\langle[\delta\psi(s+\ell)-\delta\psi(s)]^2\rangle$ [rad²]")
    ax.set_title("winding-angle structure function (solid a=0, dashed a=1)"); ax.legend(fontsize=7)
    ax = axs[1, 1]
    Lm = Ls[-1]
    for a, base in ((As[0], 0), (As[-1], 1)):
        q, conn, full, err, R = connected_polar_S(S, "cR", Lm, a)
        qO, connO, _, errO, _ = connected_polar_S(S, "cO", Lm, a)
        sc = 10.0 ** (3 * base)
        ax.plot(q, full * sc, "-", color="0.6", lw=0.8, label="R_mean (KFR Eq. 9), full" if base == 0 else None)
        ax.plot(q, conn * sc, "-", color="C0", lw=0.9, label="R_mean, connected" if base == 0 else None)
        ax.plot(qO, connO * sc, "-", color="C3", lw=0.9, label="R_out (envelope), connected" if base == 0 else None)
        m, e, _ = spsi(Lm, a)
        ax.plot(qc, m / qc ** 2 * sc, "ko", ms=4, label=r"contour: $S_\psi/q^2$" if base == 0 else None)
        ax.text(1.5e-3, full[0] * sc * 3, f"a={a}" + (" (×10³)" if base else ""), fontsize=8)
    ax.axvspan(Q_MAX_FIT, 10, color="0.93", zorder=-1)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("q"); ax.set_ylabel("S(q)")
    ax.set_title(f"polar projections vs contour, L = {Lm}"); ax.legend(fontsize=7)
    ax = axs[1, 2]
    for a, ls in ((As[0], "-"), (As[-1], "--")):
        for L in Ls:
            m = np.nanmean(S["ruler"][sel(L, a)], axis=0)
            df = 1 - np.gradient(np.log(m), np.log(EPS_RULER))
            ok = np.isfinite(df) & (EPS_RULER < 0.25 * S["R_eff"][sel(L, a)].mean())
            ax.plot(EPS_RULER[ok], df[ok], "o" + ls, color=cL[L], ms=3, lw=0.8,
                    label=f"a={a}, L={L}" if L in (Ls[0], Ls[-1]) else None)
    ax.axhline(1, color="k", lw=0.5); ax.axhline(4 / 3, color="k", lw=0.5, ls=":"); ax.axhline(7 / 4, color="k", lw=0.5, ls="--")
    ax.text(1.1, 4 / 3 + 0.01, "4/3", fontsize=7); ax.text(1.1, 7 / 4 + 0.01, "7/4 (percolation hull)", fontsize=7)
    ax.set_xscale("log"); ax.set_xlabel(r"ruler $\varepsilon$"); ax.set_ylabel(r"$d_f(\varepsilon) = 1 - d\ln P/d\ln\varepsilon$")
    ax.set_title("divider dimension (solid a=0, dashed a=1)"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(d, "contours_roughness.png"), dpi=110); plt.close(fig)

    report.append("")
    report.append(f"zeta from the contour, S_psi ~ q^(1-2 zeta); fits n >= {N_MIN_FIT} (q >= 2 pi {N_MIN_FIT}/P) to q = {Q_MAX_FIT}")
    report.append(f"  [all]   whole window;  [<q_ov] q < {q_ov:.4f} (scales > l_ov);  [>q_ov] {q_ov:.4f} < q < {Q_MAX_FIT}")
    for a in As:
        parts = []
        for L in Ls:
            z, ze = zeta_fit(L, a)
            z1, z1e = zeta_fit(L, a, qhi=q_ov)
            z2, z2e = zeta_fit(L, a, qlo=q_ov)
            parts.append(f"L={L}: {z:.2f}({ze * 100:.0f}) [<{z1:.2f}({z1e * 100:.0f}) >{z2:.2f}({z2e * 100:.0f})]")
        report.append(f"  a={a:4.2f}  " + "  ".join(parts))
    report.append(f"zeta from the polar R_mean profile (KFR Eq. 9 analogue), connected S(q), fit 12/<R> < q < {Q_MAX_FIT}:")
    for a in As:
        row = []
        for L in Ls:
            q, conn, full, err, R = connected_polar_S(S, "cR", L, a)
            mb, eb = bin_log_err(q, conn, err)
            sl, se = fit_power(qc, mb, eb, 12 / R, Q_MAX_FIT)
            row.append(f"L={L}: {-(1 + sl) / 2:.2f}({se / 2 * 100:.0f})")
        report.append(f"  a={a:4.2f}  " + "  ".join(row))
    report.append("divider dimension d_f(eps) at the largest eps < R_eff/4:")
    for a in As:
        row = []
        for L in Ls:
            m = np.nanmean(S["ruler"][sel(L, a)], axis=0)
            df = 1 - np.gradient(np.log(m), np.log(EPS_RULER))
            ok = np.isfinite(df) & (EPS_RULER < 0.25 * S["R_eff"][sel(L, a)].mean())
            row.append(f"L={L}: {df[ok][-1]:.3f}@{EPS_RULER[ok][-1]:.0f}")
        report.append(f"  a={a:4.2f}  " + "  ".join(row))

    # ================================================================== anisotropy
    fig, axs = plt.subplots(2, 3, figsize=(16, 9.5))
    ax = axs[0, 0]
    ph = np.rad2deg((np.arange(N_PHI_BINS) + 0.5) * np.pi / N_PHI_BINS)
    for a in As:
        m, e = mean_err(S["Pphi_8"][sel(Ls[-1], a)])
        ax.errorbar(ph, m * N_PHI_BINS, e * N_PHI_BINS, fmt="o-", color=ca[a], ms=3, label=f"a={a}")
    ax.axhline(1, color="k", lw=0.5)
    ax.set_xlabel(r"normal angle $\phi_n$ to b [deg] (0, 180: edge; 90: screw)")
    ax.set_ylabel(r"$P(\phi_n)$ / uniform"); ax.set_title(f"wall orientation distribution, σ=8, L={Ls[-1]}")
    ax.legend(fontsize=8)
    ax = axs[0, 1]
    for a in As:
        for sig, ls in zip(SIGMAS, (":", "-", "--")):
            m, e = zip(*[mean_err(S[f"E_{sig:g}"][sel(L, a)]) for L in Ls])
            ax.errorbar(Ls, m, e, fmt="o" + ls, color=ca[a], ms=3, label=f"a={a}" if sig == 8 else None)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xscale("log"); ax.set_xlabel("L"); ax.set_ylabel(r"$E=\langle\cos 2\phi_n\rangle$ (+1 edge, −1 screw)")
    ax.set_title("wall orientation order (σ = 2 dotted, 8 solid, 32 dashed)"); ax.legend(fontsize=8)
    ax = axs[0, 2]
    angs = np.rad2deg(np.linspace(0, 2 * np.pi, N_SUPPORT, endpoint=False))
    fold = lambda v: 0.25 * (v + np.roll(v[::-1], 1) + np.roll(v, -N_SUPPORT // 2) + np.roll(np.roll(v[::-1], 1), -N_SUPPORT // 2))
    i90 = int(np.argmin(abs(angs - 90)))
    for a in As:
        sup = np.array([fold(v) for v in S["support"][sel(Ls[-1], a)]])
        m, e = mean_err(sup)
        k = angs <= 90
        ax.errorbar(angs[k], m[k], e[k], fmt="o-", color=ca[a], ms=3, label=f"a={a}")
    ax.set_xlabel(r"direction to b [deg] (0: along b, 90: ⊥ b)")
    ax.set_ylabel(r"support $h(\phi) - R_0$"); ax.set_title(f"advance of the most advanced wall point, L={Ls[-1]}")
    ax.legend(fontsize=8)
    ax = axs[1, 0]

    def adv(L, a):
        sup = np.array([fold(v) for v in S["support"][sel(L, a)]])
        return sup[:, i90], sup[:, 0]

    for a in As:
        m, e = zip(*[ratio_of_means(*adv(L, a)) for L in Ls])
        ax.errorbar(Ls, m, e, fmt="o-", color=ca[a], label=f"a={a}")
    ax.axhline(1, color="k", lw=0.5)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("L")
    ax.set_ylabel(r"$\langle h(90°)-R_0\rangle/\langle h(0°)-R_0\rangle$")
    ax.set_title("advance ⊥ b over advance ∥ b"); ax.legend(fontsize=8)
    ih = HP_SIGMAS.index(HP_SIGMA)

    def hpvar(L, a, si=ih):
        V = S["psi_hp_var"][sel(L, a)][:, si, :]
        return mean_err(V)

    ax = axs[1, 1]
    for L, mk in ((Ls[-2], "o--"), (Ls[-1], "o-")):
        for ci, (c, clr) in enumerate((("edge", "C0"), ("screw", "C3"), ("mixed", "0.5"))):
            m, e = zip(*[(hpvar(L, a)[0][ci], hpvar(L, a)[1][ci]) for a in As])
            ax.errorbar(As, m, e, fmt=mk, color=clr, ms=4, label=f"{c}, L={L}")
    ax.set_xlabel("anisotropy a"); ax.set_ylabel(rf"var($\psi-\psi_\sigma$), σ={HP_SIGMA:g} [rad²]")
    ax.set_title(f"slope fluctuations below {HP_SIGMA:g} sites, by wall class (σ={CLASS_SIGMA:g})"); ax.legend(fontsize=7)
    ax = axs[1, 2]
    for si, (sh, ls) in enumerate(zip(HP_SIGMAS, (":", "-.", "-", "--"))):
        for L, clr in ((Ls[-2], cL[Ls[-2]]), (Ls[-1], cL[Ls[-1]])):
            r = []
            for a in As:
                (m, e) = hpvar(L, a, si)
                rr = m[1] / m[0]
                r.append((rr, rr * np.hypot(e[1] / m[1], e[0] / m[0])))
            m, e = zip(*r)
            ax.errorbar(As, m, e, fmt="o" + ls, color=clr, ms=3, label=f"σ={sh:g}, L={L}")
    ax.axhline(1, color="k", lw=0.5)
    ax.set_xlabel("anisotropy a"); ax.set_ylabel("var screw / var edge")
    ax.set_title("screw/edge slope variance vs high-pass σ (σ=32 includes tip curvature)", fontsize=9)
    ax.legend(fontsize=6.5, ncol=2)
    fig.tight_layout(); fig.savefig(os.path.join(d, "contours_anisotropy.png"), dpi=110); plt.close(fig)

    report.append("")
    report.append("anisotropy: E = <cos 2 phi_n> (sigma = 8); advance h(90) and h(0) - R0 of the support function:")
    for a in As:
        row = []
        for L in Ls:
            h9, h0 = adv(L, a)
            row.append(f"L={L}: E={np.mean(S['E_8'][sel(L, a)]):+.3f} h90={h9.mean():.0f} h0={h0.mean():.0f}")
        report.append(f"  a={a:4.2f}  " + "  ".join(row))
    report.append("fraction of wall length per class (normal within 30 deg of b: edge; of b-perp: screw), sigma = 32:")
    for a in As:
        row = []
        for L in Ls[-2:]:
            row.append(f"L={L}: " + " ".join(f"{c}={np.mean(S[f'len_32_{c}'][sel(L, a)]):.3f}" for c in names))
        report.append(f"  a={a:4.2f}  " + "  ".join(row))
    report.append("slope variance var(psi - psi_sigma) on screw / edge wall (ratio), by high-pass sigma:")
    for a in As:
        row = []
        for L in Ls[-2:]:
            parts = []
            for si, sh in enumerate(HP_SIGMAS):
                m, e = hpvar(L, a, si)
                parts.append(f"s{sh:g}: {m[1] / m[0]:.2f}")
            m, e = hpvar(L, a)
            row.append(f"L={L}: edge {m[0]:.4f} screw {m[1]:.4f} [" + " ".join(parts) + "]")
        report.append(f"  a={a:4.2f}  " + "  ".join(row))
    report.append("intrinsic overhang fraction (sigma = 32) on edge / screw wall:")
    for a in As:
        row = []
        for L in Ls[-3:]:
            row.append(f"L={L}: edge {np.nanmean(S['ovi_32_edge'][sel(L, a)]):.4f} screw {np.nanmean(S['ovi_32_screw'][sel(L, a)]):.4f}")
        report.append(f"  a={a:4.2f}  " + "  ".join(row))
    txt = "\n".join(report)
    print(txt)
    with open(os.path.join(d, "contour_report.txt"), "w") as f:
        f.write(txt + "\n")


# ---------------------------------------------------------------------------
# synthetic checks
# ---------------------------------------------------------------------------
def synthetic_field(L, R, zeta, amp, rng, aspect=1.0):
    """phi = R(theta) - r on an L x L grid, R(theta) self-affine with S ~ q^-(1+2 zeta)."""
    nt = 8192
    n = np.arange(1, nt // 2)
    q = n / R
    S = q ** (-(1 + 2 * zeta)) * np.exp(-(q * 0.3) ** 2)
    c = np.zeros(nt // 2 + 1, complex)
    c[1:nt // 2] = np.sqrt(S / (2 * np.pi * R) / 2) * (rng.normal(size=len(n)) + 1j * rng.normal(size=len(n)))
    c[1] = 0
    dR = np.fft.irfft(c, nt) * nt
    dR *= amp / dR.std()
    th = np.linspace(0, 2 * np.pi, nt, endpoint=False)
    Rth = R * (1 + (aspect - 1) * np.sin(th) ** 2) + dR
    y, x = np.mgrid[0:L, 0:L].astype(float)
    xc = yc = (L - 1) / 2
    t = np.mod(np.arctan2(y - yc, x - xc), 2 * np.pi)
    return np.interp(t, th, Rth, period=2 * np.pi) - np.hypot(x - xc, y - yc), dR


def selftest():
    import contourpy
    rng = np.random.default_rng(1)
    L, R = 1024, 300.0
    print("synthetic walls, L=1024, R=300 (no overhangs): zeta from S_psi and from polar R(theta)")
    qc = np.sqrt(Q_EDGES[1:] * Q_EDGES[:-1])
    for zeta, amp in ((0.5, 1.5), (1.0, 4.0), (1.25, 6.0)):
        rows = []
        for _ in range(12):
            phi, _ = synthetic_field(L, R, zeta, amp, rng)
            xy = max(contourpy.contour_generator(z=phi).lines(0.0), key=len)
            rows.append(analyse_wall(xy, L, R, 4096))
        Sp = np.array([r["Spsi_bins"] for r in rows])
        m, e = mean_err(Sp)
        P = np.mean([r["P"] for r in rows])
        sl, se = fit_power(qc, m, e, 2 * np.pi * N_MIN_FIT / P, Q_MAX_FIT)
        cR = np.array([r["cR"] for r in rows])
        Rm = np.abs(cR[:, 0]).mean()
        nn = np.arange(cR.shape[1])
        conn = 2 * np.pi * Rm * ((np.abs(cR) ** 2).mean(0) - np.abs(cR.mean(0)) ** 2)
        mb = bin_log(nn[3:] / Rm, conn[3:])
        sl2, _ = fit_power(qc, mb, mb * 0.1, 12 / Rm, Q_MAX_FIT)
        E = np.mean([r["E_8"] for r in rows])
        print(f"  input zeta={zeta}: S_psi -> {(1 - sl) / 2:.3f} +- {se / 2:.3f};  polar R(theta) -> {-(1 + sl2) / 2:.3f};"
              f"  E={E:+.3f} (expect 0);  winding={np.mean([r['winding_total'] for r in rows]) / (2 * np.pi):.3f} turns;"
              f"  multivalued={np.mean([r['frac_multi'] for r in rows]):.4f};  intrinsic ovh={np.mean([r['ovi_32'] for r in rows]):.4f}")
    # anisotropic mean shape: E and support ratio
    phi, _ = synthetic_field(L, R, 1.0, 1.0, rng, aspect=1.5)
    xy = max(contourpy.contour_generator(z=phi).lines(0.0), key=len)
    r = analyse_wall(xy, L, R, 4096)
    sup = r["support"]
    print(f"  ellipse-like shape elongated along y (aspect 1.5): E={r['E_8']:+.3f} (expect > 0: long sides have normal along x = b);"
          f"  h(90)/h(0)-like: support(90)={sup[N_SUPPORT // 4]:.1f}, support(0)={sup[0]:.1f} (expect ~150, ~0)")
    # ruler on a circle and on a Koch-like curve
    t = np.linspace(0, 2 * np.pi, 20000, endpoint=False)
    circ = np.c_[300 * np.cos(t), 300 * np.sin(t)]
    fine, _ = resample_closed(circ, 0.25)
    print(f"  ruler: circle R=300, P(eps=1)={ruler_length(fine, 1.0):.2f}, P(eps=64)={ruler_length(fine, 64.0):.2f}"
          f" (2 pi R = {2 * np.pi * 300:.2f})")
    # crossings on a curve with an overhang: a circle with a fold
    print("  ray crossings on the circle: mult all 1 ->", np.unique(ray_crossings(circ, np.zeros(2), 1024)[0]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=os.path.join(HERE, "runs", "fss_D0.04"))
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--recompute", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    fn = os.path.join(args.dir, "contour_stats.npz")
    if args.recompute or not os.path.exists(fn):
        compute(args.dir, args.workers)
    S = dict(np.load(fn))
    plot(args.dir, S)


if __name__ == "__main__":
    main()
