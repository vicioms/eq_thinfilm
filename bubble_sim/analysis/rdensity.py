"""Roughness density (Caballero & Giamarchi, "Roughness density", 2022): the roughness of a wall
without a single-valued height u(x).

The wall density along a column (flat wall: x = const; bubble: a ray at angle theta from the
centroid) is
    rho(x, y) = sum_k s_k delta(y - y_k(x)),
summed over all crossings y_k of the column with the wall, s_k = +-1 the direction in which the
wall crosses it (sum_k s_k = 1: one wall). This is the sharp-wall limit of the paper's
rho = d_y phi / (2 phi_0), so its form factor chi_0 = 1. Its Fourier transform along the column
    rho~(q, x) = sum_k s_k exp(-i q y_k(x))
gives the roughness density (paper Eq. 2, averaged over x and walls)
    chi(q, r) = < Re rho~(q, x) rho~*(q, x + r) >  =  exp(-q^2 B(r) / 2)   (Gaussian),
and B(r) is the slope of -2 ln chi vs q^2 at small q (paper Eq. 8): fitted as B q^2 + C q^4
on the NQFIT smallest q with YMIN <= -2 ln chi <= YMAX. (A fit over all of -2 ln chi <= YMAX
is pulled down by overhangs: at q > 1 / overhang size -2 ln chi saturates.)

Exactly, for q -> 0 (any distribution):
    B(r) = < [u(x + r) - u(x)]^2 > + 2 < M2 - u^2 >,   u = sum_k s_k y_k,  M2 = sum_k s_k y_k^2
with u the signed ("volume") height of the column. M2 - u^2 = 0 for a single crossing and > 0
with overhangs (2 a b for an overhang of the two lengths a, b), so B(0) measures the overhangs
and B(r) - B(0) the roughness of the signed height. Bdir is this exact limit, from the same
crossings, as a check of the fit.

For a bubble the radius is taken relative to the smooth shape (modes n < N_SMOOTH of u(theta)),
and r is the arc length R_eff * dtheta. chi is also kept per sector of directions from b (edge
flanks, screw tips): the average then runs over the pairs of rays both inside the sector.

Real-space B(r) and the exponent: B ∝ r^(2 zeta) only for zeta < 1. For zeta > 1 (qEW, 1.25)
the scaling is anomalous: B(r) ∝ r^2 L^(2 zeta - 2) for r << L, so the local slope saturates
at 2 and zeta shows in the L dependence of B at fixed r (B ∝ L^(2 zeta - 2)). For qKPZ
(zeta = 0.633) B ∝ r^1.27 and does not depend on L at fixed r. Both use B(r) - B(0): B(0) is the
overhang term (as large as B(32) for screw walls and bubbles), not roughness.
"""
from __future__ import annotations

import numpy as np

QS = np.geomspace(1e-4, 3.0, 64)
YMAX = 0.5
YMIN = 1e-9
NQFIT = 8
NR = 40
COL_OFF = 9 / 16           # column positions x = j + 9/16: never on a stored vertex (1/8 grid)


def r_grid(n):
    """Separations (in columns / rays) at which chi is kept: 0 and NR log-spaced up to n / 2."""
    return np.unique(np.r_[0, np.rint(np.geomspace(1, n // 2, NR)).astype(int)])


def _hits(a, b, va, vb, off):
    """Crossings of the segments (a -> b) in the column coordinate with the columns k + off:
    column k, interpolated value of v, sign of b - a."""
    k0 = np.ceil(np.minimum(a, b) - off).astype(np.int64)
    cnt = np.maximum(np.ceil(np.maximum(a, b) - off).astype(np.int64) - k0, 0)
    seg = np.repeat(np.arange(len(a)), cnt)
    k = np.repeat(k0, cnt) + (np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt))
    da = b[seg] - a[seg]
    f = (k + off - a[seg]) / np.where(da != 0, da, 1.0)
    return k, va[seg] + f * (vb[seg] - va[seg]), np.sign(da)


def column_crossings(curve, L):
    """Flat wall (open curve from x ~ 0 to x ~ L, periodic with shift (L, 0)): (column, y, sign)."""
    c = np.vstack([curve, curve[:1] + np.array([L, 0.0])])
    k, y, s = _hits(c[:-1, 0], c[1:, 0], c[:-1, 1], c[1:, 1], COL_OFF)
    return np.mod(k, L), y, s


def ray_crossings(pts, centre, nt):
    """Closed counter-clockwise curve: (ray, radius, sign) for nt rays from centre."""
    xy = np.vstack([pts, pts[:1]])
    d = xy - centre
    th = np.arctan2(d[:, 1], d[:, 0])
    r = np.hypot(d[:, 0], d[:, 1])
    T = np.r_[th[0], th[0] + np.cumsum((np.diff(th) + np.pi) % (2 * np.pi) - np.pi)] * nt / (2 * np.pi)
    k, rk, s = _hits(T[:-1], T[1:], r[:-1], r[1:], 0.0)
    return np.mod(k, nt), rk, s


def smooth(u, n_keep):
    c = np.fft.rfft(u)
    c[n_keep:] = 0
    return np.fft.irfft(c, len(u))


def _corr(f, g):
    """A(r) = sum_x f(x) g(x - r), periodic, real f and g."""
    return np.fft.ifft(np.fft.fft(f) * np.conj(np.fft.fft(g))).real


def chi_one(col, y, s, n, base=None, mask=None, rmax=None):
    """One wall: chi(q, r) on QS x r_grid(n), and the exact q -> 0 limit Bdir(r).
    mask (bool per column / ray): average only over the pairs with both columns in the mask
    (a sector of a bubble); r > rmax (in columns) is left nan, e.g. the separations that would
    pair the two opposite arcs of a sector."""
    if base is not None:
        y = y - base[col]
    rg = r_grid(n)
    u = np.bincount(col, s * y, n)
    o = np.bincount(col, s * y * y, n) - u * u           # overhang term M2 - u^2 per column
    m = np.ones(n) if mask is None else mask.astype(float)
    den = _corr(m, m)[rg]
    with np.errstate(invalid="ignore", divide="ignore"):
        bdir = (_corr(m * u * u, m) + _corr(m, m * u * u) - 2 * _corr(m * u, m * u)
                + _corr(m * o, m) + _corr(m, m * o))[rg] / den
    chi = np.empty((len(QS), len(rg)))
    for i, q in enumerate(QS):
        ph = q * y
        rho = m * (np.bincount(col, s * np.cos(ph), n) - 1j * np.bincount(col, s * np.sin(ph), n))
        F = np.fft.fft(rho)
        with np.errstate(invalid="ignore", divide="ignore"):
            chi[i] = np.fft.ifft(np.abs(F) ** 2).real[rg] / den
    bad = den < 0.5
    if rmax is not None:
        bad |= rg > rmax
    chi[:, bad] = np.nan
    bdir[bad] = np.nan
    return chi, bdir


def fit_B(chi, ymax=YMAX):
    """B(r): -2 ln chi = B q^2 + C q^4 on the NQFIT smallest q with YMIN <= -2 ln chi <= ymax;
    0 where -2 ln chi < YMIN at every q (single-valued walls at r = 0), nan with < 3 points."""
    with np.errstate(invalid="ignore", divide="ignore"):
        y = -2 * np.log(chi)
    B = np.full(y.shape[1], np.nan)
    for j in range(y.shape[1]):
        yj = y[:, j]
        if np.all(~np.isfinite(yj) | (yj < YMIN)) and np.isfinite(yj[0]):
            B[j] = 0.0
            continue
        i = np.nonzero(np.isfinite(yj) & (yj >= YMIN) & (yj <= ymax))[0][:NQFIT]
        if len(i) >= 3:
            q2 = QS[i] ** 2
            B[j] = np.linalg.lstsq(np.c_[q2, q2 * q2], yj[i], rcond=None)[0][0]
    return B, y


def local_slope(r, B, lo, hi):
    """d ln B / d ln r from a log-log fit over lo <= r <= hi."""
    ok = (r >= lo) & (r <= hi) & np.isfinite(B) & (B > 0)
    return np.polyfit(np.log(r[ok]), np.log(B[ok]), 1)[0] if ok.sum() >= 3 else np.nan


def local_slopes(r, B, half=np.log(2) / 2):
    """d ln B / d ln r at every r from a log-log fit over [r e^-half, r e^half] (>= 3 points)."""
    out = np.full(len(r), np.nan)
    for i, ri in enumerate(r):
        if ri > 0:
            out[i] = local_slope(r, B, ri * np.exp(-half), ri * np.exp(half))
    return out
