"""Curves: arc-length resampling, tangent angle, area, ray crossings, orientation.

Two kinds of wall:
  open    one period of a flat wall, x from 0 to L, repeating with shift (L, 0)
  closed  a bubble wall (made counter-clockwise)
"""
from __future__ import annotations

import numpy as np

DS = 0.5          # arc-length resampling step (lattice units)
TAN_HALF = 2      # tangent from the chord p[k + m] - p[k - m]: 2 m DS = 2 sites
CHORD = 8         # orientation from chords of CHORD * DS = 4 sites (the lattice steps of the contour average out)


def resample(curve, closed, L=None, ds=DS):
    """Uniform arc-length resampling. Returns (pts, P); for an open periodic wall the point
    after the last one is pts[0] + (L, 0); a closed curve is returned counter-clockwise."""
    xy = np.asarray(curve, float)
    if closed:
        if np.allclose(xy[0], xy[-1]):
            xy = xy[:-1]
        if signed_area(xy) < 0:
            xy = xy[::-1]
        c = np.vstack([xy, xy[:1]])
    else:
        c = np.vstack([xy, xy[:1] + np.array([L, 0.0])])
    s = np.r_[0, np.cumsum(np.hypot(*np.diff(c, axis=0).T))]
    P = s[-1]
    n = max(32, int(round(P / ds)))
    t = np.arange(n) * (P / n)
    return np.c_[np.interp(t, s, c[:, 0]), np.interp(t, s, c[:, 1])], P


def _shift(closed, L):
    return np.zeros(2) if closed else np.array([L, 0.0])


def chords(pts, k, closed, L=None):
    """pts[i + k] - pts[i - k] with the periodic wrap of the curve."""
    sh = _shift(closed, L)
    n = len(pts)
    ext = np.vstack([pts[-k:] - sh, pts, pts[:k] + sh])
    return ext[2 * k:2 * k + n] - ext[:n]


def tangent_angle(pts, P, closed, L=None, m=TAN_HALF):
    """Tangent angle psi(s) along the wall, unwrapped, minus its trend: the 2 pi s / P of a
    closed curve's winding, the mean of an open one."""
    d = chords(pts, m, closed, L)
    psi = np.unwrap(np.arctan2(d[:, 1], d[:, 0]))
    if closed:
        s = np.arange(len(pts)) * (P / len(pts))
        psi = psi - 2 * np.pi * s / P
    return psi - psi.mean()


def signed_area(xy):
    x, y = xy[:, 0], xy[:, 1]
    return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)


def centroid(xy):
    x, y = xy[:, 0], xy[:, 1]
    cr = x * np.roll(y, -1) - np.roll(x, -1) * y
    return np.array([np.sum((x + np.roll(x, -1)) * cr), np.sum((y + np.roll(y, -1)) * cr)]) / (6 * signed_area(xy))


def angle_to_b(pts, closed, b_angle_deg, L=None, k=CHORD // 2):
    """Angle between the wall line (chord over 2 k DS) and b, folded to [0, 90] degrees:
    0 = screw-like (line ∥ b), 90 = edge-like (line ⊥ b, normal ∥ b). One value per point,
    i.e. per arc length DS."""
    d = chords(pts, k, closed, L)
    a = np.degrees(np.arctan2(d[:, 1], d[:, 0])) - b_angle_deg
    a = np.mod(a, 180.0)
    return np.minimum(a, 180.0 - a)


def ray_crossings(xy, centre, n_theta):
    """All crossings of a closed polygon with n_theta rays from `centre` (KFR Eqs. 9, 11 in
    polar form). Per ray: multiplicity, mean, outermost radius and variance of the crossings."""
    xy = np.asarray(xy, float)
    if not np.allclose(xy[0], xy[-1]):
        xy = np.vstack([xy, xy[:1]])
    d = xy - centre
    th = np.arctan2(d[:, 1], d[:, 0])
    r = np.hypot(d[:, 0], d[:, 1])
    dth = (np.diff(th) + np.pi) % (2 * np.pi) - np.pi
    T = np.r_[th[0], th[0] + np.cumsum(dth)]
    dt = 2 * np.pi / n_theta
    a, b = T[:-1], T[1:]
    k0 = np.ceil(np.minimum(a, b) / dt).astype(np.int64)
    cnt = np.maximum(np.ceil(np.maximum(a, b) / dt).astype(np.int64) - k0, 0)
    seg = np.repeat(np.arange(len(a)), cnt)
    k = np.repeat(k0, cnt) + (np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt))
    f = (k * dt - a[seg]) / np.where(b[seg] != a[seg], b[seg] - a[seg], 1.0)
    rk = r[:-1][seg] + f * (r[1:][seg] - r[:-1][seg])
    kk = np.mod(k, n_theta)
    mult = np.bincount(kk, minlength=n_theta)
    s1 = np.bincount(kk, rk, n_theta)
    s2 = np.bincount(kk, rk * rk, n_theta)
    out = np.full(n_theta, -np.inf)
    np.maximum.at(out, kk, rk)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = s1 / mult
        var = np.maximum(s2 / mult - mean ** 2, 0)
    return mult, mean, out, var


def n_rays(L):
    """Rays per bubble: >= 2 per lattice unit of wall at R ~ L/4, a power of 2, >= 1024."""
    return max(1024, 2 ** int(np.ceil(np.log2(np.pi * L))))
