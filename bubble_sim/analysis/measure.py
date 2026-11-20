"""One wall -> its observables, the same for flat walls and bubbles where they make sense.

Per wall:
  P        length of the wall (one period for a flat wall)
  size     L for a flat wall, R_eff = sqrt(area / pi) for a bubble
  travel   how far the wall moved before its last pinning (flat: mean height - start,
           bubble: R_eff - R0)
  psi      modes of the tangent angle psi(s) along the arc length (intrinsic)
  hs       modes of the height along the arc length, h(s) (intrinsic, single-valued with overhangs):
           flat: y(s); bubble: distance |r(s) - centroid|. S_h ~ q^-(1+2 zeta), q conjugate to s
  prof     modes of the single-valued profile (KFR Eq. 9): u(x) per column (flat),
           R(theta) per ray about the centroid (bubble)
  multi    fraction of columns / rays crossing the wall more than once (overhangs)
  u_o      rms spread of the crossings of one column / ray (KFR Eq. 11)
  ang      histogram of the angle between the wall line and b, [0, 90] degrees, per arc length
  shape    (bubbles) <R(theta)> / R_eff in N_SHAPE angular bins, theta measured from b
  Pdir     (bubbles) wall length by direction: the arc length of the contour whose polar angle
           about the centroid, folded to [0, 90] degrees from b, falls in each of the N_DIR bins
           DIR_BINS (0 = edge flanks, normal ∥ b; 90 = screw tips, normal ⊥ b)
  chi      roughness density chi(q, r) (Caballero & Giamarchi; rdensity.py) on rdensity.QS x
           rdensity.r_grid(columns / rays), and Bdir(r) its exact q -> 0 limit
  chi_edge, chi_screw (bubbles) chi and Bdir averaged over the pairs of rays both in the sector
           RD_SECTORS of folded directions from b: edge flanks [0, 30), screw tips [60, 90]
  Sdir     the same for the smooth shape, R(theta) keeping only the modes n < N_SMOOTH (the
           shape modes left out of the zeta fits): Pdir / Sdir is the excess length of the
           roughness in that direction, without the elongation of the bubble
"""
from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from . import geometry as G
from . import rdensity as RD
from .spectra import modes

ANG_BINS = np.linspace(0, 90, 19)
N_SHAPE = 72
DIR_BINS = np.linspace(0, 90, 7)
N_SMOOTH = 12
RD_SECTORS = dict(edge=(0.0, 30.0), screw=(60.0, 90.01))   # folded direction from b (deg)
SCALAR = ("P", "size", "travel", "multi", "u_o", "A")


def flat_wall(curve, u, mult, L, b_angle, travel):
    pts, P = G.resample(curve, closed=False, L=L)
    psi = G.tangent_angle(pts, P, closed=False, L=L)
    ok = np.isfinite(u)
    uu = np.interp(np.arange(L), np.nonzero(ok)[0], u[ok], period=L) if ok.sum() > 1 else np.zeros(L)
    ang = np.histogram(G.angle_to_b(pts, False, b_angle, L), ANG_BINS)[0] / len(pts)
    chi, bdir = RD.chi_one(*RD.column_crossings(np.asarray(curve, float), L), L)
    return dict(P=P, size=float(L), travel=travel, multi=float(np.mean(mult > 1)), u_o=np.nan, A=np.nan,
                psi=modes(psi, P), hs=modes(pts[:, 1], P), prof=modes(uu, L), ang=ang, chi=chi, Bdir=bdir,
                rd_n=L)


def bubble(curve, L, R0, b_angle):
    pts, P = G.resample(curve, closed=True)
    A = G.signed_area(pts)
    R_eff = np.sqrt(A / np.pi)
    cen = G.centroid(pts)
    psi = G.tangent_angle(pts, P, closed=True)
    nt = G.n_rays(L)
    mult, Rm, _, var = G.ray_crossings(pts, cen, nt)
    good = mult > 0
    Rf = np.where(good, Rm, np.nanmean(Rm[good]))
    th = (2 * np.pi * np.arange(nt) / nt - np.deg2rad(b_angle)) % (2 * np.pi)
    k = (th / (2 * np.pi) * N_SHAPE).astype(int) % N_SHAPE
    shape = np.bincount(k, Rf, N_SHAPE) / np.maximum(np.bincount(k, minlength=N_SHAPE), 1) / R_eff
    ang = np.histogram(G.angle_to_b(pts, True, b_angle), ANG_BINS)[0] / len(pts)
    Pdir, Sdir = length_by_direction(pts, P, cen, Rf, b_angle)
    col, rk, sk = RD.ray_crossings(pts, cen, nt)
    base = RD.smooth(np.bincount(col, sk * rk, nt), N_SMOOTH)
    chi, bdir = RD.chi_one(col, rk, sk, nt, base)
    fold = _fold(np.degrees(2 * np.pi * np.arange(nt) / nt) - b_angle)
    sect = {}
    for name, (lo, hi) in RD_SECTORS.items():
        sect[name] = RD.chi_one(col, rk, sk, nt, base, mask=(fold >= lo) & (fold < hi),
                                rmax=nt * 2 * (hi - lo) / 360)
    return dict(P=P, size=R_eff, travel=R_eff - R0, multi=float(np.mean(mult > 1)),
                u_o=float(np.sqrt(var[good].mean())), A=A, xc=cen[0], yc=cen[1],
                psi=modes(psi, P), hs=modes(np.hypot(*(pts - cen).T), P), prof=modes(Rf, 2 * np.pi * R_eff),
                ang=ang, shape=shape,
                Pdir=Pdir, Sdir=Sdir, chi=chi, Bdir=bdir, rd_n=nt,
                **{f"chi_{k}": v[0] for k, v in sect.items()}, **{f"Bdir_{k}": v[1] for k, v in sect.items()})


def _fold(theta_deg):
    t = np.mod(theta_deg, 180.0)
    return np.minimum(t, 180.0 - t)


def length_by_direction(pts, P, cen, Rf, b_angle):
    """Arc length of the contour and of its smooth shape in each folded direction bin."""
    d = pts - cen
    f = _fold(np.degrees(np.arctan2(d[:, 1], d[:, 0])) - b_angle)
    Pdir = np.histogram(f, DIR_BINS, weights=np.full(len(pts), P / len(pts)))[0]
    nt = len(Rf)
    c = np.fft.rfft(Rf)
    c[N_SMOOTH:] = 0
    Rs = np.fft.irfft(c, nt)
    th = 2 * np.pi * np.arange(nt + 1) / nt
    Rs = np.r_[Rs, Rs[:1]]
    x, y = Rs * np.cos(th), Rs * np.sin(th)
    ds = np.hypot(np.diff(x), np.diff(y))
    fm = _fold(np.degrees(0.5 * (th[1:] + th[:-1])) - b_angle)
    Sdir = np.histogram(fm, DIR_BINS, weights=ds)[0]
    return Pdir, Sdir


def _bubble_job(args):
    return bubble(*args)


def _flat_job(args):
    return flat_wall(*args)


def measure_run(W, workers=None):
    """Adds the observables of every wall of a loaded run (data.load_flat / load_bubbles) to W."""
    n = len(W["L"])
    b = W["burgers_angle"]
    if W["geometry"] == "flat":
        jobs = [(W["curve"][i], W["u"][i], W["mult"][i], int(W["L"][i]), b, W["travel"][i]) for i in range(n)]
        fn = _flat_job
    else:
        jobs = [(W["curve"][i], int(W["L"][i]), W["R0"][i], b) for i in range(n)]
        fn = _bubble_job
    with ProcessPoolExecutor(workers or os.cpu_count()) as ex:
        rows = list(ex.map(fn, jobs, chunksize=8))
    for k in SCALAR:
        W[k] = np.array([r[k] for r in rows], float)
    for k in ("psi", "hs", "prof", "chi", "Bdir"):
        W[k] = [r[k] for r in rows]
    W["rd_n"] = np.array([r["rd_n"] for r in rows])
    W["rd_r"] = [RD.r_grid(n) for n in W["rd_n"]]
    W["ang"] = np.array([r["ang"] for r in rows])
    if W["geometry"] == "bubble":
        for k in ("shape", "Pdir", "Sdir"):
            W[k] = np.array([r[k] for r in rows])
        for k in RD_SECTORS:
            W[f"chi_{k}"] = [r[f"chi_{k}"] for r in rows]
            W[f"Bdir_{k}"] = [r[f"Bdir_{k}"] for r in rows]
    return W
