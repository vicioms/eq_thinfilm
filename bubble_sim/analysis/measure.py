"""One wall -> its observables, the same for flat walls and bubbles where they make sense.

Per wall:
  P        length of the wall (one period for a flat wall)
  size     L for a flat wall, R_eff = sqrt(area / pi) for a bubble
  travel   how far the wall moved before its last pinning (flat: mean height - start,
           bubble: R_eff - R0)
  psi      modes of the tangent angle psi(s) along the arc length (intrinsic)
  prof     modes of the single-valued profile (KFR Eq. 9): u(x) per column (flat),
           R(theta) per ray about the centroid (bubble)
  multi    fraction of columns / rays crossing the wall more than once (overhangs)
  u_o      rms spread of the crossings of one column / ray (KFR Eq. 11)
  ang      histogram of the angle between the wall line and b, [0, 90] degrees, per arc length
  shape    (bubbles) <R(theta)> / R_eff in N_SHAPE angular bins, theta measured from b
"""
from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from . import geometry as G
from .spectra import modes

ANG_BINS = np.linspace(0, 90, 19)
N_SHAPE = 72
SCALAR = ("P", "size", "travel", "multi", "u_o", "A")


def flat_wall(curve, u, mult, L, b_angle, travel):
    pts, P = G.resample(curve, closed=False, L=L)
    psi = G.tangent_angle(pts, P, closed=False, L=L)
    ok = np.isfinite(u)
    uu = np.interp(np.arange(L), np.nonzero(ok)[0], u[ok], period=L) if ok.sum() > 1 else np.zeros(L)
    ang = np.histogram(G.angle_to_b(pts, False, b_angle, L), ANG_BINS)[0] / len(pts)
    return dict(P=P, size=float(L), travel=travel, multi=float(np.mean(mult > 1)), u_o=np.nan, A=np.nan,
                psi=modes(psi, P), prof=modes(uu, L), ang=ang)


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
    return dict(P=P, size=R_eff, travel=R_eff - R0, multi=float(np.mean(mult > 1)),
                u_o=float(np.sqrt(var[good].mean())), A=A, xc=cen[0], yc=cen[1],
                psi=modes(psi, P), prof=modes(Rf, 2 * np.pi * R_eff), ang=ang, shape=shape)


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
    for k in ("psi", "prof"):
        W[k] = [r[k] for r in rows]
    W["ang"] = np.array([r["ang"] for r in rows])
    if W["geometry"] == "bubble":
        W["shape"] = np.array([r["shape"] for r in rows])
    return W
