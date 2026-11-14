#!/usr/bin/env python3
"""
Bubble depinning in a disordered 2D Ginzburg-Landau model
=========================================================

Adapted from Kolton, Ferrero & Rosso, "Depinning free of the elastic
approximation" (arXiv:2306.13415), but in a *bubble* geometry and with
optionally *anisotropic* disorder.

Model (lattice version of Eq. 3 of the paper, generalised):

    d phi_ij / dt =  sum_{nn b} J_b (phi_nb - phi_ij)
                   + eps0 [ (1 + r_ij) phi_ij - phi_ij^3 ] + h + g_ij

  * phi = +1 inside the bubble, -1 outside; h > 0 favours +1 (bubble grows).
  * --noise dislocation  (default): anisotropic correlated noise of
                 Geslin, Legros, ..., Rosso, Demery ("Roughness" letter, SM
                 Secs. 2-3): the shear stress of a random alloy resolved on a
                 Burgers vector b = (cos beta, sin beta), smeared over the core
                 width w. Gaussian field with power spectrum (SM Eq. 21)
                     F(k) ~ (k . b_hat)^2 / |k| * exp(-|k|^2 w^2),
                 i.e. correlations given by SM Eqs. 12-15. A straight wall
                 segment parallel to b ("screw") feels zero integrated force,
                 one perpendicular to b ("edge") feels the maximum (SM Eq. 17).
                 Since f_p in the letter is a pinning *force*, by default it
                 enters as a random field g_ij = Delta * f_p / std(f_p)
                 (--coupling field): a wall moving by du gains force ~ 2 g du,
                 exactly like the line in Eq. 5 of the SM. With
                 --coupling bond it is used as r_ij instead (Gaussian, std
                 Delta; |r| > 1 sites become locally single-well).
  * r_ij       : quenched "random-bond" site disorder, |r| <= Delta < 1.
                 Its spatial correlations can be anisotropic: Gaussian
                 correlation lengths xi_x, xi_y (in lattice units), optionally
                 rotated by an angle theta. xi = 0 -> uncorrelated (paper),
                 xi = inf along one axis -> columnar disorder.
                 The one-point distribution is always uniform in [-Delta, Delta]
                 (Gaussian field mapped through its CDF), so only the
                 correlations are anisotropic, not the amplitude.
  * J_b        : elastic couplings on bonds. Mean c_x on x-bonds, c_y on
                 y-bonds, with optional independent random-coupling disorder
                 J = c_a (1 + Delta_Ja * eta), eta ~ U[-1, 1], Delta_Ja < 1.
                 This is a second, "bond"-type way to make disorder
                 anisotropic in x/y.

Algorithm (as in the paper): "variant Monte Carlo" — each elementary move sets
phi_ij to a root of the local force (a cubic), namely the first root reached
in the direction of the local force (= what local gradient flow would reach).
Checkerboard Gauss-Seidel sweeps, parallelised with numba.

Protocol:
  1. Initial bubble of radius R0 (or ellipse R0x, R0y) with a tanh profile.
  2. Full (bidirectional) relaxation at h = h_start  -> first metastable state.
     (Laplace pressure ~ sigma/R makes a small bubble shrink; if it collapses
     you need a larger h_start or R0.)
  3. Quasi-static ramp h -> h + dh with *forward-only* moves (justified by the
     no-passing / forward-moving theorems: from a metastable state at h the
     configuration only moves forward when h increases).
  4. Depinning = the domain wall reaches the system boundary (no metastable
     state inside the box). h_d is then refined by bisection.

Flat-wall geometry (--geometry flat): instead of a bubble, a straight wall
along x at height y0 = y0_frac * Ly, +1 below, -1 above, periodic in x and cut
in y (the y bonds across the box edge are removed, so there is a single wall).
The box is L along the wall and Ly = aspect * L (default 2, --aspect) in the
direction of motion, so that the wall can travel far before depinning.
The wall is pushed upwards and depins when the phi > 0 domain attached to the
bottom reaches the top band (the KFR geometry, with our border criterion). The
wall direction relative to b is set by --burgers-angle: 0 -> wall line || b
(screw wall, normal _|_ b), 90 -> wall line _|_ b (edge wall, normal || b).

There is also a plain explicit-Euler "dynamics" mode at constant h (like
Fig. 8 of the paper) to look at bubble growth.

Usage examples
--------------
  python bubble_depinning.py                                  # dislocation noise, b along x, w = 3
  python bubble_depinning.py --burgers-angle 90 --w 3         # b along y
  python bubble_depinning.py --check-noise                    # compare noise with SM Eqs. 12-13
  python bubble_depinning.py --noise white                    # isotropic, uncorrelated (Kolton et al.)
  python bubble_depinning.py --noise gaussian --xi-x 8 --xi-y 1   # Gaussian anisotropic correlations
  python bubble_depinning.py --geometry flat --burgers-angle 90 --L 256   # flat edge wall, 256 x 512 box
  python bubble_depinning.py --geometry flat --L 256 --aspect 4           # 256 x 1024 box
  python bubble_depinning.py --xi-x inf --xi-y 0              # columnar disorder along x
  python bubble_depinning.py --xi-x 6 --xi-y 1 --theta 30     # rotated anisotropy
  python bubble_depinning.py --dJx 0.6 --dJy 0.0              # anisotropic bond disorder
  python bubble_depinning.py --mode dynamics --h 0.05 --tmax 400
"""

from __future__ import annotations

import argparse
import math
import os
import time
from dataclasses import dataclass, field, asdict

import numpy as np
from numba import njit, prange
from scipy.ndimage import map_coordinates, label
from scipy.special import erf


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
@dataclass
class Params:
    L: int = 256                 # linear size (must be even): Lx = L; Ly = L (bubble) or aspect * L (flat)
    eps0: float = 1.0
    cx: float = 1.0              # mean elastic coupling on x-bonds
    cy: float = 1.0              # mean elastic coupling on y-bonds
    # disorder type and strength
    noise: str = "dislocation"   # dislocation | gaussian | white
    coupling: str = "auto"       # field | bond | auto (field for dislocation, bond otherwise)
    Delta: float | None = None   # strength; default 0.08 (field) or 0.6 (bond)
    # dislocation noise (letter SM Eq. 21)
    w: float = 3.0               # core width in lattice units (letter: w = 1.5 b, grid b/2 -> 3)
    burgers_angle: float = 0.0   # direction of b w.r.t. x [deg]
    aniso: float = 1.0           # 0 = isotropic, 1 = full dislocation noise (see dislocation_spectrum)
    # gaussian noise: correlation lengths
    xi_x: float = 0.0            # correlation length along x (0: uncorrelated, inf: columnar)
    xi_y: float = 0.0            # correlation length along y
    theta: float = 0.0           # rotation of the correlation ellipse [deg]
    # bond (coupling) disorder
    dJx: float = 0.0             # J_x = cx (1 + dJx * eta), eta ~ U[-1,1]
    dJy: float = 0.0
    # geometry: "bubble" (circular domain) or "flat" (straight wall along x)
    geometry: str = "bubble"
    y0_frac: float = 0.1         # flat: initial wall height / Ly
    aspect: float = 2.0          # flat: box height Ly = aspect * L (rounded to even); L = wall length
    # initial bubble
    R0: float = 30.0
    R0x: float | None = None     # if given, elliptic initial bubble
    R0y: float | None = None
    # protocol
    h_start: float = 0.01
    dh: float = 0.001
    h_max: float = 0.5
    n_bisect: int = 8
    tol: float = 1e-6            # convergence: max |dphi| per sweep
    max_sweeps: int = 2_000_000
    margin: int = 3              # depinned when phi > 0 within `margin` of the border
    seed: int = 0
    # dynamics mode
    dt: float = 0.05
    tmax: float = 200.0
    n_frames: int = 10

    @property
    def Ly(self) -> int:
        """Box size along y (rows): L for bubbles, aspect * L rounded to even for flat walls."""
        return 2 * int(round(self.aspect * self.L / 2)) if self.geometry == "flat" else self.L


# ---------------------------------------------------------------------------
# Disorder
# ---------------------------------------------------------------------------
def _shape(L):
    """(Ly, Lx) from an int (square box) or a (Ly, Lx) tuple."""
    return (int(L), int(L)) if np.isscalar(L) else (int(L[0]), int(L[1]))


def _kgrid(Ly, Lx):
    """Wavevector grids ky, kx with array index [i=y, j=x]."""
    return np.meshgrid(2 * np.pi * np.fft.fftfreq(Ly), 2 * np.pi * np.fft.fftfreq(Lx), indexing="ij")


def anisotropic_gaussian_field(L, xi_x, xi_y, theta_deg, rng):
    """Unit-variance Gaussian field with anisotropic Gaussian correlations.

    C(dx, dy) ~ exp(-dx'^2/(2 xi_x^2) - dy'^2/(2 xi_y^2)) in rotated
    coordinates (approx. correlation length xi along each principal axis).
    xi = 0 -> white along that axis; xi = inf -> perfectly correlated (columnar).
    """
    Ly, Lx = _shape(L)
    inf_x, inf_y = math.isinf(xi_x), math.isinf(xi_y)
    if (inf_x or inf_y) and theta_deg % 180 != 0:
        raise ValueError("Infinite correlation length only supported with theta = 0 (mod 180).")
    if inf_x and inf_y:
        return np.full((Ly, Lx), rng.standard_normal())
    if inf_x:   # constant along x: each row (fixed y) is one value -> stripes along x
        col = anisotropic_gaussian_field_1d(Ly, xi_y, rng)
        return np.repeat(col[:, None], Lx, axis=1)
    if inf_y:   # constant along y
        row = anisotropic_gaussian_field_1d(Lx, xi_x, rng)
        return np.repeat(row[None, :], Ly, axis=0)

    w = rng.standard_normal((Ly, Lx))
    if xi_x == 0 and xi_y == 0:
        return w
    ky, kx = _kgrid(Ly, Lx)
    th = np.deg2rad(theta_deg)
    kpar = kx * np.cos(th) + ky * np.sin(th)      # along the rotated "x" axis
    kperp = -kx * np.sin(th) + ky * np.cos(th)
    filt = np.exp(-0.5 * (kpar**2 * xi_x**2 + kperp**2 * xi_y**2) / 2)
    g = np.real(np.fft.ifft2(np.fft.fft2(w) * filt))
    g -= g.mean()
    return g / g.std()


def anisotropic_gaussian_field_1d(L, xi, rng):
    w = rng.standard_normal(L)
    if xi == 0:
        return w
    k = 2 * np.pi * np.fft.fftfreq(L)
    g = np.real(np.fft.ifft(np.fft.fft(w) * np.exp(-0.25 * k**2 * xi**2)))
    g -= g.mean()
    return g / g.std()


def dislocation_spectrum(L, w, beta_deg, aniso=1.0):
    """Power spectrum of the pinning force, interpolating between isotropic and
    SM Eq. 21:
        S_a(k) ~ |k| exp(-k^2 w^2) [1 + a cos 2(chi - beta)]
               = [(1+a) (k.b)^2 + (1-a) (k.b_perp)^2] / |k| * exp(-k^2 w^2)
    chi = angle of k. a = 1 is SM Eq. 21 ((k.b)^2/|k|), a = 0 its angular
    average (isotropic, same radial spectrum and variance). Equivalently the
    field is sqrt((1+a)/2) f_b + sqrt((1-a)/2) f_bperp with f_b, f_bperp
    independent dislocation fields for b and b rotated by 90 degrees.
    Array index [i=y, j=x]; L is an int (square) or (Ly, Lx).
    Normalised so that the field has unit variance."""
    ky, kx = _kgrid(*_shape(L))
    b = np.deg2rad(beta_deg)
    kk = np.hypot(kx, ky)
    kk[0, 0] = 1.0
    kb = kx * np.cos(b) + ky * np.sin(b)
    kp = -kx * np.sin(b) + ky * np.cos(b)
    S = ((1 + aniso) * kb**2 + (1 - aniso) * kp**2) / kk * np.exp(-(kk * w) ** 2)
    S[0, 0] = 0.0
    return S / S.mean()          # var = mean(S) for the construction below


def dislocation_noise(L, w, beta_deg, rng, aniso=1.0):
    """Gaussian field with spectrum SM Eq. 21 (unit variance on average).
    FFT of real white noise gives Hermitian Gaussian amplitudes with <|W|^2> = N;
    multiplying by sqrt(S) is equivalent to the letter's SM Sec. 3 recipe."""
    S = dislocation_spectrum(L, w, beta_deg, aniso)
    W = np.fft.fft2(rng.standard_normal(_shape(L)))
    return np.real(np.fft.ifft2(W * np.sqrt(S)))


def make_disorder(p: Params, rng):
    """Returns r (random bond), g (random field), Jx, Jy."""
    sh = (p.Ly, p.L)                                              # rows (y), columns (x)
    if p.noise == "dislocation":
        f = dislocation_noise(sh, p.w, p.burgers_angle, rng, p.aniso)
    elif p.noise == "gaussian":
        f = anisotropic_gaussian_field(sh, p.xi_x, p.xi_y, p.theta, rng)
    else:
        f = rng.uniform(-np.sqrt(3), np.sqrt(3), size=sh)         # unit variance
    zero = np.zeros(sh)
    if p.coupling == "field":
        r, g = zero, p.Delta * f
    elif p.noise == "white":
        r, g = p.Delta / np.sqrt(3) * f, zero                     # uniform in [-Delta, Delta]
    elif p.noise == "gaussian":
        r, g = p.Delta * erf(f / np.sqrt(2)), zero                # uniform marginal, as before
    else:
        r, g = p.Delta * f, zero                                  # Gaussian, std Delta
    Jx = p.cx * (1 + p.dJx * rng.uniform(-1, 1, size=sh))
    Jy = p.cy * (1 + p.dJy * rng.uniform(-1, 1, size=sh))
    if p.geometry == "flat":
        Jy[-1, :] = 0.0          # cut the bonds between row L-1 and row 0: open in y, one wall
        g = g.copy()
        g[:2, :] += 1.0          # pin the bottom rows to +1 so the wall cannot leave through the bottom
    return r.astype(np.float64), g.astype(np.float64), Jx.astype(np.float64), Jy.astype(np.float64)


# ---------------------------------------------------------------------------
# Initial condition
# ---------------------------------------------------------------------------
def bubble(p: Params):
    L = p.L
    delta = math.sqrt(0.5 * (p.cx + p.cy) / p.eps0)       # wall width
    y, x = np.mgrid[0:p.Ly, 0:L].astype(float)
    if p.geometry == "flat":
        return np.tanh((p.y0_frac * p.Ly - y) / (math.sqrt(2) * delta))
    xc = yc = (L - 1) / 2
    Rx = p.R0x if p.R0x is not None else p.R0
    Ry = p.R0y if p.R0y is not None else p.R0
    # signed distance (approximate for ellipses)
    rho = np.sqrt(((x - xc) / Rx) ** 2 + ((y - yc) / Ry) ** 2)
    d = (1 - rho) * math.sqrt(Rx * Ry)
    return np.tanh(d / (math.sqrt(2) * delta))


# ---------------------------------------------------------------------------
# Numba kernels
# ---------------------------------------------------------------------------
@njit(cache=True, inline="always")
def _select_root(x, f, eps0, p, q):
    """Root of eps0*t^3 + p*t + q = 0 first reached from x moving along sign(f)."""
    P = p / eps0
    Q = q / eps0
    D = 0.25 * Q * Q + (P * P * P) / 27.0
    if D > 0.0:
        sD = math.sqrt(D)
        return np.cbrt(-0.5 * Q + sD) + np.cbrt(-0.5 * Q - sD)
    if P == 0.0:
        return 0.0
    m = 2.0 * math.sqrt(-P / 3.0)
    arg = (3.0 * Q / (2.0 * P)) * math.sqrt(-3.0 / P)
    if arg > 1.0:
        arg = 1.0
    elif arg < -1.0:
        arg = -1.0
    phi0 = math.acos(arg) / 3.0
    best = x
    bestd = 1e300
    for k in range(3):
        t = m * math.cos(phi0 - 2.0 * math.pi * k / 3.0)
        dd = (t - x) if f > 0 else (x - t)
        if dd >= -1e-12 and dd < bestd:
            bestd = dd
            best = t
    return best


@njit(parallel=True, cache=True)
def _sweep(phi, r, g, Jx, Jy, h, eps0, forward_only, color):
    Ly, Lx = phi.shape
    rowmax = np.zeros(Ly)
    for i in prange(Ly):
        ip = i + 1 if i + 1 < Ly else 0
        im = i - 1 if i > 0 else Ly - 1
        m = 0.0
        start = (i + color) & 1
        for j in range(start, Lx, 2):
            jp = j + 1 if j + 1 < Lx else 0
            jm = j - 1 if j > 0 else Lx - 1
            jr = Jx[i, j]
            jl = Jx[i, jm]
            ju = Jy[i, j]
            jd = Jy[im, j]
            S = jr * phi[i, jp] + jl * phi[i, jm] + ju * phi[ip, j] + jd * phi[im, j]
            Js = jr + jl + ju + jd
            x = phi[i, j]
            a = eps0 * (1.0 + r[i, j])
            hl = h + g[i, j]
            f = S + hl - Js * x + a * x - eps0 * x * x * x
            if f == 0.0 or (forward_only and f < 0.0):
                continue
            new = _select_root(x, f, eps0, Js - a, -(S + hl))
            dd = abs(new - x)
            phi[i, j] = new
            if dd > m:
                m = dd
        rowmax[i] = m
    return rowmax.max()


@njit(parallel=True, cache=True)
def _force(phi, r, g, Jx, Jy, h, eps0):
    Ly, Lx = phi.shape
    out = np.empty_like(phi)
    for i in prange(Ly):
        ip = i + 1 if i + 1 < Ly else 0
        im = i - 1 if i > 0 else Ly - 1
        for j in range(Lx):
            jp = j + 1 if j + 1 < Lx else 0
            jm = j - 1 if j > 0 else Lx - 1
            x = phi[i, j]
            out[i, j] = (Jx[i, j] * (phi[i, jp] - x) + Jx[i, jm] * (phi[i, jm] - x)
                         + Jy[i, j] * (phi[ip, j] - x) + Jy[im, j] * (phi[im, j] - x)
                         + eps0 * ((1.0 + r[i, j]) * x - x * x * x) + h + g[i, j])
    return out


def bubble_mask(phi):
    """The phi > 0 domain containing the centre of the box (the bubble), so that
    domains nucleated elsewhere by a strong random field are ignored.
    Falls back to the largest domain if the centre is not inside phi > 0."""
    lab, n = label(phi > 0)
    if n == 0:
        return np.zeros(phi.shape, bool)
    if n == 1:
        return lab == 1
    c = lab[phi.shape[0] // 2, phi.shape[1] // 2]
    if c > 0:
        return lab == c
    sizes = np.bincount(lab.ravel())
    sizes[0] = 0
    return lab == sizes.argmax()


def wall_mask(phi):
    """Flat geometry: the phi > 0 domain attached to the bottom row."""
    lab, n = label(phi > 0)
    if n == 0:
        return np.zeros(phi.shape, bool)
    ids = np.unique(lab[0][lab[0] > 0])
    return np.isin(lab, ids)


def domain_mask(phi, geometry="bubble"):
    return wall_mask(phi) if geometry == "flat" else bubble_mask(phi)


def touches_border(phi, margin, geometry="bubble"):
    if geometry == "flat":                       # the wall reached the top band
        return bool(wall_mask(phi)[-margin:].any())
    m = bubble_mask(phi)
    return bool(m[:margin].any() or m[-margin:].any() or m[:, :margin].any() or m[:, -margin:].any())


# ---------------------------------------------------------------------------
# Relaxation to a metastable state
# ---------------------------------------------------------------------------
def relax(phi, r, g, Jx, Jy, h, p: Params, forward_only, check_every=50):
    """Relax phi in place. Returns 'pinned', 'depinned' (wall hit the border),
    'collapsed' (bubble vanished) or 'maxiter'."""
    for s in range(p.max_sweeps):
        m0 = _sweep(phi, r, g, Jx, Jy, h, p.eps0, forward_only, 0)
        m1 = _sweep(phi, r, g, Jx, Jy, h, p.eps0, forward_only, 1)
        if max(m0, m1) < p.tol:
            if phi.max() < 0:
                return "collapsed", s + 1
            if touches_border(phi, p.margin, p.geometry):
                return "depinned", s + 1
            return "pinned", s + 1
        if s % check_every == 0:
            if touches_border(phi, p.margin, p.geometry):
                return "depinned", s + 1
            if phi.max() < 0:
                return "collapsed", s + 1
    return "maxiter", p.max_sweeps


# ---------------------------------------------------------------------------
# Geometry analysis (polar analogue of Eqs. 9 and 11 of the paper)
# ---------------------------------------------------------------------------
def analyse(phi, n_theta=1024):
    L = phi.shape[0]
    inside = bubble_mask(phi)
    phi = np.where(inside, phi, -np.abs(phi))     # hide other domains from the ray analysis
    A = inside.sum()
    if A == 0:
        return None
    y, x = np.nonzero(inside)
    xc, yc = x.mean(), y.mean()
    # gyration tensor -> aspect ratio and orientation of the bubble
    cxx = ((x - xc) ** 2).mean(); cyy = ((y - yc) ** 2).mean(); cxy = ((x - xc) * (y - yc)).mean()
    ev, evec = np.linalg.eigh(np.array([[cxx, cxy], [cxy, cyy]]))
    aspect = math.sqrt(ev[1] / max(ev[0], 1e-12))
    orient = math.degrees(math.atan2(evec[1, 1], evec[0, 1])) % 180

    # rays from the centre: count phi=0 crossings
    rmax = min(xc, yc, L - 1 - xc, L - 1 - yc)
    rs = np.arange(0, rmax, 0.25)
    th = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
    R_mean = np.full(n_theta, np.nan)
    mult = np.zeros(n_theta, dtype=int)
    var = np.zeros(n_theta)
    for c0 in range(0, n_theta, 512):               # in chunks: L = 4096 has ~10^8 ray points
        t = th[c0:c0 + 512]
        X = xc + np.outer(np.cos(t), rs)
        Y = yc + np.outer(np.sin(t), rs)
        vals = map_coordinates(phi, [Y.ravel(), X.ravel()], order=1, mode="wrap").reshape(X.shape)
        s = np.sign(vals)
        for kk in range(len(t)):
            idx = np.nonzero(s[kk, :-1] * s[kk, 1:] < 0)[0]
            if len(idx) == 0:
                continue
            v0, v1 = vals[kk, idx], vals[kk, idx + 1]
            rc = rs[idx] + (rs[idx + 1] - rs[idx]) * v0 / (v0 - v1)    # linear interpolation
            k = c0 + kk
            mult[k] = len(rc)
            R_mean[k] = rc.mean()                                      # Eq. (9) analogue
            var[k] = rc.var()                                          # Eq. (11) analogue
    ok = mult > 0
    Rm = np.where(ok, R_mean, np.nanmean(R_mean))
    dR = Rm - Rm.mean()
    Sq = np.abs(np.fft.rfft(dR)) ** 2                               # structure factor vs angular mode
    a2 = 2 * np.abs(np.fft.rfft(Rm)[2]) / n_theta                   # elliptic (n=2) Fourier amplitude
    return dict(area=int(A), xc=xc, yc=yc, R_eff=math.sqrt(A / np.pi),
                aspect=aspect, orientation_deg=orient, a2=a2,
                frac_multivalued=float((mult > 1).mean()),
                u_o=float(math.sqrt(var.mean())),
                theta=th, R_theta=R_mean, mult=mult, S_n=Sq)


def analyse_flat(phi):
    """Flat geometry: wall height profile u(x) = mean of the phi = 0 crossings in
    each column (KFR Eq. 9) of the domain attached to the bottom, with the
    crossing multiplicity (> 1: overhang) and u_o (KFR Eq. 11)."""
    L = phi.shape[1]
    m = wall_mask(phi)
    if not m.any():
        return None
    f = np.where(m, phi, -np.abs(phi))
    s = np.sign(f)
    u = np.full(L, np.nan)
    mult = np.zeros(L, int)
    var = np.zeros(L)
    for j in range(L):
        col = f[:, j]
        idx = np.nonzero(s[:-1, j] * s[1:, j] < 0)[0]
        if len(idx) == 0:
            continue
        v0, v1 = col[idx], col[idx + 1]
        yc = idx + v0 / (v0 - v1)
        mult[j] = len(yc)
        u[j] = yc.mean()
        var[j] = yc.var()
    ok = mult > 0
    return dict(area=int(m.sum()), mean_height=float(np.nanmean(u)), width=float(np.nanstd(u)),
                frac_multivalued=float((mult > 1).mean()), u_o=float(np.sqrt(var[ok].mean())) if ok.any() else np.nan,
                u=u, mult=mult)


# ---------------------------------------------------------------------------
# Backends: where phi lives and how it is relaxed
# ---------------------------------------------------------------------------
class CPUBackend:
    """numba (multi-threaded) CPU backend; phi arrays are numpy arrays."""

    def __init__(self, r, g, Jx, Jy, p):
        self.r, self.g, self.Jx, self.Jy, self.p = r, g, Jx, Jy, p

    def put(self, a):
        return np.array(a, dtype=np.float64)

    def get(self, a):
        return a.copy()

    def copy(self, a):
        return a.copy()

    def mean_diff(self, a, b):
        return float((a - b).mean())

    def relax(self, phi, h, forward_only):
        return relax(phi, self.r, self.g, self.Jx, self.Jy, h, self.p, forward_only)


def make_backend(name, r, g, Jx, Jy, p, device=None):
    if name == "cpu":
        return CPUBackend(r, g, Jx, Jy, p)
    if name == "cuda":
        from gpu_backend import CUDABackend
        return CUDABackend(r, g, Jx, Jy, p, device=device)
    raise ValueError(f"unknown backend {name!r}")


# ---------------------------------------------------------------------------
# Quasi-static protocol
# ---------------------------------------------------------------------------
def run_depinning(p: Params, verbose=True, backend="cpu", keep_snapshots=True, device=None):
    rng = np.random.default_rng(p.seed)
    r, g, Jx, Jy = make_disorder(p, rng)
    be = make_backend(backend, r, g, Jx, Jy, p, device=device)
    phi = be.put(bubble(p))

    t0 = time.time()
    status, ns = be.relax(phi, p.h_start, forward_only=False)
    phi_h = be.get(phi)
    if verbose:
        print(f"[init] h={p.h_start:.5f}  {status} after {ns} sweeps  area={domain_mask(phi_h, p.geometry).sum()}")
    if status != "pinned":
        return dict(params=p, status_init=status, r=r, g=g, Jx=Jx, Jy=Jy, phi=phi_h,
                    h=np.array([]), area=np.array([]), dphi=np.array([]),
                    snapshots=[], h_d=np.nan, phi_c=phi_h)

    hs, areas, dphis = [p.h_start], [domain_mask(phi_h, p.geometry).sum()], [0.0]
    snaps = [(p.h_start, phi_h)] if keep_snapshots else []
    h = p.h_start
    h_fail = None
    while h < p.h_max:
        h_new = h + p.dh
        trial = be.copy(phi)
        status, ns = be.relax(trial, h_new, forward_only=True)
        if status == "pinned":
            dphis.append(be.mean_diff(trial, phi))     # Eq. (4)
            phi = trial; h = h_new
            phi_h = be.get(phi)
            hs.append(h); areas.append(domain_mask(phi_h, p.geometry).sum())
            if keep_snapshots:
                snaps.append((h, phi_h))
            if verbose:
                print(f"  h={h:.5f}  sweeps={ns:6d}  area={areas[-1]:7d}  <dphi>={dphis[-1]:.3e}")
        else:
            h_fail = h_new
            if verbose:
                print(f"  h={h_new:.5f}  -> {status}")
            break

    # bisection between last pinned h and h_fail
    if h_fail is not None:
        lo, hi, phi_lo = h, h_fail, be.copy(phi)
        for _ in range(p.n_bisect):
            mid = 0.5 * (lo + hi)
            trial = be.copy(phi_lo)
            status, _ = be.relax(trial, mid, forward_only=True)
            if status == "pinned":
                lo, phi_lo = mid, trial
            else:
                hi = mid
        h_d = 0.5 * (lo + hi)
        phi_c = be.get(phi_lo)
        if verbose:
            print(f"[depinning] h_d = {h_d:.6f} +- {0.5*(hi-lo):.1e}   ({time.time()-t0:.1f}s)")
    else:
        h_d, phi_c = np.nan, be.get(phi)
        if verbose:
            print(f"[warning] no depinning up to h_max={p.h_max}")

    return dict(params=p, status_init="pinned", r=r, g=g, Jx=Jx, Jy=Jy,
                h=np.array(hs), area=np.array(areas), dphi=np.array(dphis),
                snapshots=snaps, h_d=h_d, phi_c=phi_c, phi=be.get(phi))


# ---------------------------------------------------------------------------
# Direct time integration (explicit Euler) at constant h
# ---------------------------------------------------------------------------
def run_dynamics(p: Params, h, verbose=True):
    rng = np.random.default_rng(p.seed)
    r, g, Jx, Jy = make_disorder(p, rng)
    phi = bubble(p)
    dt_max = 1.0 / (2 * (Jx.max() + Jy.max()) + 2 * p.eps0)
    dt = min(p.dt, 0.9 * dt_max)
    nsteps = int(p.tmax / dt)
    every = max(1, nsteps // p.n_frames)
    frames = [(0.0, phi.copy())]
    for n in range(1, nsteps + 1):
        phi += dt * _force(phi, r, g, Jx, Jy, h, p.eps0)
        if n % every == 0:
            frames.append((n * dt, phi.copy()))
            if verbose:
                print(f"  t={n*dt:8.2f}  area={domain_mask(phi, p.geometry).sum()}")
            if touches_border(phi, p.margin, p.geometry):
                if verbose:
                    print("  wall reached the border, stopping")
                break
    return dict(params=p, r=r, g=g, Jx=Jx, Jy=Jy, frames=frames, h=h, dt=dt)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def active_disorder(res):
    p = res["params"]
    return (res["g"], "g (random field)") if p.coupling == "field" else (res["r"], "r (random bond)")


def noise_label(p):
    if p.noise == "dislocation":
        return f"dislocation, w={p.w}, b at {p.burgers_angle}°, a={p.aniso}, Δ={p.Delta}"
    if p.noise == "gaussian":
        return f"gaussian, ξx={p.xi_x}, ξy={p.xi_y}, θ={p.theta}°, Δ={p.Delta}"
    return f"white, Δ={p.Delta}"


def check_noise(p: Params, out_png, n_samples=20):
    """Reproduce Fig. 1 of the letter's SM: sample correlation map and its cuts
    along / across b, compared with SM Eqs. 12-13."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.special import ive
    rng = np.random.default_rng(p.seed)
    C = np.zeros((p.L, p.L))
    for _ in range(n_samples):
        f = dislocation_noise(p.L, p.w, p.burgers_angle, rng, p.aniso)
        F = np.fft.fft2(f)
        C += np.real(np.fft.ifft2(np.abs(F) ** 2)) / f.size
    C /= n_samples
    C0 = C[0, 0]
    C = np.fft.fftshift(C) / C0
    c = p.L // 2
    n = int(min(c - 1, 6 * p.w))
    s = np.arange(-n, n + 1)
    b = np.deg2rad(p.burgers_angle)
    # cuts along b (longitudinal, Eq. 12) and perpendicular to b (transverse, Eq. 13)
    lon = map_coordinates(C, [c + s * np.sin(b), c + s * np.cos(b)], order=1)
    tra = map_coordinates(C, [c + s * np.cos(b), c - s * np.sin(b)], order=1)
    z = s**2 / (8 * p.w**2)
    u = s**2 / (2 * p.w**2)
    E12 = (1 - u) * ive(0, z) + (1 + u) * ive(1, z)     # SM Eq. 12
    E13 = ive(0, z) - ive(1, z)                         # SM Eq. 13
    wa, wb = (1 + p.aniso) / 2, (1 - p.aniso) / 2       # mixture weights (aniso = 1 -> pure SM)
    F12 = wa * E12 + wb * E13
    F13 = wa * E13 + wb * E12
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
    im = axs[0].imshow(C[c - n:c + n + 1, c - n:c + n + 1], origin="lower", cmap="coolwarm",
                       extent=[-n, n, -n, n])
    plt.colorbar(im, ax=axs[0]); axs[0].set_xlabel("x"); axs[0].set_ylabel("y")
    axs[0].set_title(r"$\mathcal{F}(x,y)/\langle f_p^2\rangle$ (sample)")
    axs[1].plot(s, lon, "ko", ms=3, label="sample"); axs[1].plot(s, F12, "r-", label="SM Eq. 12" if p.aniso == 1 else "SM Eqs. 12-13 mix")
    axs[1].set_title("along b (longitudinal)"); axs[1].legend()
    axs[2].plot(s, tra, "ko", ms=3, label="sample"); axs[2].plot(s, F13, "b-", label="SM Eq. 13" if p.aniso == 1 else "SM Eqs. 12-13 mix")
    axs[2].set_title("perpendicular to b (transverse)"); axs[2].legend()
    fig.tight_layout(); fig.savefig(out_png, dpi=130); plt.close(fig)
    return float(np.max(np.abs(lon - F12))), float(np.max(np.abs(tra - F13)))


def plot_depinning(res, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm

    p = res["params"]
    fig = plt.figure(figsize=(15, 9.5))
    gs = fig.add_gridspec(2, 3)

    ax = fig.add_subplot(gs[0, 0])
    dis, name = active_disorder(res)
    ax.imshow(dis, cmap="RdBu_r", origin="lower", vmin=-2 * p.Delta, vmax=2 * p.Delta)
    ax.set_title(f"disorder {name}: {noise_label(p)}", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])

    ax = fig.add_subplot(gs[0, 1])
    ax.imshow(res["phi_c"], cmap="Greens", origin="lower", vmin=-1.5, vmax=1.2, alpha=0.5)
    snaps = res["snapshots"]
    sel = snaps[:: max(1, len(snaps) // 12)] + [(res["h_d"], res["phi_c"])]
    colors = cm.viridis(np.linspace(0, 1, len(sel)))
    for (hh, ph), col in zip(sel, colors):
        ax.contour(ph, levels=[0], colors=[col], linewidths=0.8, origin="lower")
    ax.set_title(rf"metastable walls, $h$ from {p.h_start} to $h_d$={res['h_d']:.4f}")
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])

    ax = fig.add_subplot(gs[0, 2])
    ax.imshow(res["phi_c"], cmap="RdBu_r", origin="lower", vmin=-1, vmax=1)
    ax.contour(res["phi_c"], levels=[0], colors="k", linewidths=0.8, origin="lower")
    ax.set_title("critical configuration (last pinned state)")
    ax.set_xticks([]); ax.set_yticks([])

    ax = fig.add_subplot(gs[1, 0])
    ax.plot(res["h"], res["area"], ".-")
    ax.axvline(res["h_d"], color="r", ls="--", label=r"$h_d$")
    ax.set_xlabel("h"); ax.set_ylabel("bubble area"); ax.legend()
    ax2 = ax.twinx()
    ax2.plot(res["h"][1:], res["dphi"][1:], "C1.", ms=3)
    ax2.set_ylabel(r"$\langle\Delta\phi\rangle$ per step", color="C1")

    if p.geometry == "flat":
        a = analyse_flat(res["phi_c"])
        if a is not None:
            ax = fig.add_subplot(gs[1, 1])
            ax.plot(a["u"], lw=0.8)
            mv = a["mult"] > 1
            ax.plot(np.nonzero(mv)[0], a["u"][mv], "r.", ms=2, label="multivalued")
            ax.set_xlabel("x"); ax.set_ylabel("u(x) (mean crossing)")
            ax.set_title(f"wall profile: width={a['width']:.1f}, multivalued frac={a['frac_multivalued']:.3f}", fontsize=9)
            ax.legend(fontsize=7)
            ax = fig.add_subplot(gs[1, 2])
            uu = np.nan_to_num(a["u"] - np.nanmean(a["u"]))
            Sq = np.abs(np.fft.rfft(uu)) ** 2 / len(uu)
            q = 2 * np.pi * np.arange(len(Sq)) / len(uu)
            ax.loglog(q[1:], Sq[1:], ".", ms=3)
            ax.set_xlabel("q"); ax.set_ylabel("S(q)"); ax.set_title("structure factor of u(x)")
        fig.tight_layout()
        fig.savefig(out_png, dpi=130)
        plt.close(fig)
        return a
    a = analyse(res["phi_c"])
    if a is not None:
        ax = fig.add_subplot(gs[1, 1], projection="polar")
        ax.plot(a["theta"], a["R_theta"], lw=0.8)
        mv = a["mult"] > 1
        ax.plot(a["theta"][mv], a["R_theta"][mv], "r.", ms=2, label="multivalued")
        ax.set_title(f"R(θ): aspect={a['aspect']:.2f}, orient={a['orientation_deg']:.0f}°,\n"
                     f"multivalued frac={a['frac_multivalued']:.2f}, u_o={a['u_o']:.2f}", fontsize=9, pad=28)
        ax.legend(loc="lower left", fontsize=7)

        ax = fig.add_subplot(gs[1, 2])
        n = np.arange(len(a["S_n"]))
        ax.loglog(n[1:], a["S_n"][1:], ".", ms=3)
        ax.set_xlabel("angular mode n"); ax.set_ylabel(r"$S(n)=|\tilde R_n|^2$")
        ax.set_title("structure factor of R(θ)")

    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    return a


def plot_dynamics(res, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm

    p = res["params"]
    fig, axs = plt.subplots(1, 2, figsize=(11, 5.5))
    dis, name = active_disorder(res)
    axs[0].imshow(dis, cmap="RdBu_r", origin="lower", vmin=-2 * p.Delta, vmax=2 * p.Delta)
    axs[0].set_title(f"disorder {name}: {noise_label(p)}", fontsize=9)
    frames = res["frames"]
    colors = cm.viridis(np.linspace(0, 1, len(frames)))
    for (t, ph), col in zip(frames, colors):
        axs[1].contour(ph, levels=[0], colors=[col], linewidths=0.8, origin="lower")
    axs[1].set_xlim(0, p.L); axs[1].set_ylim(0, p.Ly); axs[1].set_aspect("equal")
    axs[1].set_title(f"bubble growth at h={res['h']} (dt={res['dt']:.3f}, equally spaced times)")
    for a in axs:
        a.set_xticks([]); a.set_yticks([])
    fig.tight_layout(); fig.savefig(out_png, dpi=130); plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    d = Params()
    ap.add_argument("--mode", choices=["depinning", "dynamics"], default="depinning")
    for name, val in asdict(d).items():
        flag = "--" + name.replace("_", "-")
        if val is None:
            ap.add_argument(flag, type=float, default=None)
        else:
            ap.add_argument(flag, type=type(val), default=val)
    ap.add_argument("--h", type=float, default=0.05, help="field for --mode dynamics")
    ap.add_argument("--out", default="runs/bubble_run",
                    help="output prefix (folders are created; default: runs/bubble_run)")
    ap.add_argument("--backend", choices=["cpu", "cuda"], default="cpu",
                    help="cpu: numba multi-threaded; cuda: NVIDIA GPU via numba.cuda (depinning mode)")
    ap.add_argument("--gpu", type=int, default=None, help="CUDA device index (default: current)")
    ap.add_argument("--check-noise", action="store_true",
                    help="only generate dislocation noise and compare its correlations with SM Eqs. 12-13")
    args = ap.parse_args()

    p = Params(**{k: getattr(args, k) for k in asdict(d)})
    if p.L % 2:
        raise SystemExit("L must be even (checkerboard update).")
    if p.noise not in ("dislocation", "gaussian", "white"):
        raise SystemExit("--noise must be dislocation, gaussian or white.")
    if p.coupling == "auto":
        p.coupling = "field" if p.noise == "dislocation" else "bond"
    if p.coupling not in ("field", "bond"):
        raise SystemExit("--coupling must be field, bond or auto.")
    if p.Delta is None:
        p.Delta = 0.08 if p.coupling == "field" else 0.6
    if p.Delta < 0 or (p.coupling == "bond" and p.noise != "dislocation" and p.Delta >= 1):
        raise SystemExit("Need 0 <= Delta (< 1 for bounded random-bond disorder).")
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    if args.check_noise:
        e12, e13 = check_noise(p, args.out + "_noise_check.png")
        print(f"max deviation from SM Eq. 12: {e12:.3f}, Eq. 13: {e13:.3f}")
        print("saved", args.out + "_noise_check.png")
        return
    if p.dJx >= 1 or p.dJy >= 1:
        raise SystemExit("Need dJx, dJy < 1 so that couplings stay positive.")
    if p.geometry not in ("bubble", "flat"):
        raise SystemExit("--geometry must be bubble or flat.")

    print("Parameters:", asdict(p))
    if args.mode == "dynamics":
        res = run_dynamics(p, args.h)
        plot_dynamics(res, args.out + "_dynamics.png")
        np.savez_compressed(args.out + "_dynamics.npz", r=res["r"], g=res["g"],
                            t=np.array([f[0] for f in res["frames"]]),
                            phi=np.array([f[1] for f in res["frames"]], dtype=np.float32))
        print("saved", args.out + "_dynamics.png")
        return

    res = run_depinning(p, backend=args.backend, device=args.gpu)
    if res["status_init"] != "pinned":
        print(f"Initial relaxation ended with '{res['status_init']}'. "
              "If 'collapsed', increase --h-start or --R0; if 'depinned', lower --h-start.")
        return
    a = plot_depinning(res, args.out + ".png")
    np.savez_compressed(args.out + ".npz", r=res["r"], g=res["g"], Jx=res["Jx"], Jy=res["Jy"],
                        h=res["h"], area=res["area"], dphi=res["dphi"], h_d=res["h_d"],
                        phi_c=res["phi_c"].astype(np.float32),
                        snap_h=np.array([s[0] for s in res["snapshots"]]),
                        snap_phi=np.array([s[1] for s in res["snapshots"]], dtype=np.float32),
                        params=str(asdict(p)))
    if a is not None and p.geometry == "flat":
        print(f"critical wall: mean height={a['mean_height']:.1f}  width={a['width']:.2f}  "
              f"multivalued frac={a['frac_multivalued']:.3f}  u_o={a['u_o']:.2f}")
    elif a is not None:
        print(f"critical bubble: R_eff={a['R_eff']:.1f}  aspect={a['aspect']:.3f}  "
              f"orientation={a['orientation_deg']:.1f} deg  multivalued frac={a['frac_multivalued']:.3f}  "
              f"u_o={a['u_o']:.2f}")
    print("saved", args.out + ".png", args.out + ".npz")


if __name__ == "__main__":
    main()
