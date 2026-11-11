"""
CUDA backend for bubble_depinning.py (numba.cuda).

Same algorithm as the CPU kernel: checkerboard Gauss-Seidel, each site jumps
to the root of its local cubic force first reached along the force. All sites
of one colour are independent, so one colour = one kernel launch with one
thread per site (this is how Kolton, Ferrero & Rosso run it on GPUs).

Host <-> device traffic during a relaxation is limited to two tiny reads every
`check_every` sweeps (max |dphi| of the last sweep, and border/positive flags);
the full field is copied back only when the wall touches the border band (to
check, with connected-component labelling on the host, that it is really the
bubble and not an island nucleated elsewhere) and once per pinned field step.

Requirements on the cluster: an NVIDIA GPU, CUDA toolkit, and
    pip install numba-cuda        (or numba <= 0.61 with its built-in CUDA target)
Double precision is used throughout (fine on A100/H100/V100; on consumer
GeForce cards fp64 is much slower).

To test without a GPU:  NUMBA_ENABLE_CUDASIM=1 python test_gpu_backend.py
"""
from __future__ import annotations

import math

import numpy as np
from numba import cuda


@cuda.jit(device=True, inline=True)
def _cbrt(v):
    return math.copysign(abs(v) ** (1.0 / 3.0), v)


@cuda.jit(device=True, inline=True)
def _select_root(x, f, eps0, p, q):
    """Root of eps0*t^3 + p*t + q = 0 first reached from x moving along sign(f)."""
    P = p / eps0
    Q = q / eps0
    D = 0.25 * Q * Q + (P * P * P) / 27.0
    if D > 0.0:
        sD = math.sqrt(D)
        return _cbrt(-0.5 * Q + sD) + _cbrt(-0.5 * Q - sD)
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


@cuda.jit
def _sweep_kernel(phi, r, g, Jx, Jy, h, eps0, forward_only, color, measure, maxd):
    jj, i = cuda.grid(2)             # jj: half-column index, i: row (y)
    Ly, Lx = phi.shape
    if i >= Ly or jj >= Lx // 2:
        return
    j = 2 * jj + ((i + color) & 1)
    ip = i + 1 if i + 1 < Ly else 0
    im = i - 1 if i > 0 else Ly - 1
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
        return
    new = _select_root(x, f, eps0, Js - a, -(S + hl))
    phi[i, j] = new
    if measure:
        cuda.atomic.max(maxd, 0, abs(new - x))


@cuda.jit
def _flags_kernel(phi, margin, flags):
    """flags[0] = 1 if phi > 0 in the border band, flags[1] = 1 if phi > 0 anywhere."""
    j, i = cuda.grid(2)
    Ly, Lx = phi.shape
    if i >= Ly or j >= Lx:
        return
    if phi[i, j] > 0.0:
        flags[1] = 1
        if i < margin or i >= Ly - margin or j < margin or j >= Lx - margin:
            flags[0] = 1


class CUDABackend:
    """Device-resident disorder; phi arrays handled by the backend are device arrays."""

    def __init__(self, r, g, Jx, Jy, p, block=(32, 8), check_every=32, device=None):
        if device is not None:
            cuda.select_device(device)
        L = p.L
        self.p = p
        self.check_every = check_every
        self.r, self.g, self.Jx, self.Jy = (cuda.to_device(np.ascontiguousarray(a, dtype=np.float64))
                                            for a in (r, g, Jx, Jy))
        self.maxd = cuda.to_device(np.zeros(1))
        self.flags = cuda.to_device(np.zeros(2, np.int32))
        self._z1 = np.zeros(1)
        self._z2 = np.zeros(2, np.int32)
        self.block = block
        self.grid_half = ((L // 2 + block[0] - 1) // block[0], (L + block[1] - 1) // block[1])
        self.grid_full = ((L + block[0] - 1) // block[0], (L + block[1] - 1) // block[1])

    # array helpers -------------------------------------------------------
    def put(self, a):
        return cuda.to_device(np.ascontiguousarray(a, dtype=np.float64))

    def get(self, a):
        return a.copy_to_host()

    def copy(self, a):
        b = cuda.device_array_like(a)
        b.copy_to_device(a)
        return b

    def mean_diff(self, a, b):
        return float((a.copy_to_host() - b.copy_to_host()).mean())

    # relaxation ----------------------------------------------------------
    def _sweep(self, phi, h, forward_only, measure):
        for color in (0, 1):
            _sweep_kernel[self.grid_half, self.block](phi, self.r, self.g, self.Jx, self.Jy, h,
                                                     self.p.eps0, forward_only, color, measure, self.maxd)

    def _flags(self, phi):
        self.flags.copy_to_device(self._z2)
        _flags_kernel[self.grid_full, self.block](phi, self.p.margin, self.flags)
        f = self.flags.copy_to_host()
        return bool(f[0]), bool(f[1])

    def relax(self, phi, h, forward_only):
        from bubble_depinning import touches_border
        p, ce = self.p, self.check_every
        h = float(h)
        fo = bool(forward_only)
        for s in range(p.max_sweeps):
            measure = (s % ce == ce - 1)
            if measure:
                self.maxd.copy_to_device(self._z1)
            self._sweep(phi, h, fo, measure)
            if not measure:
                continue
            m = self.maxd.copy_to_host()[0]
            border, positive = self._flags(phi)
            if not positive:
                return "collapsed", s + 1
            if border and touches_border(phi.copy_to_host(), p.margin):
                return "depinned", s + 1
            if m < p.tol:
                return "pinned", s + 1
        return "maxiter", p.max_sweeps
