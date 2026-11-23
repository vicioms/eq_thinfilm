"""Diagnostic for sq.py: single-site force autocorrelation C(t) and displacement spread
B(t) measured in the hot sliding state (random field, m = 1024, f = 0.12), compared with
(i) the bare C = Delta(v t), (ii) the Gaussian smearing with the measured B(t), and
(iii) the closure B = T b(t).  Usage: python3 sq_diag.py [m f]. Writes data/sq_diag_m{m}_f{f}.npz.
"""
import numpy as np
from protocols import _run, default_dt
from numba import njit
from sim import _g


@njit(cache=True)
def site_force(U, seeds, sites, A):
    S = U.shape[0]
    out = np.zeros((S, sites.shape[0]))
    for s in range(S):
        for k in range(sites.shape[0]):
            i = sites[k]
            n = int(np.floor(U[s, i])); x = U[s, i] - n
            out[s, k] = _g(seeds[s], i, n, A, 0) * (1 - x) + _g(seeds[s], i, n + 1, A, 0) * x
    return out


from sq import C_of_t

import sys
m = float(sys.argv[1]) if len(sys.argv) > 1 else 1024.0
f = float(sys.argv[2]) if len(sys.argv) > 2 else 0.12
A, L = 1.0, 2048
seeds = np.arange(4) + 52000
S = len(seeds); dt = default_dt(m)
U = np.zeros((S, L)); V = np.full((S, L), f)
_run(U, V, seeds, np.full(S, f), m, int(3e5 / dt), 10, A, 1.0, 1.0, dt, 0.0, False)
sites = np.arange(0, L, 8)
nrec, every = 12000, 5            # snapshots every 5 steps = 2 time units
Ur = np.zeros((nrec, S, len(sites))); Fr = np.zeros_like(Ur)
for k in range(nrec):
    _run(U, V, seeds, np.full(S, f), m, every, every, A, 1.0, 1.0, dt, 0.0, False)
    Ur[k] = (U - U.mean(1, keepdims=True))[:, sites]
    Fr[k] = site_force(U, seeds.astype(np.int64), sites, A)
tau = dt * every
v = None
# mean velocity from the record of the full line is not stored; use displacement of the mean
# (recomputed from a short extra run)
u0 = U.mean(); _run(U, V, seeds, np.full(S, f), m, 50000, 10, A, 1.0, 1.0, dt, 0.0, False)
v = (U.mean() - u0) / (50000 * dt)
nl = 1500
F0 = Fr - Fr.mean()
Fk = np.fft.rfft(F0, n=2 * nrec, axis=0)
ac = np.fft.irfft(np.abs(Fk) ** 2, axis=0)[:nl].mean(axis=(1, 2)) / (nrec - np.arange(nl))
B = np.array([((Ur[l:] - Ur[:nrec - l]) ** 2).mean() for l in range(nl)])
t = tau * np.arange(nl)
CG = C_of_t(t, v, np.sqrt(B), A)
Cb = C_of_t(t, v, np.zeros_like(t), A)
np.savez(f'data/sq_diag_m{m:g}_f{f:g}.npz', t=t, C=ac, B=B, CG=CG, Cb=Cb, v=v, Fmean=Fr.mean())
I = lambda c: 2 * np.trapezoid(c, t)
print('v', round(v, 4), '<F>', round(Fr.mean(), 4))
print('S_xi(0) = 2 int C:  measured', round(I(ac), 3), ' Gaussian(meas B)', round(I(CG), 3), ' bare', round(I(Cb), 3))
for l in [0, 5, 10, 20, 40, 80, 160, 320, 640, 1200]:
    print(f'  t={t[l]:7.1f}  C {ac[l]:.4f}  CG {CG[l]:.4f}  Cb {Cb[l]:.4f}  B {B[l]:.3f}')
