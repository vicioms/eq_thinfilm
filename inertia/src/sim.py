"""Driven elastic line with inertia in a random-field landscape.

    m u_i'' + eta u_i' = kappa (u_{i+1} + u_{i-1} - 2 u_i) + F_i(u_i) + f + k0 (w - u_i)

i = 0..L-1 periodic.  F_i(u) is the linear interpolation of i.i.d. values
g_i(n), uniform on [-A, A], at integer n (r_f = 1).  The values come from a
counter-based hash of (seed, i, n mod Pu), so the landscape is infinite (or
Pu-periodic in u) and costs no memory.  Statistics of the interpolated force:
<F^2> = 2 A^2 / 9 at a random point, correlator Fourier transform
Delta(k) = (A^2/3) sinc^4(k/2).  A < 0 selects random-bond disorder of amplitude
|A| (node values are differences of i.i.d. values): Delta(k) = (A^2/3)(2 - 2 cos k)
sinc^4(k/2), <F^2> = A^2/3.

Time step (first order, uniformly valid in m):
    v <- v e^{-eta dt/m} + (F_tot/eta)(1 - e^{-eta dt/m}),   u <- u + dt v
which is the exact damped update at frozen force; at m = 0 it is Euler for the
overdamped line.  Many samples run in parallel (one numba thread per sample).
"""
import numpy as np
from numba import njit, prange

_M1 = np.uint64(0x9E3779B97F4A7C15)
_M2 = np.uint64(0xBF58476D1CE4E5B9)
_M3 = np.uint64(0x94D049BB133111EB)


@njit(inline='always')
def _hash(seed, i, n):
    x = np.uint64(seed) * _M1 ^ (np.uint64(i) + np.uint64(0x632BE59BD9B4E019)) * _M2 \
        ^ np.uint64(n) * _M3
    x ^= x >> np.uint64(30)
    x *= _M2
    x ^= x >> np.uint64(27)
    x *= _M3
    x ^= x >> np.uint64(31)
    return (x >> np.uint64(11)) * (1.0 / 9007199254740992.0)


@njit(inline='always')
def _g0(seed, i, n, A, Pu):
    if Pu > 0:
        n = n % Pu
        if n < 0:
            n += Pu
    return A * (2.0 * _hash(seed, i, n + 1073741824) - 1.0)


@njit(inline='always')
def _g(seed, i, n, A, Pu):
    """Node value of the force at integer u = n.

    A >= 0: random field, i.i.d. uniform on [-A, A].
    A < 0:  random bond of amplitude |A|: node values g(n+1) - g(n), so the
            potential is bounded (force correlator integrates to zero),
            <F^2> = A^2/3 at a random point.
    """
    if A >= 0:
        return _g0(seed, i, n, A, Pu)
    return _g0(seed, i, n + 1, -A, Pu) - _g0(seed, i, n, -A, Pu)


@njit(parallel=True, cache=True)
def run(U, V, seeds, fs, m, eta, kappa, A, Pu, k0, w, wdot, dt, nsteps, nrec,
        vtol, stop_on_rest):
    """Advance S samples (rows of U, V) for nsteps.

    fs[s]: constant force of sample s; w[s]: well position, moved by wdot*dt
    per step (k0 = 0: no well).  Records every nrec steps the mean velocity
    and the mean of v^2.  With stop_on_rest, a sample stops as soon as
    max|v| < vtol (checked at record times) and its row of `stopped` is set.

    Returns vbar[S, nr], v2[S, nr], ubar[S, nr], stopped[S], t_done[S], nr = nsteps//nrec.
    """
    S, L = U.shape
    nr = nsteps // nrec
    vbar = np.zeros((S, nr))
    v2 = np.zeros((S, nr))
    ubar = np.zeros((S, nr))
    stopped = np.zeros(S, np.bool_)
    tdone = np.zeros(S, np.int64)
    if m > 0:
        e = np.exp(-eta * dt / m)
    else:
        e = 0.0
    for s in prange(S):
        u = U[s]
        v = V[s]
        seed = seeds[s]
        f = fs[s]
        ww = w[s]
        cn = np.empty(L, np.int64)
        g0 = np.empty(L)
        g1 = np.empty(L)
        for i in range(L):
            n = int(np.floor(u[i]))
            cn[i] = n
            g0[i] = _g(seed, i, n, A, Pu)
            g1[i] = _g(seed, i, n + 1, A, Pu)
        Ft = np.empty(L)
        r = 0
        for step in range(nsteps):
            for i in range(L):
                ui = u[i]
                n = int(np.floor(ui))
                if n != cn[i]:
                    if n == cn[i] + 1:
                        g0[i] = g1[i]
                        g1[i] = _g(seed, i, n + 1, A, Pu)
                    elif n == cn[i] - 1:
                        g1[i] = g0[i]
                        g0[i] = _g(seed, i, n, A, Pu)
                    else:
                        g0[i] = _g(seed, i, n, A, Pu)
                        g1[i] = _g(seed, i, n + 1, A, Pu)
                    cn[i] = n
                t = ui - n
                im = i - 1 if i > 0 else L - 1
                ip = i + 1 if i < L - 1 else 0
                Ft[i] = kappa * (u[ip] + u[im] - 2.0 * ui) + g0[i] * (1.0 - t) \
                    + g1[i] * t + f + k0 * (ww - ui)
            sv = 0.0
            su = 0.0
            sv2 = 0.0
            vmax = 0.0
            for i in range(L):
                vi = v[i] * e + (Ft[i] / eta) * (1.0 - e)
                v[i] = vi
                u[i] += dt * vi
                su += u[i]
                sv += vi
                sv2 += vi * vi
                av = abs(vi)
                if av > vmax:
                    vmax = av
            ww += wdot * dt
            if (step + 1) % nrec == 0:
                vbar[s, r] = sv / L
                v2[s, r] = sv2 / L
                ubar[s, r] = su / L
                r += 1
                if stop_on_rest and vmax < vtol:
                    stopped[s] = True
                    tdone[s] = step + 1
                    break
        if not stopped[s]:
            tdone[s] = nsteps
        w[s] = ww
    return vbar, v2, ubar, stopped, tdone


@njit(parallel=True, cache=True)
def pinning_force(U, seeds, A, Pu):
    """Mean disorder force sum_i F_i(u_i)/L for each sample."""
    S, L = U.shape
    out = np.zeros(S)
    for s in prange(S):
        acc = 0.0
        for i in range(L):
            n = int(np.floor(U[s, i]))
            t = U[s, i] - n
            acc += _g(seeds[s], i, n, A, Pu) * (1 - t) + _g(seeds[s], i, n + 1, A, Pu) * t
        out[s] = acc / L
    return out
