"""Exact sample depinning threshold for long chains (no-passing relaxation).

Free-end chain, kappa = 1, triangle-wave pinning (amplitude A, period P,
phase phi_i):  F_i(u) = (Delta_N^2 u)_i + f - V_i'(u_i).

No-passing (Middleton): forward-only relaxation from any start u^0 either
stops in a stable configuration, or reaches u >= u^0 + P at every site; by
the shift symmetry u -> u + P it then moves forever.  So "pinned at f" is
decided exactly, and f_c(sample) follows by bisection.  Each local move
solves F_i = 0 exactly on the piecewise-linear force (first zero ahead).
This is the polynomial-time feasibility oracle of the exact formula
f_c = A - min over stable (lift, branches) of max_i L_i.

    python3 relax.py   -> relax.json
"""
import json
import numpy as np
from numba import njit
from multiprocessing import Pool

P = 4.0


@njit(cache=True)
def _tri(u, ph, A, P):
    half = 0.5 * P
    x = (u - ph) / half
    n = np.floor(x)
    t = x - n
    if (int(n) & 1) == 0:
        return -A + 2.0 * A * t, 1.0, ph + half * (n + 1.0)
    return A - 2.0 * A * t, -1.0, ph + half * (n + 1.0)


@njit(cache=True)
def _local_zero(u, ph, A, P, k, a, b, f):
    """First zero ahead of F(v) = a - b v + f - g(v), b = 2 (1 at ends),
    given F(u) > 0.  Walks piece by piece."""
    v = u
    for _ in range(1000000):
        g, sg, vnext = _tri(v, ph, A, P)
        F = a - b * v + f - g
        slope = -b - sg * k
        if slope < 0.0:
            vz = v - F / slope
            if vz <= vnext:
                return vz
        v = vnext + 1e-15 * (1.0 + abs(vnext))
    return v


@njit(cache=True)
def pinned(phase, A, P, f, max_sweeps):
    """True if forward relaxation from u = 0 stops before advancing P everywhere."""
    K = phase.shape[0]
    k = 4.0 * A / P
    u = np.zeros(K)
    active = np.ones(K, np.bool_)
    tol = 1e-11
    for sweep in range(max_sweeps):
        moved = False
        for i in range(K):
            if not active[i]:
                continue
            active[i] = False
            if i == 0:
                a = u[1]; b = 1.0
            elif i == K - 1:
                a = u[K - 2]; b = 1.0
            else:
                a = u[i - 1] + u[i + 1]; b = 2.0
            g, sg, vn = _tri(u[i], phase[i], A, P)
            F = a - b * u[i] + f - g
            if F > tol:
                u[i] = _local_zero(u[i], phase[i], A, P, k, a, b, f)
                moved = True
                if i > 0:
                    active[i - 1] = True
                if i < K - 1:
                    active[i + 1] = True
        if not moved:
            return True, sweep
        if u.min() >= P:
            return False, sweep
    return False, max_sweeps


def fc_relax(phase, A, tol=1e-7, max_sweeps=10_000_000):
    lo, hi = 0.0, A
    while hi - lo > tol * A:
        mid = 0.5 * (lo + hi)
        ok, _ = pinned(phase, A, P, mid, max_sweeps)
        if ok:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def work(args):
    seed, K, A = args
    ph = np.random.default_rng(seed).uniform(0, P, K)
    return fc_relax(ph, A)


if __name__ == "__main__":
    from annealed import gle, f_ann
    out = {}
    with Pool(9) as pool:
        for A in (0.3, 0.6, 1.0, 2.0, 3.0, 6.0):
            rows = []
            for K in (16, 32, 64, 128, 256, 512, 1024):
                n = 600 if K <= 256 else 300
                v = np.array(pool.map(work, [(71_000_000 + 1000 * K + s, K, A) for s in range(n)]))
                rows.append(dict(K=K, mean=float(v.mean()), err=float(v.std() / np.sqrt(n)), sd=float(v.std())))
                print(f"A={A} K={K}: <f_c>={v.mean():.5f}+-{v.std()/np.sqrt(n):.5f} sd*sqrt(K)={v.std()*np.sqrt(K):.3f}", flush=True)
            out[str(A)] = dict(A=A, k=4 * A / P, rows=rows, f_ann=f_ann(A, P, gle(1.0, 4 * A / P)))
            json.dump(out, open("relax.json", "w"), indent=1)
