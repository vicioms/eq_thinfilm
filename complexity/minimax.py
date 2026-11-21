"""Strong-pinning threshold: minimax curvature of integer lifts.

For V'' = k -> infinity every site sits at a force peak q_i = p_i + n_i P
(p_i = phi_i + P/2 for the triangle wave), up to O(P/k), and force balance
(Delta_N^2 u)_i = V_i'(u_i) - f <= A - f gives

    f_c = A - M*,   M* = min_{n in Z^K} max_i (Delta_N^2 q)_i .

With slopes w_j = q_{j+1} - q_j in d_j + P Z (d_j = p_{j+1} - p_j), the
free-end conditions are w_0 <= M, w_j - w_{j-1} <= M, -w_{K-2} <= M.  The
greedy choice of the largest admissible slope at each step is optimal, so
feasibility of M is a single pass, and M* follows by bisection.  The greedy
increments are uniform on (M - P, M], so M* -> P/2 as K -> infinity, with
O(K^{-1/2}) corrections.

    python3 minimax.py   -> minimax.json
"""
import json
import numpy as np
from numba import njit

P = 4.0


@njit(cache=True)
def feasible(d, M):
    w = d[0] + P * np.floor((M - d[0]) / P)
    for j in range(1, d.shape[0]):
        w = d[j] + P * np.floor((w + M - d[j]) / P)
    return w >= -M


@njit(cache=True)
def mstar(d):
    lo, hi = 0.0, 2.0 * P
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if feasible(d, mid):
            hi = mid
        else:
            lo = mid
    return hi


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    out = {}
    for K in [4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 32, 64, 128, 256, 1024, 4096]:
        ms = []
        for s in range(20000 if K <= 256 else 2000):
            ph = rng.uniform(0, P, K)
            p = ph + P / 2
            d = np.diff(p)
            ms.append(mstar(d))
        ms = np.array(ms)
        out[K] = dict(mean=float(ms.mean()), err=float(ms.std() / np.sqrt(len(ms))), sd=float(ms.std()))
        print(f"K={K:5d}  <M*> = {ms.mean():.4f} +- {ms.std()/np.sqrt(len(ms)):.4f}   sd = {ms.std():.4f}   "
              f"sqrt(K)(<M*>-P/2) = {np.sqrt(K)*(ms.mean()-P/2):.3f}", flush=True)
    json.dump(out, open("minimax.json", "w"), indent=1)
