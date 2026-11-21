"""Open point: how marginal are freshly stopped regions?

Compares a deeply relaxed pinned line (overdamped relaxation from flat at f)
with a line that has just stopped after sliding (hot line quenched to f until
it stops, then settled with m = 0).  For every site, the local stability
margin x_i is the extra force needed to push site i over its next force
barrier with its neighbours held fixed; also the lowest Hessian eigenvalue.
Then repeats the half-hot / half-pinned edge experiment with the freshly
stopped line as the cold half.  Writes data/marg_*.npz, prints a summary.
Run with a thread cap, e.g. NUMBA_NUM_THREADS=3.
"""
import sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from sim import _g
from protocols import _run, default_dt
from edges import hot_length

L = 4096


def node(seed, i, n):
    return _g(np.int64(seed), np.int64(i), np.int64(n), 1.0, np.int64(0))


def margins_slope(u, seed):
    Lc = len(u)
    slope = np.empty(Lc)
    for i in range(Lc):
        n0 = int(np.floor(u[i]))
        slope[i] = node(seed, i, n0 + 1) - node(seed, i, n0)
    return None, slope


def margins(u, seed, f):
    """x_i = -max of g(u') = kappa(u_{i-1}+u_{i+1}-2u') + F_i(u') + f over the first barrier ahead."""
    Lc = len(u)
    x = np.empty(Lc); slope = np.empty(Lc)
    for i in range(Lc):
        s = u[i - 1] + u[(i + 1) % Lc]
        n0 = int(np.floor(u[i]))
        vals = [node(seed, i, n0 + k) for k in range(0, 8)]
        slope[i] = vals[1] - vals[0]                      # F_i' on the current cell
        g = lambda uu, k: (s - 2 * uu) + vals[k] + f       # at node n0+k
        best = None
        prev = (s - 2 * u[i]) + vals[0] * (1 - (u[i] - n0)) + vals[1] * (u[i] - n0) + f
        rising = False
        for k in range(1, 8):
            gk = g(n0 + k, k)
            if not rising and gk > prev:                   # passed the local minimum of g
                rising = True
            elif rising and gk < prev:                     # g turned down again: barrier top at previous node
                best = prev
                break
            prev = gk
        x[i] = -(best if best is not None else prev)
    return x, slope


def low_eigs(slope, k=30):
    Lc = len(slope)
    d = 2.0 - slope
    H = sp.diags([d, -np.ones(Lc - 1), -np.ones(Lc - 1)], [0, 1, -1], format='lil')
    H[0, Lc - 1] = -1; H[Lc - 1, 0] = -1
    return np.sort(eigsh(H.tocsr(), k=k, which='SA', return_eigenvectors=False))


def trigger(U, seed, f, df=2e-5, fmax=None, jump=0.05):
    """Extra force (m = 0, quasi-static steps df) before the first avalanche: some site jumps by > jump."""
    U = U.copy(); V = np.zeros_like(U)
    g = f
    fmax = f + 0.01 if fmax is None else fmax
    while g < fmax:
        g += df
        before = U[0].copy()
        for _ in range(50):
            _, _, _, st, _ = _run(U, V, [seed], [g], 0.0, 20000, 200, 1.0, 1.0, 1.0, 0.05, 1e-7, True)
            if st[0] or np.max(U[0] - before) > jump:
                break
        if np.max(U[0] - before) > jump:
            return g - f
    return np.nan


def relax_m0(U, V, seed, f):
    for _ in range(500):
        _, _, _, st, _ = _run(U, V, [seed], [f], 0.0, 20000, 200, 1.0, 1.0, 1.0, 0.05, 1e-7, True)
        if st[0]:
            return True
    return False


def relaxed(seed, f):
    U = np.zeros((1, L)); V = np.zeros((1, L))
    relax_m0(U, V, seed, f)
    return U, V


def fresh(seed, m, f, fhot):
    dt = default_dt(m)
    U = np.zeros((1, L)); V = np.full((1, L), fhot)
    n = int(max(10 * m, 4000) / dt)
    _run(U, V, [seed], [fhot], m, n, n, 1.0, 1.0, 1.0, dt, 0.0, False)
    ns = int(max(m, 200) / dt)
    for _ in range(4000):
        u0 = U[0].mean()
        _run(U, V, [seed], [f], m, ns, ns, 1.0, 1.0, 1.0, dt, 0.0, False)
        if U[0].mean() - u0 < 0.01:
            break
    V[:] = 0.0
    relax_m0(U, V, seed, f)
    return U, V


def summary(tag, U, seed, f):
    _, slope = margins_slope(U[0], seed)
    ev = low_eigs(slope)
    dft = trigger(U, seed, f)
    print(f'  {tag:8s} lowest eigenvalues {ev[0]:.2e} {ev[4]:.2e} {ev[-1]:.2e};  trigger extra force {dft:.4f}', flush=True)
    return ev, dft


if __name__ == '__main__':
    cases = [(64.0, 0.185, 0.21), (1024.0, 0.088, 0.11)]
    for m, f, fhot in cases:
        t0 = time.time()
        print(f'm={m:g} f={f}', flush=True)
        X = {}
        for tag, maker in [('relaxed', lambda s: relaxed(s, f)), ('fresh', lambda s: fresh(s, m, f, fhot))]:
            evs, dfs = [], []
            for seed in range(4):
                U, V = maker(500 + seed)
                np.save(f'data/marg_state_{tag}_m{m:g}_s{seed}.npy', U[0])
                ev, dft = summary(tag, U, 500 + seed, f)
                evs.append(ev); dfs.append(dft)
            X[tag] = (np.array(evs), np.array(dfs))
            print(f'  {tag}: mean trigger {np.nanmean(X[tag][1]):.4f}', flush=True)
        np.savez(f'data/marg_m{m:g}_f{f:g}.npz', ev_relaxed=X['relaxed'][0], ev_fresh=X['fresh'][0],
                 trig_relaxed=X['relaxed'][1], trig_fresh=X['fresh'][1])
        print(f'  ({time.time() - t0:.0f}s)', flush=True)
