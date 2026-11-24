"""Open point: speed of the edge between hot sliding and pinned line.

A line of L sites starts half hot (taken from a sliding run) and half pinned
(overdamped relaxation at f), joined without a displacement step.  The length
of the hot region l_h(t) (blocks of 16 sites with |v| > 0.02 on average) gives
the edge speed s = (1/2) dl_h/dt (two edges).  Writes data/edge_m{m}_f{f}.npz.
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
from config import V_ACT
import sys, time
import numpy as np
from protocols import _run, default_dt

L, B = 8192, 16


def hot_length(V):
    a = (np.abs(V).reshape(-1, B).mean(1) > V_ACT)
    return a.sum() * B


def one(m, f, seed, fhot, T_run, nrec=200):
    dt = default_dt(m)
    # hot line
    Uh = np.zeros((1, L)); Vh = np.full((1, L), fhot)
    n = int(max(10 * m, 4000) / dt)
    _run(Uh, Vh, [seed], [fhot], m, n, n, 1.0, 1.0, 1.0, dt, 0.0, False)
    # pinned line at f (overdamped)
    Uc = np.zeros((1, L)); Vc = np.zeros((1, L))
    for _ in range(200):
        _, _, _, st, _ = _run(Uc, Vc, [seed], [f], 0.0, 20000, 200, 1.0, 1.0, 1.0, 0.05, 1e-6, True)
        if st[0]:
            break
    U = Uc.copy(); V = np.zeros((1, L))
    h = slice(0, L // 2)
    shift = Uc[0, h].mean() - Uh[0, h].mean()
    U[0, h] = Uh[0, h] + shift          # same seed: site disorder differs after the shift; fine for a hot line
    V[0, h] = Vh[0, h]
    # smooth the two junctions over 64 sites
    for j0 in [0, L // 2]:
        for k in range(-32, 32):
            i = (j0 + k) % L
            w = 0.5 * (1 + np.tanh(k / 8.0))
            if j0 == L // 2:            # hot region is [0, L/2): hot weight -> 1 inside it
                w = 1 - w
            U[0, i] = w * (Uh[0, i] + shift) + (1 - w) * Uc[0, i]
    ns = int(T_run / dt / nrec)
    t, lh = [], []
    for k in range(nrec):
        _run(U, V, [seed], [f], m, ns, ns, 1.0, 1.0, 1.0, dt, 0.0, False)
        t.append((k + 1) * ns * dt); lh.append(hot_length(V[0]))
    return np.array(t), np.array(lh)


if __name__ == '__main__':
    jobs = {1024.0: ([0.080, 0.086, 0.090, 0.094, 0.098, 0.105, 0.120], 0.11),
            256.0: ([0.120, 0.128, 0.134, 0.140, 0.146, 0.155, 0.170], 0.16),
            64.0: ([0.180, 0.186, 0.190, 0.194, 0.198, 0.204, 0.215], 0.24)}
    ms = [float(x) for x in sys.argv[1:]] or list(jobs)
    for m in ms:
        fs, fhot = jobs[m]
        for f in fs:
            t0 = time.time()
            T, LH = [], []
            for seed in range(4):
                t, lh = one(m, f, 300 + seed, max(fhot, f), T_run=60 * m)
                T.append(t); LH.append(lh)
            LH = np.array(LH)
            np.savez(f'data/edge_m{m:g}_f{f:g}.npz', m=m, f=f, L=L, t=T[0], lh=LH)
            i1 = len(T[0]) // 4
            s = 0.5 * (LH[:, i1] - LH[:, 0]).mean() / (T[0][i1] - T[0][0])
            print(time.strftime('%H:%M:%S'), 'edge', m, f, 'l_h start', LH[:, 0].mean(), 'end', LH[:, -1].mean(),
                  'early s', round(s, 4), 'c', round(1 / np.sqrt(m), 4), f'{time.time() - t0:.0f}s', flush=True)
