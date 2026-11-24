"""Spreading from a single hot seed at the end of sliding (m = 64).

Substrate: the natural absorbing state.  A line slides at f_prep = 0.25, is quenched to
f_dead = 0.1955 (below f_inf) and left until it is at rest (max|v| < 1e-3).  Each dead
substrate is replicated NK times, raised to the test force f, relaxed, and kicked: sites
x0 .. x0 + w get velocity v_kick, at a different x0 in each replica.  Recorded every m:
number of active sites N (|v| > 0.02), their mean squared distance R2 from the seed
(minimal image), and the half-extent of the active region.

At a directed-percolation critical point (1+1 d):
    P_surv ~ t^-delta, delta = 0.159;  <N> ~ t^theta, theta = 0.313;  R2 ~ t^(2/z), z = 1.581.

Usage: [SPREAD_L=32768 SPREAD_NSUB=32 SPREAD_NK=16] python3 spread.py sub                 (prepare substrates -> data/spread_sub.npz)
       python3 spread.py run f [f ...] [--T T_over_m] [--w w] [--v v_kick]
Writes data/spread_f{f}.npz (data/spread_f{f}_L{L}.npz for L != 16384).
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import os
import sys
import time
import numpy as np
from protocols import _run, default_dt
from config import V_ACT, V_REST

m = 64.0
# defaults: the first production run; SPREAD_L / SPREAD_NSUB / SPREAD_NK select a larger set
L = int(os.environ.get('SPREAD_L', 16384))
NSUB = int(os.environ.get('SPREAD_NSUB', 12)); NK = int(os.environ.get('SPREAD_NK', 8))
TAG = '' if L == 16384 else f'_L{L}'
SUB_FILE = f'data/spread_sub{TAG}.npz'
SUB_SEEDS = np.arange(NSUB) + (33000 if L == 16384 else 34000)
f_prep, f_dead, vact = 0.25, 0.1955, V_ACT
dt = default_dt(m)


def substrates():
    S = NSUB
    U = np.zeros((S, L)); V = np.full((S, L), f_prep)
    n0 = int(round(4000.0 / dt))
    _run(U, V, SUB_SEEDS, np.full(S, f_prep), m, n0, n0 // 10, 1.0, 1.0, 1.0, dt, 0.0, False)
    ns = int(round(50 * m / dt))
    dead = np.zeros(S, bool); t = 0.0; t_dead = np.full(S, np.inf)
    while not dead.all() and t < 12000 * m:
        idx = np.where(~dead)[0]
        Ua, Va = U[idx].copy(), V[idx].copy()
        _run(Ua, Va, SUB_SEEDS[idx], np.full(len(idx), f_dead), m, ns, ns // 10, 1.0, 1.0, 1.0, dt,
             V_REST, True)
        U[idx], V[idx] = Ua, Va
        t += ns * dt
        for j, s in enumerate(idx):
            if np.abs(Va[j]).max() < V_REST:
                dead[s] = True; t_dead[s] = t
        print(time.strftime('%H:%M:%S'), f'sub t/m {t / m:.0f} dead {dead.sum()}/{S}', flush=True)
    np.savez(SUB_FILE, U=U[dead], V=V[dead], seeds=SUB_SEEDS[dead], t_dead=t_dead[dead])


def spread(f, T_over_m=3000.0, w=64, v_kick=0.3):
    d = np.load(SUB_FILE)
    U0, V0, sd = d['U'], d['V'], d['seeds']
    U = np.repeat(U0, NK, 0); V = np.repeat(V0, NK, 0); seeds = np.repeat(sd, NK)
    S = len(seeds)
    fs = np.full(S, f)
    # relax at the test force without a seed: raising f from f_dead may itself trigger events
    nrel = int(round(20 * m / dt))
    _run(U, V, seeds, fs, m, nrel, nrel // 10, 1.0, 1.0, 1.0, dt, 0.0, False)
    pre_act = (np.abs(V) > vact).sum(1)
    x0 = ((np.arange(S) % NK) * (L // NK) + 17 * (np.arange(S) // NK)) % L
    for k in range(S):
        V[k, (x0[k] + np.arange(w)) % L] = v_kick
    x = np.arange(L)
    ns = int(round(m / dt)); nt = int(T_over_m)
    N = np.zeros((S, nt), np.int32); R2 = np.zeros((S, nt)); ext = np.zeros((S, nt))
    vmax = np.zeros((S, nt))
    alive = np.ones(S, bool)
    t0 = time.time()
    for it in range(nt):
        idx = np.where(alive)[0]
        if len(idx) == 0:
            break
        Ua, Va = U[idx].copy(), V[idx].copy()
        _run(Ua, Va, seeds[idx], fs[idx], m, ns, ns, 1.0, 1.0, 1.0, dt, 0.0, False)
        U[idx], V[idx] = Ua, Va
        for j, k in enumerate(idx):
            a = np.abs(Va[j]) > vact
            vmax[k, it] = np.abs(Va[j]).max()
            n = a.sum(); N[k, it] = n
            if n:
                dx = (x[a] - (x0[k] + w // 2) + L // 2) % L - L // 2
                R2[k, it] = (dx.astype(float) ** 2).mean(); ext[k, it] = np.abs(dx).max()
            if vmax[k, it] < V_REST:
                alive[k] = False
        if (it + 1) % 250 == 0:
            print(time.strftime('%H:%M:%S'), f'f {f} t/m {it + 1} alive {alive.sum()}/{S} '
                  f'<N> {N[:, it].mean():.1f} max ext {ext[:, it].max():.0f}', f'{time.time() - t0:.0f}s',
                  flush=True)
    np.savez(f'data/spread_f{f:g}{TAG}.npz', N=N, R2=R2, ext=ext, vmax=vmax, x0=x0, pre_act=pre_act,
             seeds=seeds, f=f, m=m, L=L, w=w, v_kick=v_kick)
    return N, ext, pre_act


if __name__ == '__main__':
    if sys.argv[1] == 'sub':
        substrates()
    else:
        args = sys.argv[2:]
        opt = dict(T_over_m=3000.0, w=64, v_kick=0.3)
        fl = []
        i = 0
        while i < len(args):
            if args[i] == '--T':
                opt['T_over_m'] = float(args[i + 1]); i += 2
            elif args[i] == '--w':
                opt['w'] = int(args[i + 1]); i += 2
            elif args[i] == '--v':
                opt['v_kick'] = float(args[i + 1]); i += 2
            else:
                fl.append(float(args[i])); i += 1
        for f in fl:
            N, ext, pre = spread(f, **opt)
            surv = (N > 0)
            print(time.strftime('%H:%M:%S'), 'done f', f, 'pre-activity rows', int((pre > 0).sum()),
                  'P_surv at t/m 10,100,1000:', [round(float(surv[:, k].mean()), 3) for k in [9, 99, 999]
                                                  if k < N.shape[1]], flush=True)
