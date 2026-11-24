"""Short runs that fix the analysis thresholds from the data instead of by hand.

Modes (one at a time; each takes minutes with NUMBA_NUM_THREADS=6):

  vhist       single-site speed distribution |du/dt| in partial and uniform sliding
              (m = 64 and 1024), to place the active-site threshold in the gap between
              pinned and moving sites                          -> data/check_vhist.npz
  subs F      substrates (pinned lines) quenched at force F instead of 0.1955
                                                               -> data/check_sub_F{F}.npz
  spread TAG  spreading at f = 0.1981 on the L = 16384 substrates, recording N, R2 and
              survival for several activity thresholds at once; TAG selects the kick or
              the substrate:  base (64 sites, v 0.3), narrow (16 = l_m sites, v 0.1),
              wide (256 sites, v 0.3), sub0.1908 (substrates quenched at 0.95 f_min)
                                                               -> data/check_spread_{TAG}.npz
  aval C M..  avalanches under the standard well with detection threshold vbar > C wdot
                                                               -> data/check_aval_c{C}_m{M}.npz
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import time
import numpy as np
from protocols import _run, default_dt

VACTS = np.array([0.0025, 0.005, 0.01, 0.02, 0.05, 0.1])   # activity thresholds tested
from config import V_REST as V_DEAD                         # a run ends when max|v| < V_REST


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def vhist():
    edges = np.linspace(-6, 1, 141)                        # log10 |v|
    out = {}
    for m, f in [(64.0, 0.198), (64.0, 0.25), (1024.0, 0.093), (1024.0, 0.12)]:
        L, S, dt = 4096, 4, default_dt(m)
        seeds = np.arange(S) + 71000
        U = np.zeros((S, L)); V = np.full((S, L), f)
        n = int(round(500 * m / dt))
        _run(U, V, seeds, np.full(S, f), m, n, n // 10, 1.0, 1.0, 1.0, dt, 0.0, False)
        h = np.zeros(len(edges) - 1); rho = []
        ns = int(round(m / dt))
        for _ in range(500):
            _run(U, V, seeds, np.full(S, f), m, ns, ns, 1.0, 1.0, 1.0, dt, 0.0, False)
            a = np.abs(V).ravel()
            h += np.histogram(np.log10(np.maximum(a, 1e-12)), edges)[0]
            rho.append([(np.abs(V) > v).mean() for v in VACTS])
        out[f'm{m:g}_f{f:g}'] = h
        out[f'rho_m{m:g}_f{f:g}'] = np.array(rho)
        log('vhist', m, f, 'rho(0.02)', np.array(rho)[:, 3].mean())
    np.savez('data/check_vhist.npz', edges=edges, vacts=VACTS, **out)


M, L_SP, F_TEST = 64.0, 16384, 0.1981


def substrates(f_dead, nsub=12, f_prep=0.25):
    dt = default_dt(M)
    seeds = np.arange(nsub) + 33000                        # same disorder as the L = 16384 set
    U = np.zeros((nsub, L_SP)); V = np.full((nsub, L_SP), f_prep)
    n0 = int(round(4000.0 / dt))
    _run(U, V, seeds, np.full(nsub, f_prep), M, n0, n0 // 10, 1.0, 1.0, 1.0, dt, 0.0, False)
    ns = int(round(50 * M / dt)); t = 0.0
    dead = np.zeros(nsub, bool)
    while not dead.all() and t < 12000 * M:
        idx = np.where(~dead)[0]
        Ua, Va = U[idx].copy(), V[idx].copy()
        _run(Ua, Va, seeds[idx], np.full(len(idx), f_dead), M, ns, ns // 10, 1.0, 1.0, 1.0, dt,
             V_DEAD, True)
        U[idx], V[idx] = Ua, Va; t += ns * dt
        dead[idx[np.abs(Va).max(1) < V_DEAD]] = True
    log('subs', f_dead, 'dead', dead.sum(), 'of', nsub, f't/m {t / M:.0f}')
    np.savez(f'data/check_sub_F{f_dead:g}.npz', U=U[dead], V=V[dead], seeds=seeds[dead])


def spread(tag, T_over_m=2000, nk=8):
    w, v_kick, sub = {'base': (64, 0.3, 'data/spread_sub.npz'),
                      'narrow': (16, 0.1, 'data/spread_sub.npz'),
                      'wide': (256, 0.3, 'data/spread_sub.npz'),
                      'sub0.1908': (64, 0.3, 'data/check_sub_F0.1908.npz')}[tag]
    dt = default_dt(M)
    d = np.load(sub)
    U = np.repeat(d['U'], nk, 0); V = np.repeat(d['V'], nk, 0); seeds = np.repeat(d['seeds'], nk)
    S, L = U.shape
    fs = np.full(S, F_TEST)
    nrel = int(round(20 * M / dt))
    _run(U, V, seeds, fs, M, nrel, nrel // 10, 1.0, 1.0, 1.0, dt, 0.0, False)
    pre_act = (np.abs(V) > 0.02).sum(1)
    x0 = ((np.arange(S) % nk) * (L // nk) + 17 * (np.arange(S) // nk)) % L
    for k in range(S):
        V[k, (x0[k] + np.arange(w)) % L] = v_kick
    x = np.arange(L); ns = int(round(M / dt)); nt = int(T_over_m); nv = len(VACTS)
    N = np.zeros((nv, S, nt), np.int32); R2 = np.zeros((nv, S, nt))
    alive = np.ones(S, bool)
    for it in range(nt):
        idx = np.where(alive)[0]
        if len(idx) == 0:
            break
        Ua, Va = U[idx].copy(), V[idx].copy()
        _run(Ua, Va, seeds[idx], fs[idx], M, ns, ns, 1.0, 1.0, 1.0, dt, 0.0, False)
        U[idx], V[idx] = Ua, Va
        for j, k in enumerate(idx):
            av = np.abs(Va[j])
            for iv, vt in enumerate(VACTS):
                a = av > vt; n = a.sum(); N[iv, k, it] = n
                if n:
                    dx = (x[a] - (x0[k] + w // 2) + L // 2) % L - L // 2
                    R2[iv, k, it] = (dx.astype(float) ** 2).mean()
            if av.max() < V_DEAD:
                alive[k] = False
        if (it + 1) % 500 == 0:
            log('spread', tag, f't/m {it + 1} alive {alive.sum()}/{S}')
    np.savez(f'data/check_spread_{tag}.npz', N=N, R2=R2, vacts=VACTS, seeds=seeds, pre_act=pre_act,
             f=F_TEST, m=M, L=L, w=w, v_kick=v_kick, sub=sub)


if __name__ == '__main__':
    what = sys.argv[1]
    t0 = time.time()
    if what == 'vhist':
        vhist()
    elif what == 'subs':
        substrates(float(sys.argv[2]))
    elif what == 'spread':
        spread(sys.argv[2])
    elif what == 'aval':
        from aval import avalanches
        c = float(sys.argv[2])
        for m in [float(x) for x in sys.argv[3:]]:
            r = avalanches(m, W_total=70.0, W_skip=2.0, c_th=c)
            np.savez(f'data/check_aval_c{c:g}_m{m:g}.npz', **r)
            log('aval c', c, 'm', m, 'events', len(r['s']))
    log('done', what, f'{time.time() - t0:.0f}s')
