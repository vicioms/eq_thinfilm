"""Driving protocols for the inertial line (sim.run does the integration).

flow:       constant force, start sliding (u = 0, v = f/eta); after a transient,
            mean velocity, velocity variance and active fraction.
ramp_down:  start sliding at f0 and lower f by df every dwell; f_down is the force at
            which the line comes to rest.
ramp_up:    relax at f0 with m = 0 (a pinned state), then raise f by df each time the
            line is at rest; f_up is the force at which it travels more than D.
well:       parabolic well k0 (w - u_i), w advanced at constant rate wdot; returns the
            force k0 (w - ubar) and the velocity as time series.
All functions run S samples in parallel, one numba thread per sample.
"""
import numpy as np
from sim import run


def default_dt(m):
    if m <= 0:
        return 0.05
    return float(min(0.4, max(0.05, 0.05 * np.sqrt(m))))


def _run(U, V, seeds, f, m, nsteps, nrec, A, eta, kappa, dt, vtol, stop,
         k0=0.0, w=None, wdot=0.0):
    S = U.shape[0]
    w = np.zeros(S) if w is None else w
    return run(U, V, np.asarray(seeds, np.int64), np.asarray(f, float).copy(), m, eta,
               kappa, A, 0, k0, w, wdot, dt, nsteps, nrec, vtol, stop)


def flow(m, L, seeds, fs, T_relax, T_meas, A=1.0, eta=1.0, kappa=1.0, dt=None,
         vact=0.02, nsnap=200):
    """One sample per (seed, f) pair: arrays of shape (len(fs), len(seeds))."""
    dt = dt or default_dt(m)
    F, Sd = np.meshgrid(fs, seeds, indexing='ij')
    F, Sd = F.ravel(), Sd.ravel()
    S = len(F)
    U = np.zeros((S, L))
    V = np.repeat(F[:, None] / eta, L, axis=1)
    n0 = int(round(T_relax / dt))
    _run(U, V, Sd, F, m, n0, max(1, n0 // 10), A, eta, kappa, dt, 0.0, False)
    ns = max(1, int(round(T_meas / dt / nsnap)))
    vbar = np.zeros(S); dv2 = np.zeros(S); act = np.zeros(S)
    u0 = U.mean(1)
    for k in range(nsnap):
        vb, v2, ub, _, _ = _run(U, V, Sd, F, m, ns, ns, A, eta, kappa, dt, 0.0, False)
        dv2 += v2[:, -1] - vb[:, -1] ** 2
        act += (np.abs(V) > vact).mean(1)
    vbar = (U.mean(1) - u0) / (ns * nsnap * dt)
    sh = (len(fs), len(seeds))
    return dict(f=F.reshape(sh), v=vbar.reshape(sh), dv2=(dv2 / nsnap).reshape(sh),
                act=(act / nsnap).reshape(sh))


def ramp_down(m, L, seeds, f0, df, dwell, A=1.0, eta=1.0, kappa=1.0, fmin=0.0, dt=None,
              dstop=0.02):
    dt = dt or default_dt(m)
    S = len(seeds)
    U = np.zeros((S, L))
    V = np.full((S, L), f0 / eta)
    f = np.full(S, float(f0))
    seeds = np.asarray(seeds, np.int64)
    nsteps = int(round(dwell / dt))
    nrec = max(1, nsteps // 200)
    active = np.ones(S, bool)
    fdown = np.full(S, np.nan)
    rec_f, rec_v = [], []
    while active.any():
        idx = np.where(active)[0]
        Ua, Va = U[idx].copy(), V[idx].copy()
        vb, v2, ub, st, _ = _run(Ua, Va, seeds[idx], f[idx], m, nsteps, nrec, A, eta, kappa,
                                 dt, 1e-5, True)
        U[idx], V[idx] = Ua, Va
        nr = vb.shape[1]
        vrow = np.full(S, np.nan)
        for j, s in enumerate(idx):
            disp = ub[j, -1] - ub[j, nr // 2] if not st[j] else 0.0
            if st[j] or disp < dstop:
                fdown[s] = f[s]
                active[s] = False
                vrow[s] = 0.0
            else:
                vrow[s] = disp / ((nr - 1 - nr // 2) * nrec * dt)
        rec_f.append(f.copy()); rec_v.append(vrow)
        f[active] -= df
        low = active & (f < fmin)
        active &= ~low
    return dict(fdown=fdown, f=np.array(rec_f), v=np.array(rec_v))


def ramp_up(m, L, seeds, f0, df, D=30.0, A=1.0, eta=1.0, kappa=1.0, fmax=1.0, dt=None,
            vtol=1e-4):
    S = len(seeds)
    seeds = np.asarray(seeds, np.int64)
    U = np.zeros((S, L))
    V = np.zeros((S, L))
    f = np.full(S, float(f0))
    # pinned initial state: overdamped relaxation at f0
    for _ in range(1000):
        _, _, _, st, _ = _run(U, V, seeds, f, 0.0, 20000, 200, A, eta, kappa, 0.05, 1e-6, True)
        if st.all():
            break
    V[:] = 0.0
    dt = dt or default_dt(m)
    chunk = max(2000, int(round(max(200.0, 4 * m) / dt)))
    nrec = max(1, chunk // 50)
    ustart = U.mean(1)
    active = np.ones(S, bool)
    fup = np.full(S, np.nan)
    while active.any():
        idx = np.where(active)[0]
        Ua, Va = U[idx].copy(), V[idx].copy()
        _, _, ub, st, _ = _run(Ua, Va, seeds[idx], f[idx], m, chunk, nrec, A, eta, kappa,
                               dt, vtol, True)
        U[idx], V[idx] = Ua, Va
        for j, s in enumerate(idx):
            if U[s].mean() - ustart[s] > D:
                fup[s] = f[s]
                active[s] = False
            elif st[j]:
                f[s] += df
                ustart[s] = U[s].mean()
                V[s] = 0.0
                if f[s] > fmax:
                    active[s] = False
    return dict(fup=fup)


def well(m, L, seeds, k0, wdot, T, w0=0.0, A=1.0, eta=1.0, kappa=1.0, dt=None, nrec=None):
    dt = dt or default_dt(m)
    S = len(seeds)
    U = np.zeros((S, L)); V = np.zeros((S, L))
    w = np.full(S, float(w0))
    nsteps = int(round(T / dt))
    nrec = nrec or max(1, nsteps // 20000)
    vb, v2, ub, _, _ = _run(U, V, seeds, np.zeros(S), m, nsteps, nrec, A, eta, kappa, dt, 0.0,
                            False, k0=k0, w=w, wdot=wdot)
    t = (np.arange(vb.shape[1]) + 1) * nrec * dt
    wt = w0 + wdot * t
    return dict(t=t, force=k0 * (wt[None, :] - ub), v=vb)


def lifetimes(m, L, seeds, f, f_prep, T_prep, T_max, A=1.0, eta=1.0, kappa=1.0, dt=None,
              chunk=None, vact=0.02, dstop=0.01, T_skip=None, trace=False):
    """Quench a sliding line from f_prep to f and record when it stops.

    Returns t_stop[S] (inf if still sliding at T_max), and the time-averaged mean
    velocity and active fraction over the surviving part after T_skip.
    """
    dt = dt or default_dt(m)
    S = len(seeds)
    seeds = np.asarray(seeds, np.int64)
    U = np.zeros((S, L)); V = np.full((S, L), f_prep / eta)
    n0 = int(round(T_prep / dt))
    _run(U, V, seeds, np.full(S, f_prep), m, n0, max(1, n0 // 10), A, eta, kappa, dt, 0.0, False)
    chunk = chunk or max(200.0, float(m))
    ns = int(round(chunk / dt))
    T_skip = chunk * 5 if T_skip is None else T_skip
    tstop = np.full(S, np.inf)
    alive = np.ones(S, bool)
    vsum = np.zeros(S); asum = np.zeros(S); nacc = np.zeros(S)
    tr_t, tr_a = [], []   # optional activity trace: per chunk, active fraction of every line (0 if stopped)
    t = 0.0
    while alive.any() and t < T_max:
        idx = np.where(alive)[0]
        Ua, Va = U[idx].copy(), V[idx].copy()
        _, _, ub, _, _ = _run(Ua, Va, seeds[idx], np.full(len(idx), f), m, ns, ns // 4, A, eta,
                              kappa, dt, 0.0, False)
        U[idx], V[idx] = Ua, Va
        t += ns * dt
        if trace:
            tr_t.append(t); tr_a.append((np.abs(V) > vact).mean(axis=1) * alive)
        for j, s in enumerate(idx):
            disp = ub[j, -1] - ub[j, 0]
            if disp < dstop:
                tstop[s] = t
                alive[s] = False
            elif t > T_skip:
                vsum[s] += disp / (ub.shape[1] - 1) / (ns // 4) / dt
                asum[s] += (np.abs(V[s]) > vact).mean()
                nacc[s] += 1
    with np.errstate(invalid='ignore', divide='ignore'):
        out = dict(tstop=tstop, v=vsum / nacc, act=asum / nacc, nacc=nacc)
    if trace:
        out.update(trace_t=np.array(tr_t), trace_act=np.array(tr_a))
    return out
