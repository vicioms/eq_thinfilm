"""Production runs for the inertial-depinning note.  Results in data/*.npz.

    python3 runs.py flow        # flow curves, L = 4096, A = 1
    python3 runs.py thresholds  # f_down, f_up vs m and L, A = 1 (and A = 0.5 at L = 1024)
    python3 runs.py well        # parabolic-well driving, time series
    python3 runs.py kymo        # space-time velocity maps
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
from config import RAMP_GO
import sys, time
import numpy as np
from protocols import flow, ramp_down, ramp_up, well, default_dt, _run

MS = [0.0, 0.25, 1.0, 4.0, 16.0, 64.0, 256.0, 1024.0, 4096.0]
OUT = 'data/'


def f_est(m, A):
    """Rough stopping force, used only to choose where ramps start."""
    if m == 0:
        return 0.30 * A ** (4 / 3)
    return min(0.30 * A ** (4 / 3), 0.55 * A * m ** -0.25)


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def do_flow():
    L = 4096
    seeds = np.arange(4) + 11
    for m in MS:
        fe = f_est(m, 1.0)
        fs = np.round(np.linspace(0.6 * fe, max(0.45, fe + 0.15), 26), 4)
        t0 = time.time()
        r = flow(m, L, seeds, fs, T_relax=max(2000, 20 * m), T_meas=max(4000, 20 * m))
        np.savez(f'{OUT}flow_m{m:g}.npz', m=m, L=L, A=1.0, **r)
        log('flow', m, f'{time.time() - t0:.0f}s')


def do_thresholds(cases=((1.0, (256, 1024, 4096)), (0.5, (1024,)))):
    for A, Ls in cases:
        tscale = A ** (-4 / 3)
        for L in Ls:
            seeds = np.arange(12) + 1000 + L
            for m in MS:
                fe = f_est(m, A)
                df = 0.0025 * A ** (4 / 3)
                t0 = time.time()
                rd = ramp_down(m, L, seeds, f0=1.25 * fe + 0.02 * A ** (4 / 3), df=df,
                               dwell=tscale * max(1000, 10 * m), A=A)
                t1 = time.time()
                ru = ramp_up(m, L, seeds, f0=0.5 * fe, df=df, D=RAMP_GO, A=A)
                np.savez(f'{OUT}thr_A{A:g}_L{L}_m{m:g}.npz', m=m, L=L, A=A,
                         fdown=rd['fdown'], fup=ru['fup'], ramp_f=rd['f'], ramp_v=rd['v'])
                log('thr', A, L, m, 'down', np.nanmean(rd['fdown']).round(4),
                    'up', np.nanmean(ru['fup']).round(4), f'{t1 - t0:.0f}+{time.time() - t1:.0f}s')


def do_well():
    L = 1024
    seeds = np.arange(3) + 77
    k0, wdot, T = 0.01, 1e-3, 4e5
    for m in [0.0, 64.0, 1024.0]:
        t0 = time.time()
        r = well(m, L, seeds, k0, wdot, T, nrec=max(1, int(round(2.0 / default_dt(m)))))
        np.savez(f'{OUT}well_m{m:g}.npz', m=m, L=L, k0=k0, wdot=wdot, **r)
        log('well', m, f'{time.time() - t0:.0f}s')


def do_kymo():
    L = 4096
    for m, f in [(0.0, 0.32), (1024.0, 0.09), (1024.0, 0.12)]:
        dt = default_dt(m)
        U = np.zeros((1, L)); V = np.full((1, L), f)
        seeds = [7]
        n0 = int(max(2000, 30 * m) / dt)
        _run(U, V, seeds, [f], m, n0, max(1, n0 // 10), 1.0, 1.0, 1.0, dt, 0.0, False)
        nch, T1 = 400, max(4000, 20 * m)
        ns = int(T1 / dt / nch)
        K = np.zeros((nch, L // 4), np.float32)
        for k in range(nch):
            _run(U, V, seeds, [f], m, ns, ns, 1.0, 1.0, 1.0, dt, 0.0, False)
            K[k] = V[0].reshape(-1, 4).mean(1)
        np.savez(f'{OUT}kymo_m{m:g}_f{f:g}.npz', m=m, f=f, L=L, T=T1, K=K)
        log('kymo', m, f)


if __name__ == '__main__':
    for job in sys.argv[1:]:
        {'flow': do_flow, 'thresholds': do_thresholds, 'well': do_well, 'kymo': do_kymo,
         'thr05': lambda: do_thresholds(((0.5, (1024,)),))}[job]()
