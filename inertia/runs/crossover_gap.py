"""Second round of open points.  Results in data/*.npz.

    python3 runs4.py cross   # crossover 1 < M < 5: paired thresholds on a fine m grid (RF and RB)
    python3 runs4.py gap     # long lifetimes just below f_min at m = 64, 256, up to L = 16384
    python3 runs4.py heavy   # lifetimes at m = 16384 around f_min = 0.0412
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
from config import RAMP_GO
import sys, time
import numpy as np
from protocols import ramp_down, ramp_up, lifetimes


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def do_cross():
    L = 2048
    seeds = np.arange(24) + 7000          # same seeds at every m: paired differences
    for A, fc in [(1.0, 0.308), (-1.0, 0.514)]:
        for m in [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 11.0, 16.0]:
            t0 = time.time()
            rd = ramp_down(m, L, seeds, f0=fc + 0.04, df=0.001, dwell=max(2000.0, 20 * m), A=A)
            ru = ramp_up(m, L, seeds, f0=fc - 0.08, df=0.001, D=RAMP_GO, A=A)
            np.savez(f'data/cross_A{A:g}_m{m:g}.npz', m=m, A=A, L=L, fdown=rd['fdown'], fup=ru['fup'])
            log('cross', A, m, 'down', np.nanmean(rd['fdown']).round(4), 'up', np.nanmean(ru['fup']).round(4),
                f'{time.time() - t0:.0f}s')


def do_gap():
    for m, fs, fprep in [(256.0, [0.128, 0.132, 0.134, 0.136, 0.138], 0.17),
                         (64.0, [0.180, 0.188, 0.194, 0.198], 0.25)]:
        for L in [1024, 4096, 16384]:
            seeds = np.arange(16) + 9000 + L
            for f in fs:
                t0 = time.time()
                r = lifetimes(m, L, seeds, f, fprep, T_prep=max(4000.0, 4 * m), T_max=1000 * m)
                np.savez(f'data/gap_m{m:g}_L{L}_f{f:g}.npz', m=m, L=L, f=f, **r)
                log('gap', m, L, f, 'alive', np.isinf(r['tstop']).mean(), 'median t/m',
                    np.median(r['tstop']) / m, f'{time.time() - t0:.0f}s')


def do_heavy():
    m = 16384.0
    for L in [1024, 4096]:
        seeds = np.arange(16) + 11000 + L
        for f in [0.034, 0.037, 0.040, 0.042, 0.044]:
            t0 = time.time()
            r = lifetimes(m, L, seeds, f, 0.06, T_prep=4 * m, T_max=200 * m, dt=1.0)
            np.savez(f'data/life_m{m:g}_L{L}_f{f:g}.npz', m=m, L=L, f=f, **r)
            log('heavy', m, L, f, 'alive', np.isinf(r['tstop']).mean(), 'median t/m',
                np.median(r['tstop']) / m, 'v', np.nanmean(r['v']).round(4), f'{time.time() - t0:.0f}s')


if __name__ == '__main__':
    for job in sys.argv[1:]:
        {'cross': do_cross, 'gap': do_gap, 'heavy': do_heavy}[job]()
