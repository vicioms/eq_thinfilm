"""Runs for the open points of the note.  Results in data/*.npz.

    python3 runs2.py life      # quench lifetimes near f_down (m = 256, 1024; several L)
    python3 runs2.py heavy     # thresholds at m = 16384, 65536 (L = 1024, A = 1)
    python3 runs2.py rb        # random-bond thresholds (A = -1 in sim.py), L = 1024
"""
import sys, time
import numpy as np
from protocols import lifetimes, ramp_down, ramp_up

OUT = 'data/'


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def do_life():
    for m, fs, fprep in [(1024.0, [0.072, 0.076, 0.080, 0.084, 0.088, 0.092], 0.12),
                         (256.0, [0.108, 0.112, 0.116, 0.120, 0.124, 0.128], 0.17)]:
        for L in [512, 1024, 2048, 4096]:
            seeds = np.arange(16) + 5000 + L
            for f in fs:
                t0 = time.time()
                r = lifetimes(m, L, seeds, f, fprep, T_prep=max(4000.0, 4 * m), T_max=200 * m)
                np.savez(f'{OUT}life_m{m:g}_L{L}_f{f:g}.npz', m=m, L=L, f=f, **r)
                alive = np.isinf(r['tstop']).mean()
                log('life', m, L, f, 'alive', alive, 'median t', np.median(r['tstop']),
                    f'{time.time() - t0:.0f}s')


def do_heavy():
    L = 1024
    seeds = np.arange(12) + 1000 + L
    for m in [16384.0, 65536.0]:
        fe = 0.55 * m ** -0.25
        t0 = time.time()
        rd = ramp_down(m, L, seeds, f0=1.25 * fe + 0.02, df=0.0025, dwell=10 * m, dt=1.0)
        t1 = time.time()
        ru = ramp_up(m, L, seeds, f0=0.5 * fe, df=0.0025, D=30.0, dt=1.0)
        np.savez(f'{OUT}thr_A1_L{L}_m{m:g}.npz', m=m, L=L, A=1.0, fdown=rd['fdown'], fup=ru['fup'],
                 ramp_f=rd['f'], ramp_v=rd['v'])
        log('heavy', m, 'down', np.nanmean(rd['fdown']).round(4), 'up', np.nanmean(ru['fup']).round(4),
            f'{t1 - t0:.0f}+{time.time() - t1:.0f}s')


def do_rb():
    L = 1024
    seeds = np.arange(12) + 3000 + L
    for m in [0.0, 1.0, 4.0, 16.0, 64.0, 256.0, 1024.0, 4096.0]:
        fe = 0.51 if m == 0 else min(0.51, 0.5 * m ** -0.25)
        t0 = time.time()
        rd = ramp_down(m, L, seeds, f0=1.25 * fe + 0.03, df=0.0025, dwell=max(1000, 10 * m), A=-1.0)
        t1 = time.time()
        ru = ramp_up(m, L, seeds, f0=0.5 * fe, df=0.0025, D=30.0, A=-1.0)
        np.savez(f'{OUT}thr_A-1_L{L}_m{m:g}.npz', m=m, L=L, A=-1.0, fdown=rd['fdown'], fup=ru['fup'],
                 ramp_f=rd['f'], ramp_v=rd['v'])
        log('rb', m, 'down', np.nanmean(rd['fdown']).round(4), 'up', np.nanmean(ru['fup']).round(4),
            f'{t1 - t0:.0f}+{time.time() - t1:.0f}s')


if __name__ == '__main__':
    for job in sys.argv[1:]:
        {'life': do_life, 'heavy': do_heavy, 'rb': do_rb}[job]()
