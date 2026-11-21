"""Patchy sliding phase near the crossover: stationary velocity and activity vs f (m = 64).

Quench a sliding line (prepared at f = 0.25) to f, run up to T_max, record whether it
survives and, for survivors, the time-averaged velocity and active fraction.
Writes data/patchy_m64_L{L}_f{f}.npz.  Run with NUMBA_NUM_THREADS capped.
"""
import sys, time
import numpy as np
from protocols import lifetimes

m = 64.0
L = int(sys.argv[1]) if len(sys.argv) > 1 else 16384
fs = [0.189, 0.190, 0.192, 0.193, 0.195, 0.196, 0.200, 0.205]
for f in fs:
    t0 = time.time()
    seeds = np.arange(16) + 13000 + L
    r = lifetimes(m, L, seeds, f, 0.25, T_prep=4000.0, T_max=2000 * m, T_skip=200 * m)
    np.savez(f'data/patchy_m64_L{L}_f{f:g}.npz', m=m, L=L, f=f, **r)
    a = np.isinf(r['tstop'])
    print(time.strftime('%H:%M:%S'), 'patchy', L, f, 'alive', a.mean(),
          'v', np.round(np.nanmean(r['v'][a]), 4) if a.any() else '-',
          'act', np.round(np.nanmean(r['act'][a]), 3) if a.any() else '-', f'{time.time() - t0:.0f}s', flush=True)
