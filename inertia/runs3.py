"""Lifetimes straddling the closure f_min (m = 256: f_min = 0.1396; m = 1024: 0.0943)."""
import time
import numpy as np
from protocols import lifetimes

for m, fs, fprep in [(1024.0, [0.094, 0.096, 0.098], 0.12), (256.0, [0.132, 0.136, 0.140, 0.144], 0.17)]:
    for L in [1024, 4096]:
        seeds = np.arange(16) + 5000 + L
        for f in fs:
            t0 = time.time()
            r = lifetimes(m, L, seeds, f, fprep, T_prep=max(4000.0, 4 * m), T_max=200 * m)
            np.savez(f'data/life_m{m:g}_L{L}_f{f:g}.npz', m=m, L=L, f=f, **r)
            print(time.strftime('%H:%M:%S'), 'life', m, L, f, 'alive', np.isinf(r['tstop']).mean(),
                  'median t', np.median(r['tstop']), 'v', np.nanmean(r['v']).round(4),
                  'act', np.nanmean(r['act']).round(3), f'{time.time() - t0:.0f}s', flush=True)
