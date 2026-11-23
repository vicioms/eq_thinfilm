"""Size dependence of the edge state at m = 64 (f = 0.1955, 0.196): weak jump vs critical point.

A first-order end keeps the survivors' active fraction ~0.1 at every L; an
absorbing-state critical point gives activity ~ L^(-beta/nu_perp) (DP: 0.25) and
lifetimes growing as a power of L.  Writes data/endslide_m64_L{L}_f{f}.npz.
"""
import time
import numpy as np
from protocols import lifetimes

m = 64.0
for L in [4096, 65536]:
    for f in [0.1955, 0.196]:
        t0 = time.time()
        seeds = np.arange(12) + 19000 + L
        r = lifetimes(m, L, seeds, f, 0.25, T_prep=4000.0, T_max=4000 * m, T_skip=500 * m)
        np.savez(f'data/endslide_m64_L{L}_f{f:g}.npz', m=m, L=L, f=f, **r)
        a = np.isinf(r['tstop'])
        print(time.strftime('%H:%M:%S'), 'endL', L, f, 'alive', round(a.mean(), 3),
              'v', np.round(np.nanmean(r['v'][a]), 4) if a.any() else '-',
              'act', np.round(np.nanmean(r['act'][a]), 3) if a.any() else '-',
              'median stop/m', np.round(np.median(r['tstop'][~a]) / m, 0) if (~a).any() else '-',
              f'{time.time() - t0:.0f}s', flush=True)
