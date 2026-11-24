"""End of sliding at m = 64: is the drop of the velocity at f_inf continuous or a jump?

Fine force grid just above the long-time threshold (f_inf ~ 0.195 = 0.97 f_min),
runs of 4000m after a quench from sliding at f = 0.25, L = 16384, 12 samples.
For survivors: time-averaged velocity and active fraction after 500m; for all: stop times.
Writes data/endslide_m64_f{f}.npz.  Run with NUMBA_NUM_THREADS capped.
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import time
import numpy as np
from protocols import lifetimes

m, L = 64.0, 16384
for f in [0.1945, 0.1950, 0.1955, 0.1960, 0.1970, 0.1985]:
    t0 = time.time()
    seeds = np.arange(12) + 17000
    r = lifetimes(m, L, seeds, f, 0.25, T_prep=4000.0, T_max=4000 * m, T_skip=500 * m)
    np.savez(f'data/endslide_m64_f{f:g}.npz', m=m, L=L, f=f, **r)
    a = np.isinf(r['tstop'])
    print(time.strftime('%H:%M:%S'), 'endslide', f, 'alive', round(a.mean(), 3),
          'v', np.round(np.nanmean(r['v'][a]), 4) if a.any() else '-',
          'sd_v', np.round(np.nanstd(r['v'][a]), 4) if a.sum() > 1 else '-',
          'act', np.round(np.nanmean(r['act'][a]), 3) if a.any() else '-',
          'median stop/m', np.round(np.median(r['tstop'][~a]) / m, 0) if (~a).any() else '-',
          f'{time.time() - t0:.0f}s', flush=True)
