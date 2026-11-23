"""Time decay of activity at the end of sliding (m = 64, L = 65536).

Fixed-time survivor averages mix finite-size conditioning into both scenarios. The
infinite-system decay rho(t) separates them: an absorbing-state critical point gives
rho ~ t^(-delta) (DP: delta = 0.159); a weak first-order end gives a plateau
(sliding branch) or exponential decay (pinned branch). Writes data/endslide_t_f{f}.npz.
"""
import time
import numpy as np
from protocols import lifetimes

m, L = 64.0, 65536
import sys
FS = [float(x) for x in sys.argv[1:]] or [0.196, 0.1955]
NS = 6 if len(sys.argv) < 2 else 4
for f in FS:
    t0 = time.time()
    seeds = np.arange(NS) + 31000
    r = lifetimes(m, L, seeds, f, 0.25, T_prep=4000.0, T_max=12000 * m, T_skip=500 * m,
                  trace=True)
    np.savez(f'data/endslide_t_f{f:g}.npz', m=m, L=L, f=f, **r)
    a = r['trace_act'].mean(axis=1); tt = r['trace_t'] / m
    pick = [np.searchsorted(tt, x) for x in [250, 1000, 3000, 6000, 11900] if x < tt[-1]]
    print(time.strftime('%H:%M:%S'), 'endt', f, 'alive', np.isfinite(r['tstop']).__invert__().mean(),
          'rho(t/m)', ' '.join(f'{tt[i]:.0f}:{a[i]:.3f}' for i in pick),
          f'{time.time() - t0:.0f}s', flush=True)
