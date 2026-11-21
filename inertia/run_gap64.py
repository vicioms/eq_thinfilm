import numpy as np, time
from protocols import lifetimes
m = 64.0
for f in [0.188, 0.191, 0.194]:
    seeds = np.arange(16) + 9000 + 32768
    r = lifetimes(m, 32768, seeds, f, 0.25, T_prep=4000.0, T_max=1000 * m)
    np.savez(f'data/gap_m64_L32768_f{f:g}.npz', m=m, L=32768, f=f, **r)
    print(time.strftime('%H:%M:%S'), 'gap', m, 32768, f, 'alive', np.isinf(r['tstop']).mean(), 'median t/m', np.median(r['tstop']) / m, flush=True)
