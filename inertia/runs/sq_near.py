"""S(q) of the hot line approaching f_min (m = 1024, random field): tests
T_lw = S_xi(0) / (2 df/dv), which diverges as df/dv -> 0 at v_min.
Writes data/sq_rf1024_f{f}.npz; S_xi(0) from sq_diag.py at the same f.
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import sys, time
import numpy as np
from sq import measure

for f in [float(x) for x in sys.argv[1:]]:
    t0 = time.time()
    r = measure(1024.0, f, 1.0)
    np.savez(f'data/sq_rf1024_f{f:g}.npz', **r)
    lam = 2 * (1 - np.cos(r['q']))
    print(time.strftime('%H:%M:%S'), 'near', f, 'v', round(r['v'], 4), 'Tkin', round(r['Tkin'], 3),
          'T_lw', round((lam * r['S'])[1:6].mean(), 3), f'{time.time() - t0:.0f}s', flush=True)
