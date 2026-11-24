"""Self-consistent cavity solutions at the four hot-line points measured by sq_diag.py."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import numpy as np
from cavity import solve

CASES = [(1024.0, 0.0891), (1024.0, 0.0708), (1024.0, 0.0512), (256.0, 0.1169)]
for m, v in CASES:
    r, hist, C, tC, gamma, f0 = solve(m, v, it=6, T=1.5e5, T0=1.5e4)
    np.savez(f'data/cavity_m{m:g}_v{v:g}.npz', hist=np.array([list(h.values()) for h in hist]),
             hist_keys=list(hist[0].keys()), C_final=C, tC=tC, gamma=gamma, f=f0, **r)
print('done', flush=True)
