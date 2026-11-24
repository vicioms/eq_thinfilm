"""Closure for the open points: random bond (A = -1) and heavier lines (A = 1)."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import numpy as np
from sc import Line

jobs = [('A-1', -1.0, [16.0, 64.0, 256.0, 1024.0, 4096.0], np.geomspace(0.002, 0.5, 48)),
        ('A1_heavy', 1.0, [16384.0, 65536.0], np.geomspace(0.001, 0.2, 44))]
for name, A, MS, vs in jobs:
    out = {'vs': vs, 'ms': np.array(MS)}
    for m in MS:
        dtau = 0.05 if m <= 4096 else 0.2
        Lm = Line(m, A=A, dtau=dtau, tau_max=max(20000.0, 12 * m))
        f, T, Fp = Lm.flow(vs)
        out[f'f_{m:g}'] = f; out[f'T_{m:g}'] = T
        i = np.argmin(f)
        print(name, m, 'f_min', round(f[i], 4), 'v', round(vs[i], 4), 'T', round(T[i], 3), flush=True)
    np.savez(f'data/sc_{name}.npz', **out)
