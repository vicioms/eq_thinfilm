"""Flow curves f(v) and their minimum from the self-energy Gaussian DMFT (gdmft2.py), m = 64, 256, 1024.
Compare: closure f_min 0.2008, 0.1396, 0.0943; measured f_inf ~ 0.198, 0.137, 0.092.
Writes data/gdmft2_flow.npz."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import numpy as np
from gdmft2 import solve

out = {}
for m, vs in [(1024.0, [0.2, 0.12, 0.09, 0.07, 0.06, 0.05, 0.045, 0.04, 0.035, 0.03, 0.025, 0.02, 0.015, 0.01, 0.007, 0.005]),
              (256.0, [0.3, 0.16, 0.12, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 0.006]),
              (64.0, [0.4, 0.25, 0.2, 0.17, 0.15, 0.13, 0.12, 0.11, 0.1, 0.09, 0.08, 0.06, 0.04, 0.03, 0.02, 0.01])]:
    rows = []
    for v in vs:
        r = solve(m, v, verbose=False, it=120)
        rows.append([v, r['f'], r['Sxi0'], r['eta_eff'], r['T_lw'], r['Tkin'], r['T_pb'], r['stab'],
                     float(r['converged'])])
        print(f'm {m:g} v {v}: f {r["f"]:.4f} S_xi0 {r["Sxi0"]:.3f} eta_eff {r["eta_eff"]:.3f} '
              f'T_kin {r["Tkin"]:.3f} (pb {r["T_pb"]:.3f}) min Im z/w {r["stab"]:.3f} conv {r["converged"]}',
              flush=True)
    out[f'm{m:g}'] = np.array(rows)
np.savez('data/gdmft2_flow.npz', cols=['v', 'f', 'Sxi0', 'eta_eff', 'T_lw', 'Tkin', 'T_pb', 'stab', 'converged'],
         **out)
print('done', flush=True)
