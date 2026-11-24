"""Flow curves f(v) and f_min from the Gaussian DMFT (gdmft.py) for m = 1024, 256, 64.
Writes data/gdmft_flow.npz.  Closure minima for comparison: 0.0943, 0.1396, 0.2008."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import numpy as np
from gdmft import solve

out = {}
for m, vs in [(1024.0, [0.4, 0.2, 0.12, 0.08, 0.06, 0.05, 0.046, 0.044, 0.042, 0.04, 0.035]),
              (256.0, [0.5, 0.25, 0.16, 0.12, 0.1, 0.09, 0.08, 0.075, 0.07, 0.065, 0.06]),
              (64.0, [0.5, 0.3, 0.22, 0.18, 0.16, 0.14, 0.13, 0.12, 0.11, 0.1])]:
    rows = []
    for v in vs:
        r = solve(m, v, verbose=False)
        rows.append([v, r['f'], r['Sxi0'], r['eta_eff'], r['T_lw'], r['Tkin'], float(r['converged']),
                     float(r['no_solution'])])
        print(f'm {m:g} v {v}: f {r["f"]:.4f} S_xi0 {r["Sxi0"]:.3f} (bare {1 / (3 * v):.3f}) eta_eff {r["eta_eff"]:.3f} '
              f'T_lw {r["T_lw"]:.3f} T_kin {r["Tkin"]:.3f} conv {r["converged"]} no-solution {r["no_solution"]}', flush=True)
    out[f'm{m:g}'] = np.array(rows)
np.savez('data/gdmft_flow.npz', cols=['v', 'f', 'Sxi0', 'eta_eff', 'T_lw', 'Tkin', 'converged', 'no_solution'], **out)
print('done', flush=True)
