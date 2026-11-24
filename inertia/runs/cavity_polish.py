"""Polish the self-consistent cavity solutions (cavity_runs.py): bigger ring, more rings,
central-difference eta_eff.  Writes data/cavity_polish_m{m}_v{v}.npz."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import time
import numpy as np
from cavity import polish

CASES = [(1024.0, 0.0891), (1024.0, 0.0708), (1024.0, 0.0512), (256.0, 0.1169)]
for m, v in CASES:
    t0 = time.time()
    c = np.load(f'data/cavity_m{m:g}_v{v:g}.npz')
    keys = list(c['hist_keys']); H = c['hist']
    gamma = float(np.mean(H[-3:, keys.index('gamma')]))
    f0 = float(np.mean(H[-3:, keys.index('f')]))
    r = polish(m, v, c['C_final'], c['tC'], gamma, f0)
    np.savez(f'data/cavity_polish_m{m:g}_v{v:g}.npz', **r)
    print(time.strftime('%H:%M:%S'), f'm {m:g} v {v}: f {r["f"]:.4f} S_xi(0) {r["Sxi0"]:.3f} '
          f'eta_eff {r["eta_eff"]:.3f} (gamma in {gamma:.3f}) T_lw {r["Sxi0"] / (2 * r["eta_eff"]):.3f} '
          f'Tkin {r["Tkin"]:.3f}  {time.time() - t0:.0f}s', flush=True)
print('done', flush=True)
