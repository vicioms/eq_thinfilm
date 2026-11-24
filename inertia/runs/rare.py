"""Are the long-lived active patches below f_inf quenched rare regions? (m = 64, L = 65536)

endslide_t.py found, at 0.974 and 0.976 f_min, one sample in six keeping ~1000 active
sites for thousands of m after the rest had stopped.  A Griffiths / infinite-randomness
scenario needs such patches to be quenched rare regions, sitting at a fixed x set by the
disorder.  Here each of those two samples (same disorder seed) is run twice, from lines
displaced by a constant offset before the preparation, so the configurations differ but
the disorder is the same.  The activity is stored in blocks of 64 sites every chunk.

Writes data/rare.npz.
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
from config import V_ACT
import time
import numpy as np
from protocols import _run, default_dt

m, L, B = 64.0, 65536, 64
seeds = np.array([31003, 31003, 31005, 31005], np.int64)
fs = np.array([0.1955, 0.1955, 0.196, 0.196])
offs = np.array([0.0, 37.3, 0.0, 37.3])
f_prep, T_prep, T_max, chunk, vact = 0.25, 4000.0, 12000 * m, 200.0, V_ACT
dt = default_dt(m)
S = len(seeds)
U = np.zeros((S, L)) + offs[:, None]; V = np.full((S, L), f_prep)
n0 = int(round(T_prep / dt))
_run(U, V, seeds, np.full(S, f_prep), m, n0, max(1, n0 // 10), 1.0, 1.0, 1.0, dt, 0.0, False)
ns = int(round(chunk / dt))
amap, umap, tt = [], [], []
t = 0.0; t0 = time.time()
while t < T_max:
    _run(U, V, seeds, fs, m, ns, ns // 4, 1.0, 1.0, 1.0, dt, 0.0, False)
    t += ns * dt
    tt.append(t)
    amap.append((np.abs(V) > vact).reshape(S, L // B, B).mean(2).astype(np.float32))
    umap.append(U.reshape(S, L // B, B).mean(2).astype(np.float32))
    if len(tt) % 400 == 0:
        print(time.strftime('%H:%M:%S'), f't/m {t / m:.0f}', 'rho', np.round(amap[-1].mean(1), 4),
              f'{time.time() - t0:.0f}s', flush=True)
np.savez('data/rare.npz', t=np.array(tt), act=np.array(amap), u=np.array(umap), seeds=seeds, f=fs,
         offs=offs, m=m, L=L, B=B)
print('done', flush=True)
