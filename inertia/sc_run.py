"""Self-consistent flow curves (sc.py) for the masses and amplitudes of the runs -> data/sc_A{A}.npz."""
import numpy as np
from sc import Line

MS = [4.0, 16.0, 64.0, 256.0, 1024.0, 4096.0]
for A in [1.0, 0.5]:
    vs = np.geomspace(0.002, 0.5, 48) * A
    out = {'vs': vs, 'ms': np.array(MS)}
    for m in MS:
        Lm = Line(m, A=A, tau_max=max(20000.0, 12 * m) / A)
        f, T, Fp = Lm.flow(vs)
        out[f'f_{m:g}'] = f; out[f'T_{m:g}'] = T
        i = np.argmin(f)
        print(A, m, 'f_min', round(f[i], 4), 'v', round(vs[i], 4), 'T', round(T[i], 3), flush=True)
    np.savez(f'data/sc_A{A:g}.npz', **out)
