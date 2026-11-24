"""Figure 10: end of sliding at m = 64 (endslide.py, endslide_L.py, endslide_t.py)."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fmin = 0.2008
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
cols = {0.1955: 'C0', 0.196: 'C1', 0.1975: 'C2', 0.1985: 'C3'}
for f, c in cols.items():
    d = np.load(f'data/endslide_t_f{f:g}.npz')
    t = d['trace_t'] / 64; a = d['trace_act'].mean(axis=1)
    edges = np.geomspace(20, t[-1] * 1.001, 30)
    tc, ac = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (t >= lo) & (t < hi)
        if s.any() and a[s].mean() > 0:
            tc.append(np.sqrt(lo * hi)); ac.append(a[s].mean())
    ax[0].loglog(tc, ac, 'o-', ms=3, color=c, label=f'$f/f_{{\\min}}={f / fmin:.3f}$')
tt = np.geomspace(300, 12000, 10)
ax[0].loglog(tt, 0.75 * (tt / 300) ** -0.159, 'k:', label=r'DP $t^{-0.159}$')
ax[0].set_xlabel('$t/m$'); ax[0].set_ylabel('active fraction (all lines)')
ax[0].set_ylim(1e-3, 1); ax[0].legend(fontsize=8, loc='lower left')
ax[0].set_title('(a) $m=64$, $L=65536$')

for f, mk in [(0.1955, 'o'), (0.196, 's')]:
    Ls, act = [], []
    for L in [4096, 16384, 65536]:
        fn = f'data/endslide_m64_L{L}_f{f:g}.npz' if L != 16384 else f'data/endslide_m64_f{f:g}.npz'
        try:
            d = np.load(fn)
        except FileNotFoundError:
            continue
        alive = ~np.isfinite(d['tstop'])
        if alive.any():
            Ls.append(L); act.append(np.nanmean(d['act'][alive]))
    ax[1].loglog(Ls, act, mk + '-', label=f'$f/f_{{\\min}}={f / fmin:.3f}$')
LL = np.array([4096, 65536])
ax[1].loglog(LL, 0.2 * (LL / 16384) ** -0.25, 'k:', label=r'$L^{-0.25}$')
ax[1].set_xlabel('$L$'); ax[1].set_ylabel('survivors\' active fraction at $4000m$')
from matplotlib.ticker import ScalarFormatter, NullFormatter
ax[1].yaxis.set_major_formatter(ScalarFormatter()); ax[1].yaxis.set_minor_formatter(NullFormatter())
ax[1].set_yticks([0.1, 0.15, 0.2, 0.25])
ax[1].legend(fontsize=8); ax[1].set_title('(b) fixed-time survivor average')
fig.tight_layout(); fig.savefig('figures/end_of_sliding.pdf'); fig.savefig('figures/end_of_sliding.png', dpi=120)
print('ok')
