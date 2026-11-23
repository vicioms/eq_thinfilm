"""Figure 12 and statistics6.json: avalanches under slow well driving (aval.py)."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from aval_merge import merge

L, k0 = 512, 3e-3
ms = [0, 1, 4, 16, 64, 256, 1024]
fmin = {64: 0.2008, 256: 0.1396, 1024: 0.0943}
cols = plt.cm.viridis(np.linspace(0, 0.9, len(ms)))
fig, ax = plt.subplots(1, 3, figsize=(14, 4))
bins = np.geomspace(0.3, 1e5, 22); bc = np.sqrt(bins[1:] * bins[:-1])
rows = []
for m, c in zip(ms, cols):
    r = np.load(f'data/aval_m{m}.npz'); s, T, fo = merge(r)
    h, _ = np.histogram(s, bins)
    dens = h / np.diff(bins) / len(s)
    ok = h > 0
    ax[0].loglog(bc[ok], dens[ok] * bc[ok], 'o-', ms=3, color=c, label=f'$m={m}$ ($M={0.37 * m:.2g}$)')
    big = s > 4 * L
    fs = fo[big] - k0 * s[big] / L
    rows.append(dict(m=m, M=0.37 * m, n=int(len(s)), n_slips=int(big.sum()),
                     slip_share=float(s[big].sum() / s.sum()),
                     slip_per_site=float(np.median(s[big]) / L) if big.any() else None,
                     cold_p90=float(np.percentile(s[~big], 90)),
                     f_on=float(fo[big].mean()), f_stop=float(fs.mean()),
                     f_k=float(((fo[big] + fs) / 2).mean()), mean_force=float(r['force'].mean()),
                     f_min=fmin.get(m)))
ax[0].set_xlabel('$s$'); ax[0].set_ylabel('$s\\,P(s)$'); ax[0].set_ylim(1e-4, 1)
ax[0].axvline(4 * L, color='0.6', ls=':', lw=0.8)
ax[0].legend(fontsize=7, loc='lower left'); ax[0].set_title('(a) avalanche sizes (fragments merged)')

for m, c, off in [(0, cols[0], 0.0), (64, cols[4], 0.0), (1024, cols[6], 0.0)]:
    r = np.load(f'data/aval_m{m}.npz')
    t = r['force_t'] * r['wdot']; F = r['force'][0]
    sel = (t > 20) & (t < 45)
    ax[1].plot(t[sel], F[sel], '-', lw=0.7, color=c, label=f'$m={m}$')
for m, f in fmin.items():
    ax[1].axhline(f, color='0.6', ls=':', lw=0.7)
ax[1].set_xlabel('$w$'); ax[1].set_ylabel('$k_0(w-\\bar u)$'); ax[1].legend(fontsize=7)
ax[1].set_title('(b) force, one sample; dotted: $f_{\\min}$ ($m=64,256,1024$)')

mm = [r_['m'] for r_ in rows if r_['f_min']]
for key, mk, lab in [('f_on', '^', 'slip onset'), ('f_k', 'o', '$(f_{\\rm on}+f_{\\rm stop})/2$'),
                     ('f_stop', 'v', 'slip arrest')]:
    ax[2].semilogx(mm, [r_[key] / r_['f_min'] for r_ in rows if r_['f_min']], mk + '-', label=lab)
ax[2].axhline(1, color='k', lw=0.6); ax[2].axhline(0.986, color='0.6', ls=':', lw=0.7)
ax[2].set_xlabel('$m$'); ax[2].set_ylabel('force $/f_{\\min}$'); ax[2].legend(fontsize=8)
ax[2].set_title('(c) slips as a spring block on the hot branch')
fig.tight_layout(); fig.savefig('fig12_avalanches.pdf'); fig.savefig('fig12_avalanches.png', dpi=110)
json.dump(rows, open('statistics6.json', 'w'), indent=1)
print('ok')
