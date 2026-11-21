"""Figures and numbers for the open points: lifetimes, heavy lines, random bond, local balance."""
import glob, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 9, 'axes.labelsize': 9, 'legend.fontsize': 7.5, 'savefig.bbox': 'tight'})
st = {}


def save(fig, name):
    fig.savefig(name + '.pdf'); fig.savefig(name + '.png', dpi=150); plt.close(fig)


def scmin(name, m):
    d = np.load(f'data/sc_{name}.npz')
    f = d[f'f_{m:g}']; i = np.argmin(f)
    return float(f[i]), float(d['vs'][i]), float(d[f'T_{m:g}'][i])


def fmin_of(A, m):
    for name in ([f'A{A:g}', 'A1_heavy'] if A == 1 else [f'A{A:g}']):
        try:
            return scmin(name, m)[0]
        except (KeyError, FileNotFoundError):
            pass
    return np.nan


# ------------------------------------------------------------- lifetimes
life = {}
for fn in glob.glob('data/life_*.npz'):
    d = np.load(fn)
    life[(float(d['m']), int(d['L']), float(d['f']))] = d
if life:
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.3))
    cols = {512: 'C0', 1024: 'C1', 2048: 'C2', 4096: 'C3'}
    st['lifetimes'] = {}
    for k, (m, mk) in enumerate([(1024.0, 'o'), (256.0, 's')]):
        fm = fmin_of(1.0, m)
        for L in [512, 1024, 2048, 4096]:
            fs = sorted(f for (mm, LL, f) in life if mm == m and LL == L)
            if not fs:
                continue
            med = np.array([np.median(life[(m, L, f)]['tstop']) for f in fs])
            alive = np.array([np.isinf(life[(m, L, f)]['tstop']).mean() for f in fs])
            x = fm - np.array(fs)
            ok = np.isfinite(med) & (x > 0)
            ax[k].loglog(x[ok], med[ok] / m, mk + '-', ms=4, color=cols[L], label=f'L = {L}')
            st['lifetimes'][f'm{m:g}_L{L}'] = {'f': fs, 'median_t_over_m': list(med / m), 'alive': list(alive)}
            if L == 4096 and ok.sum() >= 3:
                p = np.polyfit(np.log(x[ok]), np.log(med[ok]), 1)
                st['lifetimes'][f'm{m:g}_exponent'] = float(-p[0])
        xx = np.geomspace(2e-3, 3e-2, 10)
        ax[k].loglog(xx, 30.0 * (xx / 0.01) ** -1, 'k:', lw=0.8, label=r'$(f_{\min}-f)^{-1}$')
        ax[k].set(xlabel=r'$f_{\min}-f$  ($f_{\min}$ from the closure)', ylabel=r'median lifetime $/\,m$')
        ax[k].set_title(f'({"ab"[k]}) m = {m:g}, quench from sliding, $f_{{\\min}}$ = {fm:.4f}', fontsize=8.5)
        ax[k].legend(frameon=False)
    fig.tight_layout(); save(fig, 'fig5_lifetimes')

# ------------------------------------------------------------- thresholds incl. heavy and random bond
thr = {}
for fn in glob.glob('data/thr_*.npz'):
    d = np.load(fn)
    thr[(float(d['A']), int(d['L']), float(d['m']))] = d
rows = []
for (A, L, m), d in sorted(thr.items()):
    if L != 1024 or m < 16:
        continue
    F2 = A ** 2 / 3 if A < 0 else 2 * A ** 2 / 9
    fstar = np.sqrt(2 * F2) * m ** -0.25
    fm = fmin_of(A, m)
    rows.append(dict(A=A, m=m, fdown=float(np.nanmean(d['fdown'])), fup=float(np.nanmean(d['fup'])),
                     fstar=float(fstar), fmin=fm))
st['thresholds_L1024'] = rows
if rows:
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.3))
    for A, mk, lab in [(1.0, 'o', 'random field'), (-1.0, 's', 'random bond')]:
        r = [x for x in rows if x['A'] == A]
        if not r:
            continue
        m = np.array([x['m'] for x in r])
        ax[0].loglog(m, [x['fdown'] for x in r], mk, mfc='none', color='C0', label=f'$f_\\downarrow$ (dwell 10m), {lab}')
        ax[0].loglog(m, [x['fup'] for x in r], mk, color='C3', ms=3.5, label=f'$f_\\uparrow$, {lab}')
        ok = np.isfinite([x['fmin'] for x in r])
        ax[0].loglog(m[ok], np.array([x['fmin'] for x in r])[ok], '-' if A > 0 else '--', color='k', lw=0.9,
                     label=f'closure $f_{{\\min}}$, {lab}')
        ax[1].semilogx(m[ok], (np.array([x['fdown'] for x in r]) / np.array([x['fmin'] for x in r]))[ok], mk + '-',
                       mfc='none', color='C0', label=f'$f_\\downarrow/f_{{\\min}}$, {lab}')
        ax[1].semilogx(m, [x['fmin'] / x['fstar'] for x in r], mk + ':', color='k', ms=3.5, label=f'$f_{{\\min}}/f_*$, {lab}')
    mm = np.geomspace(16, 7e4, 10)
    ax[0].loglog(mm, 0.667 * mm ** -0.25, ':', color='0.5', lw=0.8, label=r'$f_*$ (random field)')
    ax[0].set(xlabel='m', ylabel='force'); ax[0].legend(frameon=False, fontsize=6.5)
    ax[0].set_title('(a) thresholds, L = 1024, A = 1', fontsize=8.5)
    ax[1].set(xlabel='m', ylabel='ratio', ylim=(0.3, 1.1)); ax[1].legend(frameon=False, fontsize=6.5)
    ax[1].set_title('(b) ramp-down $f_\\downarrow$ and closure $f_{\\min}$', fontsize=8.5)
    fig.tight_layout(); save(fig, 'fig6_heavy_rb')

json.dump(st, open('statistics2.json', 'w'), indent=1, default=float)
for r in rows:
    print(f"A={r['A']:g} m={r['m']:7g}  f_down={r['fdown']:.4f} f_up={r['fup']:.4f} f_min={r['fmin']:.4f} f*={r['fstar']:.4f}"
          f"  f_down/f_min={r['fdown'] / r['fmin']:.3f}  f_min/f*={r['fmin'] / r['fstar']:.3f}")
print({k: v for k, v in st.get('lifetimes', {}).items() if 'exponent' in k})
