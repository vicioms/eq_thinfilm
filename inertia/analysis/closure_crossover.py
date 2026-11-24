"""Figures and numbers for the second round of open points (fig7-fig9, data/stats/closure_crossover.json)."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import glob, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 9, 'axes.labelsize': 9, 'legend.fontsize': 7.5, 'savefig.bbox': 'tight'})
st = {}
FMIN = {64.0: 0.2008, 256.0: 0.1396, 1024.0: 0.0943, 16384.0: 0.0412}


def save(fig, name):
    fig.savefig(name + '.pdf'); fig.savefig(name + '.png', dpi=150); plt.close(fig)


# ---------------------------------------------------------------- fig7: closure check + edge speeds
fig, ax = plt.subplots(1, 3, figsize=(12, 3.3))
cc = sorted(glob.glob('data/closure_check_*.npz'), key=lambda s: float(s.split('_m')[1].split('_')[0]))
st['closure_check'] = {}
for k, fn in enumerate(cc):
    d = np.load(fn); m = float(d['m'])
    c = plt.get_cmap('viridis')(k / max(1, len(cc) - 1))
    x = d['taus'] * float(d['v'])
    ax[0].semilogx(x, d['B'] / (float(d['Tcl']) * d['bcl']), 'o-', ms=2.5, lw=0.8, color=c, label=f'm = {m:g}')
    ax[1].semilogx(x, d['kurt'], 'o-', ms=2.5, lw=0.8, color=c)
    i = np.argmin(np.abs(x - 1.0))
    st['closure_check'][str(m)] = dict(f=float(d['f']), v=float(d['v']), T_kin=float(d['Tkin']), T_closure=float(d['Tcl']),
                                       F_meas=float(d['F_meas']), F_closure=float(d['F_cl']),
                                       B_ratio_at_vtau_1=float(d['B'][i] / (d['Tcl'] * d['bcl'][i])))
ax[0].axhline(1, color='k', lw=0.5, ls='--'); ax[0].axvline(1, color='k', lw=0.5, ls=':')
ax[0].set(xlabel=r'$v\tau/r_f$', ylabel=r'$B(\tau)\,/\,T\,b(\tau)$'); ax[0].legend(frameon=False)
ax[0].set_title('(a) measured / closure displacement structure', fontsize=8.5)
ax[1].axhline(0, color='k', lw=0.5)
ax[1].set(xlabel=r'$v\tau/r_f$', ylabel='excess kurtosis of increments', ylim=(-0.5, 0.5))
ax[1].set_title('(b) Gaussianity of the fluctuations', fontsize=8.5)
st['edges'] = {}
for m, mk in [(1024.0, 'o'), (256.0, 's')]:
    fs = sorted(float(fn.split('_f')[1][:-4]) for fn in glob.glob(f'data/edge_m{m:g}_f*.npz'))
    if not fs:
        continue
    S = []
    for f in fs:
        d = np.load(f'data/edge_m{m:g}_f{f:g}.npz'); t = d['t']; lh = d['lh']
        i1 = len(t) // 4
        s = 0.5 * (lh[:, i1] - lh[:, 0]) / (t[i1] - t[0])
        S.append((s.mean(), s.std(ddof=1) / np.sqrt(len(s))))
    S = np.array(S); c = 1 / np.sqrt(m)
    ax[2].errorbar(np.array(fs) / FMIN[m], S[:, 0] / c, S[:, 1] / c, fmt=mk + '-', ms=4, label=f'm = {m:g}')
    st['edges'][str(m)] = dict(f=fs, s=list(S[:, 0]), s_err=list(S[:, 1]), c=c)
    # zero crossing by linear interpolation
    j = np.where(np.diff(np.sign(S[:, 0])) > 0)[0]
    if len(j):
        j = j[0]; f0 = fs[j] - S[j, 0] * (fs[j + 1] - fs[j]) / (S[j + 1, 0] - S[j, 0])
        st['edges'][str(m)]['f0'] = f0; st['edges'][str(m)]['f0_over_fmin'] = f0 / FMIN[m]
ax[2].axhline(0, color='k', lw=0.5); ax[2].axvline(1, color='k', lw=0.5, ls='--')
ax[2].axhline(1, color='k', lw=0.5, ls=':')
ax[2].set(xlabel=r'$f/f_{\min}$', ylabel=r'edge speed $s/c$', ylim=(-1.2, 1.2))
ax[2].legend(frameon=False); ax[2].set_title('(c) hot/pinned edge speed (positive: hot invades)', fontsize=8.5)
fig.tight_layout(); save(fig, 'figures/closure_edges')

# ---------------------------------------------------------------- fig8: crossover
cr = {}
for fn in glob.glob('data/cross_*.npz'):
    d = np.load(fn); cr[(float(d['A']), float(d['m']))] = d
if cr:
    fig, ax = plt.subplots(1, 1, figsize=(5.2, 3.4)); ax = [ax]
    st['crossover'] = {}
    for A, mk, lab, Lc2 in [(1.0, 'o', 'random field', 1.65 ** 2), (-1.0, 's', 'random bond', None)]:
        ms = sorted(m for (a, m) in cr if a == A)
        if 0.0 not in ms:
            continue
        d0, u0 = cr[(A, 0.0)]['fdown'], cr[(A, 0.0)]['fup']
        rows = []
        for m in ms:
            if m == 0:
                continue
            dd = cr[(A, m)]['fdown'] - d0; du = cr[(A, m)]['fup'] - u0
            fc = 0.5 * (np.nanmean(d0) + np.nanmean(u0))
            rows.append((m, -np.nanmean(dd) / fc, np.nanstd(dd) / np.sqrt(np.isfinite(dd).sum()) / fc,
                         -np.nanmean(du) / fc, np.nanstd(du) / np.sqrt(np.isfinite(du).sum()) / fc))
        r = np.array(rows)
        ax[0].errorbar(r[:, 0], r[:, 1], r[:, 2], fmt=mk, mfc='none', color='C0', label=f'$1-f_\\downarrow/f_c$, {lab}')
        ok = np.array([not (A < 0 and x >= 6) for x in r[:, 0]])   # RB ramp-up censored at its start force for m >= 6
        ax[0].errorbar(r[ok, 0], r[ok, 3], r[ok, 4], fmt=mk, color='C3', ms=3.5, label=f'$1-f_\\uparrow/f_c$, {lab}')
        st['crossover'][str(A)] = [dict(m=x[0], shift_down=x[1], err_down=x[2], shift_up=x[3], err_up=x[4]) for x in rows]
        pos = r[:, 1] > 2 * r[:, 2]
        if pos.sum() >= 3:
            p = np.polyfit(np.log(r[pos, 0]), np.log(r[pos, 1]), 1)
            st['crossover'][f'exponent_down_{A:g}'] = float(p[0])
    ax[0].set(xscale='log', yscale='symlog', xlabel='m', ylabel='relative threshold shift')
    ax[0].set_yscale('symlog', linthresh=1e-3)
    ax[0].legend(frameon=False, fontsize=6.5); ax[0].set_title('(a) shift below $f_c$ (paired samples, L = 2048)', fontsize=8.5)
    fig.tight_layout(); save(fig, 'figures/crossover')

# ---------------------------------------------------------------- fig9: gap and heavy lifetimes
fig, ax = plt.subplots(1, 3, figsize=(12, 3.3))
st['gap'] = {}
for k, m in enumerate([256.0, 64.0]):
    for L, c in [(1024, 'C1'), (4096, 'C3'), (16384, 'C4')]:
        fs = sorted(float(fn.split('_f')[1][:-4]) for fn in glob.glob(f'data/gap_m{m:g}_L{L}_f*.npz'))
        if not fs:
            continue
        al = [np.isinf(np.load(f'data/gap_m{m:g}_L{L}_f{f:g}.npz')['tstop']).mean() for f in fs]
        ax[k].plot(np.array(fs) / FMIN[m], al, 'o-', color=c, label=f'L = {L}')
        st['gap'][f'm{m:g}_L{L}'] = dict(f=fs, alive=al)
    ax[k].axvline(1, color='k', lw=0.5, ls='--')
    ax[k].set(xlabel=r'$f/f_{\min}$', ylabel='fraction sliding after $1000m$', ylim=(-0.05, 1.05))
    ax[k].legend(frameon=False); ax[k].set_title(f'({"ab"[k]}) m = {m:g}', fontsize=8.5)
st['heavy'] = {}
for L, c in [(1024, 'C1'), (4096, 'C3')]:
    fs = sorted(float(fn.split('_f')[1][:-4]) for fn in glob.glob(f'data/life_m16384_L{L}_f*.npz'))
    if not fs:
        continue
    al = [np.isinf(np.load(f'data/life_m16384_L{L}_f{f:g}.npz')['tstop']).mean() for f in fs]
    ax[2].plot(np.array(fs) / FMIN[16384.0], al, 'o-', color=c, label=f'L = {L}')
    st['heavy'][f'L{L}'] = dict(f=fs, alive=al)
ax[2].axvline(1, color='k', lw=0.5, ls='--')
ax[2].set(xlabel=r'$f/f_{\min}$', ylabel='fraction sliding after $200m$', ylim=(-0.05, 1.05))
ax[2].legend(frameon=False); ax[2].set_title('(c) m = 16384', fontsize=8.5)
fig.tight_layout(); save(fig, 'figures/long_runs')

json.dump(st, open('data/stats/closure_crossover.json', 'w'), indent=1, default=float)
print(json.dumps(st, indent=1, default=float)[:4000])
