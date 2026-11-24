"""Figures and numbers for the inertial-depinning note (reads data/*.npz)."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import glob, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from linresp import flow_curve, power

plt.rcParams.update({'font.size': 9, 'axes.labelsize': 9, 'legend.fontsize': 7.5,
                     'figure.dpi': 150, 'savefig.bbox': 'tight'})
MS = [0.0, 0.25, 1.0, 4.0, 16.0, 64.0, 256.0, 1024.0, 4096.0]
cmap = plt.get_cmap('viridis')
col = {m: cmap(i / (len(MS) - 1)) for i, m in enumerate(MS)}
F2 = lambda A: 2 * A ** 2 / 9                 # <F^2> of the interpolated force
Lc = lambda A: (1.0 / np.sqrt(F2(A))) ** (2 / 3)  # Larkin length, kappa = r_f = 1
fstar = lambda m, A: np.sqrt(2 * F2(A)) * m ** -0.25  # radiation minimum, eta = kappa = 1
vstar = lambda m, A: np.sqrt(F2(A) / 2) * m ** -0.25
stats = {}
SC = {}
for A in [1.0, 0.5]:
    try:
        d = np.load(f'data/sc_A{A:g}.npz')
        SC[A] = {float(m): (d['vs'], d[f'f_{m:g}'], d[f'T_{m:g}']) for m in d['ms'] if f'f_{m:g}' in d}
    except FileNotFoundError:
        pass
def sc_min(A, m):
    vs, f, T = SC[A][m]
    i = np.argmin(f)
    return f[i], vs[i], T[i]


def save(fig, name):
    fig.savefig(name + '.pdf'); fig.savefig(name + '.png', dpi=150); plt.close(fig)


def load_thr():
    T = {}
    for fn in glob.glob('data/thr_*.npz'):
        d = np.load(fn)
        T[(float(d['A']), int(d['L']), float(d['m']))] = d
    return T


def mean_err(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) == 0:
        return np.nan, np.nan
    return x.mean(), x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan


# ---------------------------------------------------------------- thresholds
T = load_thr()
fc = {}
for A in [1.0, 0.5]:
    Ls = sorted({L for (a, L, m) in T if a == A})
    if not Ls:
        continue
    d = T.get((A, Ls[-1], 0.0))
    if d is not None:
        fc[A] = 0.5 * (np.nanmean(d['fdown']) + np.nanmean(d['fup']))
stats['f_c'] = {str(A): fc[A] for A in fc}
stats['L_c'] = {str(A): Lc(A) for A in fc}
thr_rows = []
for (A, L, m), d in sorted(T.items()):
    fd, efd = mean_err(d['fdown']); fu, efu = mean_err(d['fup'])
    thr_rows.append(dict(A=A, L=L, m=m, M=m / Lc(A) ** 2, fdown=fd, efdown=efd, fup=fu, efup=efu,
                         fdown_over_fstar=fd / fstar(m, A) if m > 0 else np.nan))
stats['thresholds'] = thr_rows

# ---------------------------------------------------------------- fig1: flow curves
flows = {}
for fn in glob.glob('data/flow_m*.npz'):
    d = np.load(fn); flows[float(d['m'])] = d
if flows:
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.3))
    vv = np.geomspace(3e-3, 0.6, 120)
    lr_min = {}
    for m in MS:
        if m not in flows:
            continue
        d = flows[m]
        f = d['f'][:, 0]; v = d['v'].mean(1); ev = d['v'].std(1) / np.sqrt(d['v'].shape[1])
        ax[0].errorbar(f, v, ev, fmt='o', ms=2.5, color=col[m], label=f'm = {m:g}')
        fl, P = flow_curve(vv, m, 1.0, 1.0)
        ax[0].plot(fl, vv, ':', lw=0.7, color=col[m])
        if 1.0 in SC and m in SC[1.0]:
            vs_, f_, _ = SC[1.0][m]
            ax[0].plot(f_, vs_, '-', lw=0.9, color=col[m])
        i = np.argmin(fl); lr_min[m] = (fl[i], vv[i])
        ok = v > 2e-3
        pin = f[ok] - v[ok]
        shake = d['dv2'].mean(1)[ok] / v[ok]
        ax[1].loglog(v[ok], pin, 'o', ms=2.5, color=col[m])
        ax[1].loglog(v[ok], shake, 'x', ms=3, color=col[m])
        if 1.0 in SC and m in SC[1.0]:
            vs_, f_, _ = SC[1.0][m]
            ax[1].loglog(vs_, f_ - vs_, '-', lw=0.8, color=col[m])
        else:
            ax[1].loglog(vv, fl - vv, '-', lw=0.7, color=col[m])
        if m > 0:
            ax[1].loglog(vv, F2(1) / (2 * np.sqrt(m) * vv), ':', lw=0.6, color=col[m])
        # reduced: (f - eta v) v / P_rad versus v / v*
        if m >= 16:
            x = v[ok] / vstar(m, 1.0); y = pin * v[ok] / (F2(1) / (2 * np.sqrt(m)))
            ax[2].semilogx(x, y, 'o', ms=2.5, color=col[m])
            ax[2].semilogx(vv / vstar(m, 1.0), (fl - vv) * vv / (F2(1) / (2 * np.sqrt(m))), ':',
                           lw=0.7, color=col[m])
            if 1.0 in SC and m in SC[1.0]:
                vs_, f_, _ = SC[1.0][m]
                ax[2].semilogx(vs_ / vstar(m, 1.0), (f_ - vs_) * vs_ / (F2(1) / (2 * np.sqrt(m))), '-',
                               lw=0.9, color=col[m])
    ax[0].set(xlabel='$f$', ylabel=r'$\bar v$', xlim=(0, 0.46), ylim=(0, 0.42))
    ax[0].legend(ncol=2, frameon=False)
    ax[0].set_title('(a) flow curves: simulation (markers), closure (solid), 2nd order (dotted)', fontsize=7.5)
    ax[1].set(xlabel=r'$\bar v$', ylabel=r'$f-\eta\bar v$  (o),  $\eta\langle\delta\dot u^2\rangle/\bar v$  (x)',
              xlim=(3e-3, 0.6), ylim=(5e-3, 1))
    ax[1].set_title('(b) pinning force = shaking; lines: closure, dotted: $\\langle F^2\\rangle/2\\sqrt{\\kappa m}\\,v$', fontsize=7.5)
    ax[2].axhline(1, color='k', lw=0.5, ls='--')
    ax[2].set(xlabel=r'$\bar v / v_*$', ylabel=r'$P/P_{\rm rad}$,  $P=\bar v(f-\eta\bar v)$', ylim=(0, 1.3))
    ax[2].set_title('(c) absorbed power, $m\\geq16$ (solid: closure)', fontsize=8.5)
    fig.tight_layout(); save(fig, 'figures/flow_curves')
    stats['lr_minimum'] = {str(m): {'f': lr_min[m][0], 'v': lr_min[m][1]} for m in lr_min}
    # identity check
    rel = []
    for m, d in flows.items():
        v = d['v'].ravel(); f = d['f'].ravel(); s = d['dv2'].ravel()
        ok = v > 0.02
        rel += list(((f - v)[ok] - s[ok] / v[ok]) / (f - v)[ok])
    stats['energy_identity_rel_dev'] = {'mean': float(np.mean(rel)), 'rms': float(np.sqrt(np.mean(np.square(rel))))}
    # second-order accuracy at large v
    lr = {}
    for m, d in flows.items():
        f = d['f'][:, 0]; v = d['v'].mean(1)
        sel = v > 0.25
        if sel.any():
            pr = np.array([power(x, m, 1.0, 1.0) / x for x in v[sel]])
            lr[str(m)] = float(np.mean((f[sel] - v[sel]) / pr - 1))
    stats['second_order_rel_dev_v_gt_0.25'] = lr

    # ------------------------------------------------------------ fig2b: active fraction
    fig2b = {}
    for m in [0.0, 64.0, 256.0, 1024.0, 4096.0]:
        if m in flows:
            d = flows[m]
            fig2b[m] = (d['f'][:, 0], d['act'].mean(1), d['v'].mean(1))

# ---------------------------------------------------------------- fig2: kymographs + active fraction
ky = sorted(glob.glob('data/kymo_*.npz'))
if ky:
    fig = plt.figure(figsize=(11, 3.4))
    gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1.15])
    for k, fn in enumerate(ky[:3]):
        d = np.load(fn)
        ax = fig.add_subplot(gs[k])
        K = np.asarray(d['K']); vmax = np.percentile(np.abs(K), 99.5)
        ax.imshow(K.T, aspect='auto', origin='lower', cmap='magma', vmin=0, vmax=vmax,
                  extent=[0, float(d['T']), 0, int(d['L'])], interpolation='nearest')
        ax.set_title(f"({'abc'[k]}) m = {float(d['m']):g}, f = {float(d['f']):g}", fontsize=8.5)
        ax.set_xlabel('t'); ax.set_ylabel('i' if k == 0 else '')
    if flows:
        ax = fig.add_subplot(gs[3])
        for m, (f, act, v) in fig2b.items():
            ax.plot(f, act, 'o-', ms=2.5, lw=0.7, color=col[m], label=f'm = {m:g}')
        ax.set(xlabel='$f$', ylabel=r'active fraction $\phi$ ($|\dot u_i|>0.02$)', ylim=(0, 1.02))
        ax.legend(frameon=False); ax.set_title('(d) sliding fraction', fontsize=8.5)
    fig.tight_layout(); save(fig, 'figures/coexistence')

# ---------------------------------------------------------------- fig3: thresholds
if T:
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.3))
    mk = {256: 's', 1024: 'o', 4096: 'D'}
    for L in [256, 1024, 4096]:
        ms_ = [m for m in MS if (1.0, L, m) in T and m > 0]
        if not ms_:
            continue
        fd = [mean_err(T[(1.0, L, m)]['fdown']) for m in ms_]
        fu = [mean_err(T[(1.0, L, m)]['fup']) for m in ms_]
        ax[0].errorbar(ms_, [x[0] for x in fd], [x[1] for x in fd], fmt=mk[L], mfc='none', ms=4,
                       color='C0', label=f'$f_\\downarrow$, L = {L}')
        ax[0].errorbar(ms_, [x[0] for x in fu], [x[1] for x in fu], fmt=mk[L], ms=3.5,
                       color='C3', label=f'$f_\\uparrow$, L = {L}')
    mm = np.geomspace(0.2, 6000, 100)
    if 1.0 in fc:
        ax[0].axhline(fc[1.0], color='k', lw=0.6, ls='--', label='$f_c$ (m = 0)')
    ax[0].plot(mm, fstar(mm, 1.0), 'k:', lw=0.8, label=r'$f_*=(2\eta\langle F^2\rangle)^{1/2}(\kappa m)^{-1/4}$')
    if 1.0 in SC:
        msc = sorted(SC[1.0])
        ax[0].plot(msc, [sc_min(1.0, m)[0] for m in msc], 'k-', lw=1.0, label='$f_{\\min}$ (closure)')
    ax[0].set(xscale='log', yscale='log', xlabel='m', ylabel='force', ylim=(0.04, 0.45))
    ax[0].legend(frameon=False, fontsize=6.5, ncol=1)
    ax[0].set_title('(a) thresholds, A = 1', fontsize=8.5)
    # collapse
    for A, mkA in [(1.0, 'o'), (0.5, 's')]:
        if A not in fc:
            continue
        Ls = sorted({L for (a, L, m) in T if a == A})
        L = 1024 if 1024 in Ls else Ls[-1]
        ms_ = [m for m in MS if (A, L, m) in T and m > 0]
        M = np.array(ms_) / Lc(A) ** 2
        fd = np.array([mean_err(T[(A, L, m)]['fdown'])[0] for m in ms_])
        fu = np.array([mean_err(T[(A, L, m)]['fup'])[0] for m in ms_])
        ax[1].loglog(M, fd / fc[A], mkA, mfc='none', color='C0', label=f'$f_\\downarrow/f_c$, A = {A:g}, L = {L}')
        ax[1].loglog(M, fu / fc[A], mkA, color='C3', ms=3.5, label=f'$f_\\uparrow/f_c$, A = {A:g}, L = {L}')
        if A in SC:
            msc = sorted(SC[A])
            ax[1].loglog(np.array(msc) / Lc(A) ** 2, [sc_min(A, m)[0] / fc[A] for m in msc],
                         '-' if A == 1.0 else '--', color='k', lw=0.9, label=f'closure $f_{{\\min}}/f_c$, A = {A:g}')
    MM = np.geomspace(0.05, 3000, 100)
    if 1.0 in fc:
        ax[1].loglog(MM, fstar(MM * Lc(1.0) ** 2, 1.0) / fc[1.0], 'k:', lw=0.8, label='$f_*/f_c$ (A = 1)')
        if 0.5 in fc:
            ax[1].loglog(MM, fstar(MM * Lc(0.5) ** 2, 0.5) / fc[0.5], ':', color='0.5', lw=0.8, label='$f_*/f_c$ (A = 0.5)')
    ax[1].axhline(1, color='k', lw=0.5, ls='--')
    ax[1].set(xlabel=r'$M=\kappa m/(\eta L_c)^2$', ylabel='$f/f_c$', ylim=(0.15, 1.3))
    ax[1].legend(frameon=False, fontsize=6.5); ax[1].set_title('(b) reduced units', fontsize=8.5)
    # hysteresis vs L
    for m in [4.0, 16.0, 64.0, 256.0, 1024.0, 4096.0]:
        Ls = [L for L in [256, 1024, 4096] if (1.0, L, m) in T]
        if len(Ls) < 2:
            continue
        h0 = [np.nanmean(T[(1.0, L, 0.0)]['fup']) - np.nanmean(T[(1.0, L, 0.0)]['fdown']) for L in Ls]
        h = [np.nanmean(T[(1.0, L, m)]['fup']) - np.nanmean(T[(1.0, L, m)]['fdown']) - x for L, x in zip(Ls, h0)]
        ax[2].loglog(Ls, h, 'o-', ms=3.5, color=col[m], label=f'm = {m:g}')
        stats.setdefault('hysteresis', {})[str(m)] = dict(zip(map(str, Ls), map(float, h)))
    LL = np.array([256, 4096]); ax[2].loglog(LL, 0.07 * (LL / 256.) ** -0.5, 'k:', lw=0.7, label='$L^{-1/2}$')
    ax[2].set(xlabel='L', ylabel=r'$\Delta f_{\rm hyst}(m)-\Delta f_{\rm hyst}(0)$')
    ax[2].legend(frameon=False); ax[2].set_title('(c) hysteresis width', fontsize=8.5)
    fig.tight_layout(); save(fig, 'figures/ramp_thresholds')

# ---------------------------------------------------------------- fig4: well driving
wl = {}
for fn in glob.glob('data/well_m*.npz'):
    d = np.load(fn); wl[float(d['m'])] = d
if wl:
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.0), gridspec_kw={'width_ratios': [2.2, 1]})
    for m in sorted(wl):
        d = wl[m]
        t = d['t']; F = d['force'][0]
        ax[0].plot(t, F, lw=0.6, color=col[m], label=f'm = {m:g}')
        sel = t > 0.3 * t[-1]
        ax[1].hist(d['force'][:, sel].ravel(), bins=80, density=True, histtype='step', color=col[m])
        stats.setdefault('well', {})[str(m)] = {'mean_force': float(d['force'][:, sel].mean()),
                                               'min': float(d['force'][:, sel].min()),
                                               'max': float(d['force'][:, sel].max())}
    for x, ls in [(fc.get(1.0), '--')]:
        if x:
            ax[0].axhline(x, color='k', lw=0.5, ls=ls); ax[1].axvline(x, color='k', lw=0.5, ls=ls)
    ax[0].set(xlabel='t', ylabel=r'$k_0(w-\bar u)$', ylim=(0, 0.4)); ax[0].legend(frameon=False)
    ax[0].set_title('(a) force on the line, parabolic well ($k_0=0.01$, $\\dot w=10^{-3}$, L = 1024)', fontsize=8.5)
    ax[1].set(xlabel='force', ylabel='density'); ax[1].set_title('(b) force distribution (3 samples)', fontsize=8.5)
    fig.tight_layout(); save(fig, 'figures/well_driving')

stats['closure'] = {str(A): {str(m): dict(zip(['f_min', 'v_min', 'T_min'], map(float, sc_min(A, m))))
                             for m in sorted(SC[A])} for A in SC}
json.dump(stats, open('data/stats/sliding.json', 'w'), indent=1, default=float)
print(json.dumps({k: v for k, v in stats.items() if k != 'thresholds'}, indent=1, default=float))
for r in thr_rows:
    print(f"A={r['A']:g} L={r['L']:5d} m={r['m']:7g} M={r['M']:8.2f}  f_down={r['fdown']:.4f}+-{r['efdown']:.4f}"
          f"  f_up={r['fup']:.4f}+-{r['efup']:.4f}  f_down/f*={r['fdown_over_fstar']:.3f}")
