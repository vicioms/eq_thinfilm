"""figures/thresholds.pdf and data/stats/thresholds.json: where each analysis threshold in
src/config.py comes from, and how much the results move when it is varied.

(a) single-site speed distribution (runs/thresholds.py vhist): V_ACT sits in the gap between
    pinned and moving sites;
(b) active fraction over that of Gaussian velocities with the measured mean and variance
    (data/flow_m*.npz): uniform sliding sits at 1, partial sliding below; RHO_TOL;
(c) time gaps between avalanche fragments (data/aval_m*.npz): ringing peak, trough, waiting
    times; MERGE_GAP = 4m;
(d) avalanche sizes against the largest overdamped event (slip_cut);
(e) spreading exponents at f = 0.1981 against the activity threshold, and for other kicks and
    substrates (runs/thresholds.py spread);
(f) avalanche sizes for detection thresholds C_TH = 1.5, 3, 10 (runs/thresholds.py aval).
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
from config import V_ACT, V_REST, RHO_TOL, merge_gap, slip_cut
import json
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import erf
from aval_merge import merge

rng = np.random.default_rng(3)
stats = {}
fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
ax = ax.ravel()

# (a) single-site speeds ---------------------------------------------------------------------
if os.path.exists('data/check_vhist.npz'):
    d = np.load('data/check_vhist.npz')
    e = d['edges']; c = 0.5 * (e[1:] + e[:-1])
    st = {}
    for key, lab in [('m64_f0.198', '$m=64$, $f=0.198$ (partial)'), ('m64_f0.25', '$m=64$, $f=0.25$ (uniform)'),
                     ('m1024_f0.12', '$m=1024$, $f=0.12$ (uniform)'), ('m1024_f0.093', '$m=1024$, $f=0.093$ (stopped)')]:
        h = d[key].astype(float)
        if h.sum() == 0:
            continue
        p = h / h.sum() / np.diff(e)
        ax[0].semilogy(10 ** c, np.where(p > 0, p, np.nan), label=lab)
        rho = d['rho_' + key].mean(0)
        st[key] = dict(rho=dict(zip([f'{v:g}' for v in d['vacts']], rho.round(4).tolist())))
        if 'partial' in lab:
            # moving sites form a peak; pinned sites ring down and fill a flat tail in log|v|
            hs = np.convolve(h, np.ones(5) / 5, 'same')
            imov = int(np.argmax(hs))
            st[key].update(moving_peak=float(10 ** c[imov]),
                           fraction_per_decade_below_1em3=float(h[c < -3].sum() / h.sum() / 3))
    ax[0].set_xscale('log'); ax[0].axvline(V_ACT, color='k', ls=':', lw=0.8)
    ax[0].set_xlabel('$|\\dot u_i|$'); ax[0].set_ylabel('density in $\\log_{10}|\\dot u_i|$')
    ax[0].legend(fontsize=7); ax[0].set_title('(a) single-site speeds; dotted: $V_{\\rm act}$')
    stats['speeds'] = st

# (b) uniform versus partial ------------------------------------------------------------------
fmin = {64: 0.2008, 256: 0.1396, 1024: 0.0943, 4096: 0.063}
st = {}
for m in [64, 256, 1024, 4096]:
    d = np.load(f'data/flow_m{m}.npz')
    f, v, sd, act = d['f'][:, 0], d['v'].mean(1), np.sqrt(d['dv2'].mean(1)), d['act'].mean(1)
    rho_g = 1 - 0.5 * (erf((V_ACT - v) / (sd * np.sqrt(2))) - erf((-V_ACT - v) / (sd * np.sqrt(2))))
    ok = v > V_REST                                   # lines that stopped are left out
    ratio = act[ok] / rho_g[ok]
    ax[1].plot(f[ok] / fmin[m], ratio, 'o-', ms=3, label=f'$m={m}$')
    uni = ratio >= 1 - RHO_TOL
    st[m] = dict(min_uniform=float(ratio[uni].min()) if uni.any() else None,
                 max_partial=float(ratio[~uni].max()) if (~uni).any() else None,
                 n_uniform=int(uni.sum()), n_partial=int((~uni).sum()),
                 old_cut_disagrees=int(((act[ok] > 0.95) != uni).sum()))
ax[1].axhline(1 - RHO_TOL, color='k', ls=':', lw=0.8)
ax[1].set_xlabel('$f/f_{\\min}$'); ax[1].set_ylabel('$\\rho/\\rho_G$')
ax[1].set_xlim(0.8, 2); ax[1].set_ylim(0, 1.05); ax[1].legend(fontsize=7)
ax[1].set_title('(b) activity over that of Gaussian velocities')
stats['uniform_vs_partial'] = st


# (c) fragment gaps ---------------------------------------------------------------------------
def gaps(r):
    wd = float(r['wdot']); g = []
    for k in np.unique(r['sample']):
        sel = np.where(r['sample'] == k)[0]; sel = sel[np.argsort(r['w_on'][sel])]
        ton = r['w_on'][sel] / wd; tend = ton + r['T'][sel]
        g += list(ton[1:] - tend[:-1])
    g = np.array(g); return g[g > 0]


st = {}
e = np.arange(-0.0625, 6.5, 0.125); c = 0.5 * (e[1:] + e[:-1])
for m in [1, 4, 16, 64, 256, 1024]:
    r = np.load(f'data/aval_m{m}.npz'); g = gaps(r)
    h = np.histogram(np.log10(g), e)[0]
    ax[2].plot(10 ** c / np.sqrt(m), h / h.sum(), '-', lw=1, label=f'$m={m}$')
    ax[2].plot(merge_gap(m) / np.sqrt(m), 0.002, 'v', color=ax[2].lines[-1].get_color())
    # merged fraction for gaps from the end of the ringing peak to the start of the waiting times
    st[m] = {f'{x:g}m': float(np.mean(g < x * max(m, 1))) for x in [2, 4, 8, 16]}
ax[2].set_xscale('log'); ax[2].set_xlabel('gap$/\\sqrt{m}$'); ax[2].set_ylabel('fraction')
ax[2].legend(fontsize=7); ax[2].set_title('(c) gaps between fragments; $\\nabla$: merge gap $4m$')
stats['merge_gap_merged_fraction'] = st

# (d) slips -------------------------------------------------------------------------------------
L = 512
s0 = merge(np.load('data/aval_m0.npz'))[0]; cut = slip_cut(s0)
bins = np.geomspace(0.3, 1e5, 25); bc = np.sqrt(bins[1:] * bins[:-1])
st = dict(cut_over_L=cut / L)
for m in [0, 16, 64, 256, 1024]:
    r = np.load(f'data/aval_m{m}.npz'); s, T, fo = merge(r); k0 = float(r['k0'])
    hh = np.histogram(s, bins)[0]
    ax[3].loglog(bc, np.where(hh > 0, hh / len(s), np.nan), 'o-', ms=3, label=f'$m={m}$')
    if m:
        row = {}
        for name, x in [('4L', 4 * L), ('max_m0', cut), ('16L', 16 * L)]:
            big = s > x
            if big.any():
                fs = fo[big] - k0 * s[big] / L
                row[name] = dict(n=int(big.sum()), share=float(s[big].sum() / s.sum()),
                                 f_k=float(((fo[big] + fs) / 2).mean()))
        st[m] = row
ax[3].axvline(cut, color='k', ls=':', lw=0.8)
ax[3].set_xlabel('$s$'); ax[3].set_ylabel('fraction'); ax[3].legend(fontsize=7)
ax[3].set_title('(d) sizes; dotted: largest event at $m=0$')
stats['slip_definition'] = st


# (e) spreading -----------------------------------------------------------------------------------
def curves(N, R2):
    act = N > 0
    last = np.where(act.any(1), N.shape[1] - 1 - np.argmax(act[:, ::-1], 1), -1)
    P = (last[:, None] >= np.arange(N.shape[1])[None, :]).mean(0)
    with np.errstate(invalid='ignore', divide='ignore'):
        R = (N * R2).sum(0) / N.sum(0)
    return P, N.mean(0), R


def slopes(P, Nm, R, lo=100, hi=1000):
    t = np.arange(1, len(P) + 1, dtype=float)
    w = (t >= lo) & (t <= hi) & (P > 0) & np.isfinite(R) & (Nm > 0)
    f = lambda y: np.polyfit(np.log(t[w]), np.log(y[w]), 1)[0]
    return np.array([-f(P), f(Nm), f(R)])


st = {}
tags = [t for t in ['base', 'narrow', 'wide', 'sub0.1908'] if os.path.exists(f'data/check_spread_{t}.npz')]
for j, tag in enumerate(tags):
    d = np.load(f'data/check_spread_{tag}.npz')
    bad = np.unique(d['seeds'][d['pre_act'] > 0]); keep = ~np.isin(d['seeds'], bad)
    seeds = d['seeds'][keep]; useq = np.unique(seeds)
    row = {}
    for iv, vt in enumerate(d['vacts']):
        N, R2 = d['N'][iv][keep].astype(float), d['R2'][iv][keep]
        e0 = slopes(*curves(N, R2))
        b = []
        for _ in range(200):
            pick = rng.choice(useq, len(useq))
            idx = np.concatenate([np.where(seeds == q)[0] for q in pick])
            b.append(slopes(*curves(N[idx], R2[idx])))
        sd = np.array(b).std(0)
        row[f'{vt:g}'] = {k: [float(x), float(s)] for k, x, s in zip(['delta', 'theta', 'twoz'], e0, sd)}
    st[tag] = dict(runs=int(keep.sum()), w=int(d['w']), v_kick=float(d['v_kick']), exponents=row)
    if tag == 'base':
        vts = d['vacts']
        for k, mk, lab in [('delta', 'o', '$\\delta$'), ('theta', 's', '$\\theta$'), ('twoz', '^', '$2/z$')]:
            y = np.array([row[f'{v:g}'][k] for v in vts])
            ax[4].errorbar(vts, y[:, 0], y[:, 1], fmt=mk + '-', ms=4, capsize=2, label=lab + ' (base kick)')
    else:
        y = row[f'{V_ACT:g}']
        for k, mk in [('delta', 'o'), ('theta', 's'), ('twoz', '^')]:
            ax[4].errorbar(V_ACT * (1.25 + 0.2 * j), y[k][0], y[k][1], fmt=mk, mfc='none', ms=5, capsize=2,
                           color=f'C{j + 2}', label=tag if k == 'delta' else None)
if tags:
    ax[4].set_xscale('log'); ax[4].axvline(V_ACT, color='k', ls=':', lw=0.8)
    ax[4].set_xlabel('activity threshold'); ax[4].set_ylabel('exponent, $t/m\\in[100,1000]$')
    ax[4].legend(fontsize=6.5, ncol=2); ax[4].set_title('(e) spreading at $f=0.1981$, $m=64$')
stats['spreading'] = st

# (f) detection threshold -------------------------------------------------------------------------
st = {}
for j, m in enumerate([4, 64, 1024]):
    for cth, ls in [(1.5, ':'), (3.0, '-'), (10.0, '--')]:
        fn = f'data/aval_m{m}.npz' if cth == 3.0 else f'data/check_aval_c{cth:g}_m{m}.npz'
        if not os.path.exists(fn):
            continue
        s = merge(np.load(fn))[0]
        hh = np.histogram(s, bins)[0]
        ax[5].loglog(bc, np.where(hh > 0, hh / len(s), np.nan), ls, color=f'C{j}',
                     label=f'$m={m}$, $c={cth:g}$')
        st[f'm{m}_c{cth:g}'] = dict(n=int(len(s)), median=float(np.median(s)),
                                   p90=float(np.percentile(s, 90)), slip_share=float(s[s > cut].sum() / s.sum()))
ax[5].set_xlabel('$s$'); ax[5].set_ylabel('fraction'); ax[5].legend(fontsize=6.5)
ax[5].set_title('(f) detection $\\bar v>c\\,\\dot w$')
stats['detection'] = st

fig.tight_layout(); fig.savefig('figures/thresholds.pdf'); fig.savefig('figures/thresholds.png', dpi=110)
json.dump(stats, open('data/stats/thresholds.json', 'w'), indent=1, default=float)
print(json.dumps({k: v for k, v in stats.items() if k != 'spreading'}, indent=1, default=float)[:3000])
for tag, r in stats.get('spreading', {}).items():
    print(tag, r['runs'], {v: {k: f'{x[0]:.3f}+-{x[1]:.3f}' for k, x in e.items()} for v, e in r['exponents'].items()})
