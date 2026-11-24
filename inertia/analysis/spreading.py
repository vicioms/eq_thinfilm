"""Figure 15 (final) and data/stats/spreading.json: spreading exponents from the large run.

Large run: L = 32768, 32 substrates x 16 kicks = 512 runs at f = 0.1981 and 0.1983, 5000 m
(SPREAD_L=32768 SPREAD_NSUB=32 SPREAD_NK=16 python3 spread.py run 0.1981 0.1983 --T 5000).
Bracketing run (L = 16384, 96 runs, 2000 m) at 0.1975, 0.198, 0.1985 is overlaid.

Exponents: bootstrap over runs (resampling whole substrates, since the copies of one
substrate share disorder) in windows [100,1000], [300,3000], [1000,5000] m.  The two forces
bracket f_inf; 2/z is the least sensitive to the distance from it, so its value at both
forces is the main result.  DP: delta 0.1595, theta 0.3137, 2/z 1.265.
Conserved DP / qEW: z = 1.39-1.43, 2/z = 1.40-1.44.
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DP = dict(delta=0.1595, theta=0.3137, twoz=1.265)
CDP_twoz = (1.40, 1.44)
fmin = 0.2008
rng = np.random.default_rng(10)
WINS = [(100, 1000), (300, 3000), (1000, 5000)]


def curves(N, R2):
    act = N > 0
    last = np.where(act.any(1), N.shape[1] - 1 - np.argmax(act[:, ::-1], 1), -1)
    P = (last[:, None] >= np.arange(N.shape[1])[None, :]).mean(0)
    Nm = N.mean(0)
    with np.errstate(invalid='ignore', divide='ignore'):
        R = (N * R2).sum(0) / N.sum(0)
    return P, Nm, R


def slopes(t, P, Nm, R, lo, hi):
    w = (t >= lo) & (t <= hi) & (P > 0) & np.isfinite(R)
    f = lambda y: np.polyfit(np.log(t[w]), np.log(y[w]), 1)[0]
    return -f(P), f(Nm), f(R)


def load(f, tag):
    d = np.load(f'data/spread_f{f:g}{tag}.npz')
    # drop every copy of a substrate that was still active before the kick (raising f from the
    # dead force to the test force can trigger a spontaneous event that would mix with the seed)
    bad = np.unique(d['seeds'][d['pre_act'] > 0])
    keep = ~np.isin(d['seeds'], bad)
    if bad.size:
        print(f'  f {f}{tag}: dropping {(~keep).sum()} runs on {bad.size} substrates active before the kick')
    return d['N'][keep].astype(float), d['R2'][keep], d['seeds'][keep], d['ext'][keep].max(), int(d['L'])


def local(t, y, w=0.5):
    lt = np.log10(t); out_t, out_s = [], []
    for c in np.arange(0.5, lt[-1] - w + 1e-9, 0.1):
        sel = (lt >= c) & (lt <= c + w) & (y > 0) & np.isfinite(y)
        if sel.sum() > 3:
            out_t.append(10 ** (c + w / 2)); out_s.append(np.polyfit(lt[sel], np.log10(y[sel]), 1)[0])
    return np.array(out_t), np.array(out_s)


stats = {}
big = [(0.1981, '_L32768'), (0.1983, '_L32768')]
small = [(0.1975, ''), (0.198, ''), (0.1985, '')]
fig, ax = plt.subplots(1, 4, figsize=(17, 4))
cols = dict(zip([0.1975, 0.198, 0.1981, 0.1983, 0.1985], plt.cm.viridis([0.05, 0.3, 0.5, 0.7, 0.9])))
for f, tag in small + big:
    N, R2, seeds, emax, L = load(f, tag)
    t = np.arange(1, N.shape[1] + 1, dtype=float)
    P, Nm, R = curves(N, R2)
    ls = '-' if tag else ':'
    lab = f'$f={f}$ ($={f / fmin:.4f}f_{{\\min}}$), {len(N)} runs'
    ax[0].loglog(t, P, ls, color=cols[f], label=lab)
    ax[1].loglog(t, Nm, ls, color=cols[f]); ax[2].loglog(t, R, ls, color=cols[f])
    if not tag:
        continue
    s = dict(runs=int(len(N)), L=L, max_ext=float(emax), P_end=float(P[-1]))
    useq = np.unique(seeds)
    for lo, hi in WINS:
        if hi > t[-1]:
            continue
        e = slopes(t, P, Nm, R, lo, hi)
        b = []
        for _ in range(300):
            pick = rng.choice(useq, len(useq))
            idx = np.concatenate([np.where(seeds == q)[0] for q in pick])
            b.append(slopes(t, *curves(N[idx], R2[idx]), lo, hi))
        b = np.array(b)
        s[f'{lo}-{hi}'] = {k: [float(v), float(sd)] for k, v, sd in zip(['delta', 'theta', 'twoz'], e, b.std(0))}
    for key, y, sign, mk in [('delta', P, -1, 'o'), ('theta', Nm, 1, 's'), ('twoz', R, 1, '^')]:
        lt, sl = local(t, y)
        s['local_' + key] = [lt.tolist(), (sign * sl).tolist()]
        ax[3].semilogx(lt, sign * sl, mk + '-', ms=3, lw=0.8, color=cols[f])
    stats[str(f)] = s
    print(f, {k: v for k, v in s.items() if not k.startswith('local')})
t = np.array([10, 5000])
ax[0].loglog(t, 1.0 * (t / 10) ** -DP['delta'], 'k:', lw=0.8, label='DP slopes')
ax[1].loglog(t, 100 * (t / 10) ** DP['theta'], 'k:', lw=0.8)
ax[2].loglog(t, 2e3 * (t / 10) ** DP['twoz'], 'k:', lw=0.8)
ax[2].loglog(t, 2e3 * (t / 10) ** 1.42, 'k--', lw=0.8, label='$t^{1.42}$ (conserved DP)')
for v in DP.values():
    ax[3].axhline(v, color='k', ls=':', lw=0.7)
ax[3].axhspan(*CDP_twoz, color='0.85', zorder=0)
ax[0].set_ylabel('$P_{\\rm surv}(t)$'); ax[1].set_ylabel('$\\langle N(t)\\rangle$'); ax[2].set_ylabel('$R^2(t)$')
ax[3].set_ylabel('local exponent'); ax[3].set_ylim(-0.1, 1.7)
for a in ax:
    a.set_xlabel('$t/m$')
ax[0].legend(fontsize=6); ax[2].legend(fontsize=7)
ax[0].set_title('(a) survival'); ax[1].set_title('(b) active sites'); ax[2].set_title('(c) spread')
ax[3].set_title('(d) local $\\delta$ (o), $\\theta$ (s), $2/z$ ($\\triangle$), large run', fontsize=9)
fig.tight_layout(); fig.savefig('figures/spreading.pdf'); fig.savefig('figures/spreading.png', dpi=110)
json.dump(stats, open('data/stats/spreading.json', 'w'), indent=1)
print('ok')
