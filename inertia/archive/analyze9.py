"""Figure 15 and statistics9.json: spreading from a hot seed at the end of sliding (spread.py).

P(t): fraction of runs active at some time >= t.  N(t): mean number of active sites over
all runs.  R2(t): sum over runs of N R2 / sum of N (spread of the active sites about
the seed).  Local exponents are slopes over half-decade windows.  DP (1+1 d): delta =
0.1595, theta = 0.3137, 2/z = 1.265; hyperscaling theta + 2 delta = d/z.
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FS = [0.1975, 0.198, 0.1985]
fmin = 0.2008
DP = dict(delta=0.1595, theta=0.3137, twoz=1.265)
cols = plt.cm.viridis([0.1, 0.5, 0.85])


def obs(f):
    d = np.load(f'data/spread_f{f:g}.npz')
    N = d['N'].astype(float); R2 = d['R2']; ext = d['ext']; L = int(d['L'])
    t = np.arange(1, N.shape[1] + 1, dtype=float)
    act = N > 0
    last = np.where(act.any(1), N.shape[1] - 1 - np.argmax(act[:, ::-1], 1), -1)
    P = (last[:, None] >= np.arange(N.shape[1])[None, :]).mean(0)
    Nm = N.mean(0)
    with np.errstate(invalid='ignore', divide='ignore'):
        R = (N * R2).sum(0) / N.sum(0)
    return t, P, Nm, R, ext.max(), L, len(N)


def local(t, y, w=0.5):
    """Slope of log y vs log t over windows [t, t 10^w], at the window centres."""
    lt = np.log10(t); out_t, out_s = [], []
    for c in np.arange(0.5, lt[-1] - w + 1e-9, 0.1):
        sel = (lt >= c) & (lt <= c + w) & (y > 0) & np.isfinite(y)
        if sel.sum() > 3:
            out_t.append(10 ** (c + w / 2)); out_s.append(np.polyfit(lt[sel], np.log10(y[sel]), 1)[0])
    return np.array(out_t), np.array(out_s)


stats = {}
fig, ax = plt.subplots(1, 4, figsize=(17, 4))
for f, c in zip(FS, cols):
    t, P, Nm, R, emax, L, n = obs(f)
    lab = f'$f={f}$ ($={f / fmin:.4f}f_{{\\min}}$)'
    ax[0].loglog(t, P, color=c, label=lab)
    ax[1].loglog(t, Nm, color=c)
    ax[2].loglog(t, R, color=c)
    s = {}
    for key, y, sign in [('delta', P, -1), ('theta', Nm, 1), ('twoz', R, 1)]:
        lt, sl = local(t, y)
        s[key] = dict(t=lt.tolist(), slope=(sign * sl).tolist())
        win = (t >= 100) & (t <= t[-1])
        s[key + '_100_end'] = float(sign * np.polyfit(np.log(t[win]), np.log(y[win]), 1)[0])
    s['P_end'] = float(P[-1]); s['N_end'] = float(Nm[-1]); s['max_ext'] = float(emax); s['L'] = L
    s['runs'] = n
    s['hyperscaling_theta_plus_2delta'] = s['theta_100_end'] + 2 * s['delta_100_end']
    s['one_over_z'] = s['twoz_100_end'] / 2
    stats[str(f)] = s
    for key, mk in [('delta', 'o'), ('theta', 's'), ('twoz', '^')]:
        ax[3].semilogx(s[key]['t'], s[key]['slope'], mk + '-', ms=3, color=c, lw=0.8)
    print(f, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in s.items() if not isinstance(v, dict)})
t = np.array([10, 2000])
ax[0].loglog(t, 1.0 * (t / 10) ** -DP['delta'], 'k:', lw=0.8, label='DP slopes')
ax[1].loglog(t, 100 * (t / 10) ** DP['theta'], 'k:', lw=0.8)
ax[2].loglog(t, 2e3 * (t / 10) ** DP['twoz'], 'k:', lw=0.8)
for v in DP.values():
    ax[3].axhline(v, color='k', ls=':', lw=0.7)
ax[0].set_ylabel('$P_{\\rm surv}(t)$'); ax[1].set_ylabel('$\\langle N(t)\\rangle$'); ax[2].set_ylabel('$R^2(t)$')
ax[3].set_ylabel('local exponent'); ax[3].set_ylim(-0.1, 1.6)
ax[3].text(1.5e3, DP['delta'] + 0.03, '$\\delta$', fontsize=9); ax[3].text(1.5e3, DP['theta'] + 0.03, '$\\theta$', fontsize=9)
ax[3].text(1.5e3, DP['twoz'] + 0.03, '$2/z$', fontsize=9)
for a in ax:
    a.set_xlabel('$t/m$')
ax[0].legend(fontsize=7)
ax[0].set_title('(a) survival'); ax[1].set_title('(b) active sites'); ax[2].set_title('(c) spread')
ax[3].set_title('(d) local exponents (circles $\\delta$, squares $\\theta$, triangles $2/z$)', fontsize=9)
fig.tight_layout(); fig.savefig('archive/spreading_small_run.pdf'); fig.savefig('archive/spreading_small_run.png', dpi=110)
json.dump(stats, open('archive/spreading_small_run.json', 'w'), indent=1)
print('ok')
