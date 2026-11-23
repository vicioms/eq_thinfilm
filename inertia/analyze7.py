"""Figure 13 and statistics7.json: avalanche exponents from two well stiffnesses (aval.py scan, tau).

Small well: k0 = 3e-3 (xi ~ 18), L = 512   -> data/aval_m{m}.npz
Soft well:  k0 = 5e-4 (xi ~ 45), L = 2048  -> data/avalT_m{m}.npz

tau from three estimators (all at m = 0):
  * pure power law (maximum likelihood) on [3, 300], well below the cutoff;
  * moment scaling: <s> ~ S_m^(2 - tau), S_m = <s^2>/2<s>, between the two wells, with <s> either
    from the detected events or from the sum rule <s> = L dW / n (includes sub-threshold creep);
  * size-duration scaling s ~ T^gamma, gamma = (1 + zeta)/z, compared with qEW 1.57.
Also zeta_eff from S_m ~ k0^-(1+zeta)/2, and the inertial shift of the mean force.
"""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from aval_merge import merge
from tau_fit import fit, powerfit

rng = np.random.default_rng(7)


def load(m, soft):
    r = dict(np.load(f'data/avalT_m{m:g}.npz' if soft else f'data/aval_m{m:g}.npz'))
    s, T, fo = merge(r)
    ok = s > 0
    n_samp = len(np.unique(r['sample']))
    sum_rule = int(r['L']) * (float(r['W_total']) - float(r['W_skip'])) * n_samp / ok.sum()
    return r, s[ok], T[ok], r['sample'], sum_rule


def pl_ml(x, lo, hi):
    x = x[(x >= lo) & (x <= hi)]
    taus = np.linspace(0.5005, 2.5005, 2000)
    ll = [-(t * np.log(x)).sum() - len(x) * np.log((hi ** (1 - t) - lo ** (1 - t)) / (1 - t)) for t in taus]
    return float(taus[int(np.argmax(ll))]), int(len(x))


def pl_err(x, lo, hi, nb=200):
    return float(np.std([pl_ml(rng.choice(x, len(x)), lo, hi)[0] for _ in range(nb)]))


def moments(s):
    return s.mean(), (s ** 2).mean() / (2 * s.mean())


def boot_moments(s, nb=400):
    b = []
    for _ in range(nb):
        x = rng.choice(s, len(s)); b.append(moments(x))
    return np.log(np.array(b))


def sT(s, T):
    lb = np.geomspace(max(T.min(), 1), T.max(), 20)
    idx = np.digitize(T, lb)
    keep = [i for i in range(1, len(lb)) if (idx == i).sum() > 3]
    return (np.array([T[idx == i].mean() for i in keep]), np.array([s[idx == i].mean() for i in keep]))


k_small, k_soft = 3e-3, 5e-4
stats = dict(soft={}, pair={})
ms_soft = [0, 1, 2, 4, 8, 16]
cols = plt.cm.viridis(np.linspace(0, 0.9, len(ms_soft)))
cmap = dict(zip(ms_soft, cols))
fig, ax = plt.subplots(1, 4, figsize=(17, 4))

# (a) size distributions in the soft well
bins = np.geomspace(0.3, 2e5, 26); bc = np.sqrt(bins[1:] * bins[:-1])
for m in ms_soft:
    r, s, T, _, _ = load(m, True)
    h, _ = np.histogram(s, bins); dens = h / np.diff(bins) / len(s); ok = h > 0
    ax[0].loglog(bc[ok], dens[ok] * bc[ok] ** 1.11, 'o-', ms=3, color=cmap[m], label=f'$m={m}$')
    tau_pl, n_pl = pl_ml(s, 3, 300)
    F = fit(s, 3, nboot=100)
    Tm, sm = sT(s, T)
    stats['soft'][m] = dict(n=int(len(s)), mean_force=float(r['force'].mean()),
                            force_sd=float(r['force'].std()),
                            tau_pl_3_300=tau_pl, tau_pl_err=pl_err(s, 3, 300), n_pl=n_pl,
                            tau_cut_3=float(F['tau']), tau_cut_err=F['err'], s_c=float(F['s_c']),
                            gamma_sT=float(powerfit(Tm, sm, 10, 1000)),
                            alpha_pl_5_300=pl_ml(T, 5, 300)[0])
ax[0].set_xlabel('$s$'); ax[0].set_ylabel('$s^{1.11}P(s)$')
ax[0].legend(fontsize=7, loc='lower left')
ax[0].set_title('(a) soft well $\\xi\\simeq45$: flat = qEW $\\tau$')

# (b) <s> against S_m between the two wells
for m in [0, 1, 4, 16]:
    pts = []
    for soft in [False, True]:
        r, s, T, smp, sr = load(m, soft)
        pts.append((moments(s), boot_moments(s), sr))
    (a1, S1), B1, sr1 = pts[0]; (a2, S2), B2, sr2 = pts[1]
    sl = np.log(a2 / a1) / np.log(S2 / S1)
    sl_b = (B2[:, 0] - B1[:, 0]) / (B2[:, 1] - B1[:, 1])
    sl_sr = np.log(sr2 / sr1) / np.log(S2 / S1)
    zeta = 2 * np.log(S2 / S1) / np.log(k_small / k_soft) - 1
    stats['pair'][m] = dict(S_m=[float(S1), float(S2)], mean_s=[float(a1), float(a2)],
                            mean_s_sum_rule=[float(sr1), float(sr2)],
                            tau_moment=float(2 - sl), tau_moment_err=float(np.std(sl_b)),
                            tau_moment_sum_rule=float(2 - sl_sr), zeta_eff=float(zeta))
    c = cmap[m]
    ax[1].loglog([S1, S2], [a1, a2], 'o-', color=c, label=f'$m={m}$: $\\tau={2 - sl:.2f}\\pm{np.std(sl_b):.2f}$')
    ax[1].loglog([S1, S2], [sr1, sr2], 's--', color=c, mfc='none', lw=0.8)
x = np.array([400, 1.2e4])
ax[1].loglog(x, 230 * (x / 503) ** 0.89, 'k:', lw=0.8, label='slope $2-1.11$')
ax[1].set_xlabel('$S_m=\\langle s^2\\rangle/2\\langle s\\rangle$'); ax[1].set_ylabel('$\\langle s\\rangle$')
ax[1].legend(fontsize=7); ax[1].set_title('(b) $\\xi\\simeq18\\to45$; open: sum rule')

# (c) size against duration
for m in [0, 4, 16]:
    r, s, T, _, _ = load(m, True)
    Tm, sm = sT(s, T)
    ax[2].loglog(Tm, sm, 'o-', ms=3, color=cmap[m], label=f'$m={m}$')
x = np.array([10, 3000])
ax[2].loglog(x, 0.15 * x ** 1.57, 'k:', lw=0.8, label='$T^{1.57}$ (qEW)')
ax[2].loglog(x, 0.012 * x ** 2.0, 'k--', lw=0.8, label='$T^{2}$')
ax[2].set_xlabel('$T$'); ax[2].set_ylabel('$\\langle s\\,|\\,T\\rangle$'); ax[2].legend(fontsize=7)
ax[2].set_title('(c) size against duration, soft well')

# (d) inertial shift of the mean force, in units of the m = 0 force spread
for soft, mk, lab, msl in [(False, 'o-', '$\\xi\\simeq18$', [0, 1, 4, 16]), (True, 's-', '$\\xi\\simeq45$', ms_soft)]:
    f0 = load(0, soft)[0]['force']; sd = f0.std(); mu = f0.mean()
    y = [(mu - load(m, soft)[0]['force'].mean()) / sd for m in msl]
    ax[3].plot(msl, y, mk, label=lab)
    stats['shift_' + ('soft' if soft else 'small')] = dict(zip([str(m) for m in msl], map(float, y)))
ax[3].set_xscale('symlog', linthresh=1); ax[3].axhline(0, color='k', lw=0.6)
ax[3].set_xlabel('$m$'); ax[3].set_ylabel('$(\\langle f\\rangle_0-\\langle f\\rangle_m)/\\sigma_f(0)$')
ax[3].legend(fontsize=8); ax[3].set_title('(d) inertial force shift / critical window')

fig.tight_layout(); fig.savefig('fig13_tau.pdf'); fig.savefig('fig13_tau.png', dpi=110)
json.dump(stats, open('statistics7.json', 'w'), indent=1)
for m, d in stats['soft'].items():
    print(m, {k: round(v, 3) if isinstance(v, float) else v for k, v in d.items()})
for m, d in stats['pair'].items():
    print(m, {k: (np.round(v, 3).tolist() if isinstance(v, list) else round(v, 3)) for k, v in d.items()})
print(stats['shift_small'], stats['shift_soft'])
