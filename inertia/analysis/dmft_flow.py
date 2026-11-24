"""Figure 17: flow curves of the hot line: simulation, closure, Gaussian DMFT (constant eta_eff,
gdmft_flow.py) and self-energy Gaussian DMFT (gdmft2_flow.py), m = 64, 256, 1024."""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
from config import V_ACT, V_REST, RHO_TOL
import numpy as np
from scipy.special import erf
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sc = np.load('data/sc_flow.npz')
g1 = np.load('data/gdmft_flow.npz'); g2 = np.load('data/gdmft2_flow.npz')
finf = {64: 0.198, 256: 0.137, 1024: 0.092}        # measured long-time thresholds
fmin = {64: 0.2008, 256: 0.1396, 1024: 0.0943}     # closure minima
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
for a, m in zip(ax, [64, 256, 1024]):
    d = np.load(f'data/flow_m{m}.npz')
    fs, vs = d['f'][:, 0], d['v'].mean(1)
    act = d['act'].mean(1); sd = np.sqrt(d['dv2'].mean(1))
    # uniform sliding: active fraction of Gaussian velocities with the measured mean and variance
    rho_g = 1 - 0.5 * (erf((V_ACT - vs) / (sd * np.sqrt(2))) - erf((-V_ACT - vs) / (sd * np.sqrt(2))))
    uni = act >= (1 - RHO_TOL) * rho_g
    ok = vs > V_REST                                  # lines that stopped are not plotted
    a.plot(vs[ok & uni], fs[ok & uni], 'ko', ms=4, label='simulation (uniform)')
    a.plot(vs[ok & ~uni], fs[ok & ~uni], 'o', mfc='none', mec='k', ms=4, label='simulation (patchy)')
    a.plot(sc['vs'], sc[f'm{m}_f'], '-', color='C0', label='closure')
    r1 = g1[f'm{m}']; s1 = r1[:, -1] < 0.5
    a.plot(r1[s1, 0], r1[s1, 1], '--', color='0.5', label='Gaussian DMFT, constant $\\eta_{\\rm eff}$')
    r2 = g2[f'm{m}']
    a.plot(r2[:, 0], r2[:, 1], '-', color='C3', lw=1.5, label='Gaussian DMFT, $\\Sigma(\\omega)$')
    a.axhline(finf[m], color='C2', ls=':', label='measured $f_\\infty$')
    a.set_xscale('log'); a.set_xlim(4e-3, 0.5)
    lo = min(r2[:, 1].min(), finf[m]) * 0.9
    a.set_ylim(lo, finf[m] * 2.2)
    a.set_xlabel('$\\bar v$'); a.set_title(f'$m={m}$')
ax[0].set_ylabel('$f$'); ax[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig('figures/dmft_flow.pdf'); fig.savefig('figures/dmft_flow.png', dpi=110)
for m in [64, 256, 1024]:
    r2 = g2[f'm{m}']; i = np.argmin(r2[:, 1])
    pb = np.max(np.abs(r2[:, 5] / r2[:, 6] - 1))
    print(m, 'Sigma-DMFT min f', round(r2[i, 1], 4), 'at v', r2[i, 0], '(edge of scan)' if i in (0, len(r2) - 1) else '',
          '| closure f_min', fmin[m], '| f_inf', finf[m], '| max power-balance error', round(pb, 4))
