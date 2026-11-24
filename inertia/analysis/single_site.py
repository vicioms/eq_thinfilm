"""Figure 16 and data/stats/single_site.json: single-site cavity theory against the simulated hot line.

Uses the polished pass (cavity_polish.py: ring of 128, 32 rings, central-difference eta_eff)
when present, else the average of the last three self-consistency passes.

Cavity: data/cavity_m{m}_v{v}.npz (cavity_runs.py); Gaussian DMFT with self-energy: data/gdmft2_m{m}_v{v}.npz (gdmft2.py).  Simulation: data/sq_diag_m{m}_f{f}.npz
(C(t), B(t)) and data/stats/roughness.json (T_lw, T_kin, S_xi(0), eta_eff needed).
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CASES = [(1024, 0.0891, 0.12), (1024, 0.0708, 0.11), (1024, 0.0512, 0.10), (256, 0.1169, 0.17)]
s5 = {(r['m'], r['f']): r for r in json.load(open('data/stats/roughness.json'))}
rows = []
fig, ax = plt.subplots(2, 4, figsize=(17, 7))
for j, (m, v, f) in enumerate(CASES):
    c = np.load(f'data/cavity_m{m:g}_v{v:g}.npz')
    d = np.load(f'data/sq_diag_m{m:g}_f{f:g}.npz')
    ref = s5[(m, f)]
    keys = list(c['hist_keys']); H = c['hist']
    last = dict(zip(keys, H[-1]))
    eta_eff = float(np.mean(H[-3:, keys.index('eta_eff')]))       # average the last three passes
    Sxi0 = float(np.mean(H[-3:, keys.index('Sxi0')]))
    fpred = float(np.mean(H[-3:, keys.index('f')]))
    import os
    pol = f'data/cavity_polish_m{m:g}_v{v:g}.npz'
    if os.path.exists(pol):   # polished pass: bigger ring, central-difference eta_eff
        c = np.load(pol)
        eta_eff, Sxi0, fpred = float(c['eta_eff']), float(c['Sxi0']), float(c['f'])
        H = np.array([[float(c['Tkin'])]] * 3); keys = ['Tkin']
    row = dict(m=m, v=v, f_sim=f, f_cav=fpred, Sxi0_sim=ref['S_xi0'], Sxi0_cav=Sxi0,
               Sxi0_bare=1 / (3 * v), eta_eff_sim=ref['eta_eff_needed'], eta_eff_cav=eta_eff,
               eta_eff_closure=ref['eta_eff_closure'], T_lw_sim=ref['T_lw'],
               T_lw_cav=Sxi0 / (2 * eta_eff), T_lw_closure=ref['T_lw_pred'], T_lw_bare=ref['T_lw_bare'],
               Tkin_sim=ref['Tkin'], Tkin_cav=float(np.mean(H[-3:, keys.index('Tkin')])),
               passes=int(len(H)))
    g = np.load(f'data/gdmft2_m{m:g}_v{v:g}.npz')     # Gaussian DMFT with Sigma(w) (gdmft2.py)
    row.update(f_gd=float(g['f']), Sxi0_gd=float(g['Sxi0']), eta_eff_gd=float(g['eta_eff']),
               T_lw_gd=float(g['T_lw']), Tkin_gd=float(g['Tkin']))
    rows.append(row)
    a = ax[0, j]
    a.plot(g['t'] * v, g['C'], '--', color='C2', lw=1, label='Gaussian DMFT, $\\Sigma(\\omega)$')
    a.plot(d['t'] * v, d['C'], 'k-', lw=1.2, label='simulation')
    a.plot(c['t'] * v, c['C'], '-', color='C3', lw=1, label='cavity')
    a.plot(d['t'] * v, d['Cb'], ':', color='0.5', label='bare $\\Delta(vt)$')
    a.plot(d['t'] * v, d['CG'], '--', color='C0', lw=0.8, label='Gaussian, measured $B$')
    a.set_xlim(0, 6); a.set_xlabel('$vt$'); a.set_title(f'$m={m}$, $f={f}$ ($v={v}$)', fontsize=9)
    b = ax[1, j]
    b.loglog(d['t'][1:], d['B'][1:], 'k-', lw=1.2, label='simulation')
    b.loglog(c['tB'][1:], c['B'][1:], '-', color='C3', lw=1, label='cavity')
    b.loglog(g['t'][1:], g['B'][1:], '--', color='C2', lw=1, label='Gaussian DMFT, $\\Sigma(\\omega)$')
    b.set_xlim(1, 4000)
    b.set_xlabel('$t$')
ax[0, 0].set_ylabel('$C(t)$'); ax[1, 0].set_ylabel('$B(t)$')
ax[0, 0].legend(fontsize=7); ax[1, 0].legend(fontsize=7)
fig.tight_layout(); fig.savefig('figures/single_site.pdf'); fig.savefig('figures/single_site.png', dpi=110)
json.dump(rows, open('data/stats/single_site.json', 'w'), indent=1)
hdr = ['m', 'v', 'f_sim', 'f_cav', 'f_gd', 'Sxi0_sim', 'Sxi0_cav', 'Sxi0_gd', 'Sxi0_bare', 'eta_eff_sim',
       'eta_eff_cav', 'eta_eff_gd', 'T_lw_sim', 'T_lw_cav', 'T_lw_gd', 'Tkin_sim', 'Tkin_cav', 'Tkin_gd']
print(' '.join(f'{h:>11}' for h in hdr))
for r in rows:
    print(' '.join(f'{r[h]:11.4g}' for h in hdr))
