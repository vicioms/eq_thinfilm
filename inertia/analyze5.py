"""Figure 11 and statistics5.json: roughness of the hot sliding line (sq.py, sq_diag.py, sq_near.py)."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sc import Line


def Tpot(d):
    q = d['q']; return q, 2 * (1 - np.cos(q)) * d['S']


def logbin(q, y, nb=28):
    e = np.geomspace(q[1], q[-1] * 1.0001, nb); qc, yc = [], []
    for lo, hi in zip(e[:-1], e[1:]):
        s = (q >= lo) & (q < hi)
        if s.any(): qc.append(np.sqrt(lo * hi)); yc.append(y[s].mean())
    return np.array(qc), np.array(yc)


fig, ax = plt.subplots(1, 3, figsize=(14, 4))
# (a) mode temperatures
cases = [('data/sq_rf1024_f0.1.npz', 'RF $m=1024$, $f=0.10$', 'C3'),
         ('data/sq_rf1024_f0.11.npz', 'RF $m=1024$, $f=0.11$', 'C1'),
         ('data/sq_rf1024.npz', 'RF $m=1024$, $f=0.12$', 'C0'),
         ('data/sq_rf256.npz', 'RF $m=256$, $f=0.17$', 'C2'),
         ('data/sq_rb1024.npz', 'RB $m=1024$, $f=0.12$', 'C4')]
for fn, lab, c in cases:
    d = np.load(fn); q, T = Tpot(d); qc, Tc = logbin(q, T)
    ax[0].loglog(qc, Tc, 'o', ms=3, color=c, label=lab)
for lab, c in [('rf1024', 'C0'), ('rf256', 'C2'), ('rb1024', 'C4')]:
    t = np.load(f'data/sq_theory_{lab}.npz'); q = t['q']
    ax[0].loglog(q[1:], (2 * (1 - np.cos(q)) * t['S_b'])[1:], '-', color=c, lw=0.8)
ax[0].set_xlabel('$q$'); ax[0].set_ylabel(r'$T(q)=\lambda_q S(q)$')
ax[0].set_ylim(3e-3, 30); ax[0].legend(fontsize=7, loc='lower left')
ax[0].set_title('(a) mode temperatures; lines: bare 2nd order')

# (b) single-site force correlator
d = np.load('data/sq_diag_m1024_f0.12.npz'); t = d['t']
ax[1].semilogy(t, np.abs(d['C']), 'k.', ms=2, label='measured $|C(t)|$')
ax[1].semilogy(t, d['CG'], 'C0-', label='Gaussian smearing, measured $B(t)$')
ax[1].semilogy(t, np.maximum(d['Cb'], 1e-9), 'C3--', label=r'bare $\Delta(vt)$')
ax[1].set_xlim(0, 300); ax[1].set_ylim(1e-5, 0.3); ax[1].set_xlabel('$t$')
ax[1].legend(fontsize=7); ax[1].set_title('(b) $m=1024$, $f=0.12$')

# (c) large-scale temperature: measured vs S_xi(0)/(2 eta_eff), eta_eff = eta + dF/dv|_T
rows = []
for m, f, sqf, dg in [(1024, 0.10, 'data/sq_rf1024_f0.1.npz', 'data/sq_diag_m1024_f0.1.npz'),
                      (1024, 0.11, 'data/sq_rf1024_f0.11.npz', 'data/sq_diag_m1024_f0.11.npz'),
                      (1024, 0.12, 'data/sq_rf1024.npz', 'data/sq_diag_m1024_f0.12.npz'),
                      (256, 0.17, 'data/sq_rf256.npz', 'data/sq_diag_m256_f0.17.npz')]:
    s = np.load(sqf); g = np.load(dg)
    q, T = Tpot(s); Tlw = T[1:6].mean()
    v = float(s['v']); Tk = float(s['Tkin'])
    Sx = 2 * np.trapezoid(g['C'], g['t'])
    ln = Line(float(m)); h = 1e-3 * v
    eta_eff = 1 + (ln.F(v + h, Tk) - ln.F(v - h, Tk)) / (2 * h)
    rows.append(dict(m=m, f=f, v=v, Tkin=Tk, T_lw=Tlw, S_xi0=Sx, eta_eff_needed=Sx / (2 * Tlw),
                     eta_eff_closure=eta_eff, T_lw_pred=Sx / (2 * eta_eff), T_lw_bare=1 / (6 * v)))
for r in rows:
    mk = 'o' if r['m'] == 1024 else 's'
    ax[2].plot(r['v'], r['T_lw'], mk, color='k')
    ax[2].plot(r['v'], r['T_lw_pred'], mk, mfc='none', color='C0')
    ax[2].plot(r['v'], r['T_lw_bare'], mk, mfc='none', color='C3')
    ax[2].plot(r['v'], r['Tkin'], mk, mfc='none', color='0.6')
ax[2].plot([], [], 'ko', label='measured $T_{\\rm lw}$')
ax[2].plot([], [], 'o', mfc='none', color='C0', label=r'$S_\xi(0)/2\eta_{\rm eff}$')
ax[2].plot([], [], 'o', mfc='none', color='C3', label=r'bare $A^2/6\eta v$')
ax[2].plot([], [], 'o', mfc='none', color='0.6', label=r'$T_{\rm kin}$')
ax[2].set_xlabel('$v$'); ax[2].set_ylabel('temperature'); ax[2].legend(fontsize=7)
ax[2].set_title('(c) $T_{\\rm lw}$; circles $m=1024$, squares $256$')
fig.tight_layout(); fig.savefig('fig11_roughness.pdf'); fig.savefig('fig11_roughness.png', dpi=110)
json.dump(rows, open('statistics5.json', 'w'), indent=1)
for r in rows: print({k: round(v, 3) if isinstance(v, float) else v for k, v in r.items()})
