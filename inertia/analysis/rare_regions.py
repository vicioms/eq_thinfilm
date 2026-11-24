"""Figure 14: are long-lived active patches below f_inf quenched rare regions? (rare.py)

For each disorder sample, two copies start from displaced lines.  Quenched rare regions
would put the late activity at the same x in both copies.  Printed: overlap of the late
time-averaged activity profiles between copies (cosine similarity), against the same
quantity between different samples (the chance level), plus the drift of each patch.
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = np.load('data/rare.npz')
t = d['t'] / float(d['m']); A = d['act']; B = int(d['B']); L = int(d['L'])   # A[t, copy, block]
x = (np.arange(A.shape[2]) + 0.5) * B
rho = A.mean(2)
print('rho at t/m = 1000, 3000, 6000, 11900:')
for k in range(4):
    print(f'  seed {d["seeds"][k]} f {d["f"][k]} off {d["offs"][k]}:',
          ' '.join(f'{rho[np.searchsorted(t, s), k]:.4f}' for s in [1000, 3000, 6000, 11900]))


def late_profile(k, frac=0.5):
    """Time-averaged activity profile over the last `frac` of this copy's active life."""
    alive = np.where(rho[:, k] > 0)[0]
    if len(alive) == 0:
        return None, 0
    t_end = alive[-1]; i0 = int(t_end * (1 - frac))
    return A[i0:t_end + 1, k].mean(0), t[t_end]


def cos(a, b):
    return float(a @ b / np.sqrt((a @ a) * (b @ b)))


P = [late_profile(k) for k in range(4)]
for k in range(4):
    print(f'copy {k}: last activity at t/m = {P[k][1]:.0f}')
pairs = {'same disorder (31003)': (0, 1), 'same disorder (31005)': (2, 3),
         'different disorder': (0, 2), 'different disorder ': (1, 3)}
for lab, (i, j) in pairs.items():
    if P[i][0] is not None and P[j][0] is not None:
        print(f'  overlap {lab}: {cos(P[i][0], P[j][0]):.3f}')

# patch drift: centre and width of activity in the last active stretch, per copy
for k in range(4):
    alive = np.where(rho[:, k] > 0)[0]
    if len(alive) == 0:
        continue
    sel = alive[alive > alive[-1] * 0.5]
    w = A[sel, k]
    cen = (w * x).sum(1) / w.sum(1)
    n_act = (w > 0).sum(1)
    print(f'copy {k}: centre x from {cen[0]:.0f} to {cen[-1]:.0f} (range {cen.min():.0f}-{cen.max():.0f}); '
          f'active blocks median {np.median(n_act):.0f}')

fig, ax = plt.subplots(1, 4, figsize=(17, 4.2), sharey=True)
for k in range(4):
    ax[k].imshow(A[:, k].T, aspect='auto', origin='lower', cmap='magma', vmin=0, vmax=0.6,
                 extent=[t[0], t[-1], 0, L], interpolation='nearest')
    ax[k].set_xlabel('$t/m$')
    ax[k].set_title(f'disorder {d["seeds"][k]}, $f={d["f"][k]}$, offset {d["offs"][k]:g}', fontsize=9)
ax[0].set_ylabel('$x$')
fig.suptitle('Activity (fraction of sites with $|\\dot u|>0.02$, blocks of 64), $m=64$, $L=65536$', fontsize=10)
fig.tight_layout(); fig.savefig('figures/rare_regions.png', dpi=110); fig.savefig('figures/rare_regions.pdf')
print('ok')
