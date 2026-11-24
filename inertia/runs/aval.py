"""Avalanches under slow parabolic-well driving, for a range of masses.

The well k0 (w - u_i) advances at wdot.  An avalanche is a maximal interval in which the
mean velocity vbar exceeds vth = c_th * wdot; its size is s = L int (vbar - wdot) dt and its
duration the interval length.  Detection is done chunk by chunk (state carried over), so no
long time series is stored.  The force k0 (w - ubar) is kept every T_force time units.

Expectations (main.tex, Sec. "Exponents"): qEW tau ~ 1.11 for small M, cut off at
s ~ xi^(1+zeta), xi = k0^(-1/2); at large M, cold avalanches up to a nucleation size, then
system-spanning hot slips.

Usage: python3 aval.py pilot | scan [m ...] | dtcheck | tau [m ...]
Writes data/aval_m{m}.npz (scan) and data/avalT_m{m}.npz (tau: k0 = 5e-4, L = 2048).
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import sys
import time
import numpy as np
from protocols import _run, default_dt
from config import C_TH

L, K0, WDOT = 512, 3e-3, 1e-5
SEEDS = np.arange(12) + 61000


def avalanches(m, W_total, W_skip, L=L, k0=K0, wdot=WDOT, seeds=SEEDS, A=1.0, f0=0.26,
               c_th=C_TH, t_rec=1.0, chunk=1e5, T_force=100.0, W_prep=30.0, dt=None):
    dt = dt or default_dt(m)
    S = len(seeds)
    nrec = max(1, int(round(t_rec / dt))); t_rec = nrec * dt
    ns = int(round(chunk / t_rec)) * nrec
    U = np.zeros((S, L)); V = np.zeros((S, L))
    w = np.full(S, f0 / k0)
    # preparation: fast loading into the steady state, then slow driving
    wp = 2e-4
    _run(U, V, seeds, np.zeros(S), m, int(round(W_prep / wp / dt)), 1000, A, 1.0, 1.0, dt, 0.0,
         False, k0=k0, w=w, wdot=wp)
    vth = c_th * wdot
    # carried state per sample: inside avalanche?, accumulated size and duration, onset force
    inside = np.zeros(S, bool); acc_s = np.zeros(S); acc_T = np.zeros(S); f_on = np.zeros(S)
    out_s, out_T, out_k, out_f, out_w = [], [], [], [], []
    w0 = w.copy()
    force_t, force = [], []
    t = 0.0; T_total = W_total / wdot; t_skip = W_skip / wdot
    while t < T_total:
        vb, _, ub, _, _ = _run(U, V, seeds, np.zeros(S), m, ns, nrec, A, 1.0, 1.0, dt, 0.0,
                               False, k0=k0, w=w, wdot=wdot)
        tt = t + t_rec * (np.arange(vb.shape[1]) + 1)
        wt = w[:, None] - wdot * (tt[-1] - tt)[None, :]
        F = k0 * (wt - ub)
        kf = max(1, int(round(T_force / t_rec)))
        force_t.append(tt[::kf]); force.append(F[:, ::kf])
        for k in range(S):
            on = vb[k] > vth
            # boundaries of runs of `on`
            d = np.diff(np.concatenate([[inside[k]], on, [False]]).astype(int))
            starts = np.where(d == 1)[0]; ends = np.where(d == -1)[0]
            if inside[k]:
                starts = np.concatenate([[-1], starts])
            for a, b in zip(starts, ends):
                lo = max(a, 0)
                seg = vb[k, lo:b]
                s_add = L * (seg - wdot).sum() * t_rec
                T_add = (b - lo) * t_rec
                if a == -1:
                    s_tot, T_tot, fo = acc_s[k] + s_add, acc_T[k] + T_add, f_on[k]
                else:
                    s_tot, T_tot, fo = s_add, T_add, F[k, max(lo - 1, 0)]
                if b == len(on):   # still running at chunk end
                    inside[k] = True; acc_s[k] = s_tot; acc_T[k] = T_tot; f_on[k] = fo
                    continue
                inside[k] = False
                if tt[lo] > t_skip:
                    out_s.append(s_tot); out_T.append(T_tot); out_k.append(k); out_f.append(fo)
                    out_w.append(wt[k, lo])
            if len(starts) == 0 and not on.any():
                inside[k] = False
        t = tt[-1]
    return dict(s=np.array(out_s), T=np.array(out_T), sample=np.array(out_k),
                f_on=np.array(out_f), w_on=np.array(out_w), force_t=np.concatenate(force_t),
                force=np.concatenate(force, axis=1), m=m, L=L, k0=k0, wdot=wdot, c_th=c_th,
                W_skip=W_skip, W_total=W_total)


def summary(r):
    s = r['s']; s = s[s > 0]
    if len(s) < 5:
        return f'n {len(s)}'
    # crude tau: ML estimate for s in [s_lo, s_hi] (pure power law between the bounds)
    lo, hi = 3.0, np.percentile(s, 90)
    x = s[(s >= lo) & (s <= hi)]
    taus = np.linspace(0.5005, 2.5005, 400)
    ll = [-(t * np.log(x)).sum() - len(x) * np.log((hi ** (1 - t) - lo ** (1 - t)) / (1 - t))
          for t in taus]
    tau = taus[int(np.argmax(ll))]
    return (f'n {len(s)} median {np.median(s):.1f} p90 {np.percentile(s, 90):.0f} '
            f'max {s.max():.0f} tau_ML[3,p90] {tau:.2f} <f> {r["force"].mean():.4f}')


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'pilot'
    if what == 'pilot':
        t0 = time.time()
        r = avalanches(0.0, W_total=12.0, W_skip=1.0)
        np.savez('data/aval_pilot_m0.npz', **r)
        print(time.strftime('%H:%M:%S'), 'pilot m=0', summary(r), f'{time.time() - t0:.0f}s',
              flush=True)
    elif what == 'dtcheck':
        for dt in [0.05, 0.1]:
            t0 = time.time()
            r = avalanches(0.0, W_total=26.0, W_skip=2.0, dt=dt)
            np.savez(f'data/aval_dtcheck_{dt:g}.npz', **r)
            print(time.strftime('%H:%M:%S'), 'dtcheck', dt, summary(r), f'{time.time() - t0:.0f}s', flush=True)
    elif what == 'tau':
        # soft well: xi = k0^(-1/2) ~ 45, L = 2048; dt = 0.1 for m <= 4 (checked by 'dtcheck')
        for m in [float(x) for x in sys.argv[2:]] or [0, 1, 2, 4, 8, 16]:
            t0 = time.time()
            r = avalanches(m, W_total=50.0, W_skip=2.0, L=2048, k0=5e-4, wdot=1e-5, f0=0.28,
                           W_prep=60.0, dt=max(0.1, default_dt(m)))
            np.savez(f'data/avalT_m{m:g}.npz', **r)
            print(time.strftime('%H:%M:%S'), 'tau m', m, summary(r), f'{time.time() - t0:.0f}s',
                  flush=True)
    else:
        ms = [float(x) for x in sys.argv[2:]] or [0, 1, 4, 16, 64, 256, 1024]
        for m in ms:
            t0 = time.time()
            r = avalanches(m, W_total=70.0, W_skip=2.0)
            np.savez(f'data/aval_m{m:g}.npz', **r)
            print(time.strftime('%H:%M:%S'), 'aval m', m, summary(r), f'{time.time() - t0:.0f}s',
                  flush=True)
