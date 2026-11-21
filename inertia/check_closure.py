"""Open point: test the closure's Gaussian, equilibrium-shaped fluctuations.

For a uniformly sliding line, measure
  v (mean velocity), T_kin = m <du'^2>,
  B(tau) = <[du_i(t+tau) - du_i(t)]^2>, du_i = u_i - ubar(t),
  the excess kurtosis of the increments,
and compare with the closure: T(v) from the power balance and B = T b(tau).
Writes data/closure_check_m{m}_f{f}.npz.
"""
import sys, time
import numpy as np
from protocols import _run, default_dt
from sc import Line

cases = [(64.0, 0.25), (256.0, 0.17), (1024.0, 0.12), (4096.0, 0.09)]
if len(sys.argv) > 1:
    cases = [(float(sys.argv[1]), float(sys.argv[2]))]
L = 2048
for m, f in cases:
    t0 = time.time()
    dt = default_dt(m)
    U = np.zeros((1, L)); V = np.full((1, L), f); seeds = [21]
    n = int(max(20 * m, 4000) / dt)
    _run(U, V, seeds, [f], m, n, n, 1.0, 1.0, 1.0, dt, 0.0, False)
    dts = max(0.5, dt)                     # snapshot spacing
    T_rec = max(4000.0, 8 * m)
    ns = int(round(dts / dt)); nsnap = int(T_rec / dts)
    D = np.zeros((nsnap, L), np.float32)
    v2 = 0.0
    u0 = U[0].mean()
    for k in range(nsnap):
        _run(U, V, seeds, [f], m, ns, ns, 1.0, 1.0, 1.0, dt, 0.0, False)
        D[k] = U[0] - U[0].mean()
        v2 += ((V[0] - V[0].mean()) ** 2).mean()
    v = (U[0].mean() - u0) / (nsnap * ns * dt)
    Tkin = m * v2 / nsnap
    lags = np.unique(np.round(np.geomspace(1, nsnap // 4, 40)).astype(int))
    B = np.array([np.mean((D[l:] - D[:-l]).astype(float) ** 2) for l in lags])
    kurt = []
    for l in lags:
        x = (D[l:] - D[:-l]).astype(float).ravel()
        kurt.append(np.mean(x ** 4) / np.mean(x ** 2) ** 2 - 3)
    taus = lags * dts
    Lm = Line(m, tau_max=float(taus[-1] + 10))
    Tcl = Lm.T_of_v(v)
    bcl = np.interp(taus, Lm.taus, Lm.b)
    np.savez(f'data/closure_check_m{m:g}_f{f:g}.npz', m=m, f=f, v=v, Tkin=Tkin, Tcl=Tcl,
             taus=taus, B=B, bcl=bcl, kurt=np.array(kurt),
             F_meas=f - v, F_cl=Lm.F(v, Tcl))
    r = B / (Tcl * bcl)
    print(f'm={m:g} f={f}: v={v:.4f}  T_kin={Tkin:.3f}  T_closure={Tcl:.3f}  '
          f'F_pin meas={f - v:.4f} closure={Lm.F(v, Tcl):.4f}  '
          f'B/(T b): short {r[0]:.2f} mid {r[len(r)//2]:.2f} long {r[-1]:.2f}  '
          f'max|kurt|={np.max(np.abs(kurt)):.2f}  ({time.time() - t0:.0f}s)', flush=True)
