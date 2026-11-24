"""Gaussian DMFT of the hot sliding line: deterministic, no simulation.

The cavity solution (cavity.py) shows that the force autocorrelation C(t) is the disorder
correlator smeared by the true relative-displacement spread B(t) of a site.  Keeping only
that (Gaussian displacements) closes the single-site problem on three deterministic maps:

  B -> C:  C(t)   = <Delta(v t + x)>_{x ~ N(0, B(t))}   (real-space quadrature; Delta compact)
  C -> S:  S_xi(w) = 2 int_0^inf C(t) cos(w t) dt
  S -> B:  B(t)   = (2/pi) int_0^inf dw S_w(w) (1 - cos w t),
           S_w(w) = S_xi(w) Im g(z)/(eta_eff w),  g(z) = int dq/2pi 1/(lam_q - z) = 1/sqrt(z(z - 4 kappa)),
           z = m w^2 + i eta_eff w     (a linear chain driven by independent noises S_xi, friction eta_eff)

and the friction of the medium from the velocity dependence of the mean force at fixed bath:

  <F>     = int_0^inf R_loc(t) <Delta'(v t + x)> dt,      f(v) = eta v - <F>
  eta_eff = eta - d<F>/dv |_{B, R} = eta - int_0^inf R_loc(t) t <Delta''(v t + x)> dt,

R_loc(t) = int dq/2pi R_q(t), R_q the response of a damped mode with friction eta_eff.
Delta = (A^2/3) M4 (cubic B-spline) for the random field of sim.py.
Outputs: f(v), S_xi(0), eta_eff, T_lw = S_xi(0)/2 eta_eff, T_kin = m <dv^2>, B(t), C(t).

Usage: python3 gdmft.py m v [v ...]
"""
import sys
import numpy as np
import config  # noqa: F401  (run from the inertia folder)

KAPPA, ETA, A = 1.0, 1.0, 1.0
ETA_FLOOR = 0.03
qs = (np.arange(4096) + 0.5) * np.pi / 4096                 # q in (0, pi): even integrand
lam_q = 2 * KAPPA * (1 - np.cos(qs))


def _M4(x, d=0):
    """Cubic B-spline M4 and its first two derivatives (Delta = (A^2/3) M4)."""
    a = np.abs(x); s = np.sign(x)
    if d == 0:
        return np.where(a < 1, 2 / 3 - a ** 2 + a ** 3 / 2, np.where(a < 2, (2 - a) ** 3 / 6, 0.0))
    if d == 1:
        return np.where(a < 1, s * (-2 * a + 1.5 * a ** 2), np.where(a < 2, -s * (2 - a) ** 2 / 2, 0.0))
    return np.where(a < 1, -2 + 3 * a, np.where(a < 2, 2 - a, 0.0))


_y = np.linspace(-2, 2, 801); _dy = _y[1] - _y[0]


def smeared(t, v, B, deriv=0):
    """<Delta^(deriv)(v t + x)>, x ~ N(0, B): real-space quadrature over the support of Delta
    (exact for B -> 0 through a second-order Taylor term when sqrt(B) is below the grid)."""
    D = (A ** 2 / 3) * _M4(_y, deriv)
    out = np.empty_like(t)
    for i0 in range(0, len(t), 1024):
        x = v * t[i0:i0 + 1024]; b = B[i0:i0 + 1024]
        small = b < (3 * _dy) ** 2
        bb = np.where(small, 1.0, b)
        phi = np.exp(-(_y[None, :] - x[:, None]) ** 2 / (2 * bb[:, None])) / np.sqrt(2 * np.pi * bb[:, None])
        val = (phi * D[None, :]).sum(1) * _dy
        if small.any():
            xs = x[small]
            d0 = (A ** 2 / 3) * _M4(xs, deriv)
            d2 = (A ** 2 / 3) * _M4(xs, deriv + 2) if deriv == 0 else 0.0
            val[small] = d0 + 0.5 * b[small] * d2
        out[i0:i0 + 1024] = val
    return out


def R_loc(t, m, eta):
    """Local response int dq/2pi R_q(t) of a chain of damped modes (friction eta)."""
    g = eta / m
    W2 = lam_q / m - 0.25 * g * g
    out = np.zeros_like(t)
    for i0 in range(0, len(t), 256):
        tt = t[i0:i0 + 256, None]
        e = np.exp(-0.5 * g * tt)
        with np.errstate(invalid='ignore', over='ignore'):
            Wu = np.sqrt(np.abs(W2))[None, :]
            under = e * np.sin(Wu * tt) / (m * Wu)
            over = (np.exp((-0.5 * g + Wu) * tt) - np.exp((-0.5 * g - Wu) * tt)) / (2 * m * Wu)
        r = np.where(W2[None, :] > 0, under, over)
        out[i0:i0 + 256] = r.mean(1)          # mean over q in (0, pi) = int dq/2pi over (-pi, pi)
    return out


def img_over_imz(w, m, eta):
    z = m * w ** 2 + 1j * eta * w
    g = 1.0 / np.sqrt(z * (z - 4 * KAPPA))
    g = np.where(g.imag < 0, -g, g)       # retarded branch: Im g > 0 for Im z > 0
    return g.imag / (eta * w)


def solve(m, v, it=60, mix=0.4, dt=0.5, nt=2 ** 14, tol=1e-4, verbose=True):
    t = np.arange(nt) * dt
    w = np.geomspace(1e-7, 20.0, 8000); lw = np.log(w)
    sub = np.arange(0, nt, 8); ts = t[sub]          # B is smooth: evaluate on every 8th time
    nf = 4 * nt; wf = 2 * np.pi * np.fft.rfftfreq(nf, dt)
    tc = t < 3000                                     # C', C'' only multiply the response
    B = np.zeros(nt); eta_eff = ETA
    for n in range(it):
        C = smeared(t, v, B)
        # S_xi(w) = 2 int_0^T C cos(w t) dt (trapezoid) by FFT, interpolated onto the log grid
        Sf = 2 * dt * (np.fft.rfft(C, nf).real - 0.5 * C[0])
        Sx = np.interp(w, wf, Sf)
        Sw = Sx * img_over_imz(w, m, eta_eff)
        Bs = (2 / np.pi) * np.trapezoid(Sw[None, :] * (1 - np.cos(w[None, :] * ts[:, None])) * w[None, :],
                                        lw, axis=1)
        Bn = np.interp(t, ts, Bs)
        R = R_loc(t, m, eta_eff)
        C1 = smeared(t[tc], v, B[tc], 1); C2 = smeared(t[tc], v, B[tc], 2)
        Fm = np.trapezoid(R[tc] * C1, t[tc])
        eta_n = ETA - np.trapezoid(R[tc] * t[tc] * C2, t[tc])
        # the bare start has eta_eff < 0 below the bare minimum v_*: clip during the iteration and
        # call it 'no hot solution' only if the converged friction sits at the floor
        eta_n = max(eta_n, ETA_FLOOR)
        err = np.max(np.abs(Bn - B)[t < 2000] / (1 + Bn[t < 2000])) + abs(eta_n - eta_eff)
        B = mix * Bn + (1 - mix) * B
        eta_eff = mix * eta_n + (1 - mix) * eta_eff
        if verbose and (n % 5 == 0 or err < tol):
            print(f'  it {n}: S_xi(0) {Sx[0]:.3f} eta_eff {eta_eff:.3f} f {ETA * v - Fm:.4f} '
                  f'B(50) {np.interp(50, t, B):.2f} err {err:.2e}', flush=True)
        if err < tol:
            break
    no_sol = bool(eta_eff < 1.5 * ETA_FLOOR)   # friction pinned at the floor: uniform hot line unstable
    vel2 = (1 / np.pi) * np.trapezoid(w ** 2 * Sw * w, lw)
    out = dict(t=t, C=C, B=B, w=w, Sxi=Sx, Sxi0=float(Sx[0]), eta_eff=float(eta_eff), f=float(ETA * v - Fm),
               T_lw=float(Sx[0] / (2 * eta_eff)), Tkin=float(m * vel2), converged=bool(err < tol), m=m, v=v,
               no_solution=no_sol)
    return out


def local_time(t, v, B):
    """Local-time estimate S_xi(0)/Delta_hat(0) = 2 int_0^inf exp(-v^2 t^2 / 2B)/sqrt(2 pi B) dt."""
    with np.errstate(divide='ignore', invalid='ignore'):
        g = np.where(B > 0, np.exp(-v * v * t * t / (2 * B)) / np.sqrt(2 * np.pi * B), 0.0)
    return 2 * np.trapezoid(g, t)


if __name__ == '__main__':
    m = float(sys.argv[1])
    for v in [float(x) for x in sys.argv[2:]]:
        print(f'm {m:g} v {v}', flush=True)
        r = solve(m, v)
        np.savez(f'data/gdmft_m{m:g}_v{v:g}.npz', **r)
        print(f'  -> f {r["f"]:.4f} S_xi(0) {r["Sxi0"]:.3f} (bare {1 / (3 * v):.3f}) eta_eff {r["eta_eff"]:.3f} '
              f'T_lw {r["T_lw"]:.3f} T_kin {r["Tkin"]:.3f} converged {r["converged"]}', flush=True)
