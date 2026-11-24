"""Gaussian DMFT with the frequency-dependent disorder self-energy Sigma(w).

gdmft.py damps the medium with a constant eta_eff = eta - int t R C'' dt at all frequencies.
But the disorder friction comes from a correlator that decorrelates at rate ~ v, so it must
vanish at high frequency; a medium damped by eta_eff < eta at all frequencies runs too hot
and violates the power balance.  Here the mean disorder force responds to a perturbation of
the site's past displacement through (one-loop MSRJD, smeared correlators)

    d<F(t)> = int ds R(t-s) <Delta''(u_t - u_s)> (du(t) - du(s))
    =>  K(w) = int_0^inf dt R_loc(t) C2(t) (1 - e^{i w t}),    C2(t) = <Delta''(v t + x)>_{N(0,B)},

K(0) = 0 (no static restoring force), -dK/d(iw)|_0 = -int t R C2 dt (the old eta_eff - eta).
Every site of the medium then has  chi_q(w) = 1/(lam_q - z(w)),  z = m w^2 + i eta w + K(w), so

    B(t)     = (2/pi) int_0^inf dw S_xi(w) Im g(z)/Im z (1 - cos w t),   g(z) = 1/sqrt(z (z - 4 kappa))
    R_loc(t) = (2/pi) int_0^inf dw Im g(z(w)) sin(w t)
    <F>      = int R_loc(t) <Delta'(v t + x)> dt,   f = eta v - <F>,
    C(t)     = <Delta(v t + x)>,  S_xi(w) = 2 int_0^inf C cos(w t) dt.

Diagnostics: power balance  T_kin = m <dv^2>  against  m v (f - eta v)/eta  (exact identity);
eta_eff = eta - int t R C2 dt (zero-frequency friction); min_w Im z/w (stability of the medium).

Usage: python3 gdmft2.py m v [v ...]     (writes data/gdmft2_m{m}_v{v}.npz)
"""
import sys
import numpy as np
import config  # noqa: F401  (run from the inertia folder)
from gdmft import smeared, KAPPA, ETA

DT, NT = 0.5, 2 ** 14            # time grid for C, B
TCUT = 3000.0                    # R, C', C'' are needed below this lag
NU = 2 ** 17                     # uniform frequency grid w_k = pi k / (NU DT), k < NU
WU = np.pi * np.arange(NU) / (NU * DT); DWU = WU[1]
WL = np.geomspace(1e-7, 20.0, 8000); LWL = np.log(WL)


def g_of(z):
    g = 1.0 / np.sqrt(z * (z - 4 * KAPPA))
    return np.where(g.imag < 0, -g, g)


def solve(m, v, it=80, mix=0.4, tol=2e-4, verbose=True):
    t = np.arange(NT) * DT
    tc = t < TCUT; nc = int(tc.sum())
    sub = np.unique(np.r_[np.arange(64), np.arange(0, NT, 8)]); ts = t[sub]
    wtr = np.full(nc, DT); wtr[0] = DT / 2                       # trapezoid weights on t < TCUT
    nf = 4 * NT; wf = 2 * np.pi * np.fft.rfftfreq(nf, DT)
    B = np.zeros(NT)
    R = None; K_u = np.zeros(NU, complex); K_l = np.zeros(len(WL), complex)
    hist = []
    for n in range(it):
        C = smeared(t, v, B)
        C1 = smeared(t[tc], v, B[tc], 1); C2 = smeared(t[tc], v, B[tc], 2)
        Sf = 2 * DT * (np.fft.rfft(C, nf).real - 0.5 * C[0])
        Sx = np.interp(WL, wf, Sf)
        # self-energy from the current response (bare medium on the first pass)
        if R is not None:
            b = np.zeros(2 * NU); b[:nc] = wtr * R * C2
            Fb = np.fft.fft(b)[:NU]
            K_u = b.sum() - np.conj(Fb)                          # sum_n b_n (1 - e^{+i w_k t_n})
            slope = -1j * (t[:nc] * b[:nc]).sum()                # small-w: K ~ -i w int t R C2
            K_l = np.where(WL < DWU, slope * WL,
                           np.interp(WL, WU, K_u.real) + 1j * np.interp(WL, WU, K_u.imag))
        z_u = m * WU ** 2 + 1j * ETA * WU + K_u
        z_l = m * WL ** 2 + 1j * ETA * WL + K_l
        imz_l = z_l.imag
        stab = float(np.min(imz_l / WL))                         # effective friction at its lowest
        imz_l = np.maximum(imz_l, 1e-3 * ETA * WL)               # keep the medium causal while iterating
        z_l = z_l.real + 1j * imz_l
        z_u = z_u.real + 1j * np.maximum(z_u.imag, 1e-3 * ETA * WU)
        Sw = Sx * g_of(z_l).imag / imz_l
        Bs = (2 / np.pi) * np.trapezoid(Sw[None, :] * (1 - np.cos(WL[None, :] * ts[:, None])) * WL[None, :],
                                        LWL, axis=1)
        Bn = np.interp(t, ts, Bs)
        # local response: R(t) = (2/pi) int Im g sin(w t) dw, by FFT on the uniform grid
        a = np.zeros(2 * NU); gi = g_of(z_u[1:]).imag; a[1:NU] = gi
        Rn = -(2 / np.pi) * DWU * np.fft.fft(a).imag[:nc]
        Fm = float(np.sum(wtr * Rn * C1))
        eta_eff = ETA - float(np.sum(wtr * Rn * t[:nc] * C2))
        err = np.max(np.abs(Bn - B)[t < 2000] / (1 + Bn[t < 2000]))
        if R is not None:
            err += np.max(np.abs(Rn - R)) / (np.max(np.abs(Rn)) + 1e-12)
        B = Bn if n == 0 else mix * Bn + (1 - mix) * B
        R = Rn if R is None else mix * Rn + (1 - mix) * R
        f = ETA * v - Fm
        vel2 = (1 / np.pi) * np.trapezoid(WL ** 2 * Sw * WL, LWL)
        Tkin = m * vel2; Tpb = m * v * (f - ETA * v) / ETA
        hist.append((n, f, Sx[0], eta_eff, Tkin, Tpb, stab, err))
        if verbose and (n % 10 == 0 or err < tol):
            print(f'  it {n}: f {f:.4f} S_xi(0) {Sx[0]:.3f} eta_eff {eta_eff:.3f} T_kin {Tkin:.3f} '
                  f'(power balance {Tpb:.3f}) min Im z/w {stab:.3f} err {err:.1e}', flush=True)
        if err < tol and n > 3:
            break
    return dict(t=t, C=C, B=B, R=R, w=WL, Sxi=Sx, K=K_l, Sxi0=float(Sx[0]), eta_eff=eta_eff, f=f,
                T_lw=float(Sx[0] / (2 * eta_eff)) if eta_eff > 0 else np.nan, Tkin=Tkin, T_pb=Tpb,
                stab=stab, converged=bool(err < tol), m=m, v=v, hist=np.array(hist))


if __name__ == '__main__':
    m = float(sys.argv[1])
    for v in [float(x) for x in sys.argv[2:]]:
        print(f'm {m:g} v {v}', flush=True)
        r = solve(m, v)
        np.savez(f'data/gdmft2_m{m:g}_v{v:g}.npz', **r)
        print(f'  -> f {r["f"]:.4f} S_xi(0) {r["Sxi0"]:.3f} eta_eff {r["eta_eff"]:.3f} T_lw {r["T_lw"]:.3f} '
              f'T_kin {r["Tkin"]:.3f} vs power balance {r["T_pb"]:.3f}  min Im z/w {r["stab"]:.3f} '
              f'converged {r["converged"]}', flush=True)
