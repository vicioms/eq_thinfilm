"""Second-order (linear-response) friction of the moving line.

At velocity v, site i sees the force F_i(v t) with power spectrum
S(omega) = (1/v) Delta(omega/v), Delta(k) = sigma^2 sinc^4(k/2), sigma^2 = A^2/3.
The line responds with the local susceptibility

    chi(omega) = int dq/2pi 1/(lambda_q - m omega^2 - i eta omega),
    lambda_q = 2 kappa (1 - cos q),

and absorbs the power  P(v) = int_0^inf (d omega/pi) S(omega) omega Im chi(omega).
Energy balance f v = eta <u'^2> = eta v^2 + P gives the flow curve

    f(v) = eta v + P(v)/v,      eta <du'^2> = P(v)    (exact identity: f - eta v = eta <du'^2>/v).

Radiation limit (eta -> 0 at fixed omega in the phonon band, omega < 2 sqrt(kappa/m)):
omega Im chi -> 1/sqrt(m (4 kappa - m omega^2)) -> 1/(2 sqrt(kappa m)), so for v small enough
that the whole spectrum is in the band, P -> <F^2>/(2 sqrt(kappa m)), <F^2> = 2 A^2/9.
"""
import numpy as np
from scipy.integrate import quad


def chi_loc(omega, m, eta, kappa=1.0):
    """Exact lattice local susceptibility (complex), any m >= 0, eta > 0."""
    z = m * omega**2 + 1j * eta * omega
    a = 2 * kappa - z
    s = np.sqrt(a - 2 * kappa + 0j) * np.sqrt(a + 2 * kappa + 0j)
    # choose the branch with |a - s| < 2 kappa (decaying Green's function)
    s = np.where(np.abs(a - s) < 2 * kappa, s, -s)
    return 1.0 / s


def Delta_k(k, A):
    x = 0.5 * k
    sinc = np.where(np.abs(x) < 1e-8, 1.0, np.sin(x) / np.where(x == 0, 1, x))
    return (A**2 / 3.0) * sinc**4


def power(v, m, eta, A, kappa=1.0):
    """P(v) = int_0^inf dk/pi Delta(k) omega Im chi(omega), omega = v k."""
    def integrand(k):
        w = v * k
        return Delta_k(k, A) * w * np.imag(chi_loc(w, m, eta, kappa)) / np.pi
    # split at the band edge and at the sinc zeros for accuracy
    pts = [2 * np.pi * j for j in range(1, 12)]
    if m > 0:
        kb = 2 * np.sqrt(kappa / m) / v
        pts = sorted(set(pts + [kb]))
    tot = 0.0
    lo = 0.0
    for p in pts + [200.0]:
        if p <= lo:
            continue
        tot += quad(integrand, lo, p, limit=400)[0]
        lo = p
    return tot


def flow_curve(vs, m, eta, A, kappa=1.0):
    P = np.array([power(v, m, eta, A, kappa) for v in vs])
    return eta * vs + P / vs, P
