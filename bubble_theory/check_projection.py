"""Numerical checks of the phi^4 -> interface projection (notes_phi4_interface.md).

Checks, for eps0 = c = 1 by default:
  1. kink constants N1 = int phi0'^2 = 4/(3w), sigma = c N1, int phi0' = 2;
  2. Poschl-Teller spectrum of L = -c d^2 + eps0 (3 phi0^2 - 1): 0, 3 eps0/2, continuum >= 2 eps0;
  3. random-field overlap  int phi0'(y) phi0'(y-v) dy  vs (4/w)[b coth b - 1]/sinh^2 b;
  4. random-bond potential correlator  int sech^2 sech^2  and Delta_RB = -R'';
  5. Fourier transform of phi0': pi q w / sinh(pi q w / 2);
  6. bubble projection: residual of the clean equation on phi = tanh((R - rho)/w)
     projected on phi0' gives  eta Rdot = 2h - sigma/R  (Laplace term) up to O(w^2/R^3).

Usage: python3 check_projection.py [--eps0 1] [--c 1]
"""
import argparse
import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh_tridiagonal


def kink(eps0, c):
    w = np.sqrt(2 * c / eps0)
    phi0 = lambda s: np.tanh(s / w)
    dphi0 = lambda s: 1.0 / (w * np.cosh(s / w) ** 2)
    return w, phi0, dphi0


def overlap_exact(v, w):
    b = np.asarray(v, float) / w
    out = np.full_like(b, 4.0 / (3.0 * w))
    nz = np.abs(b) > 1e-6
    bb = b[nz]
    out[nz] = (4.0 / w) * (bb / np.tanh(bb) - 1.0) / np.sinh(bb) ** 2
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eps0", type=float, default=1.0)
    ap.add_argument("--c", type=float, default=1.0)
    a = ap.parse_args()
    eps0, c = a.eps0, a.c
    w, phi0, dphi0 = kink(eps0, c)
    lim = 40 * w

    # 1. constants
    N1 = quad(lambda s: dphi0(s) ** 2, -lim, lim)[0]
    N3 = quad(dphi0, -lim, lim)[0]
    sigma = c * N1
    print(f"[1] w = {w:.6f}")
    print(f"    N1 = {N1:.6f}   4/(3w) = {4 / (3 * w):.6f}")
    print(f"    int phi0' = {N3:.6f} (expect 2)")
    print(f"    sigma = c N1 = {sigma:.6f}   (2 sqrt2/3) sqrt(c eps0) = {2 * np.sqrt(2) / 3 * np.sqrt(c * eps0):.6f}")
    # surface tension from the energy of the kink
    e_kink = quad(lambda s: 0.5 * c * dphi0(s) ** 2 + 0.25 * eps0 * (1 - phi0(s) ** 2) ** 2, -lim, lim)[0]
    print(f"    kink energy = {e_kink:.6f} (= sigma)")

    # 2. spectrum
    n, Ls = 4000, 30 * w
    s = np.linspace(-Ls, Ls, n)
    ds = s[1] - s[0]
    diag = 2 * c / ds ** 2 + eps0 * (3 * phi0(s) ** 2 - 1)
    off = -c / ds ** 2 * np.ones(n - 1)
    ev = eigh_tridiagonal(diag, off, select="i", select_range=(0, 3), eigvals_only=True)
    print(f"[2] lowest eigenvalues of L: {np.round(ev, 5)}")
    print(f"    expect 0, {1.5 * eps0}, then continuum from {2 * eps0}")

    # 3. random-field overlap
    vs = np.array([0.0, 0.3, 1.0, 2.0, 4.0]) * w
    num = np.array([quad(lambda y: dphi0(y) * dphi0(y - v), -lim, lim + v)[0] for v in vs])
    print("[3] RF overlap C(v)/C_exact(v):", np.round(num / overlap_exact(vs, w), 6))

    # 4. random-bond: R(v) ∝ int sech^2(y/w) sech^2((y-v)/w) dy = w^2 C(v); Delta = -R''
    Rv = lambda v: quad(lambda y: (1 - phi0(y) ** 2) * (1 - phi0(y - v) ** 2), -lim, lim + abs(v))[0]
    Rnum = np.array([Rv(v) for v in vs])
    print("[4] RB  int sech^2 sech^2 / (w^2 C_exact):", np.round(Rnum / (w ** 2 * overlap_exact(vs, w)), 6))
    h = 1e-2 * w
    D0 = -(Rv(h) - 2 * Rv(0) + Rv(-h)) / h ** 2
    Rp = lambda v: (Rv(v + h) - Rv(v - h)) / (2 * h)
    Dint = -(Rp(12 * w) - Rp(-12 * w))  # int -R'' dv = -[R']
    print(f"    -R''(0) = {D0:.5f} > 0,  int -R'' dv = {Dint:.2e} (random bond: 0)")

    # 5. Fourier transform
    for q in [0.5, 1.0, 2.0]:
        ft = quad(lambda y: np.cos(q * y) * dphi0(y), -lim, lim, limit=400)[0]
        x = np.pi * q * w
        print(f"[5] q = {q}: FT = {ft:.6f}   pi q w/sinh(pi q w/2) = {x / np.sinh(x / 2):.6f}")
    q = 2 * np.pi
    x = np.pi * q * w
    print(f"    lattice harmonic q = 2pi: |FT|^2 = {(x / np.sinh(x / 2)) ** 2:.2e} (Peierls-Nabarro scale)")

    # 6. bubble projection: residual E = c lap(phi) + eps0 (phi - phi^3) on phi = tanh((R - rho)/w)
    #    eta Rdot = (1/R) int rho drho phi0'(R-rho) E / N1 ... we compare projected E with -sigma/R
    print("[6] bubble projection (h = 0):  P[E]/(N1)  vs  -sigma/R  (velocity = P/eta)")
    for R in [5.0, 10.0, 20.0, 40.0]:
        def integrand(rho):
            z = (R - rho) / w
            ph = np.tanh(z)
            d1 = -1.0 / (w * np.cosh(z) ** 2)           # d phi / d rho
            d2 = -2.0 * ph / (w ** 2 * np.cosh(z) ** 2)  # d2 phi / d rho2
            E = c * (d2 + d1 / rho) + eps0 * (ph - ph ** 3)
            return rho * (1.0 / (w * np.cosh(z) ** 2)) * E
        P = quad(integrand, max(1e-6, R - 20 * w), R + 20 * w, limit=400)[0] / R
        Rdot = P / N1  # eta = N1 with unit mobility; projection of dphi/dt = phi0' Rdot gives N1 Rdot
        print(f"    R = {R:5.1f}: eta Rdot = {N1 * Rdot:+.6f}   -sigma/R = {-sigma / R:+.6f}")


if __name__ == "__main__":
    main()
