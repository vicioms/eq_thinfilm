"""Self-consistent (hot-line) theory of the sliding state.

Closure: the fluctuation du_i = u_i - v t is Gaussian, with the displacement
structure function of a damped lattice string in equilibrium at temperature T,

    B_T(tau) = <(du(t+tau) - du(t))^2> = T b(tau),
    b(tau)   = 2 int dq/2pi (1/lambda_q) [1 - c_q(tau)],

c_q the normalised position autocorrelation of mode q (damped oscillator).
The pinning force is the one-loop (MSR) expression with a Debye-Waller factor,

    F_pin(v, T) = int_0^inf dtau R(tau) int dk/2pi Delta(k) k sin(k v tau) exp(-k^2 T b(tau)/2)
                = - int_0^inf dtau R(tau) Delta_s'(v tau),     s = T b(tau),

R(tau) = int dq/2pi R_q(tau) the local response function, Delta_s the force
correlator smeared by a Gaussian of variance s.  For the interpolated disorder,
Delta(x) = (A^2/3) M4(x) with M4 the cubic B-spline on [-2, 2], so the tau
integral has finite support.  Random bond (A < 0): Delta(x) = (A^2/3)[2 M4(x) -
M4(x - 1) - M4(x + 1)].  At T = 0 this is the second-order result
(linresp.py).  T is fixed by the exact power balance

    T = m <du'^2> = m v F_pin(v, T) / eta,

and the flow curve is f(v) = eta v + F_pin(v, T(v)).  As v -> 0, F_pin ~ v and
T ~ v^2: the closure has no self-sustained noise.
"""
import math
import numpy as np
from numba import njit, prange



@njit(parallel=True, cache=True)
def response_and_b(taus, m, eta, kappa, nq):
    nt = taus.shape[0]
    R = np.zeros(nt)
    b = np.zeros(nt)
    for j in prange(nt):
        t = taus[j]
        r = 0.0
        bb = 0.0
        for iq in range(nq):
            q = np.pi * (iq + 0.5) / nq
            lam = 2.0 * kappa * (1.0 - np.cos(q))
            if m == 0.0:
                e = np.exp(-lam * t / eta)
                r += e / eta
                c = e
            else:
                g = eta / m
                disc = lam / m - 0.25 * g * g
                if disc > 0:
                    W = np.sqrt(disc)
                    e = np.exp(-0.5 * g * t)
                    r += e * np.sin(W * t) / (m * W)
                    c = e * (np.cos(W * t) + 0.5 * g / W * np.sin(W * t))
                else:
                    K = np.sqrt(-disc)
                    if K < 1e-12:
                        e = np.exp(-0.5 * g * t)
                        r += e * t / m
                        c = e * (1.0 + 0.5 * g * t)
                    else:
                        ep = np.exp((-0.5 * g + K) * t)
                        em = np.exp((-0.5 * g - K) * t)
                        r += (ep - em) / (2.0 * m * K)
                        c = 0.5 * (ep + em) + 0.25 * g / K * (ep - em)
            bb += (1.0 - c) / lam
        R[j] = r / nq
        b[j] = 2.0 * bb / nq
    return R, b


@njit(inline='always')
def _dM4(y):
    a = abs(y)
    if a <= 1.0:
        return -2.0 * y + 1.5 * y * a
    if a <= 2.0:
        d = 2.0 - a
        return -0.5 * d * d * (1.0 if y > 0 else -1.0)
    return 0.0


_PIECES = np.array([[-2.0, -1.0, 2.0, 2.0, 0.5],
                    [-1.0, 0.0, 0.0, -2.0, -1.5],
                    [0.0, 1.0, 0.0, -2.0, 1.5],
                    [1.0, 2.0, -2.0, 2.0, -0.5]])   # a, b, p0, p1, p2 of M4'(y) on [a, b]


@njit(cache=True)
def _dDelta_s(x, s, pieces):
    """d/dx of M4 smeared by a Gaussian of variance s (exact: M4' is piecewise quadratic)."""
    if s < 1e-12:
        return _dM4(x)
    sd = np.sqrt(s)
    tot = 0.0
    for p in range(4):
        a, bb, p0, p1, p2 = pieces[p, 0], pieces[p, 1], pieces[p, 2], pieces[p, 3], pieces[p, 4]
        z1 = (x - bb) / sd
        z2 = (x - a) / sd
        if z1 > 9.0 or z2 < -9.0:
            continue
        P1 = 0.5 * (1.0 + math.erf(z1 / np.sqrt(2.0)))
        P2 = 0.5 * (1.0 + math.erf(z2 / np.sqrt(2.0)))
        f1 = np.exp(-0.5 * z1 * z1) / np.sqrt(2.0 * np.pi)
        f2 = np.exp(-0.5 * z2 * z2) / np.sqrt(2.0 * np.pi)
        I0 = P2 - P1
        I1 = f1 - f2
        I2 = I0 + z1 * f1 - z2 * f2
        c0 = p0 + p1 * x + p2 * x * x
        c1 = -p1 * sd - 2.0 * p2 * x * sd
        c2 = p2 * s
        tot += c0 * I0 + c1 * I1 + c2 * I2
    return tot


@njit(parallel=True, cache=True)
def fpin(dtau, R, b, v, T, A2, pieces, rb):
    nt = R.shape[0]
    acc = np.zeros(nt)
    for j in prange(nt):
        t = (j + 0.5) * dtau
        s = T * b[j]
        x = v * t
        if x > 2.0 + rb + 8.0 * np.sqrt(s) + 1e-9:
            continue
        if rb:
            d = 2.0 * _dDelta_s(x, s, pieces) - _dDelta_s(x - 1.0, s, pieces) \
                - _dDelta_s(x + 1.0, s, pieces)
        else:
            d = _dDelta_s(x, s, pieces)
        acc[j] = -R[j] * A2 / 3.0 * d
    return acc.sum() * dtau


class Line:
    def __init__(self, m, A=1.0, eta=1.0, kappa=1.0, nq=2048, dtau=0.05, tau_max=20000.0):
        """A < 0: random-bond disorder of amplitude |A| (as in sim.py)."""
        self.rb = 1 if A < 0 else 0
        A = abs(A)
        self.m, self.A, self.eta, self.kappa = m, A, eta, kappa
        self.dtau = dtau
        taus = (np.arange(int(tau_max / dtau)) + 0.5) * dtau
        self.taus = taus
        self.R, self.b = response_and_b(taus, m, eta, kappa, nq)

    def F(self, v, T):
        return fpin(self.dtau, self.R, self.b, v, T, self.A ** 2, _PIECES, self.rb)

    def T_of_v(self, v, Tmax=100.0, it=50):
        """Largest solution of T = m v F(v, T)/eta (the hot branch), by bisection."""
        if self.m == 0:
            return 0.0
        g = lambda T: self.m * v * self.F(v, T) / self.eta - T
        Ts = np.concatenate([[0.0], np.geomspace(1e-4, Tmax, 50)])
        gs = np.array([g(T) for T in Ts])
        idx = np.where((gs[:-1] > 0) & (gs[1:] <= 0))[0]
        if len(idx) == 0:
            return 0.0
        lo, hi = Ts[idx[-1]], Ts[idx[-1] + 1]
        for _ in range(it):
            mid = 0.5 * (lo + hi)
            if g(mid) > 0:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def flow(self, vs):
        Ts = np.array([self.T_of_v(v) for v in vs])
        Fp = np.array([self.F(v, T) for v, T in zip(vs, Ts)])
        return self.eta * vs + Fp, Ts, Fp
