"""Fractional moment + change of measure: rigorous upper bounds on the quenched
pressure phi_q(tau) and on f_c for the triangle model (P = 4, slope k = A).

Hoelder with any positive g(env):  <ln Z> <= (1/m) ln <Z^m>
    <= ((1-m)/m) ln <g^{-m/(1-m)}> + ln <g Z>.
For g = dQ/dP the best choice is m -> 0, i.e.  <ln Z>_P <= ln <Z>_Q + D(P||Q).

Environment: phases phi_i iid uniform  <=>  second differences
V_i = phi_i - 2 phi_{i-1} + phi_{i-2} iid uniform on [0, P).
Q: V_i iid with density h(v)/P,  D(P||Q) per site = -E_unif ln h.
A strand's position relative to the phase, theta_i = u_i - phi_i, obeys
    theta_{i+1} = 2 theta_i - theta_{i-1} + g(theta_i) - V_{i+1},
so <Z>_Q grows as the top eigenvalue Lambda_Q(tau) of an operator on
(theta_i, theta_{i-1}, tangent direction):
    phi_q(tau) <= B_h(tau) = Lambda_Q(tau) - E_unif ln h,
    f_c <= inf_tau B_h(tau) / tau.
h = 1 gives the annealed Lambda(tau) and f_ann.

    python3 fracbound.py           -> fracbound.json
    python3 fracbound.py --check   (h = 1 against the exact annealed rate)
"""
import json
import os
import numpy as np
from numba import njit, prange
from scipy.optimize import minimize

P = 4.0


def g_tri(x, A):
    y = np.mod(x, P) / (P / 2)
    n = np.floor(y)
    t = y - n
    even = (n % 2) == 0
    return np.where(even, -A + 2 * A * t, A - 2 * A * t), np.where(even, 1.0, -1.0)


class Op:
    def __init__(self, A, n=96, na=48):
        self.A, self.n, self.na = A, n, na
        k = A
        self.th = (np.arange(n) + 0.5) * P / n
        self.g, self.s = g_tri(self.th, A)
        # tangent (psi_i, psi_{i-1}) -> ((2 + k s) psi_i - psi_{i-1}, psi_i), angle grid on [0, pi)
        al = (np.arange(na) + 0.5) * np.pi / na
        self.aidx = np.zeros((2, na), np.int64); self.aw = np.zeros((2, na)); self.ar = np.zeros((2, na))
        for si, s in enumerate((1.0, -1.0)):
            x = (2 + k * s) * np.cos(al) - np.sin(al); y = np.cos(al)
            r = np.hypot(x, y); a2 = np.mod(np.arctan2(y, x), np.pi)
            pos = a2 / (np.pi / na) - 0.5
            i0 = np.floor(pos).astype(np.int64); w = pos - i0
            self.aidx[si] = i0; self.aw[si] = w; self.ar[si] = r
        self.sidx = np.where(self.s > 0, 0, 1).astype(np.int64)
        # c(i1, i2) = 2 th_i1 - th_i2 + g(th_i1): mean position of theta_{i+1} before subtracting V
        self.c = (2 * self.th[:, None] - self.th[None, :] + self.g[:, None])

    def kernel(self, logh):
        """H[i1, i2, j] = q(c - th_j) * dth, q = h / P with h = exp(logh(v)) normalised."""
        n = self.n
        v = np.mod(self.c[:, :, None] - self.th[None, None, :], P)
        hv = np.exp(logh(v))
        vg = (np.arange(4000) + 0.5) * P / 4000
        norm = np.exp(logh(vg)).mean()
        return hv / norm / n, float(np.log(norm) - logh(vg).mean())

    def lam(self, logh, tau, iters=400, tol=1e-10):
        H, D = self.kernel(logh)
        lam = _power(H, np.exp(tau * (self.g - self.A)), self.sidx, self.aidx, self.aw, self.ar, iters, tol)
        return lam + tau * self.A, D


@njit(parallel=True, cache=True)
def _apply(F, H, eg, sidx, aidx, aw, ar):
    n = F.shape[0]; na = F.shape[2]
    out = np.empty_like(F)
    for i1 in prange(n):
        si = sidx[i1]
        # G[j, l] = F(j, i1, angle'(l))
        G = np.empty((n, na))
        for l in range(na):
            i0 = aidx[si, l]; w = aw[si, l]
            a0 = i0 % na; a1 = (i0 + 1) % na
            for j in range(n):
                G[j, l] = (1 - w) * F[j, i1, a0] + w * F[j, i1, a1]
        for i2 in range(n):
            for l in range(na):
                acc = 0.0
                for j in range(n):
                    acc += H[i1, i2, j] * G[j, l]
                out[i1, i2, l] = eg[i1] * ar[si, l] * acc
    return out


def _power(H, eg, sidx, aidx, aw, ar, iters, tol):
    n = H.shape[0]; na = aw.shape[1]
    F = np.ones((n, n, na))
    lam = 0.0
    for it in range(iters):
        Fn = _apply(F, H, eg, sidx, aidx, aw, ar)
        s = Fn.max()
        Fn /= s
        d = np.abs(Fn - F).max()
        F = Fn
        lam = np.log(s)
        if d < tol:
            break
    return lam


def fourier(params, nmodes):
    a = params[:nmodes]; b = params[nmodes:]
    def logh(v):
        out = np.zeros_like(v, dtype=float)
        for m in range(nmodes):
            x = 2 * np.pi * (m + 1) * v / P
            out = out + a[m] * np.cos(x) + b[m] * np.sin(x)
        return out
    return logh


def corrected(op, logh, tau, cache):
    """B_h(tau) with the grid error of the flat operator removed."""
    from annealed import gle
    key = round(tau, 8)
    if key not in cache:
        cache[key] = op.lam(lambda v: np.zeros_like(v, dtype=float), tau)[0]
    lq, D = op.lam(logh, tau)
    x = tau * op.A
    return lq + D - cache[key] + gle(1.0, op.A) + (np.log(np.sinh(x) / x) if x > 0 else 0.0), D


def optimise(A, nm=4, maxiter=1500):
    op = Op(A); cache = {}
    from annealed import f_ann, gle
    fa = f_ann(A, P, gle(1.0, A))
    # threshold bound: minimise B_h(tau)/tau over (h, tau)
    def obj(z):
        tau = np.exp(z[0]); return corrected(op, fourier(z[1:], nm), tau, cache)[0] / tau
    from scipy.optimize import minimize_scalar
    t0 = minimize_scalar(lambda lt: obj(np.r_[lt, np.zeros(2 * nm)]), bounds=(np.log(0.3 / A), np.log(4 / A)), method="bounded").x
    r = minimize(obj, np.r_[t0, np.zeros(2 * nm)], method="Nelder-Mead",
                 options=dict(maxiter=maxiter, xatol=1e-4, fatol=1e-7, adaptive=True))
    out = dict(A=A, f_ann=fa, f_bound=float(r.fun), tauA=float(np.exp(r.x[0]) * A), params=r.x[1:].tolist(),
               D=corrected(op, fourier(r.x[1:], nm), np.exp(r.x[0]), cache)[1], nfev=int(r.nfev), pressure=[])
    print(f"A={A}: f_ann={fa:.4f}  f_c bound={r.fun:.4f}  at tauA={out['tauA']:.3f}  D={out['D']:.4f}", flush=True)
    for x in (1.0, 2.0, 3.0):
        tau = x / A
        rr = minimize(lambda p: corrected(op, fourier(p, nm), tau, cache)[0], np.zeros(2 * nm), method="Nelder-Mead",
                      options=dict(maxiter=maxiter, xatol=1e-4, fatol=1e-7, adaptive=True))
        lam0 = corrected(op, fourier(np.zeros(2 * nm), nm), tau, cache)[0]
        out["pressure"].append(dict(tauA=x, Lambda=lam0, B=float(rr.fun), params=rr.x.tolist()))
        print(f"   tauA={x}: Lambda={lam0:.4f}  bound={rr.fun:.4f}  (gain {lam0-rr.fun:.4f})", flush=True)
    return out


if __name__ == "__main__":
    import sys
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    if "--check" in sys.argv:
        from annealed import gle
        A = 0.6; op = Op(A)
        flat = lambda v: np.zeros_like(v, dtype=float)
        for x in (0.0, 1.0, 1.6, 3.0):
            tau = x / A
            lg, D = op.lam(flat, tau)
            ex = gle(1.0, A) + (np.log(np.sinh(x) / x) if x > 0 else 0)
            print(f"check A={A} tauA={x}: grid {lg:.5f}  exact {ex:.5f}", flush=True)
        raise SystemExit
    res = json.load(open("fracbound.json")) if os.path.exists("fracbound.json") else {}
    for A in (0.6, 0.3, 1.0):
        if str(A) in res:
            continue
        res[str(A)] = optimise(A)
        json.dump(res, open("fracbound.json", "w"), indent=1)
