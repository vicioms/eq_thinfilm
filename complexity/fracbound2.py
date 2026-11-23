"""Change-of-measure bound with a Markov tilt of the phase second differences.

Q: V_{i+1} | V_i has density h(V_i, V_{i+1}) / P, h(v, .) normalised for every v.
D(P||Q) per site = -E[ln h(V, V')], V, V' iid uniform.
V_i is fixed by the strand's last three relative positions,
    V_i = 2 theta_{i-1} - theta_{i-2} + g(theta_{i-1}) - theta_i,
so <Z>_Q grows as the top eigenvalue of an operator on
(theta_i, theta_{i-1}, theta_{i-2}, tangent direction).
phi_q(tau) <= Lambda_Q(tau) + D,  f_c <= inf_tau [Lambda_Q(tau) + D] / tau.

    python3 fracbound2.py   -> fracbound2.json
"""
import json
import os
import numpy as np
from numba import njit, prange
from scipy.optimize import minimize
from fracbound import g_tri, P


class Op2:
    def __init__(self, A, n=48, na=24):
        self.A, self.n, self.na = A, n, na
        k = A
        self.th = (np.arange(n) + 0.5) * P / n
        self.g, self.s = g_tri(self.th, A)
        al = (np.arange(na) + 0.5) * np.pi / na
        self.aidx = np.zeros((2, na), np.int64); self.aw = np.zeros((2, na)); self.ar = np.zeros((2, na))
        for si, s in enumerate((1.0, -1.0)):
            x = (2 + k * s) * np.cos(al) - np.sin(al); y = np.cos(al)
            r = np.hypot(x, y); a2 = np.mod(np.arctan2(y, x), np.pi)
            pos = a2 / (np.pi / na) - 0.5
            i0 = np.floor(pos).astype(np.int64)
            self.aidx[si] = i0; self.aw[si] = pos - i0; self.ar[si] = r
        self.sidx = np.where(self.s > 0, 0, 1).astype(np.int64)
        th = self.th
        self.c = np.mod(2 * th[:, None] - th[None, :] + self.g[:, None], P)       # c(i0, i1)
        # V_i(i0, i1, i2) = c(i1, i2) - theta_i0 ; V_{i+1}(j; i0, i1) = c(i0, i1) - theta_j
        self.Vprev = np.mod(self.c[None, :, :] - th[:, None, None], P)            # [i0, i1, i2]
        self.Vnext = np.mod(self.c[:, :, None] - th[None, None, :], P)            # [i0, i1, j]

    def lam(self, logh, tau, iters=300, tol=1e-9):
        """logh(v, v') unnormalised; returns (Lambda_Q, D)."""
        n = self.n
        vg = (np.arange(256) + 0.5) * P / 256
        L = logh(vg[:, None], vg[None, :])
        lnorm = np.log(np.exp(L).mean(1))                      # per v
        D = float(-(L - lnorm[:, None]).mean())
        lnorm_f = lambda v: np.interp(v, np.r_[vg - P, vg, vg + P], np.r_[lnorm, lnorm, lnorm])
        # H[i0, i1, i2, j] = h_norm(Vprev[i0,i1,i2], Vnext[i0,i1,j]) / n
        vp = self.Vprev[:, :, :, None]; vn = self.Vnext[:, :, None, :]
        H = np.exp(logh(vp, vn) - lnorm_f(vp)) / n
        eg = np.exp(tau * (self.g - self.A))
        F = np.ones((n, n, n, self.na)); lam = 0.0
        for _ in range(iters):
            Fn = _apply2(F, H, eg, self.sidx, self.aidx, self.aw, self.ar)
            s = Fn.max(); Fn /= s
            d = np.abs(Fn - F).max(); F = Fn; lam = np.log(s)
            if d < tol:
                break
        return lam + tau * self.A, D


@njit(parallel=True, cache=True)
def _apply2(F, H, eg, sidx, aidx, aw, ar):
    n = F.shape[0]; na = F.shape[3]
    out = np.empty_like(F)
    for i0 in prange(n):
        si = sidx[i0]
        for i1 in range(n):
            # G[j, l] = F(j, i0, i1, angle'(l))
            G = np.empty((n, na))
            for l in range(na):
                a0 = aidx[si, l] % na; a1 = (aidx[si, l] + 1) % na; w = aw[si, l]
                for j in range(n):
                    G[j, l] = (1 - w) * F[j, i0, i1, a0] + w * F[j, i0, i1, a1]
            for i2 in range(n):
                for l in range(na):
                    acc = 0.0
                    for j in range(n):
                        acc += H[i0, i1, i2, j] * G[j, l]
                    out[i0, i1, i2, l] = eg[i0] * ar[si, l] * acc
    return out


def basis(params, nm):
    """ln h(v, v') = sum over Fourier products up to order nm in each argument."""
    c = np.asarray(params).reshape(2 * nm + 1, 2 * nm + 1)
    def f(m, v):
        if m == 0:
            return np.ones_like(v)
        q = (m + 1) // 2
        x = 2 * np.pi * q * v / P
        return np.cos(x) if m % 2 else np.sin(x)
    def logh(v, w):
        out = 0.0
        for a in range(2 * nm + 1):
            fa = f(a, v)
            for b in range(1, 2 * nm + 1):
                if c[a, b] != 0.0:
                    out = out + c[a, b] * fa * f(b, w)
        return out + 0.0 * v * w
    return logh


if __name__ == "__main__":
    import sys
    from annealed import gle, f_ann
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    nm = 2
    npar = (2 * nm + 1) ** 2
    out = json.load(open("fracbound2.json")) if os.path.exists("fracbound2.json") else {}
    for A in (0.6, 0.3, 1.0):
        if str(A) in out:
            continue
        op = Op2(A)
        flat = basis(np.zeros(npar), nm)
        cache = {}
        def B(p, tau):
            key = round(tau, 8)
            if key not in cache:
                cache[key] = op.lam(flat, tau)[0]
            lq, D = op.lam(basis(p, nm), tau)
            x = tau * A
            return lq + D - cache[key] + gle(1.0, A) + np.log(np.sinh(x) / x)
        fa = f_ann(A, P, gle(1.0, A))
        prev = json.load(open("fracbound.json")).get(str(A)) if os.path.exists("fracbound.json") else None
        tau = (prev["tauA"] if prev else 1.6) / A
        p0 = np.zeros(npar)
        if prev:  # start from the optimal iid tilt: c[0, b] terms
            c0 = np.zeros((2 * nm + 1, 2 * nm + 1)); pr = prev["params"]; m4 = len(pr) // 2
            for q in range(1, nm + 1):
                if q <= m4:
                    c0[0, 2 * q - 1] = pr[q - 1]; c0[0, 2 * q] = pr[m4 + q - 1]
            p0 = c0.ravel()
        print(f"A={A}: start {B(p0, tau)/tau:.4f} (f_ann {fa:.4f})", flush=True)
        r = minimize(lambda p: B(p, tau) / tau, p0, method="Powell",
                     options=dict(maxiter=4000, xtol=1e-3, ftol=1e-7))
        out[str(A)] = dict(A=A, f_ann=fa, tauA=tau * A, f_bound=float(r.fun), params=r.x.tolist(), nfev=int(r.nfev))
        print(f"A={A}: Markov-tilt f_c bound {r.fun:.4f} at tauA={tau*A:.3f}  (nfev {r.nfev})", flush=True)
        json.dump(out, open("fracbound2.json", "w"), indent=1)
