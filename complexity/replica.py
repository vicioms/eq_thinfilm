"""Replicated Kac-Rice: moments <Z_K(tau)^n> of the counting partition function
Z_K(tau) = int_0^P da e^{tau W_K(a)} |W_K'(a)| = int dx N_K(x) e^{tau x}
for the triangle model (P = 4, amplitude A, slope k = 4A/P = A).

n replicas a_1..a_n shoot through the SAME environment (phases phi_i).  With
theta_i = u_i - phi_i uniform and fresh at every site, the forces seen by the
replicas depend only on their offsets s^b = u^b - u^1 mod P, so

    <Z^n>  =  P  E[ prod_i  prod_b e^{tau g(theta_i + s_i^b)} |psi^b growth| ]

is a Feynman-Kac average over the Markov chain of the offsets (s^b, d^b),
d = s_i - s_{i-1}, driven by theta_i, and the n tangent vectors psi^b
(psi_{i+1} = (2 + V''(u_i)) psi_i - psi_{i-1}, |det H| = |psi_K|).
Lambda_n(tau) = lim K^-1 ln <Z^n> is its top eigenvalue, estimated here by
population dynamics (systematic resampling, N walkers).

n = 1 reproduces the exact annealed Lambda(tau) = Lambda(1) + ln(sinh(tau A)/(tau A)).

    python3 replica.py   -> replica.json
"""
import json
import os
import numpy as np
from numba import njit, prange

P = 4.0


@njit(cache=True, inline="always")
def _g(x, A):
    """Triangle force V'(x) (period P, values in [-A, A]) and slope sign."""
    h = 0.5 * P
    y = x / h
    n = np.floor(y)
    t = y - n
    if (int(n) & 1) == 0:
        return -A + 2.0 * A * t, 1.0
    return A - 2.0 * A * t, -1.0


@njit(parallel=True, cache=True)
def _step(s, d, p, q, A, tau, k, nrep, theta, w):
    """One site for all walkers; p, q: tangent (psi_i, psi_{i-1}) per replica."""
    N = s.shape[0]
    for j in prange(N):
        lw = 0.0
        g0 = 0.0
        for b in range(nrep):
            off = 0.0 if b == 0 else s[j, b]
            gb, sg = _g(theta[j] + off, A)
            if b == 0:
                g0 = gb
            lw += tau * (gb - A)  # shifted: weight <= 1 per replica
            a_ = (2.0 + k * sg) * p[j, b] - q[j, b]
            q[j, b] = p[j, b]
            p[j, b] = a_
            nrm = np.sqrt(a_ * a_ + q[j, b] * q[j, b])
            p[j, b] /= nrm
            q[j, b] /= nrm
            lw += np.log(nrm)
            if b > 0:
                # offset dynamics: s_{i+1} = 2 s_i - s_{i-1} + (U^b - U^1)
                d[j, b] = (d[j, b] + gb - g0) % P
                s[j, b] = (s[j, b] + d[j, b]) % P
        w[j] = lw


def lambda_n(A, tau, nrep, N=200_000, T=3000, burn=600, seed=0):
    rng = np.random.default_rng(seed)
    k = 4.0 * A / P
    s = rng.uniform(0, P, (N, nrep))
    d = np.zeros((N, nrep))
    ang = rng.uniform(0, np.pi, (N, nrep))
    p = np.cos(ang); q = np.sin(ang)
    w = np.empty(N)
    acc = []
    for t in range(T):
        theta = rng.uniform(0, P, N)
        _step(s, d, p, q, A, tau, k, nrep, theta, w)
        m = w.max()
        ww = np.exp(w - m)
        mean = ww.mean()
        if t >= burn:
            acc.append(np.log(mean) + m + nrep * tau * A)
        c = np.cumsum(ww); c /= c[-1]
        idx = np.minimum(np.searchsorted(c, (rng.random() + np.arange(N)) / N), N - 1)
        s, d, p, q = s[idx], d[idx], p[idx], q[idx]
    acc = np.array(acc)
    blocks = acc[: len(acc) // 20 * 20].reshape(20, -1).mean(1)
    return float(acc.mean()), float(blocks.std() / np.sqrt(20))


def lam_exact(A, tau):
    from annealed import gle
    x = tau * A
    return gle(1.0, A) + (np.log(np.sinh(x) / x) if x > 1e-8 else 0.0)


if __name__ == "__main__":
    import sys
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    out = json.load(open("replica.json")) if os.path.exists("replica.json") else {}
    if "--check" in sys.argv:
        for A in (0.6,):
            for tau in (0.0, 1.67, 5.0):
                l1 = lambda_n(A, tau, 1, N=100_000, T=1500, burn=300)
                print(f"A={A} tau={tau}: pop {l1[0]:.4f}+-{l1[1]:.4f}  exact {lam_exact(A, tau):.4f}", flush=True)
        raise SystemExit
    for A in (0.3, 0.6, 1.0):
        rows = out.get(str(A), [])
        done = {round(r["tauA"], 3) for r in rows}
        for tauA in ((0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0) if A == 0.3 else (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0)):
            if round(tauA, 3) in done:
                continue
            tau = tauA / A
            l2, e2 = lambda_n(A, tau, 2)
            L1 = lam_exact(A, tau); L2t = lam_exact(A, 2 * tau)
            rows.append(dict(tauA=tauA, tau=tau, Lambda2=l2, Lambda2_err=e2, Lambda=L1, Lambda_2tau=L2t))
            print(f"A={A} tauA={tauA}: Lambda2={l2:.4f}+-{e2:.4f}  2Lambda={2*L1:.4f}  Lambda(2tau)={L2t:.4f}  "
                  f"excess={(l2-2*L1):.4f}", flush=True)
            out[str(A)] = rows
            json.dump(out, open("replica.json", "w"), indent=1)
