"""Dynamical cavity (single-site) theory of the hot sliding line.

The disorder vertex of the MSRJD action is local in x, so one site can be solved exactly
in a self-consistent Gaussian bath formed by the rest of the line.  Site 0 carries its own
random-field force F_0(u) (the same landscape as sim.py).  Its neighbours form a linear
"medium": a ring of N - 1 sites, each with friction eta_m = eta + gamma and an independent
Gaussian force noise xi_j(t) of autocorrelation C(t), drifting at the imposed velocity v:

    m u_j'' + eta_m u_j' = kappa Lap u_j + eta_m v + xi_j(t)              (medium, j != 0)
    m u_0'' + eta   u_0' = kappa Lap u_0 + F_0(u_0) + f_0                 (cavity site)

Self-consistency (fixed point of the map):
  * C(t)  = autocorrelation of the force F_0(u_0(t)) felt by site 0 (the medium noise is the
            disorder force of the other sites, all statistically alike);
  * gamma = d f_0 / d v - eta at fixed bath: the medium's disorder friction is the slope of
            the single-site pinning force with velocity, at fixed noise and friction (eta_eff);
  * f_0   = eta v - <F_0> keeps site 0 moving at v; it is the predicted force, f(v).
Outputs: f(v), S_xi(0) = int C dt (two-sided, = 2 int_0^inf C), eta_eff, T_lw = S_xi(0)/2 eta_eff,
C(t) and B(t) of site 0 (B relative to the ring mean, as in sq_diag.py), T_kin of site 0.

Noise: each stream is a sum of overlapping blocks, each block an exact circulant-embedding
Gaussian sample of C, cross-faded with sin/cos windows so the variance is constant.

Usage: python3 cavity.py m v [v ...]     (writes data/cavity_m{m}_v{v}.npz)
"""
import sys
import time
import numpy as np
import config  # noqa: F401  (run from the inertia folder)
from numba import njit, prange
from sim import _g


@njit(parallel=True, cache=True)
def ring_chunk(U, V, cn, g0, g1, seeds, f0, v, m, eta, etam, kappa, A, dt, xi, Frec, Urec, Mrec):
    S, N = U.shape
    H = xi.shape[1]
    e = np.exp(-eta * dt / m); em = np.exp(-etam * dt / m)
    for s in prange(S):
        u = U[s]; w = V[s]
        Ft = np.empty(N)
        for step in range(H):
            u0 = u[0]
            n = int(np.floor(u0))
            if n != cn[s]:
                if n == cn[s] + 1:
                    g0[s] = g1[s]; g1[s] = _g(seeds[s], 0, n + 1, A, 0)
                elif n == cn[s] - 1:
                    g1[s] = g0[s]; g0[s] = _g(seeds[s], 0, n, A, 0)
                else:
                    g0[s] = _g(seeds[s], 0, n, A, 0); g1[s] = _g(seeds[s], 0, n + 1, A, 0)
                cn[s] = n
            x = u0 - n
            Fd = g0[s] * (1.0 - x) + g1[s] * x
            for i in range(N):
                im = i - 1 if i > 0 else N - 1
                ip = i + 1 if i < N - 1 else 0
                lap = kappa * (u[ip] + u[im] - 2.0 * u[i])
                if i == 0:
                    Ft[i] = lap + Fd + f0[s]
                else:
                    Ft[i] = lap + etam * v + xi[s, step, i - 1]
            su = 0.0
            for i in range(N):
                if i == 0:
                    w[i] = w[i] * e + (Ft[i] / eta) * (1.0 - e)
                else:
                    w[i] = w[i] * em + (Ft[i] / etam) * (1.0 - em)
                u[i] += dt * w[i]
                su += u[i]
            Frec[s, step] = Fd
            Urec[s, step] = u[0]
            Mrec[s, step] = su / N


class Noise:
    """Stationary Gaussian streams with autocorrelation C (sampled at the step dt)."""

    def __init__(self, c, nstream, H, rng):
        # c: covariance at lags 0..H (steps); block length 2H, circulant embedding
        B = 2 * H
        emb = np.concatenate([c[:H + 1], c[1:H][::-1]])
        lam = np.fft.fft(emb).real
        self.sq = np.sqrt(np.clip(lam, 0, None) / B)
        self.H, self.n, self.rng = H, nstream, rng
        self.prev = self._block()
        self.win = np.sin(0.5 * np.pi * np.arange(H) / H)

    def _block(self):
        Z = self.rng.standard_normal((self.n, 2 * self.H)) + 1j * self.rng.standard_normal((self.n, 2 * self.H))
        return np.fft.fft(self.sq[None, :] * Z, axis=1).real

    def next(self):
        cur = self._block()
        out = self.win * cur[:, :self.H] + self.win[::-1] * self.prev[:, self.H:]
        self.prev = cur
        return out


def run(m, v, C_t, tC, gamma, f0, S=16, N=64, T=2e5, T0=2e4, dt=0.4, H=2048, seed=0, A=1.0,
        eta=1.0, kappa=1.0, adapt=True):
    """One pass at fixed bath (C, gamma).  Returns the record and the adapted f0."""
    rng = np.random.default_rng(seed)
    c = np.interp(np.arange(H + 1) * dt, tC, C_t, right=0.0)
    noise = Noise(c, S * (N - 1), H, rng)
    U = np.zeros((S, N)) + v * 0.0; V = np.full((S, N), v)
    U[:, 0] = 0.0
    cn = np.full(S, -10 ** 9, np.int64); g0 = np.zeros(S); g1 = np.zeros(S)
    seeds = np.arange(S, dtype=np.int64) + 77000
    f0s = np.full(S, f0)
    nch = int(np.ceil(T / (H * dt))); n0 = int(np.ceil(T0 / (H * dt)))
    Fr, Ur, Mr = [], [], []
    Fsum, nF = 0.0, 0
    for k in range(nch):
        xi = noise.next().reshape(S, N - 1, H).transpose(0, 2, 1).copy()
        Fc = np.zeros((S, H)); Uc = np.zeros((S, H)); Mc = np.zeros((S, H))
        ring_chunk(U, V, cn, g0, g1, seeds, f0s, v, m, eta, eta + gamma, kappa, A, dt, xi, Fc, Uc, Mc)
        if adapt:   # keep site 0 at velocity v: f0 = eta v - <F_0> (running mean, after a short start)
            Fsum += Fc.sum(); nF += Fc.size
            f0s[:] = eta * v - Fsum / nF
        if k >= n0:
            Fr.append(Fc); Ur.append(Uc); Mr.append(Mc)
    Fr = np.concatenate(Fr, 1); Ur = np.concatenate(Ur, 1); Mr = np.concatenate(Mr, 1)
    return Fr, Ur, Mr, float(f0s[0])


def analyse(Fr, Ur, Mr, dt, v, m, nl=1500, every=5):
    S, n = Fr.shape
    F0 = Fr - Fr.mean()
    Fk = np.fft.rfft(F0, n=2 * n, axis=1)
    ac = np.fft.irfft(np.abs(Fk) ** 2, axis=1)[:, :nl * every].mean(0) / (n - np.arange(nl * every))
    tt = dt * np.arange(nl * every)
    Sxi0 = 2 * np.trapezoid(ac, tt)
    du = Ur - Mr
    lags = np.arange(0, nl * every, every)
    B = np.array([((du[:, l:] - du[:, :n - l]) ** 2).mean() for l in lags])
    vel = np.diff(Ur, axis=1) / dt
    vsite = (Ur[:, -1] - Ur[:, 0]).mean() / ((n - 1) * dt)
    Tkin = m * vel.var()
    return dict(t=tt, C=ac, tB=lags * dt, B=B, Sxi0=float(Sxi0), Fmean=float(Fr.mean()),
                vsite=float(vsite), Tkin=float(Tkin))


def solve(m, v, it=8, gamma0=-0.35, eps=0.15, mix=0.6, **kw):
    """Iterate C and gamma to self-consistency at fixed v."""
    dt = kw.get('dt', 0.4)
    tC = np.arange(0, 4000, dt)
    x = v * tC   # bare start: C = Delta(v t) = (A^2/3) M4(v t)
    xa = np.abs(x)
    C = np.where(xa < 1, 2 / 3 - xa ** 2 + xa ** 3 / 2, np.where(xa < 2, (2 - xa) ** 3 / 6, 0.0)) / 3
    gamma = gamma0; f0 = 1.0 * v + 0.04
    hist = []
    for k in range(it):
        t0 = time.time()
        Fr, Ur, Mr, f0 = run(m, v, C, tC, gamma, f0, seed=k, **kw)
        r = analyse(Fr, Ur, Mr, dt, v, m)
        # eta_eff at fixed bath: same noise and disorder, velocity v (1 + eps)
        _, _, _, f1 = run(m, v * (1 + eps), C, tC, gamma, f0 + eps * v * (1 + gamma), seed=k, **kw)
        eta_eff = (f1 - f0) / (eps * v)
        Cn = np.interp(tC, r['t'], r['C'], right=0.0)
        hist.append(dict(it=k, f=f0, Sxi0=r['Sxi0'], eta_eff=eta_eff, gamma=gamma, Tkin=r['Tkin'],
                         vsite=r['vsite']))
        print(time.strftime('%H:%M:%S'), f'm {m:g} v {v} it {k}: f {f0:.4f} S_xi(0) {r["Sxi0"]:.3f} '
              f'eta_eff {eta_eff:.3f} gamma {gamma:.3f} T_lw {r["Sxi0"] / (2 * eta_eff):.3f} '
              f'Tkin {r["Tkin"]:.3f} v_site {r["vsite"]:.4f}  {time.time() - t0:.0f}s', flush=True)
        C = mix * Cn + (1 - mix) * C
        gamma = max(-0.8, mix * (eta_eff - 1.0) + (1 - mix) * gamma)
    return r, hist, C, tC, gamma, f0


if __name__ == '__main__':
    m = float(sys.argv[1])
    for v in [float(x) for x in sys.argv[2:]]:
        r, hist, C, tC, gamma, f0 = solve(m, v)
        np.savez(f'data/cavity_m{m:g}_v{v:g}.npz', hist=np.array([list(h.values()) for h in hist]),
                 hist_keys=list(hist[0].keys()), C_final=C, tC=tC, gamma=gamma, f=f0, **r)


def polish(m, v, C, tC, gamma, f0, eps=0.15, S=32, N=128, T=2e5, T0=2e4, seed=99, dt=0.4):
    """Final pass at the converged bath: f, S_xi(0), C, B, T_kin at v, and eta_eff by a
    central difference f(v(1 +- eps)) with common random numbers."""
    kw = dict(S=S, N=N, T=T, T0=T0, dt=dt, seed=seed)
    Fr, Ur, Mr, fc = run(m, v, C, tC, gamma, f0, **kw)
    r = analyse(Fr, Ur, Mr, dt, v, m)
    del Fr, Ur, Mr
    fp = run(m, v * (1 + eps), C, tC, gamma, f0 + eps * v * (1 + gamma), **kw)[3]
    fm = run(m, v * (1 - eps), C, tC, gamma, f0 - eps * v * (1 + gamma), **kw)[3]
    r.update(f=fc, f_plus=fp, f_minus=fm, eta_eff=(fp - fm) / (2 * eps * v), gamma=gamma)
    return r
