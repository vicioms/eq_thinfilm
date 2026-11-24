"""Structure factor of the hot sliding line: measurement and second-order prediction.

Prediction (no free parameter): each site feels xi_j(t) = F(u_j(t), j), independent
between sites, with autocorrelation C(t) = Delta_s(t)(v t), the disorder correlator
Gaussian-smeared by the relative displacement of the site, s(t)^2 = T b(t) (closure
b from sc.response_and_b, T the kinetic temperature).  Each lattice mode responds
linearly, D_q(w) = lam_q - m w^2 - i eta w, lam_q = 2 kappa (1 - cos q), so

    S(q) = <|u_q|^2>/L = int dw/2pi  S_xi(w) / |D_q(w)|^2,
    K(q) = <|du_q/dt|^2>/L = int dw/2pi w^2 S_xi(w) / |D_q(w)|^2,

with S_xi(w) = 2 int_0^inf C(t) cos(w t) dt.  Mode temperatures: T_pot(q) = lam_q S(q),
T_kin(q) = m K(q).  For q -> 0 (overdamped): T(q) -> S_xi(0)/2 eta = Delta_hat(0)/(2 eta v)
(random field: A^2/(6 eta v); random bond: 0); underdamped modes: T(q) ~ S_xi(c q)/2 eta.
Hence zeta = 1/2 (EW) at large scales, with a scale-dependent temperature.

Usage:  python3 sq.py run    (simulations; writes data/sq_*.npz)
        python3 sq.py theory (predictions; writes data/sq_theory_*.npz)
"""
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
import config  # noqa: F401  (src/ on the path; run from the inertia folder)
import sys
import time
import numpy as np
from protocols import _run, default_dt

CASES = [  # (label, m, f, A)
    ('rf1024', 1024.0, 0.12, 1.0),
    ('rf256', 256.0, 0.17, 1.0),
    ('rb1024', 1024.0, 0.12, -1.0),
]
L = 2048
SEEDS = np.arange(8) + 52000


def measure(m, f, A, T_relax=3e5, T_meas=6e5, dt_snap=100.0):
    dt = default_dt(m)
    S = len(SEEDS)
    U = np.zeros((S, L)); V = np.full((S, L), f)
    n0 = int(round(T_relax / dt))
    _run(U, V, SEEDS, np.full(S, f), m, n0, 10, A, 1.0, 1.0, dt, 0.0, False)
    ns = int(round(dt_snap / dt)); nsnap = int(T_meas / dt_snap)
    Su = np.zeros(L // 2 + 1); Sv = np.zeros(L // 2 + 1)
    u0 = U.mean(); vk = 0.0
    for _ in range(nsnap):
        _run(U, V, SEEDS, np.full(S, f), m, ns, 2, A, 1.0, 1.0, dt, 0.0, False)
        Su += (np.abs(np.fft.rfft(U - U.mean(1, keepdims=True), axis=1)) ** 2).mean(0) / L
        Sv += (np.abs(np.fft.rfft(V - V.mean(1, keepdims=True), axis=1)) ** 2).mean(0) / L
        vk += ((V - V.mean(1, keepdims=True)) ** 2).mean()
    v = (U.mean() - u0) / (ns * nsnap * dt)
    q = 2 * np.pi * np.arange(L // 2 + 1) / L
    return dict(q=q, S=Su / nsnap, K=Sv / nsnap, v=v, Tkin=m * vk / nsnap, m=m, f=f, A=A, L=L)


# ---------------- theory ----------------
def _M4(x):
    x = np.abs(x)
    return np.where(x < 1, 2 / 3 - x ** 2 + x ** 3 / 2, np.where(x < 2, (2 - x) ** 3 / 6, 0.0))


_GH_x, _GH_w = np.polynomial.hermite_e.hermegauss(40)
_GH_w = _GH_w / _GH_w.sum()


def C_of_t(t, v, s, A):
    """Smeared single-site force correlator at lag t (vectorised in t)."""
    x = v * t[:, None] + s[:, None] * _GH_x[None, :]
    if A > 0:
        D = (A ** 2 / 3) * _M4(x)
    else:  # random bond: F = g(n+1) - g(n) interpolated -> Delta = (A^2/3)(2 M4(x) - M4(x+1) - M4(x-1))
        D = (A ** 2 / 3) * (2 * _M4(x) - _M4(x + 1) - _M4(x - 1))
    return (D * _GH_w[None, :]).sum(1)


def predict(m, v, T, A, q, smear=True, eta=1.0, kappa=1.0):
    from sc import response_and_b
    tmax = (2.0 + 10 * np.sqrt(max(T, 1e-9))) / v * 4
    dt = 0.02 / v
    t = (np.arange(int(tmax / dt)) + 0.5) * dt
    if smear and T > 0:
        _, b = response_and_b(t, m, eta, kappa, 2048)
        s = np.sqrt(T * b)
    else:
        s = np.zeros_like(t)
    C = C_of_t(t, v, s, abs(A) if A > 0 else -abs(A))
    lam = 2 * kappa * (1 - np.cos(q))
    wmax = 4.0 * max(np.sqrt(4 * kappa / m), 2 * np.pi * v)
    w = np.concatenate([np.linspace(0, wmax, 200001)])
    # S_xi(w) = 2 int C cos(w t) dt, by direct quadrature on the w grid in blocks
    Sxi = np.empty_like(w)
    for i0 in range(0, len(w), 2000):
        ww = w[i0:i0 + 2000]
        Sxi[i0:i0 + 2000] = 2 * (C[None, :] * np.cos(ww[:, None] * t[None, :])).sum(1) * dt
    Sq = np.empty_like(q); Kq = np.empty_like(q)
    dw = w[1] - w[0]
    for i, l in enumerate(lam):
        if l == 0:
            Sq[i] = Kq[i] = np.nan
            continue
        den = (l - m * w ** 2) ** 2 + (eta * w) ** 2
        # exact: int dw/2pi 1/|D|^2 = 1/(2 eta lam), int dw/2pi w^2/|D|^2 = 1/(2 eta m);
        # integrate only the smooth remainder S_xi(w) - S_xi(0) numerically
        g = (Sxi - Sxi[0]) / den
        Sq[i] = Sxi[0] / (2 * eta * l) + 2 * np.trapezoid(g, dx=dw) / (2 * np.pi)
        Kq[i] = Sxi[0] / (2 * eta * m) + 2 * np.trapezoid(w ** 2 * g, dx=dw) / (2 * np.pi)
    return Sq, Kq, w, Sxi


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'run'
    if what == 'run':
        for lab, m, f, A in CASES:
            t0 = time.time()
            r = measure(m, f, A)
            np.savez(f'data/sq_{lab}.npz', **r)
            print(time.strftime('%H:%M:%S'), 'sq', lab, 'v', round(r['v'], 4), 'Tkin', round(r['Tkin'], 3),
                  f'{time.time() - t0:.0f}s', flush=True)
    else:
        for lab, m, f, A in CASES:
            d = np.load(f'data/sq_{lab}.npz')
            q = d['q']
            out = {}
            for sm in [False, True]:
                Sq, Kq, w, Sxi = predict(m, float(d['v']), float(d['Tkin']), A, q, smear=sm)
                out['S_s' if sm else 'S_b'] = Sq; out['K_s' if sm else 'K_b'] = Kq
            np.savez(f'data/sq_theory_{lab}.npz', q=q, **out)
            print('theory', lab, 'done', flush=True)


# ---------------- self-consistent mode temperatures ----------------
def _mode_corr(q, t, m, eta=1.0, kappa=1.0):
    """Normalised equal-site correlation c_q(t) of a free damped mode (rows q, cols t)."""
    lam = 2 * kappa * (1 - np.cos(q))[:, None]
    g = eta / m
    disc = lam / m - 0.25 * g * g
    W = np.sqrt(np.abs(disc))
    e = np.exp(-0.5 * g * t[None, :])
    with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
        under = e * (np.cos(W * t) + 0.5 * g / W * np.sin(W * t))
        ep = np.exp((-0.5 * g + W) * t); em = np.exp((-0.5 * g - W) * t)
        over = 0.5 * (ep + em) + 0.25 * g / W * (ep - em)
    return np.where(disc > 0, under, over), lam[:, 0]


def predict_sc(m, v, A, q, Tq0, it=30, mix=0.5, eta=1.0, kappa=1.0, Lsum=None):
    """Iterate T(q) -> B(t) = (2/N) sum_q T(q)(1 - c_q(t))/lam_q -> C(t) -> S_xi -> T(q).

    The sum runs over the lattice modes q = 2 pi k / Lsum, k = 1..Lsum/2 (the simulated line).
    """
    Lsum = Lsum or 2 * (len(q) - 1)
    qs = 2 * np.pi * np.arange(1, Lsum // 2 + 1) / Lsum
    tmax = 4000.0 / v ** 0 if False else 60.0 / v
    dt = 0.02 / v
    t = (np.arange(int(tmax / dt)) + 0.5) * dt
    cq, lam_s = _mode_corr(qs, t, m, eta, kappa)
    Tq = np.interp(qs, q[1:], Tq0[1:])
    lam = 2 * kappa * (1 - np.cos(q))
    wmax = 4.0 * max(np.sqrt(4 * kappa / m), 2 * np.pi * v)
    w = np.linspace(0, wmax, 100001); dw = w[1] - w[0]
    hist = []
    for k in range(it):
        B = (2.0 / Lsum) * ((Tq / lam_s)[:, None] * (1 - cq)).sum(0)
        s = np.sqrt(np.maximum(B, 0))
        C = C_of_t(t, v, s, A)
        Sxi = np.empty_like(w)
        for i0 in range(0, len(w), 2000):
            ww = w[i0:i0 + 2000]
            Sxi[i0:i0 + 2000] = 2 * (C[None, :] * np.cos(ww[:, None] * t[None, :])).sum(1) * dt
        Tn = np.empty_like(q); Tk = np.empty_like(q)
        for i, l in enumerate(lam):
            if l == 0:
                Tn[i] = Tk[i] = np.nan; continue
            den = (l - m * w ** 2) ** 2 + (eta * w) ** 2
            g = (Sxi - Sxi[0]) / den
            Tn[i] = l * (Sxi[0] / (2 * eta * l) + 2 * np.trapezoid(g, dx=dw) / (2 * np.pi))
            Tk[i] = m * (Sxi[0] / (2 * eta * m) + 2 * np.trapezoid(w ** 2 * g, dx=dw) / (2 * np.pi))
        Tnew = np.interp(qs, q[1:], Tn[1:])
        hist.append(Tn[1])
        if np.max(np.abs(Tnew - Tq) / Tq) < 1e-3:
            Tq = Tnew; break
        Tq = mix * Tnew + (1 - mix) * Tq
    return Tn, Tk, np.array(hist), t, B
