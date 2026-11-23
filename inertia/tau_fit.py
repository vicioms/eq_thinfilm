"""Maximum-likelihood fit of P(s) ~ s^-tau exp(-s/s_c) for s >= s_lo, with bootstrap errors."""
import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaincc, gamma


def _norm(tau, sc, slo):
    # int_slo^inf s^-tau e^{-s/sc} ds = sc^(1-tau) Gamma(1-tau, slo/sc); for tau > 1 use recursion
    a = 1.0 - tau
    x = slo / sc
    if a > 0:
        return sc ** a * gamma(a) * gammaincc(a, x)
    # Gamma(a, x) for a <= 0 via Gamma(a, x) = (Gamma(a+1, x) - x^a e^-x) / a
    G1 = gamma(a + 1) * gammaincc(a + 1, x)
    return sc ** a * (G1 - x ** a * np.exp(-x)) / a


def nll(p, s, slo):
    tau, lsc = p
    sc = np.exp(lsc)
    if not (0.0 < tau < 2.5) or abs(1 - tau) < 1e-6:
        return 1e30
    Z = _norm(tau, sc, slo)
    if not np.isfinite(Z) or Z <= 0:
        return 1e30
    return tau * np.log(s).sum() + s.sum() / sc + len(s) * np.log(Z)


def fit(s, slo, nboot=200, seed=0):
    s = np.asarray(s, float); s = s[s >= slo]
    best = None
    for t0 in [0.9, 1.1, 1.3]:
        r = minimize(nll, [t0, np.log(np.percentile(s, 95))], args=(s, slo), method='Nelder-Mead',
                     options=dict(xatol=1e-4, fatol=1e-6, maxiter=4000))
        if best is None or r.fun < best.fun:
            best = r
    tau, sc = best.x[0], np.exp(best.x[1])
    rng = np.random.default_rng(seed); bt = []
    for _ in range(nboot):
        b = rng.choice(s, len(s))
        rb = minimize(nll, best.x, args=(b, slo), method='Nelder-Mead',
                      options=dict(xatol=1e-4, fatol=1e-6, maxiter=2000))
        bt.append(rb.x[0])
    return dict(tau=tau, s_c=sc, err=float(np.std(bt)), n=int(len(s)))


def powerfit(x, y, lo, hi):
    """Slope of log y against log x for lo <= x <= hi (least squares)."""
    sel = (x >= lo) & (x <= hi) & (y > 0)
    return np.polyfit(np.log(x[sel]), np.log(y[sel]), 1)[0]
