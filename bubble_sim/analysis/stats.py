"""Bootstrap, power laws and the depinning threshold.

Threshold models for <h_d>(L):
  flat     h_c + A L^(-1/nu)
  bubble   h_c + B / L + A L^(-1/nu)       (B / L: Laplace pressure of a bubble with R0 ∝ L)
with nu fixed to qEW (4/3) or qKPZ (1.733), or free; nu also from var(h_d) ~ L^(-2/nu).
"""
from __future__ import annotations

import numpy as np

RNG = np.random.default_rng(11)
B = 300
NU_REF = {"qEW": 4 / 3, "qKPZ": 1.733}


def sem(x):
    x = np.asarray(x, float)
    return x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan


def loglog_slope(x, y, w=None):
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    if ok.sum() < 2:
        return np.nan
    return np.polyfit(np.log(x[ok]), np.log(y[ok]), 1, w=None if w is None else w[ok])[0]


def groups(size, min_n=8, max_unique=12, n_bins=14):
    """Indices grouped by size: by value when there are few distinct sizes (flat walls: L),
    else by log bins with at least min_n members (bubbles: R_eff)."""
    u = np.unique(size)
    if len(u) <= max_unique:
        return [np.nonzero(size == s)[0] for s in u]
    e = np.geomspace(size.min() * 0.999, size.max() * 1.001, n_bins + 1)
    k = np.digitize(size, e) - 1
    return [np.nonzero(k == i)[0] for i in range(n_bins) if (k == i).sum() >= min_n]


def length_scaling(size, P, min_n=8):
    """P ∝ size^alpha: alpha from all walls (bootstrap) and local alpha between successive
    size groups; group table (mean size, mean P, its error, n)."""
    gs = [g for g in groups(size, min_n) if len(g) >= 2]
    tab = np.array([(np.exp(np.log(size[g]).mean()), P[g].mean(), sem(P[g]), len(g)) for g in gs])
    al = loglog_slope(size, P)
    ab = []
    for _ in range(B):
        i = RNG.integers(0, len(P), len(P))
        ab.append(loglog_slope(size[i], P[i]))
    loc = np.log(tab[1:, 1] / tab[:-1, 1]) / np.log(tab[1:, 0] / tab[:-1, 0]) if len(tab) > 1 else np.zeros(0)
    locb = []
    for _ in range(B):
        m = np.array([P[g][RNG.integers(0, len(g), len(g))].mean() for g in gs])
        locb.append(np.log(m[1:] / m[:-1]) / np.log(tab[1:, 0] / tab[:-1, 0]))
    return (al, np.std(ab)), list(zip(loc, np.std(locb, 0) if len(loc) else [])), tab


def _wls_nu(Ls, v):
    return -2 / loglog_slope(np.asarray(Ls, float), np.asarray(v, float))


def threshold(L, h, model, lmin=0):
    """Per-L table and fits of the depinning field. Returns (rows, fits)."""
    from scipy.optimize import curve_fit
    Ls = np.unique(L)
    rows = []
    for l in Ls:
        x = h[L == l]
        vb = [np.var(x[RNG.integers(0, len(x), len(x))], ddof=1) for _ in range(B)] if len(x) > 2 else [np.nan]
        rows.append(dict(L=int(l), n=len(x), mean=x.mean(), err=sem(x), var=np.var(x, ddof=1) if len(x) > 1 else np.nan,
                         var_err=np.std(vb)))
    fits = {}
    R = [r for r in rows if r["L"] >= lmin and r["n"] > 2]
    if len(R) >= 2:
        Lf = np.array([r["L"] for r in R], float)
        nu = _wls_nu(Lf, [r["var"] for r in R])
        nb = []
        for _ in range(B):
            vb = [np.var(h[L == l][RNG.integers(0, (L == l).sum(), (L == l).sum())], ddof=1) for l in Lf]
            if min(vb) > 0:
                nb.append(_wls_nu(Lf, vb))
        fits["nu_var"] = (nu, np.std(nb))
        m = np.array([r["mean"] for r in R])
        s = np.array([r["err"] for r in R])
        bub = model == "bubble"

        for name, nu_ in list(NU_REF.items()) + [("free", None)]:
            npar = (2 if nu_ else 3) + bub
            if len(R) < npar + 1:
                continue
            p0 = [m[-1], (m[0] - m[-1]) * Lf[0] ** (1 / (nu_ or 1.6))] + ([] if nu_ else [1 / 1.6]) + ([0.0] if bub else [])
            try:
                p, c = curve_fit(lambda Lx, *q: model_h(Lx, q, nu_, bub), Lf, m, p0=p0, sigma=s,
                                 absolute_sigma=True, maxfev=50000)
            except (RuntimeError, ValueError):
                continue
            x, xe = (1 / nu_, 0.0) if nu_ else (p[2], np.sqrt(c[2, 2]))
            chi2 = np.sum(((m - model_h(Lf, p, nu_, bub)) / s) ** 2) / max(1, len(R) - npar)
            fits[f"hc_{name}"] = dict(hc=p[0], hc_err=np.sqrt(c[0, 0]), inv_nu=x, inv_nu_err=xe, chi2=chi2,
                                      p=p, nu=nu_, bubble=bub)
    return rows, fits


def model_h(Lx, q, nu, bubble):
    """<h_d>(L) = h_c + A L^(-1/nu) [+ B / L]; q = (h_c, A[, 1/nu][, B])."""
    x = 1 / nu if nu else q[2]
    Bc = q[-1] if bubble else 0.0
    return q[0] + q[1] * np.asarray(Lx, float) ** (-x) + Bc / np.asarray(Lx, float)


def eval_fit(fit, Lx):
    return model_h(Lx, fit["p"], fit["nu"], fit["bubble"])
