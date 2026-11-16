"""Pooled mode spectra and roughness exponents.

Each wall gives modes (q_n, S_n, n), n >= 1, with S_n = P |c_n|^2 (P = period: arc length for
the tangent angle psi(s), L or 2 pi <R> for a Cartesian / polar profile). For a set of walls
(fixed L and a) the modes n >= n_min of all walls are pooled in one array and binned in log q;
adjacent fine bins are merged until each holds at least one mode per wall, so at low q, where
the mode spacing 2 pi / P exceeds a fine bin, no bin is left to the few walls whose P happens
to put a mode there. Bin position = mean log q of its modes. Fits are weighted by
sqrt(modes per bin) (the relative error of a mean of |c_n|^2), errors by bootstrap over walls
with the bins fixed.

    S_psi(q) ~ q^(1 - 2 zeta)        tangent angle along the arc length (kind "psi")
    S_u(q)   ~ q^-(1 + 2 zeta)       height / radius profile (kind "u")
"""
from __future__ import annotations

import numpy as np

FINE = np.geomspace(1e-3, 3.0, 241)
RNG = np.random.default_rng(5)
B = 300


def modes(f, P):
    """(q_n, P |c_n|^2, n), n >= 1, of a periodic signal sampled uniformly over period P."""
    c = np.fft.rfft(np.asarray(f, float) - np.mean(f)) / len(f)
    n = np.arange(1, len(c))
    return 2 * np.pi * n / P, P * np.abs(c[1:]) ** 2, n


class Pooled:
    def __init__(self, spectra, n_min):
        parts = [(q[k >= n_min], S[k >= n_min]) for q, S, k in spectra]
        self.n = len(parts)
        if self.n == 0:
            self.qc = np.zeros(0)
            self.S = self.c = np.zeros((0, 0))
            return
        q = np.concatenate([p[0] for p in parts])
        cnt = np.histogram(q, FINE)[0]
        edges, acc = [FINE[0]], 0
        for i, ci in enumerate(cnt):
            acc += ci
            if acc >= max(self.n, 8):
                edges.append(FINE[i + 1])
                acc = 0
        if acc and len(edges) > 2:
            edges[-1] = FINE[-1]
        self.edges = np.array(edges)
        nb = len(self.edges) - 1
        # per-wall bin sums: a bootstrap resample is a sum of rows
        self.S, self.c, lq = np.zeros((self.n, nb)), np.zeros((self.n, nb)), np.zeros(nb)
        for w, (qw, Sw) in enumerate(parts):
            i = np.digitize(qw, self.edges) - 1
            ok = (i >= 0) & (i < nb)
            self.S[w] = np.bincount(i[ok], Sw[ok], nb)
            self.c[w] = np.bincount(i[ok], minlength=nb)
            lq += np.bincount(i[ok], np.log(qw[ok]), nb)
        with np.errstate(invalid="ignore", divide="ignore"):
            self.qc = np.exp(lq / self.c.sum(0))

    def mean(self, idx=None):
        S, c = (self.S, self.c) if idx is None else (self.S[idx], self.c[idx])
        S, c = S.sum(0), c.sum(0)
        with np.errstate(invalid="ignore", divide="ignore"):
            return S / c, c

    def boot(self, nb=B):
        return [self.mean(RNG.integers(0, self.n, self.n)) for _ in range(nb)]

    def mean_err(self, nb=100):
        m = self.mean()[0]
        return m, np.nanstd([S for S, _ in self.boot(nb)], 0)


def _slope(qc, S, c, ok):
    ok = ok & (c > 0) & np.isfinite(S) & (S > 0)
    if ok.sum() < 3:
        return np.nan
    return np.polyfit(np.log(qc[ok]), np.log(S[ok]), 1, w=np.sqrt(c[ok]))[0]


def zeta_of_slope(sl, kind):
    return (1 - sl) / 2 if kind == "psi" else -(1 + sl) / 2


def fit(pool, qhi, kind, qlo=0.0):
    """zeta from a power law over qlo <= q <= qhi (modes n >= n_min of the pool)."""
    if pool.n < 4:
        return np.nan, np.nan
    win = (pool.qc >= qlo) & (pool.qc <= qhi)
    z0 = zeta_of_slope(_slope(pool.qc, *pool.mean(), win), kind)
    zb = [zeta_of_slope(_slope(pool.qc, S, c, win), kind) for S, c in pool.boot()]
    return z0, np.nanstd(zb, ddof=1)


def local_zeta(pool, kind, half=np.log(2) / 2, nb=100):
    """zeta_eff(q) from weighted fits over a sliding window of a factor 2 in q."""
    lq = np.log(pool.qc)

    def one(S, c):
        return np.array([zeta_of_slope(_slope(pool.qc, S, c, np.abs(lq - l) <= half), kind) for l in lq])
    return one(*pool.mean()), np.nanstd([one(S, c) for S, c in pool.boot(nb)], 0)
