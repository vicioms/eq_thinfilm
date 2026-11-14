#!/usr/bin/env python3
"""
Is the isotropic case (a = 0) in the quenched Edwards-Wilkinson (qEW) depinning class?

qEW in d = 1 (Rosso & Krauth; Duemmer & Krauth; Ferrero, Bustingorry & Kolton):
    zeta = 1.25,  nu = 1 / (2 - zeta) = 4/3  (statistical tilt symmetry)
The competing large-scale scenario of Kolton, Ferrero & Rosso (arXiv:2306.13415)
beyond the overhang length is invasion-percolation-like: also nu = 4/3, but
self-similar walls (zeta = 1, fractal hull). So zeta, not nu, separates them.

Tests, all with bootstrap errors over disorder samples:
  1. zeta from the turning-angle spectrum of the wall contour, S_psi ~ q^(1 - 2 zeta)
     (fss_contours.py), compensated spectra, local slope, and a scan of the fit
     window (q_lo, q_hi); plus a fit with an explicit short-scale cut-off
     A q^(1-2 zeta) exp(-(q/q_c)^2) over a wider window.
  2. nu from var(h_d) ~ L^(-2/nu) (L >= lmin), and from successive pairs of L.
  3. scaling relation nu = 1/(2 - zeta).
  4. shape of the h_d distribution: collapse of (h_d - <h_d>)/std over L,
     skewness and excess kurtosis vs L (one fixed point -> L-independent shape).

    python qew_a0.py --dir runs/fss_D0.04
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fss_contours import Q_EDGES, N_MIN_FIT, Q_MAX_FIT  # noqa: E402

QC = np.sqrt(Q_EDGES[1:] * Q_EDGES[:-1])
B = 400
RNG = np.random.default_rng(12345)


def fit_slope(q, S, w=None):
    x, y = np.log(q), np.log(S)
    if w is None:
        w = np.ones_like(x)
    A = np.c_[np.ones_like(x), x]
    c = np.linalg.lstsq(A * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)[0]
    return c[1]


def zeta_boot(X, qlo, qhi):
    """zeta from the mean spectrum of samples X (n, bins); bootstrap over samples.
    Weights: inverse relative variance of the bin means (from the full sample)."""
    ok = (QC >= qlo) & (QC <= qhi) & np.all(np.isfinite(X), axis=0)
    if ok.sum() < 3:
        return np.nan, np.nan, None
    m = X[:, ok].mean(0)
    rel = X[:, ok].std(0, ddof=1) / np.sqrt(len(X)) / m
    w = 1 / rel ** 2
    z0 = (1 - fit_slope(QC[ok], m, w)) / 2
    zb = np.empty(B)
    for b in range(B):
        idx = RNG.integers(0, len(X), len(X))
        zb[b] = (1 - fit_slope(QC[ok], X[idx][:, ok].mean(0), w)) / 2
    return z0, zb.std(ddof=1), zb


def cutoff_fit(X, qlo, qhi):
    """Fit log S = log A + (1 - 2 zeta) log q - (q/q_c)^2 (linear in log A, 1-2zeta, 1/q_c^2)."""
    ok = (QC >= qlo) & (QC <= qhi) & np.all(np.isfinite(X), axis=0)

    def one(Y):
        m = Y[:, ok].mean(0)
        A = np.c_[np.ones(ok.sum()), np.log(QC[ok]), -QC[ok] ** 2]
        c = np.linalg.lstsq(A, np.log(m), rcond=None)[0]
        return (1 - c[1]) / 2, 1 / np.sqrt(c[2]) if c[2] > 0 else np.inf
    z0, qc0 = one(X)
    zb = np.array([one(X[RNG.integers(0, len(X), len(X))])[0] for _ in range(B)])
    return z0, zb.std(ddof=1), qc0


def load_hd(d):
    rows = [json.loads(l) for fn in glob.glob(os.path.join(d, "results*.jsonl")) for l in open(fn)]
    return [r for r in rows if r["status"] == "pinned" and r.get("h_d") is not None and np.isfinite(r["h_d"])]


def nu_fit(Ls, hd, lmin):
    """var(h_d) ~ L^(-2/nu): log-log fit over L >= lmin; bootstrap within each L."""
    Lf = np.array([L for L in Ls if L >= lmin], float)

    def one(samples):
        v = np.array([np.var(s, ddof=1) for s in samples])
        return -2 / np.polyfit(np.log(Lf), np.log(v), 1, w=None)[0]
    base = [hd[L] for L in Ls if L >= lmin]
    n0 = one(base)
    nb = np.array([one([s[RNG.integers(0, len(s), len(s))] for s in base]) for _ in range(B)])
    return n0, nb.std(ddof=1)


def moments(x):
    z = (x - x.mean()) / x.std(ddof=1)
    return np.mean(z ** 3), np.mean(z ** 4) - 3


ELL = np.geomspace(1, 4096, 37)


def height_stats(xy):
    """Transverse displacement along the wall from the integrated slope,
    u(s) = int (dpsi - <dpsi>) ds (well defined through overhangs; for a
    single-valued wall it is the usual height). Returns per-sample
    (q_n, S_u(q_n)), B(l) = <[u(s+l)-u(s)]^2> and the local width
    w^2(l) = <<u^2>_l - <u>_l^2> over all windows of arc length l."""
    from fss_contours import resample_closed, tangent_angle
    pts, P = resample_closed(xy)
    n = len(pts)
    ds = P / n
    psi, _ = tangent_angle(pts)
    s = np.arange(n) * ds
    dpsi = psi - psi[0] - 2 * np.pi * s / P
    dpsi -= dpsi.mean()
    u = np.cumsum(dpsi) * ds
    u -= u.mean()
    c = np.fft.rfft(u) / n
    q = 2 * np.pi * np.arange(len(c)) / P
    Su = P * np.abs(c) ** 2
    F = np.fft.rfft(u)
    Cc = np.fft.irfft(np.abs(F) ** 2, n) / n
    lag = np.rint(ELL / ds).astype(int)
    okl = lag < n // 4
    Bv = np.where(okl, 2 * (Cc[0] - Cc[np.minimum(lag, n - 1)]), np.nan)
    # local width with periodic windows via cumulative sums
    uu = np.r_[u, u]
    c1 = np.r_[0, np.cumsum(uu)]
    c2 = np.r_[0, np.cumsum(uu ** 2)]
    wv = np.full(len(ELL), np.nan)
    for i, m in enumerate(lag):
        if not okl[i] or m < 2:
            continue
        st = np.arange(0, n, max(1, m // 4))
        mu = (c1[st + m] - c1[st]) / m
        wv[i] = np.mean((c2[st + m] - c2[st]) / m - mu ** 2)
    return q[1:], Su[1:], Bv, wv


def compensated_plots(d, a, Ls_plot):
    import fss_shapes
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from fss_contours import bin_log
    C = fss_shapes.load_contours(d)
    L_all, a_all = np.array(C["L"]), np.array(C["aniso"])
    cL = {L: cm.plasma(i / max(1, len(Ls_plot) - 1) * 0.85) for i, L in enumerate(Ls_plot)}
    fig, axs = plt.subplots(2, 2, figsize=(12, 9.5))
    summary = []
    for L in Ls_plot:
        idx = np.nonzero((L_all == L) & (a_all == a))[0]
        Sub, Spb, Bs, Ws, Ps = [], [], [], [], []
        for i in idx:
            q, Su, Bv, wv = height_stats(C["contours"][i][C["main"][i]])
            Sub.append(bin_log(q, Su)); Spb.append(bin_log(q, Su * q ** 2)); Bs.append(Bv); Ws.append(wv)
            Ps.append(2 * np.pi / q[0])
        Pm = np.mean(Ps)
        okq = (QC > 2 * np.pi * 6 / Pm)
        for ax, X, pw, lab in ((axs[0, 0], np.array(Sub), 3.5, r"$q^{3.5}\,S_u(q)$"),
                               (axs[0, 1], np.array(Spb), 1.5, r"$q^{1.5}\,S_\psi(q)$")):
            m = np.nanmean(X, 0); e = np.nanstd(X, 0, ddof=1) / np.sqrt(len(X))
            ax.errorbar(QC[okq], (m * QC ** pw)[okq], (e * QC ** pw)[okq], fmt="o-", ms=3, lw=0.8, color=cL[L],
                        label=f"L={L} (n={len(idx)})")
            ax.set_ylabel(lab + "   (flat for ζ = 1.25)")
        for ax, X, pw, lab in ((axs[1, 0], np.array(Bs), -2.5, r"$B(\ell)\,\ell^{-2.5}$"),
                               (axs[1, 1], np.array(Ws), -2.5, r"$w^2(\ell)\,\ell^{-2.5}$")):
            m = np.nanmean(X, 0); e = np.nanstd(X, 0, ddof=1) / np.sqrt(len(X))
            ok = np.isfinite(m)
            ax.errorbar(ELL[ok], (m * ELL ** pw)[ok], (e * ELL ** pw)[ok], fmt="o-", ms=3, lw=0.8, color=cL[L],
                        label=f"L={L}")
            ax.set_ylabel(lab + "   (flat for ζ = 1.25 if ζ < 1 held)")
        mB = np.nanmean(np.array(Bs), 0)
        okB = np.isfinite(mB) & (ELL >= 10) & (ELL <= 100)
        sB = np.polyfit(np.log(ELL[okB]), np.log(mB[okB]), 1)[0]
        mW = np.nanmean(np.array(Ws), 0)
        okW = np.isfinite(mW) & (ELL >= 10) & (ELL <= 100)
        sW = np.polyfit(np.log(ELL[okW]), np.log(mW[okW]), 1)[0]
        summary.append(f"  L={L}: local exponent 10 <= l <= 100: B -> zeta_loc = {sB / 2:.3f}, w -> zeta_loc = {sW / 2:.3f}")
    for ax in axs[0]:
        ax.axvspan(Q_MAX_FIT, 1, color="0.93", zorder=-1)
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("q (arc length)"); ax.legend(fontsize=8)
    axs[0, 0].set_title("height spectrum S_u ~ q^-(1+2ζ), compensated")
    axs[0, 1].set_title("slope (turning-angle) spectrum S_ψ ~ q^(1-2ζ), compensated")
    for ax in axs[1]:
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(r"arc length $\ell$"); ax.legend(fontsize=8)
        x = np.geomspace(3, 300, 10)
        ax.plot(x, ax.get_lines()[-1].get_ydata()[0] * 0 + np.nan * x, alpha=0)
    axs[1, 0].set_title(r"$B(\ell)=\langle[u(s+\ell)-u(s)]^2\rangle \sim \ell^{2\zeta}$, compensated by $\ell^{-2.5}$")
    axs[1, 1].set_title(r"local width $w^2(\ell) \sim \ell^{2\zeta}$, compensated by $\ell^{-2.5}$")
    for ax in axs[1]:
        yl = ax.get_ylim()
        x = np.geomspace(2, 1000, 20)
        xm = 30.0
        ym = np.exp(np.mean(np.log(yl)))
        ax.plot(x, ym * (x / xm) ** -0.5, "k:", lw=0.8, label=r"$\zeta_{loc}=1$ ($\ell^2$)")
        ax.legend(fontsize=8)
    fig.suptitle(f"a = {a}: the same walls, four compensations (u = integrated slope along the contour)")
    fig.tight_layout()
    fn = os.path.join(d, f"qew_a{a:g}_compensated.png")
    fig.savefig(fn, dpi=110)
    return fn, summary


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=os.path.join(HERE, "runs", "fss_D0.04"))
    ap.add_argument("--a", type=float, default=0.0)
    args = ap.parse_args()
    d, a = args.dir, args.a
    S = dict(np.load(os.path.join(d, "contour_stats.npz")))
    sel = S["aniso"] == a
    Ls = sorted(set(S["L"][sel].tolist()))
    out = []

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm
    cL = {L: cm.plasma(i / max(1, len(Ls) - 1) * 0.85) for i, L in enumerate(Ls)}
    fig, axs = plt.subplots(2, 3, figsize=(16, 9.5))

    # ---------------------------------------------------------------- 1. zeta
    Lz = [L for L in Ls if L >= 1024]
    X = {L: S["Spsi_bins"][sel & (S["L"] == L)] for L in Ls}
    P = {L: S["P"][sel & (S["L"] == L)].mean() for L in Ls}
    big = S["L"] >= 2048
    l_ov = S["P"][big].sum() / S["n_ov"][big].sum()
    q_ov = 2 * np.pi / l_ov

    ax = axs[0, 0]
    for L in Ls[2:]:
        m = np.nanmean(X[L], 0)
        e = np.nanstd(X[L], 0, ddof=1) / np.sqrt(len(X[L]))
        ok = QC > 2 * np.pi * 4 / P[L]
        ax.errorbar(QC[ok], (m * QC ** 1.5)[ok], (e * QC ** 1.5)[ok], fmt="o-", ms=3, lw=0.8, color=cL[L], label=f"L={L}")
    ax.axvspan(Q_MAX_FIT, 1, color="0.93", zorder=-1)
    ax.axvline(q_ov, color="C2", ls="--", lw=0.8, label=r"$q_{ov}$")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("q")
    ax.set_ylabel(r"$q^{1.5}\,S_\psi(q)$  (flat for ζ = 1.25)")
    ax.set_title(f"compensated spectrum, a = {a}"); ax.legend(fontsize=7)
    ax = axs[0, 1]
    for L in Ls[2:]:
        m = np.nanmean(X[L], 0)
        e = np.nanstd(X[L], 0, ddof=1) / np.sqrt(len(X[L]))
        ok = QC > 2 * np.pi * 4 / P[L]
        ax.errorbar(QC[ok], (m * QC)[ok], (e * QC)[ok], fmt="o-", ms=3, lw=0.8, color=cL[L], label=f"L={L}")
    ax.axvspan(Q_MAX_FIT, 1, color="0.93", zorder=-1)
    ax.axvline(q_ov, color="C2", ls="--", lw=0.8)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("q")
    ax.set_ylabel(r"$q\,S_\psi(q)$  (flat for ζ = 1, self-similar)")
    ax.set_title("compensated for ζ = 1 (invasion-percolation-like)"); ax.legend(fontsize=7)

    # fit-window scan
    ax = axs[0, 2]
    qhis = [0.05, 0.07, 0.1, 0.14]
    out.append(f"a = {a}: zeta from S_psi, bootstrap over samples (B = {B})")
    out.append(f"  l_ov = {l_ov:.0f} -> q_ov = {q_ov:.4f}")
    for L, mk in zip(Lz, ("s", "^", "o")):
        qlo0 = 2 * np.pi * N_MIN_FIT / P[L]
        qlos = [q for q in (qlo0, 0.01, 0.015, 0.02, 0.03) if q >= qlo0]
        for qhi, ls in zip(qhis, (":", "-.", "-", "--")):
            zz = [zeta_boot(X[L], ql, qhi)[:2] for ql in qlos]
            z, ze = np.array(zz).T
            ax.errorbar(np.array(qlos) * (1 + 0.02 * qhis.index(qhi)), z, ze, fmt=mk + ls, ms=4, lw=0.8, color=cL[L],
                        label=f"L={L}, q_hi={qhi}")
            out.append(f"  L={L} q_hi={qhi}: " + "  ".join(f"q_lo={ql:.4f}: {zv:.3f}±{zev:.3f}" for ql, zv, zev in zip(qlos, z, ze)))
    ax.axhline(1.25, color="k", lw=0.6, ls="--", label="qEW 1.25")
    ax.axhline(1.0, color="k", lw=0.6, ls=":", label="ζ = 1")
    ax.set_xscale("log"); ax.set_xlabel(r"lower fit limit $q_{lo}$"); ax.set_ylabel("ζ")
    ax.set_title("fit-window scan (marker: L, line: q_hi)"); ax.legend(fontsize=6, ncol=2)

    # combined L >= 2048, small-scale and large-scale windows; cut-off model
    out.append("  windows (bootstrap):")
    res = {}
    for L in Lz:
        qlo0 = 2 * np.pi * N_MIN_FIT / P[L]
        for name, ql, qh in (("main n>=12..0.1", qlo0, Q_MAX_FIT), ("q_ov..0.1", max(q_ov, qlo0), Q_MAX_FIT),
                             ("n>=12..q_ov", qlo0, q_ov), ("0.02..0.07", 0.02, 0.07)):
            z, ze, _ = zeta_boot(X[L], ql, qh)
            res[(L, name)] = (z, ze)
            out.append(f"    L={L:5d} {name:16s} zeta = {z:.3f} ± {ze:.3f}")
        z, ze, qcut = cutoff_fit(X[L], qlo0, 0.3)
        res[(L, "cutoff")] = (z, ze)
        out.append(f"    L={L:5d} cut-off model n>=12..0.3: zeta = {z:.3f} ± {ze:.3f}  (q_c = {qcut:.3f})")

    # ---------------------------------------------------------------- 2. nu
    rows = load_hd(d)
    hd = {L: np.array([r["h_d"] for r in rows if r["L"] == L and r["aniso"] == a]) for L in Ls}
    ax = axs[1, 0]
    v = np.array([np.var(hd[L], ddof=1) for L in Ls])
    ve = np.array([np.std([np.var(hd[L][RNG.integers(0, len(hd[L]), len(hd[L]))], ddof=1) for _ in range(B)]) for L in Ls])
    ax.errorbar(Ls, v, ve, fmt="ko", label="var(h_d)")
    out.append("")
    out.append("nu from var(h_d) ~ L^(-2/nu), bootstrap within each L:")
    for lmin, ls in ((256, ":"), (512, "-"), (1024, "--")):
        n0, ne = nu_fit(Ls, hd, lmin)
        Lf = np.array([L for L in Ls if L >= lmin], float)
        vf = np.array([np.var(hd[int(L)], ddof=1) for L in Lf])
        c = np.polyfit(np.log(Lf), np.log(vf), 1)
        ax.plot(Lf, np.exp(np.polyval(c, np.log(Lf))), "C3" + ls, lw=0.9, label=f"L ≥ {lmin}: ν = {n0:.2f} ± {ne:.2f}")
        out.append(f"  L >= {lmin}: nu = {n0:.3f} ± {ne:.3f}")
    pair = []
    for L1, L2 in zip(Ls[:-1], Ls[1:]):
        r = np.log(v[Ls.index(L2)] / v[Ls.index(L1)]) / np.log(L2 / L1)
        rb = []
        for _ in range(B):
            s1 = hd[L1][RNG.integers(0, len(hd[L1]), len(hd[L1]))]
            s2 = hd[L2][RNG.integers(0, len(hd[L2]), len(hd[L2]))]
            rb.append(-2 / (np.log(np.var(s2, ddof=1) / np.var(s1, ddof=1)) / np.log(L2 / L1)))
        pair.append((np.sqrt(L1 * L2), -2 / r, np.std(rb)))
        out.append(f"  pair {L1}-{L2}: nu_eff = {-2 / r:.3f} ± {np.std(rb):.3f}")
    Lx = np.array(Ls, float)
    ax.plot(Lx, v[-1] * (Lx / Lx[-1]) ** (-1.5), "k:", lw=0.6, label="ν = 4/3")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("L"); ax.set_ylabel("var(h_d)")
    ax.set_title("sample-to-sample fluctuations of h_d"); ax.legend(fontsize=7)
    ax2 = ax.inset_axes([0.6, 0.62, 0.37, 0.33])
    px, pn, pe = np.array(pair).T
    ax2.errorbar(px, pn, pe, fmt="o-", ms=3, color="C3")
    ax2.axhline(4 / 3, color="k", ls=":", lw=0.6)
    ax2.set_xscale("log"); ax2.set_ylim(0.5, 3); ax2.set_title(r"$\nu_{eff}$ from pairs", fontsize=7)
    ax2.tick_params(labelsize=6)

    # ---------------------------------------------------------------- 3. scaling relation
    ax = axs[1, 1]
    for L in Lz:
        z, ze = res[(L, "main n>=12..0.1")]
        nu_sts, nu_e = 1 / (2 - z), ze / (2 - z) ** 2
        ax.errorbar([z], [nu_sts], [nu_e], xerr=[ze], fmt="o", color=cL[L], label=f"ν = 1/(2−ζ), ζ from L={L}")
    n512, e512 = nu_fit(Ls, hd, 512)
    ax.axhspan(n512 - e512, n512 + e512, color="C3", alpha=0.2, label=f"ν from var(h_d), L ≥ 512: {n512:.2f}±{e512:.2f}")
    ax.plot([1.25], [4 / 3], "k*", ms=12, label="qEW (1.25, 4/3)")
    ax.plot([1.0], [4 / 3], "kx", ms=10, label="IP-like (1, 4/3)")
    ax.set_xlabel("ζ"); ax.set_ylabel("ν"); ax.set_xlim(0.9, 1.4); ax.set_ylim(0.9, 2.0)
    ax.set_title("statistical tilt symmetry ν = 1/(2−ζ)"); ax.legend(fontsize=7)
    out.append("")
    out.append("scaling relation nu = 1/(2 - zeta):")
    for L in Lz:
        z, ze = res[(L, "main n>=12..0.1")]
        out.append(f"  zeta(L={L}) = {z:.3f} ± {ze:.3f}  ->  1/(2-zeta) = {1 / (2 - z):.3f} ± {ze / (2 - z) ** 2:.3f};"
                   f"  measured nu (L>=512) = {n512:.3f} ± {e512:.3f}")

    # ---------------------------------------------------------------- 4. h_d distribution
    ax = axs[1, 2]
    out.append("")
    out.append("shape of the h_d distribution: skewness, excess kurtosis (bootstrap errors)")
    for L in Ls:
        x = hd[L]
        z = np.sort((x - x.mean()) / x.std(ddof=1))
        ax.plot(z, (np.arange(len(z)) + 0.5) / len(z), "-", color=cL[L], lw=1, label=f"L={L} (n={len(z)})")
        sk, ku = moments(x)
        bs = np.array([moments(x[RNG.integers(0, len(x), len(x))]) for _ in range(B)])
        out.append(f"  L={L:5d} n={len(x):4d}  skew = {sk:+.3f} ± {bs[:, 0].std():.3f}   ex.kurt = {ku:+.3f} ± {bs[:, 1].std():.3f}")
    from math import erf
    zz = np.linspace(-3.5, 3.5, 200)
    ax.plot(zz, 0.5 * (1 + np.vectorize(erf)(zz / np.sqrt(2))), "k:", lw=1, label="Gaussian")
    ax.set_xlabel(r"$(h_d-\langle h_d\rangle)/\sigma$"); ax.set_ylabel("CDF")
    ax.set_title("collapse of the h_d distribution"); ax.legend(fontsize=7)
    fig.tight_layout()
    fn = os.path.join(d, f"qew_a{a:g}.png")
    fig.savefig(fn, dpi=110)
    txt = "\n".join(out)
    print(txt)
    with open(os.path.join(d, f"qew_a{a:g}.txt"), "w") as f:
        f.write(txt + "\n")
    print("saved", fn)
    fn2, summ = compensated_plots(d, a, [L for L in Ls if L >= 512])
    print("\n".join(summ))
    with open(os.path.join(d, f"qew_a{a:g}.txt"), "a") as f:
        f.write("\nreal-space local exponents of the integrated-slope height u(s):\n" + "\n".join(summ) + "\n")
    print("saved", fn2)


if __name__ == "__main__":
    main()
