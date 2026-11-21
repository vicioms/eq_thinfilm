"""Quenched vs annealed complexity from the level runs in data_levels/.

    python3 analyze.py

Writes statistics.json, table_rates.tex, table_threshold.tex and fig1-fig4
(PDF + PNG).

Estimator.  With
    G_K(f) = ln<N_K(f)>_ann - <ln max(N_K(f), 1)>
the rate difference is Sigma_ann(f) - Sigma_q(f) = lim G_K / K.  G_K is
fitted to a + b K on the NFIT largest sizes and Sigma_q = Sigma_ann - b.
For the triangle model ln<N_K>_ann is exact at every K (Kac-Rice).  For the
interpolated Gaussian model it is K Sigma_ann(f) - (1/2) ln K, exact up to a
K-independent constant that does not affect b.
Above f_c, N = 0 in almost every sample, so Sigma_q -> 0 there.
"""
import glob
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq, curve_fit
from annealed import (gle, gle_stable, lyapunov, cramer_uniform, f_ann,
                      mean_abs_det, irwin_hall_density, _mats,
                      pl_cos_site_law, gle_tilted)

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
NU = 4.0 / 3.0  # nu = 1/(2 - zeta) for 1+1 depinning
NFIT = 4
NUB = 32
MODELS = [  # tag, label, colour, marker
    ("A0.3", r"tri, $A=0.3$", "#2a78d6", "o"),
    ("A0.6", r"tri, $A=0.6$", "#eb6834", "s"),
    ("A1", r"tri, $A=1$", "#1baf7a", "^"),
    ("cos_s0.6_M16", r"Gauss, $\sigma=0.6$", "#eda100", "D"),
]
plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "legend.frameon": False,
                     "axes.edgecolor": "#555555", "axes.labelcolor": "#222222",
                     "xtick.color": "#555555", "ytick.color": "#555555"})


def load(tag):
    out = {}
    for fn in glob.glob(f"data_levels/{tag}_K*.npz"):
        d = np.load(fn)
        out[int(d["K"])] = d
    return dict(sorted(out.items()))


def log_mean_det_stable(K, k, K0=20):
    """ln <|D_K| 1(all leading minors > 0)> for the triangle model:
    enumeration up to K0, then the exact asymptotic ratio exp(Lambda_st)."""
    Kc = min(K, K0)
    v = np.array([[1.0], [0.0]])
    ok = np.array([True])
    Ms = _mats(k)
    for j in range(Kc):
        v = np.concatenate([M @ v for M in Ms], axis=1)
        ok = np.concatenate([ok, ok])
        if j < Kc - 1:
            ok &= v[0] > 0
    ok &= v[1] > 0
    val = np.log(np.mean(np.abs(v[1]) * ok))
    if K > Kc:
        val += (K - Kc) * gle_stable(k)
    return val


class Tri:
    def __init__(self, A):
        self.A, self.P = A, 4.0
        self.k = 4 * A / self.P
        self.Lc = self.k ** (-2 / 3)
        self.lam1 = gle(1.0, self.k)
        self.lamst = gle_stable(self.k)
        h = 1e-3
        self.dlam1 = (gle(1 + h, self.k) - gle(1 - h, self.k)) / (2 * h)
        self.gamma = lyapunov(self.k, 400_000)
        self.fann = f_ann(A, self.P, self.lam1)
        self.fann_st = brentq(lambda f: self.lamst - cramer_uniform(f, A), 1e-9, A * (1 - 1e-8))
        self.chi = A ** 2 / 3
        self._lD = {}

    def S(self, f):
        return self.lam1 - cramer_uniform(f, self.A)

    def Sst(self, f):
        return self.lamst - cramer_uniform(f, self.A)

    def logEN(self, K, f):
        lD = mean_abs_det(K, self.k, 1.0, 20000)[1]
        lDs = log_mean_det_stable(K, self.k)
        out, outs = np.empty(len(f)), np.empty(len(f))
        for j, x in enumerate(f):
            p = irwin_hall_density(K, self.A, K * x)
            lp = np.log(self.P * p) if p > 0 else -np.inf
            out[j], outs[j] = lp + lD, lp + lDs
        return out, outs


class Gauss:
    """PL-interpolated FLRT disorder; annealed rate from the tilted operator."""

    def __init__(self, sigma, M):
        self.A, self.M, self.P = sigma, M, 2 * np.pi
        self.Lc = sigma ** (-2 / 3)
        self.law = pl_cos_site_law(sigma, M)
        self.taus = np.linspace(0.0, 14.0 / sigma, 57)
        self.Lt = np.array([gle_tilted(self.law, t) for t in self.taus])
        self.lam1 = self.Lt[0]
        self.lamst = np.nan
        self.dlam1 = np.nan
        self.gamma = np.nan
        self.chi = sigma ** 2
        self.fann = brentq(lambda f: self.S(f), 1e-6, 3 * sigma)
        self.fann_st = np.nan

    def S(self, f):
        # Legendre transform of the (convex, smooth) Lambda(tau), splined
        from scipy.interpolate import CubicSpline
        from scipy.optimize import minimize_scalar
        if not hasattr(self, "_sp"):
            self._sp = CubicSpline(self.taus, self.Lt)
        r = minimize_scalar(lambda t: self._sp(t) - t * f, bounds=(0.0, self.taus[-1]),
                            method="bounded", options={"xatol": 1e-6})
        return float(min(r.fun, self.Lt[0]))

    def Sst(self, f):
        return np.nan

    def logEN(self, K, f):
        a = np.array([K * self.S(x) - 0.5 * np.log(K) for x in f])
        return a, np.full(len(f), np.nan)


def analyse(tag):
    data = load(tag)
    if not data:
        return None
    d0 = next(iter(data.values()))
    is_tri = str(d0["model"]) == "tri" if "model" in d0 else True
    model = Tri(float(d0["A"])) if is_tri else Gauss(float(d0["A"]), int(d0["M"]))
    Ks = np.array(list(data.keys()))
    f = d0["f"]
    nf = len(f)
    G = np.zeros((len(Ks), nf)); Ge = np.zeros_like(G)
    Gm = np.full_like(G, np.nan); Gme = np.full_like(G, np.nan)
    P0 = np.zeros_like(G); loc = np.full_like(G, np.nan)
    ld = np.zeros(len(Ks)); lde = np.zeros(len(Ks))
    fc = np.zeros(len(Ks)); fce = np.zeros(len(Ks)); fcs = np.zeros(len(Ks))
    nsamp = []
    for a, K in enumerate(Ks):
        d = data[K]
        n = len(d["N"])
        nsamp.append(n)
        lEN, lENm = model.logEN(K, f)
        L = np.log(np.maximum(d["N"], 1))
        Lm = np.log(np.maximum(d["Nmin"], 1))
        G[a] = lEN - L.mean(0)
        Ge[a] = L.std(0) / np.sqrt(n)
        if np.all(np.isfinite(lENm)):
            Gm[a] = lENm - Lm.mean(0)
            Gme[a] = Lm.std(0) / np.sqrt(n)
        P0[a] = (d["N"] == 0).mean(0)
        with np.errstate(invalid="ignore", divide="ignore"):
            ex = (d["Y2"] - 1 / NUB - (1 - 1 / NUB) / d["N"]) / (1 - 1 / NUB)
            loc[a] = np.nanmean(np.where(d["N"] > 1, ex, np.nan), 0)
        x = d["ldet"][:, 0]
        ld[a], lde[a] = np.nanmean(x), np.nanstd(x) / np.sqrt(np.sum(~np.isnan(x)))
        fc[a], fce[a], fcs[a] = d["fc"].mean(), d["fc"].std() / np.sqrt(n), d["fc"].std()
    sel = np.arange(len(Ks))[-NFIT:]

    def slope(GG, EE):
        b = np.full(nf, np.nan); be = np.full(nf, np.nan)
        for j in range(nf):
            y, e = GG[sel, j], EE[sel, j]
            if not np.all(np.isfinite(y)):
                continue
            c, cov = np.polyfit(Ks[sel], y, 1, w=1 / np.maximum(e, 1e-4), cov="unscaled")
            b[j], be[j] = c[0], np.sqrt(cov[0, 0])
        return b, be

    b, be = slope(G, Ge)
    bm, bme = slope(Gm, Gme)
    Sann = np.array([model.S(x) for x in f])
    Sann_st = np.array([model.Sst(x) for x in f])
    Sq = np.clip(Sann - b, 0, None)
    Sqm = np.clip(Sann_st - bm, 0, None)
    cc, cov = np.polyfit(1 / Ks[sel], ld[sel], 1, w=1 / lde[sel], cov="unscaled")
    ld_inf, ld_inf_e = cc[1], np.sqrt(cov[1, 1])
    fit = lambda K, fi, c: fi + c * K ** (-1 / NU)
    msk = Ks >= 12
    pf, pcov = curve_fit(fit, Ks[msk], fc[msk], sigma=fce[msk], absolute_sigma=True, p0=[fc[-1], 0.1])
    fit3 = lambda K, fi, c, nu: fi + c * K ** (-1 / nu)
    try:
        pf3, pcov3 = curve_fit(fit3, Ks[msk], fc[msk], sigma=fce[msk], absolute_sigma=True,
                               p0=[pf[0], pf[1], NU], maxfev=20000)
        e3 = np.sqrt(np.diag(pcov3))
    except Exception:
        pf3, e3 = [np.nan] * 3, [np.nan] * 3
    # branching-front form: f_c(K) = f_c + a ln K / K + b / K
    fitB = lambda K, fi, a_, b_: fi + a_ * np.log(K) / K + b_ / K
    try:
        pB, pcovB = curve_fit(fitB, Ks[msk], fc[msk], sigma=fce[msk], absolute_sigma=True,
                              p0=[pf[0], 0.1, 0.0], maxfev=20000)
        eB = np.sqrt(np.diag(pcovB))
    except Exception:
        pB, eB = [np.nan] * 3, [np.nan] * 3
    # sample-to-sample spread of f_c(K) ~ K^{-omega}
    cs = np.polyfit(np.log(Ks[msk]), np.log(fcs[msk]), 1)
    omega = -cs[0]
    fcs_all = [pf[0], pf3[0], pB[0]]
    Ssa = None
    if not is_tri and os.path.exists("semi_annealed_cos.json"):
        from scipy.interpolate import CubicSpline
        from scipy.optimize import minimize_scalar
        sa = json.load(open("semi_annealed_cos.json"))
        if abs(sa["sigma"] - model.A) < 1e-9 and sa["M"] == model.M:
            sp = CubicSpline(sa["tau"], sa["Lambda_q"])
            Ssa = []
            for x in f:
                rr = minimize_scalar(lambda t: sp(t) - t * x, bounds=(0.0, sa["tau"][-1]), method="bounded")
                Ssa.append(float(min(rr.fun, sa["Lambda_q"][0])) if rr.x < sa["tau"][-1] - 1e-3 else np.nan)
            Ssa = np.array(Ssa)
    return dict(tag=tag, A=model.A,
                Sigma_sa=None if Ssa is None else Ssa.tolist(),
                sa_gap0=None if Ssa is None else float(Sann[0] - Ssa[0]), Lc=model.Lc, Lambda1=model.lam1, Lambda_st=model.lamst,
                dLambda1=model.dlam1, gamma=model.gamma, C=model.lam1 * model.Lc,
                C_st=model.lamst * model.Lc, chi=model.chi, f_ann=model.fann, f_ann_st=model.fann_st,
                K=Ks.tolist(), nsamp=nsamp, f=f.tolist(), fit_sizes=Ks[sel].tolist(),
                G=G.tolist(), G_err=Ge.tolist(), Gmin=Gm.tolist(), Gmin_err=Gme.tolist(),
                P_N0=P0.tolist(), loc_excess=loc.tolist(),
                rate_gap=b.tolist(), rate_gap_err=be.tolist(),
                rate_gap_min=bm.tolist(), rate_gap_min_err=bme.tolist(),
                Sigma_ann=Sann.tolist(), Sigma_ann_min=Sann_st.tolist(),
                Sigma_q=Sq.tolist(), Sigma_q_min=Sqm.tolist(),
                ldet_typ=ld.tolist(), ldet_typ_err=lde.tolist(),
                ldet_typ_inf=ld_inf, ldet_typ_inf_err=ld_inf_e,
                fc_mean=fc.tolist(), fc_err=fce.tolist(), fc_std=fcs.tolist(),
                fc_inf=float(pf[0]), fc_inf_err=float(np.sqrt(pcov[0, 0])), fc_c=float(pf[1]),
                fc_inf_freenu=float(pf3[0]), fc_inf_freenu_err=float(e3[0]),
                nu_free=float(pf3[2]), nu_free_err=float(e3[2]),
                fc_inf_front=float(pB[0]), fc_inf_front_err=float(eB[0]), fc_front_a=float(pB[1]),
                fc_front_b=float(pB[2]), fc_spread_exponent=float(omega),
                fc_range=[float(np.nanmin(fcs_all)), float(np.nanmax(fcs_all))],
                f_ann_over_fc_range=[float(model.fann / np.nanmax(fcs_all)), float(model.fann / np.nanmin(fcs_all))],
                f_ann_over_fc=model.fann / pf[0], f_ann_over_fc_freenu=model.fann / pf3[0],
                f_ann_st_over_fc=model.fann_st / pf[0])


def save(fig, name):
    fig.savefig(name + ".pdf", bbox_inches="tight")
    fig.savefig(name + ".png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def figures(R, prof):
    sty = {t: (lab, c, m) for t, lab, c, m in MODELS}
    tags = [t for t, *_ in MODELS if t in R]
    # ---- fig1: zero force ---------------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
    for t in tags:
        r = R[t]; lab, c, m = sty[t]
        K = np.array(r["K"]); G = np.array(r["G"])[:, 0]; Ge = np.array(r["G_err"])[:, 0]
        ax[0].errorbar(K, G - G[-1], Ge, color=c, marker=m, ms=4.5, lw=1.5, capsize=0, label=lab)
        Gm = np.array(r["Gmin"])[:, 0]
        if np.all(np.isfinite(Gm)):
            ax[0].errorbar(K, Gm - Gm[-1], np.array(r["Gmin_err"])[:, 0], color=c, marker=m, mfc="white",
                           ms=4.5, lw=1, ls="--", capsize=0)
        if t.startswith("A"):
            ax[1].errorbar(1 / K, np.array(r["ldet_typ"]) * r["Lc"], np.array(r["ldet_typ_err"]) * r["Lc"],
                           color=c, marker=m, ms=4.5, lw=1.5, capsize=0, label=lab)
            ax[1].plot([0, 0.13], [r["dLambda1"] * r["Lc"]] * 2, color=c, ls="--", lw=1)
            ax[1].plot([0, 0.13], [r["gamma"] * r["Lc"]] * 2, color=c, ls=":", lw=1)
    ax[0].axhline(0, color="#999999", lw=0.6)
    ax[0].set_xlabel("$K$")
    ax[0].set_ylabel(r"$G_K(0)-G_{K_{\max}}(0)$")
    ax[0].set_ylim(-0.35, 0.35)
    ax[0].legend(fontsize=7, loc="lower right")
    ax[0].text(0.03, 0.04, "filled: all; open: minima\nflat in $K$ $\\Leftrightarrow$ $\\Sigma_q=\\Sigma_{\\rm ann}$",
               fontsize=7, color="#555555", transform=ax[0].transAxes)
    ax[0].set_title("(a) $f=0$: annealed minus quenched", fontsize=9, loc="left")
    ax[1].set_xlabel("$1/K$")
    ax[1].set_ylabel(r"$L_c\,\overline{\ln|\det H|}/K$ over equilibria")
    ax[1].set_xlim(0, 0.13)
    ax[1].text(0.003, 0.655, r"$L_c\Lambda'(1)$ (dashed)", fontsize=7, color="#555555")
    ax[1].text(0.003, 0.31, r"$L_c\gamma$ (dotted)", fontsize=7, color="#555555")
    ax[1].legend(fontsize=7, loc="lower left", bbox_to_anchor=(0.02, 0.3))
    ax[1].set_title("(b) stability of typical equilibria", fontsize=9, loc="left")
    fig.tight_layout()
    save(fig, "fig1_zero_force")

    # ---- fig2: complexity vs force --------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
    for t in tags:
        r = R[t]; lab, c, m = sty[t]
        f = np.array(r["f"]); fa = r["f_ann"]; Lc = r["Lc"]
        idx = np.arange(len(f))
        for a_, key, ekey, S0, fend in [(ax[0], "Sigma_q", "rate_gap_err", "Sigma_ann", fa),
                                        (ax[1], "Sigma_q_min", "rate_gap_min_err", "Sigma_ann_min", r["f_ann_st"])]:
            if not np.isfinite(fend):
                continue
            Sa = np.array(r[S0])
            okc = (f <= fend) & np.isfinite(Sa)
            a_.plot(f[okc] / fa, Lc * Sa[okc], color=c, lw=1.5)
            if a_ is ax[0] and r.get("Sigma_sa") is not None:
                Ss = np.array(r["Sigma_sa"])
                oks = (f <= fend) & np.isfinite(Ss) & (Ss > 0)
                a_.plot(f[oks] / fa, Lc * Ss[oks], color=c, lw=1.2, ls="--")
            y = np.array(r[key]); e = np.array(r[ekey])
            ok = np.isfinite(y) & (f <= r["fc_inf"]) & (idx % 4 == 0)
            ab = np.isfinite(y) & (f > r["fc_inf"]) & (f <= 1.05 * fa) & (idx % 4 == 0)
            a_.errorbar(f[ok] / fa, Lc * y[ok], Lc * e[ok], color=c, marker=m, ms=4.5, ls="none",
                        capsize=0, label=lab)
            a_.plot(f[ab] / fa, Lc * y[ab], color=c, marker=m, mfc="white", ms=3.5, ls="none", alpha=0.6)
            a_.axvline(r["fc_inf"] / fa, color=c, lw=0.8, ls=":")
    for a_, t in zip(ax, ["(a) all equilibria", "(b) local minima (triangle)"]):
        a_.set_xlabel(r"$f/f_c^{\rm ann}$")
        a_.set_xlim(0, 1.05)
        a_.set_ylim(-0.01, 0.5)
        a_.set_title(t, fontsize=9, loc="left")
        a_.axhline(0, color="#999999", lw=0.6)
    ax[0].set_ylabel(r"$L_c\,\Sigma(f)$")
    ax[0].legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    save(fig, "fig2_force")

    # ---- fig3: G_K growth and threshold ---------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8))
    r = R["A0.6"]; K = np.array(r["K"]); G = np.array(r["G"]); f = np.array(r["f"])
    shades = ["#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#104281"]
    for s, fv in zip(shades, [0.0, 0.06, 0.10, 0.13, 0.16]):
        j = int(np.argmin(abs(f - fv)))
        ax[0].errorbar(K, G[:, j], np.array(r["G_err"])[:, j], color=s, marker="o", ms=4, lw=1.5,
                       capsize=0, label=f"$f={f[j]:.2f}$")
    ax[0].set_xlabel("$K$")
    ax[0].set_ylabel("$G_K(f)$")
    ax[0].set_title(r"(a) tri $A=0.6$: slope of $G_K$ = $\Sigma_{\rm ann}-\Sigma_q$", fontsize=9, loc="left")
    ax[0].legend(fontsize=7, loc="upper left")
    for t in tags:
        r = R[t]; lab, c, m = sty[t]
        K = np.array(r["K"]); fa = r["f_ann"]; x = K ** (-1 / NU)
        ax[1].errorbar(x, np.array(r["fc_mean"]) / fa, np.array(r["fc_err"]) / fa, color=c, marker=m,
                       ms=4.5, ls="none", capsize=0, label=lab)
        xx = np.linspace(0, x.max(), 50)
        ax[1].plot(xx, (r["fc_inf"] + r["fc_c"] * xx) / fa, color=c, lw=1)
        if np.isfinite(r["f_ann_st"]):
            ax[1].plot([0], [r["f_ann_st"] / fa], color=c, marker=m, mfc="white", ms=6, ls="none")
    ax[1].axhline(1, color="#999999", lw=0.8, ls="--")
    ax[1].axhline(1 / 1.64, color="#999999", lw=0.8, ls=":")
    ax[1].set_xlabel(r"$K^{-1/\nu}$, $\nu=4/3$")
    ax[1].set_ylabel(r"$\langle f_c(K)\rangle / f_c^{\rm ann}$")
    ax[1].set_xlim(0, None)
    ax[1].set_ylim(0.45, 1.05)
    ax[1].text(0.005, 1.01, "annealed bound", fontsize=7, color="#555555")
    ax[1].text(0.005, 0.585, "FLRT dynamics, 1/1.64", fontsize=7, color="#555555")
    ax[1].text(0.005, 0.92, "open: minima bound", fontsize=7, color="#555555")
    ax[1].legend(fontsize=7, loc="lower right")
    ax[1].set_title("(b) depinning threshold $f_c=\\max_a W_K/K$", fontsize=9, loc="left")
    fig.tight_layout()
    save(fig, "fig3_threshold")

    # ---- fig4: where the annealed picture fails --------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8))
    r = R["A0.6"]; K = np.array(r["K"]); f = np.array(r["f"]); loc = np.array(r["loc_excess"])
    for s, Kv in zip(shades, [12, 20, 28, 36, 44]):
        if Kv not in K:
            continue
        a = int(np.where(K == Kv)[0][0])
        ok = np.isfinite(loc[a])
        ax[0].plot(f[ok] / r["fc_inf"], loc[a, ok], color=s, lw=1.5, label=f"$K={Kv}$")
    ax[0].axvline(1, color="#999999", lw=0.8, ls=":")
    ax[0].set_xlim(0, 1.3)
    ax[0].set_ylim(0, None)
    ax[0].set_xlabel(r"$f/f_c$")
    ax[0].set_ylabel("excess coincidence of $u_{K-1}$ mod $P$")
    ax[0].legend(fontsize=7, loc="upper left")
    ax[0].set_title("(a) tri $A=0.6$: equilibria localise near $f_c$", fontsize=9, loc="left")
    if prof:
        Kp = np.array(prof["K"]); cc = np.array(prof["chi_c"]); ce = np.array(prof["chi_c_err"])
        cq = np.array(prof["chi_q"])
        ax[1].errorbar(Kp, Kp * cc / prof["chi"], Kp * ce / prof["chi"], color="#eb6834", marker="s", ms=4.5,
                       lw=1.5, capsize=0, label=r"$\langle c_s^2\rangle$ (not $\propto K$)")
        ax[1].plot(Kp, (cq + cc) / prof["chi"], color="#1baf7a", marker="^", ms=4.5, lw=1.5,
                   label=r"$\chi_q+\chi_c$")
        ax[1].plot(Kp, (Kp[0] * cc[0] / prof["chi"]) * Kp / Kp[0], color="#eb6834", lw=0.8, ls=":",
                   label=r"extensive reference $\propto K$")
        ax[1].axhline(1, color="#999999", lw=0.6, ls="--")
        ax[1].set_xlabel("$K$")
        ax[1].set_ylim(0, 4.2)
        ax[1].legend(fontsize=7, loc="upper left")
        ax[1].set_ylabel(r"in units of $\langle V'^2\rangle$")
        ax[1].set_title("(b) $f=0$ profile: centre spread sub-extensive", fontsize=9, loc="left")
    fig.tight_layout()
    save(fig, "fig4_localisation")


def table(R):
    def fmt(x, e=None, p=3):
        if x is None or not np.isfinite(x):
            return "--"
        return f"${x:.{p}f}\\pm{e:.{p}f}$" if e is not None else f"{x:.{p}f}"

    r1, r2 = [], []
    for t, lab, *_ in MODELS:
        if t not in R:
            continue
        r = R[t]; Lc = r["Lc"]
        name = lab.replace("tri, ", "triangle ").replace("Gauss, ", "Gaussian ")
        gmin = (fmt(r['rate_gap_min'][0] * Lc, r['rate_gap_min_err'][0] * Lc, 4)
                if np.isfinite(r['rate_gap_min'][0]) else "--")
        ldt = fmt(r['ldet_typ_inf'] * Lc, r['ldet_typ_inf_err'] * Lc) if t.startswith("A") else "--"
        r1.append(" & ".join([name, f"{Lc:.2f}", str(max(r['K'])), fmt(r['C']), fmt(r['C_st']),
                              fmt(r['rate_gap'][0] * Lc, r['rate_gap_err'][0] * Lc, 4), gmin,
                              fmt(r['gamma'] * Lc), fmt(r['dLambda1'] * Lc), ldt]) + r" \\")
        lo, hi = r['f_ann_over_fc_range']
        r2.append(" & ".join([name, fmt(r['f_ann'], p=4), fmt(r['f_ann_st'], p=4),
                              fmt(r['fc_inf'], r['fc_inf_err'], 4),
                              fmt(r['fc_inf_front'], r['fc_inf_front_err'], 4),
                              fmt(r['fc_inf_freenu'], r['fc_inf_freenu_err'], 4),
                              fmt(r['nu_free'], r['nu_free_err'], 2),
                              fmt(r['fc_spread_exponent'], p=2),
                              f"{lo:.2f}--{hi:.2f}"]) + r" \\")
    open("table_rates.tex", "w").write("\n".join(r1) + "\n")
    open("table_threshold.tex", "w").write("\n".join(r2) + "\n")


def main():
    R = {}
    for t, *_ in MODELS:
        r = analyse(t)
        if r is not None:
            R[t] = r
    prof = json.load(open("profile_center_A0.6.json")) if os.path.exists("profile_center_A0.6.json") else None
    json.dump({"models": R, "profile_center_A0.6": prof}, open("statistics.json", "w"), indent=1)
    for t, r in R.items():
        print(f"{t}: Lc={r['Lc']:.3f} K<={max(r['K'])} C={r['C']:.4f} C_st={r['C_st']:.4f} "
              f"gap(0)={r['rate_gap'][0]:+.4f}+-{r['rate_gap_err'][0]:.4f} "
              f"gap_min(0)={r['rate_gap_min'][0]:+.4f}+-{r['rate_gap_min_err'][0]:.4f} "
              f"ldet_inf={r['ldet_typ_inf']:.4f}+-{r['ldet_typ_inf_err']:.4f} Lambda'(1)={r['dLambda1']:.4f}")
        print(f"   f_ann={r['f_ann']:.4f} f_ann_st={r['f_ann_st']:.4f} fc={r['fc_inf']:.4f}+-{r['fc_inf_err']:.4f} "
              f"(free nu {r['fc_inf_freenu']:.4f}, nu={r['nu_free']:.2f}+-{r['nu_free_err']:.2f}) "
              f"front fit {r['fc_inf_front']:.4f}+-{r['fc_inf_front_err']:.4f} spread~K^-{r['fc_spread_exponent']:.2f} "
              f"f_ann/fc in {r['f_ann_over_fc_range'][0]:.3f}..{r['f_ann_over_fc_range'][1]:.3f}")
        f = np.array(r["f"]); fa = r["f_ann"]
        for fr in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]:
            j = int(np.argmin(abs(f / fa - fr)))
            print(f"   f/f_ann={f[j]/fa:.2f}: Sann={r['Sigma_ann'][j]:.4f} Sq={r['Sigma_q'][j]:.4f} "
                  f"gap={r['rate_gap'][j]:+.4f}+-{r['rate_gap_err'][j]:.4f} "
                  f"gap/Sann={r['rate_gap'][j]/r['Sigma_ann'][j]:+.3f}")
    figures(R, prof)
    table(R)


if __name__ == "__main__":
    main()
