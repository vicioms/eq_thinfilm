"""Theory checks: quenched pressure, frozen asymptote, near-threshold counts,
amplitude-quenched threshold bound.

    python3 theory.py   -> theory.json, table_theory.tex, fig5_theory.{pdf,png}

1. Quenched pressure (coarea):  Z_K(tau) = int dx N_K(x) e^{tau x}
   = int da e^{tau W_K} |det H|.  phi_q(tau) = lim K^-1 <ln Z_K>, estimated as
   Lambda(tau) - slope_K [ln E Z_K - <ln Z_K>].  ln E Z_K is exact for the
   triangle model, P (sinh(tau A)/(tau A))^K <|D_K|>, and K Lambda(tau) +
   const for the Gaussian model (tilted operator).  The x < 0 half of the
   level grid is mirrored (same law; O(1) in ln Z).
2. Frozen asymptote: phi_q(tau) -> c + f_c tau for large tau.  Its slope is
   the threshold, its intercept c = Sigma_q(f_c^-) by Legendre duality.
3. Near threshold, directly: <ln N> at an absolute depth delta below each
   sample's maximum max_a W_K (force delta/K below the sample threshold),
   for all equilibria and for minima; its slope in K is c (resp. c_min).
4. Gaussian model: amplitude-quenched bound f_sa = inf_tau Lambda_q(tau)/tau
   from semi_annealed_cos.json (rigorous: f_c <= f_sa <= f_ann).
"""
import glob
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar
from annealed import gle, gle_stable, mean_abs_det, pl_cos_site_law, gle_tilted

os.chdir(os.path.dirname(os.path.abspath(__file__)))
NFIT = 4
MODELS = [("A0.3", "tri", 0.3, r"tri, $A=0.3$", "#2a78d6", "o"),
          ("A0.6", "tri", 0.6, r"tri, $A=0.6$", "#eb6834", "s"),
          ("A1", "tri", 1.0, r"tri, $A=1$", "#1baf7a", "^"),
          ("cos_s0.6_M16", "cos", 0.6, r"Gauss, $\sigma=0.6$", "#eda100", "D")]
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "axes.edgecolor": "#555555",
                     "axes.labelcolor": "#222222", "xtick.color": "#555555", "ytick.color": "#555555"})


def lnZ(N, x, tau):
    """ln sum_j N_j (e^{tau x_j} + e^{-tau x_j}) w_j, trapezoid weights w_j."""
    dx = x[1] - x[0]
    w = np.full(len(x), dx)
    w[0] = dx / 2
    lc = np.logaddexp(tau * x, -tau * x)
    with np.errstate(divide="ignore"):
        L = np.log(N) + np.log(w)[None, :] + lc[None, :]
    m = L.max(1, keepdims=True)
    return m[:, 0] + np.log(np.exp(L - m).sum(1))


def load(tag):
    out = {}
    for fn in glob.glob(f"data_levels/{tag}_K*.npz"):
        d = np.load(fn)
        out[int(d["K"])] = d
    return dict(sorted(out.items()))


def linfit(K, y, e=None):
    c, cov = np.polyfit(K, y, 1, w=None if e is None else 1 / e, cov="unscaled" if e is not None else True)
    return c[0], float(np.sqrt(cov[0, 0]))


def analyse(tag, model, A):
    data = load(tag)
    Ks = np.array(list(data))
    sel = slice(-NFIT, None)
    if model == "tri":
        k = A
        lam1 = gle(1.0, k)
        lamst = gle_stable(k)
        Lam = lambda t: lam1 + (np.log(np.sinh(t * A) / (t * A)) if t > 0 else 0.0)
        lDs = {K: mean_abs_det(K, k, 1.0, 20000)[1] for K in Ks}
        lnEZ = lambda K, t: np.log(4.0) + K * (Lam(t) - lam1) + lDs[K]
    else:
        law = pl_cos_site_law(A, 16)
        tg = np.linspace(0, 40 / A, 81)
        sp = CubicSpline(tg, [gle_tilted(law, t) for t in tg])
        lam1 = float(sp(0.0))
        lamst = np.nan
        Lam = lambda t: float(sp(t))
        lnEZ = lambda K, t: K * Lam(t)
    taus = np.linspace(0, 24, 97) / A
    H = np.zeros((len(Ks), len(taus))); He = np.zeros_like(H)
    top = {"N": {}, "Nmin": {}}
    for a, K in enumerate(Ks):
        d = data[K]
        x = d["f"] * K
        Nf = d["N"].astype(float)
        for b, t in enumerate(taus):
            lz = lnZ(Nf, x, t)
            H[a, b] = lnEZ(K, t) - lz.mean()
            He[a, b] = lz.std() / np.sqrt(len(lz))
        Wm = d["fc"] * K
        for key in top:
            for dl in (0.5, 1.0, 2.0):
                j = np.searchsorted(x, Wm - dl * A) - 1
                ok = j >= 0
                n = d[key][np.arange(len(Wm))[ok], j[ok]]
                top[key].setdefault(dl, []).append(np.mean(np.log(np.maximum(n, 1))))
    gap = np.zeros(len(taus)); gape = np.zeros(len(taus))
    for b in range(len(taus)):
        gap[b], gape[b] = linfit(Ks[sel], H[sel, b], He[sel, b])
    lam = np.array([Lam(t) for t in taus])
    phi = lam - gap
    # frozen asymptote from the last quarter of the tau range
    m = taus * A >= 18
    cf = np.polyfit(taus[m], phi[m], 1)
    rates = {key: {dl: linfit(Ks[sel], np.array(v)[sel]) for dl, v in top[key].items()} for key in top}
    out = dict(tag=tag, model=model, A=A, K=Ks.tolist(), Lambda1=lam1, Lambda_st=lamst,
               tau=taus.tolist(), Lambda=lam.tolist(), phi_q=phi.tolist(), gap=gap.tolist(),
               gap_err=gape.tolist(), frozen_slope=float(cf[0]), frozen_intercept=float(cf[1]),
               near_top={key: {str(dl): list(v) for dl, v in r.items()} for key, r in rates.items()},
               near_top_curves={key: {str(dl): list(map(float, v)) for dl, v in tv.items()} for key, tv in top.items()})
    if model == "cos" and os.path.exists("semi_annealed_cos.json"):
        sa = json.load(open("semi_annealed_cos.json"))
        t = np.array(sa["tau"]); Lq = np.array(sa["Lambda_q"]); E = np.array(sa["Lambda_q_err"])
        fs = []
        for s in (0, 1, -1):
            spq = CubicSpline(t, Lq + s * E)
            fs.append(minimize_scalar(lambda z: spq(z) / z, bounds=(0.2, t[-1]), method="bounded").fun)
        out["f_sa"] = float(fs[0]); out["f_sa_err"] = float(abs(fs[1] - fs[2]) / 2)
        out["Lambda_q0"] = float(Lq[0]); out["Lambda_q0_err"] = float(E[0])
    return out


def figure(R, stats):
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
    anchor = None
    for tag, model, A, lab, c, mk in MODELS:
        if tag not in R:
            continue
        r = R[tag]
        K = np.array(r["K"])
        if model == "tri":
            tA = np.array(r["tau"]) * A
            g = np.array(r["gap"]) - r["gap"][0]
            e = np.array(r["gap_err"])
            ok = (tA >= 0.25) & (tA <= 2.0) & (g > 0)
            ax[0].errorbar(tA[ok], g[ok] / r["Lambda1"], e[ok] / r["Lambda1"], color=c, marker=mk, ms=4,
                           ls="-", lw=1.2, capsize=0, label=lab)
            if tag == "A0.6":
                j1 = int(np.argmin(abs(tA - 1.0)))
                anchor = g[j1] / r["Lambda1"]
        for key, mfc in (("N", c), ("Nmin", "white")):
            y = np.array(r["near_top_curves"][key]["0.5"])
            ax[1].plot(K, y, color=c, marker=mk, ms=4, mfc=mfc, lw=1.2 if key == "N" else 0.8,
                       ls="-" if key == "N" else "--", label=lab if key == "N" else None)
    xx = np.array([0.3, 2.0])
    ax[0].plot(xx, anchor * xx ** 2, color="#999999", lw=0.8, ls=":")
    ax[0].plot(xx, anchor * xx ** 4, color="#999999", lw=0.8, ls="--")
    ax[0].text(0.31, anchor * 0.31 ** 2 * 1.35, r"$\tau^2$", fontsize=8, color="#777777")
    ax[0].text(0.36, anchor * 0.36 ** 4 * 0.3, r"$\tau^4$", fontsize=8, color="#777777")
    ax[0].set_xscale("log"); ax[0].set_yscale("log")
    from matplotlib.ticker import FixedLocator, NullLocator, FixedFormatter
    ax[0].xaxis.set_major_locator(FixedLocator([0.25, 0.5, 1.0, 2.0]))
    ax[0].xaxis.set_major_formatter(FixedFormatter(["0.25", "0.5", "1", "2"]))
    ax[0].xaxis.set_minor_locator(NullLocator())
    ax[0].set_xlabel(r"$\tau A$")
    ax[0].set_ylabel(r"$[\Lambda-\varphi_q](\tau)-[\Lambda-\varphi_q](0)$, units $\Lambda(1)$")
    ax[0].legend(fontsize=7, loc="upper left")
    ax[0].set_title("(a) triangle: pressure gap is not $O(\\tau^2)$", fontsize=9, loc="left")
    ax[1].set_xlabel("$K$")
    ax[1].set_ylabel(r"$\langle\ln N\rangle$ at $\max_aW_K-A/2$")
    ax[1].set_title("(b) exponentially many states at threshold", fontsize=9, loc="left")
    ax[1].text(0.03, 0.9, "filled: all equilibria; open: minima", fontsize=7, color="#555555",
               transform=ax[1].transAxes)
    ax[1].legend(fontsize=7, loc="center left", bbox_to_anchor=(0.0, 0.62))
    fig.tight_layout()
    fig.savefig("fig5_theory.pdf", bbox_inches="tight")
    fig.savefig("fig5_theory.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def table(R, stats):
    rows = []
    for tag, model, A, lab, *_ in MODELS:
        if tag not in R:
            continue
        r = R[tag]; s = stats[tag]
        name = lab.replace("tri, ", "triangle ").replace("Gauss, ", "Gaussian ")
        ca, cae = r["near_top"]["N"]["0.5"]
        cm, cme = r["near_top"]["Nmin"]["0.5"]
        fsa = f"${r['f_sa']:.3f}\\pm{r['f_sa_err']:.3f}$" if "f_sa" in r else "$=f_c^{\\rm ann}$"
        rows.append(" & ".join([
            name, f"{s['f_ann']:.3f}", fsa,
            f"{s['fc_range'][0]:.3f}--{s['fc_range'][1]:.3f}", f"{r['frozen_slope']:.3f}",
            f"${ca:.3f}\\pm{cae:.3f}$", f"{ca / r['Lambda1']:.2f}",
            f"${cm:.3f}\\pm{cme:.3f}$" if model == "tri" else "--",
            f"{cm / r['Lambda_st']:.2f}" if model == "tri" else "--"]) + r" \\")
    open("table_theory.tex", "w").write("\n".join(rows) + "\n")


def main():
    stats = json.load(open("statistics.json"))["models"]
    R = {}
    for tag, model, A, *_ in MODELS:
        if load(tag):
            R[tag] = analyse(tag, model, A)
            r = R[tag]
            print(f"{tag}: frozen slope {r['frozen_slope']:.4f} (f_c fits {stats[tag]['fc_range']}), "
                  f"intercept {r['frozen_intercept']:.4f}; near-top rate all "
                  f"{r['near_top']['N']['0.5'][0]:.4f} (/Lambda1 = {r['near_top']['N']['0.5'][0]/r['Lambda1']:.2f}), "
                  f"minima {r['near_top']['Nmin']['0.5'][0]:.4f}"
                  + (f"; f_sa = {r['f_sa']:.4f}+-{r['f_sa_err']:.4f}" if "f_sa" in r else ""))
    json.dump(R, open("theory.json", "w"), indent=1)
    figure(R, stats)
    table(R, stats)




def figure_torus():
    """fig6: torus transfer-operator results (torus.json, torus_decomposition.json)."""
    if not (os.path.exists("torus.json") and os.path.exists("torus_decomposition.json")):
        return
    T = json.load(open("torus.json")); D = json.load(open("torus_decomposition.json"))
    E = json.load(open("theory.json"))
    fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.6))
    for A, tag, c, mk in ((0.6, "A0.6", "#eb6834", "s"), (1.0, "A1", "#1baf7a", "^")):
        sc = [x for x in T["scan"] if x["A"] == A and x["tauA"] <= 3.0]
        tA = np.array([x["tauA"] for x in sc]); L = np.array([x["Lambda"] for x in sc])
        ph = np.array([x["phi"] for x in sc])
        e = E[tag]; te = np.array(e["tau"]) * A; pe = np.array(e["phi_q"]); m = te <= 3.0
        ax[0].plot(te[m], np.array(e["Lambda"])[m] - pe[m], color=c, lw=1.5, label=f"enumeration, $A={A:g}$")
        ax[0].plot(tA, L - ph, color=c, marker=mk, ms=4.5, ls="none", label=f"torus, $A={A:g}$")
        dd = sorted([x for x in D if x["A"] == A], key=lambda x: x["tauA"])
        t2 = np.array([x["tauA"] for x in dd])
        ok = t2 > 0.3
        ax[1].plot(t2[ok], np.array([x["jensen"] + x["memory"] for x in dd])[ok], color=c, marker=mk, ms=4.5,
                   lw=1.5, label=f"Jensen+memory, $A={A:g}$")
        ax[1].plot(t2[ok], np.array([x["jensen2"] for x in dd])[ok], color=c, lw=1, ls="--")
        ax[1].plot(t2[ok], np.array([x["memory"] for x in dd])[ok], color=c, lw=1, ls=":")
    ax[0].set_xlabel(r"$\tau A$"); ax[0].set_ylabel(r"$\Lambda(\tau)-\varphi_q(\tau)$")
    ax[0].set_yscale("log"); ax[0].set_ylim(5e-4, 1)
    ax[0].legend(fontsize=6.5, loc="upper left")
    ax[0].set_title("(a) quenched gap: torus vs exact", fontsize=8.5, loc="left")
    ax[1].set_xscale("log"); ax[1].set_yscale("log")
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
    ax[1].xaxis.set_major_locator(FixedLocator([0.5, 1, 2, 3]))
    ax[1].xaxis.set_major_formatter(FixedFormatter(["0.5", "1", "2", "3"]))
    ax[1].xaxis.set_minor_locator(NullLocator())
    ax[1].set_xlabel(r"$\tau A$"); ax[1].set_ylabel("per-site gap")
    ax[1].legend(fontsize=6.5, loc="upper left")
    ax[1].set_title("(b) gap = one-step Jensen term", fontsize=8.5, loc="left")
    cv = T["convergence"]
    h = np.array([x["h"] for x in cv]); sl = np.array([x["slope"] for x in cv]); se = np.array([x["slope_err"] for x in cv])
    ic = np.array([x["intercept"] for x in cv])
    ax[2].errorbar(h, sl, se, color="#2a78d6", marker="o", ms=4.5, lw=1.5, capsize=0, label=r"slope $\to f_c$")
    ax[2].plot(h, ic, color="#e87ba4", marker="D", ms=4.5, lw=1.5, label=r"intercept $\to c$")
    ax[2].axhspan(0.171, 0.175, color="#2a78d6", alpha=0.15, lw=0)
    ax[2].axhline(E["A0.6"]["near_top"]["N"]["0.5"][0], color="#e87ba4", lw=0.8, ls="--")
    ax[2].set_xscale("log"); ax[2].set_xlabel("grid spacing $h$ ($P=4$)")
    ax[2].xaxis.set_major_locator(FixedLocator([0.004, 0.008, 0.016, 0.031]))
    ax[2].xaxis.set_major_formatter(FixedFormatter(["0.004", "0.008", "0.016", "0.031"]))
    ax[2].xaxis.set_minor_locator(NullLocator())
    ax[2].set_ylim(-0.06, 0.25)
    ax[2].legend(fontsize=6.5, loc="lower left")
    ax[2].set_title("(c) $A=0.6$, large $\\tau$: $h\\to0$", fontsize=8.5, loc="left")
    fig.tight_layout()
    fig.savefig("fig6_torus.pdf", bbox_inches="tight")
    fig.savefig("fig6_torus.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def figure_fractal():
    """fig7: dimension of the near-critical set (fractal_top.json, lyap_leader.json)."""
    if not (os.path.exists("fractal_top.json") and os.path.exists("lyap_leader.json")):
        return
    F = json.load(open("fractal_top.json")); L = json.load(open("lyap_leader.json"))
    E = json.load(open("theory.json"))
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.8))
    shades = {24: "#9ec5f4", 32: "#3987e5", 40: "#104281"}
    for K in (24, 32, 40):
        r = F[f"A0.6_K{K}_d0.5"]
        ax[0].plot(r["eps"], r["slope"], color=shades[K], lw=1.6, label=f"$K={K}$")
    lam = L["0.6"]["slopes"]["crit"]; ln = L["0.6"]["slopes"]["near"]; c = E["A0.6"]["near_top"]["N"]["0.5"][0]
    ax[0].axhline(c / lam, color="#e87ba4", lw=1, ls="--")
    ax[0].axhline(c / ln, color="#999999", lw=1, ls=":")
    ax[0].text(2e-14, c / lam - 0.06, r"$c/\lambda^*$", fontsize=8, color="#e87ba4")
    ax[0].text(2e-14, c / ln + 0.03, r"$c/\lambda_{\rm near}$", fontsize=8, color="#777777")
    ax[0].set_xscale("log"); ax[0].set_ylim(0, 1.05); ax[0].invert_xaxis()
    ax[0].set_xlabel(r"box size $\varepsilon$ in the shooting parameter $a$")
    ax[0].set_ylabel(r"$-d\ln N(\varepsilon)/d\ln\varepsilon$")
    ax[0].legend(fontsize=7, loc="upper left")
    ax[0].set_title(r"(a) tri $A=0.6$: plateau = dimension $d$", fontsize=9, loc="left")
    xs, d_meas, d_lam, d_near, labels = [], [], [], [], []
    for x, (A, key, tag, K) in enumerate((("0.3", "A0.3_K48_d0.5", "A0.3", 48),
                                          ("0.6", "A0.6_K40_d0.5", "A0.6", 40),
                                          ("1.0", "A1_K28_d0.5", "A1", 28))):
        r = F[key]; eps = np.array(r["eps"]); sl = np.array(r["slope"])
        m = (eps <= 1e-2) & (eps >= 1e-4)
        c = E[tag]["near_top"]["N"]["0.5"][0]
        lam = L[A]["slopes"]["crit"]; ln = L[A]["slopes"]["near"]
        d0 = float(np.min(sl[m])); d1 = float(np.median(sl[m & (sl < d0 + 0.05)]))
        xs.append(x); d_meas.append((d0, d1)); d_lam.append(c / lam); d_near.append(c / ln)
        labels.append(f"$A={float(A):g}$\n$K={K}$")
    lo = [a for a, b in d_meas]; hi = [max(a, b) for a, b in d_meas]
    ax[1].bar(np.array(xs) - 0.2, np.array(hi), 0.35, color="#2a78d6", label="box-counting plateau $d$")
    ax[1].bar(np.array(xs) + 0.2, d_lam, 0.35, color="#e87ba4", label=r"$c/\lambda^*$ (critical orbit)")
    ax[1].plot(np.array(xs) + 0.2, d_near, ls="none", marker="_", ms=14, mew=2, color="#555555",
               label=r"$c/\lambda_{\rm near}$ (typical near-top)")
    ax[1].set_xticks(xs); ax[1].set_xticklabels(labels, fontsize=7)
    ax[1].set_ylim(0, 0.5)
    ax[1].set_ylabel("dimension")
    ax[1].legend(fontsize=6.5, loc="upper left")
    ax[1].set_title(r"(b) $c=d\,\lambda^*$ (Ledrappier--Young)", fontsize=9, loc="left")
    fig.tight_layout()
    fig.savefig("fig7_fractal.pdf", bbox_inches="tight")
    fig.savefig("fig7_fractal.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def analyse_strong():
    """Jump rates c, c_min and thresholds across pinning stiffness k = A
    (data_levels/A*_K*.npz, incl. the strong-pinning runs of queue_strong.sh)
    -> strong_pinning.json."""
    from annealed import f_ann as _fann
    M = json.load(open("minimax.json")) if os.path.exists("minimax.json") else {}
    out = {}
    for tag, A in (("A0.3", 0.3), ("A0.6", 0.6), ("A1", 1.0), ("A2", 2.0), ("A3", 3.0), ("A6", 6.0)):
        data = load(tag)
        if len(data) < 4:
            continue
        Ks = np.array(list(data))
        y = {"N": [], "Nmin": []}
        fc, fcsd = [], []
        for K in Ks:
            d = data[K]; x = d["f"] * K; Wm = d["fc"] * K
            j = np.searchsorted(x, Wm - 0.5 * A) - 1
            ok = j >= 0
            for key in y:
                n = d[key][np.arange(len(Wm))[ok], j[ok]]
                y[key].append(float(np.mean(np.log(np.maximum(n, 1)))))
            fc.append(float(d["fc"].mean())); fcsd.append(float(d["fc"].std()))
        c, cov = np.polyfit(Ks[-4:], y["N"][-4:], 1, cov=True)
        cm, covm = np.polyfit(Ks[-4:], y["Nmin"][-4:], 1, cov=True)
        Kl = int(Ks[-1])
        L1 = gle(1.0, A)
        out[tag] = dict(A=A, k=A, Lc=A ** (-2 / 3), K=Ks.tolist(), c=float(c[0]), c_err=float(np.sqrt(cov[0, 0])),
                        cmin=float(cm[0]), cmin_err=float(np.sqrt(covm[0, 0])), Lambda1=L1,
                        f_ann=_fann(A, 4.0, L1), fc_K=fc, fc_sd_K=fcsd,
                        fc_lead=(A - M[str(Kl)]["mean"]) if str(Kl) in M else None, near_top=y)
    json.dump(out, open("strong_pinning.json", "w"), indent=1)


def figure_strong():
    """fig8: from weak to strong pinning (strong_pinning.json, minimax.json)."""
    if not (os.path.exists("strong_pinning.json") and os.path.exists("minimax.json")):
        return
    S = json.load(open("strong_pinning.json")); M = json.load(open("minimax.json"))
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
    tags = [t for t in ("A0.3", "A0.6", "A1", "A2", "A3", "A6") if t in S]
    k = np.array([S[t]["k"] for t in tags]); c = np.array([S[t]["c"] for t in tags])
    ce = np.array([S[t]["c_err"] for t in tags]); cm = np.array([S[t]["cmin"] for t in tags])
    ax[0].errorbar(k, c, ce, color="#2a78d6", marker="o", ms=5, lw=1.6, capsize=0, label="$c$ (all equilibria)")
    ax[0].plot(k, cm, color="#2a78d6", marker="o", mfc="white", ms=5, lw=1, ls="--", label=r"$c_{\min}$ (minima)")
    kk = np.logspace(np.log10(0.2), np.log10(10), 100)
    ax[0].plot(kk, 0.22 * kk ** (2 / 3), color="#999999", lw=1, ls=":")
    ax[0].axhline(np.log(2), color="#e87ba4", lw=1, ls="--")
    ax[0].text(0.22, np.log(2) + 0.02, r"$\ln 2$ (strong pinning)", fontsize=7.5, color="#e87ba4")
    ax[0].text(3.2, 0.22 * 3.2 ** (2 / 3) + 0.08, r"$0.22/L_c$", fontsize=7.5, color="#777777")
    ax[0].set_xscale("log"); ax[0].set_ylim(0, 0.8)
    ax[0].set_xlabel(r"pinning stiffness $k=|V''|$ ($L_c=k^{-2/3}$)")
    ax[0].set_ylabel("rate per site")
    ax[0].legend(fontsize=7, loc="upper left", bbox_to_anchor=(0.0, 0.88))
    ax[0].set_title("(a) the jump: from $0.22/L_c$ to $\\ln 2$", fontsize=9, loc="left")
    Ks = np.array(sorted(int(x) for x in M)); m = np.array([M[str(x)]["mean"] for x in Ks])
    sd = np.array([M[str(x)]["sd"] for x in Ks])
    ax[1].plot(Ks, 2.0 - m, color="#2a78d6", marker="o", ms=4, lw=1.5, label=r"$P/2-\langle M^*\rangle$")
    ax[1].plot(Ks, sd, color="#eb6834", marker="s", ms=4, lw=1.5, label=r"sd of $M^*$ (sample threshold)")
    ax[1].plot(Ks, 4 / np.sqrt(12 * Ks), color="#eb6834", lw=0.8, ls="--")
    ax[1].plot(Ks, 2.0 / Ks, color="#2a78d6", lw=0.8, ls=":")
    ax[1].text(300, 4 / np.sqrt(12 * 300) * 1.5, r"$P/\sqrt{12K}$", fontsize=7.5, color="#eb6834")
    ax[1].text(300, 2.0 / 300 * 0.35, r"$2/K$", fontsize=7.5, color="#2a78d6")
    ax[1].set_xscale("log"); ax[1].set_yscale("log")
    ax[1].set_xlabel("$K$")
    ax[1].legend(fontsize=7, loc="lower left")
    ax[1].set_title(r"(b) $k\to\infty$: $f_c=A-M^*$, $M^*\to P/2$", fontsize=9, loc="left")
    fig.tight_layout()
    fig.savefig("fig8_strong.pdf", bbox_inches="tight")
    fig.savefig("fig8_strong.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def analyse_relax():
    """Exact long-chain thresholds (relax.json, K <= 1024): f_c(inf) from
    f_c(K) = f_c + a K^-x over K >= 64, with x = 3/4 central and x = 1/2, 1
    as the systematic band  -> relax_summary.json."""
    from scipy.optimize import curve_fit
    if not os.path.exists("relax.json"):
        return None
    R = json.load(open("relax.json"))
    out = {}
    for key, d in R.items():
        K = np.array([r["K"] for r in d["rows"]], float)
        m = np.array([r["mean"] for r in d["rows"]]); e = np.array([r["err"] for r in d["rows"]])
        s = K >= 64
        fits = {}
        for x in (0.5, 0.75, 1.0):
            p, c = curve_fit(lambda K, f, a: f + a * K ** -x, K[s], m[s], p0=(m[-1], 0.1),
                             sigma=e[s], absolute_sigma=True)
            fits[str(x)] = (float(p[0]), float(np.sqrt(c[0, 0])))
        f0 = fits["0.75"][0]
        lo = min(v[0] for v in fits.values()); hi = max(v[0] for v in fits.values())
        stat = fits["0.75"][1]
        out[key] = dict(A=d["A"], k=d["k"], K=K.tolist(), fc_K=m.tolist(), fc_K_err=e.tolist(),
                        sd_sqrtK=[r["sd"] * np.sqrt(r["K"]) for r in d["rows"]],
                        fc_inf=f0, fc_inf_lo=lo - stat, fc_inf_hi=hi + stat, fits=fits,
                        f_ann=d["f_ann"], ratio=d["f_ann"] / f0,
                        ratio_lo=d["f_ann"] / (hi + stat), ratio_hi=d["f_ann"] / (lo - stat))
    json.dump(out, open("relax_summary.json", "w"), indent=1)
    return out


def figure_relax():
    """fig9: exact thresholds of long chains across pinning strength
    (relax_summary.json, statistics.json, strong_scan.json)."""
    S = analyse_relax()
    if not S:
        return
    st = json.load(open("statistics.json"))["models"] if os.path.exists("statistics.json") else {}
    sc = json.load(open("strong_scan.json")) if os.path.exists("strong_scan.json") else {}
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
    cols = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8a5cd6", "#e87ba4"]
    for (key, d), col in zip(S.items(), cols):
        K = np.array(d["K"]); m = np.array(d["fc_K"]); e = np.array(d["fc_K_err"])
        ax[0].errorbar(K ** -0.75, m / d["fc_inf"], e / d["fc_inf"], color=col, marker="o", ms=3.5,
                       lw=1.2, capsize=0, label=f"$k={d['k']:g}$")
    ax[0].axhline(1, color="#999999", lw=0.8, ls=":")
    ax[0].set_xlim(0, 16 ** -0.75 * 1.05)
    ax[0].set_xlabel(r"$K^{-3/4}$")
    ax[0].set_ylabel(r"$\langle f_c(K)\rangle/f_c(\infty)$")
    ax[0].legend(fontsize=7, ncol=2, loc="upper left")
    ax[0].set_title("(a) exact thresholds, $K\\leq1024$", fontsize=9, loc="left")
    k = np.array([d["k"] for d in S.values()]); r = np.array([d["ratio"] for d in S.values()])
    rl = np.array([d["ratio_lo"] for d in S.values()]); rh = np.array([d["ratio_hi"] for d in S.values()])
    ax[1].errorbar(k, r, [r - rl, rh - r], color="#2a78d6", marker="o", ms=5, lw=1.6, capsize=0,
                   label="relaxation, $K\\to\\infty$", zorder=3)
    te = [(st[t]["A"], st[t]["f_ann_over_fc_freenu"], st[t]["f_ann_over_fc"]) for t in ("A0.3", "A0.6", "A1") if t in st]
    if te:
        a, r1, r2 = map(np.array, zip(*te))
        ax[1].vlines(a * 1.04, np.minimum(r1, r2), np.maximum(r1, r2), color="#eb6834", lw=2.5,
                     label="enumeration fits, $K\\leq64$")
    if sc:
        ks = np.array([v["k"] for v in sc.values()]); rs = np.array([v["f_ann"] / v["fc_inf"] for v in sc.values()])
        ax[1].plot(ks, rs, color="#1baf7a", marker="s", ms=4, lw=1, ls="--", label="strong-pinning minimax")
    kk = np.logspace(np.log10(8), np.log10(60), 60)
    from annealed import f_ann as _fann
    fa = np.array([_fann(x, 4.0, gle(1.0, x)) for x in kk])
    ax[1].plot(kk, fa / (kk - 2.0 + 3.6 / kk), color="#999999", lw=0.9, ls=":", label=r"$f_c^{\rm ann}/(A-\frac{P}{2}+\frac{0.90P}{k})$")
    ax[1].axhline(1.64, color="#e87ba4", lw=1, ls="--")
    ax[1].text(0.32, 1.655, "FLRT dynamics 1.64", fontsize=7.5, color="#e87ba4")
    ax[1].set_xscale("log"); ax[1].set_ylim(0.95, 1.75)
    ax[1].set_xlabel(r"pinning stiffness $k=|V''|$")
    ax[1].set_ylabel(r"$f_c^{\rm ann}/f_c$")
    ax[1].legend(fontsize=6.5, loc="lower left")
    ax[1].set_title("(b) annealed overshoot across strengths", fontsize=9, loc="left")
    fig.tight_layout()
    fig.savefig("fig9_relax.pdf", bbox_inches="tight")
    fig.savefig("fig9_relax.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def figure_replica():
    """fig10: replicated Kac-Rice, n = 2 (replica.json from replica.py, theory.json,
    level counts) -> replica_summary.json."""
    from scipy.optimize import brentq
    if not os.path.exists("replica.json"):
        return
    R = json.load(open("replica.json")); th = json.load(open("theory.json"))
    tags = {"0.3": "A0.3", "0.6": "A0.6", "1.0": "A1"}
    cols = {"0.3": "#2a78d6", "0.6": "#eb6834", "1.0": "#1baf7a"}
    summ = {}
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9))
    for key, rows in R.items():
        A = float(key); tag = tags[key]; col = cols[key]
        x = np.array([r["tauA"] for r in rows]); o = np.argsort(x); x = x[o]
        l2 = np.array([r["Lambda2"] for r in rows])[o]; e2 = np.array([r["Lambda2_err"] for r in rows])[o]
        L = np.array([r["Lambda"] for r in rows])[o]; L2t = np.array([r["Lambda_2tau"] for r in rows])[o]
        L1 = gle(1.0, A)
        xx = brentq(lambda y: y / np.tanh(y) - np.exp(L1), 1e-3, 20)
        ax[0].errorbar(x, l2 - 2 * L, e2, color=col, marker="o", ms=4, lw=1.4, capsize=0, label=f"$A={A:g}$")
        ax[0].plot(x, np.maximum(L2t - 2 * L, 0), color=col, lw=0.9, ls="--")
        ax[0].axvline(xx, color=col, lw=0.6, ls=":")
        # Gaussian closure: Var ln Z / K from the level counts
        data = load(tag); Ks = np.array(list(data)); s2 = []
        for xv in x:
            tau = max(xv, 1e-9) / A
            var = [lnZ(data[K]["N"].astype(float), data[K]["f"] * K, tau).var() for K in Ks]
            s2.append(float(np.polyfit(Ks[-4:], var[-4:], 1)[0]))
        s2 = np.array(s2)
        d = th[tag]; tq = np.array(d["tau"]) * A
        gap = np.interp(x, tq, d["gap"]); gerr = np.interp(x, tq, d["gap_err"])
        ax[1].errorbar(x, gap, gerr, color=col, marker="o", ms=4, lw=0, elinewidth=1, capsize=0)
        ax[1].plot(x, s2 / 2, color=col, lw=1.2, ls="-")
        ax[1].plot(x, (l2 - 2 * L) / 2, color=col, lw=1, ls="--")
        summ[key] = dict(A=A, tauA=x.tolist(), Lambda2=l2.tolist(), Lambda2_err=e2.tolist(), Lambda=L.tolist(),
                         Lambda_2tau=L2t.tolist(), excess=(l2 - 2 * L).tolist(), sigma2=s2.tolist(),
                         gap=gap.tolist(), gap_err=gerr.tolist(), x_L2=xx)
    ax[0].set_yscale("symlog", linthresh=1e-3); ax[0].set_ylim(0, None)
    ax[0].set_xlabel(r"$\tau A$"); ax[0].set_ylabel(r"$\Lambda_2(\tau)-2\Lambda(\tau)$")
    ax[0].legend(fontsize=7, loc="upper left")
    ax[0].set_title("(a) pair moment: excess over $\\langle Z\\rangle^2$", fontsize=9, loc="left")
    ax[1].set_xlim(0, 2.6); ax[1].set_ylim(-0.005, 0.25)
    ax[1].set_xlabel(r"$\tau A$"); ax[1].set_ylabel(r"$\Lambda-\varphi_q$")
    ax[1].plot([], [], color="#555555", marker="o", lw=0, label="measured gap")
    ax[1].plot([], [], color="#555555", lw=1.2, label=r"${\rm Var}\ln Z/2K$")
    ax[1].plot([], [], color="#555555", lw=1, ls="--", label=r"$(\Lambda_2-2\Lambda)/2$")
    ax[1].legend(fontsize=7, loc="upper left")
    ax[1].set_title("(b) replica Gaussian closure", fontsize=9, loc="left")
    fig.tight_layout()
    fig.savefig("fig10_replica.pdf", bbox_inches="tight")
    fig.savefig("fig10_replica.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    json.dump(summ, open("replica_summary.json", "w"), indent=1)


if __name__ == "__main__":
    import sys
    if "--torus-only" in sys.argv:
        figure_torus()
        figure_fractal()
        analyse_strong()
        figure_strong()
        figure_relax()
        figure_replica()
    elif "--relax-only" in sys.argv:
        figure_relax()
    elif "--replica-only" in sys.argv:
        figure_replica()
    else:
        main()
        figure_torus()
        figure_fractal()
        analyse_strong()
        figure_strong()
        figure_relax()
        figure_replica()
