"""Figures and text report shared by analyze_flat.py and analyze_bubbles.py.

A "set" is one run after data.load_* and measure.measure_run: a dict of per-wall arrays with
'label' and 'geometry' ('flat' | 'bubble'). Every figure has one row per set.

    configs.png     example critical walls at the largest L (+ b), the distribution of the angle
                    between the wall line and b for every L, and (bubbles) the mean shape
    roughness.png   pooled S_psi at the largest L, local zeta_eff(q), zeta vs L (intrinsic and
                    profile)
    threshold.png   <h_d>(L) with the fits, var(h_d)(L) -> nu, extrapolated h_c vs a
    length.png      wall length P vs size, P / size, fraction of columns / rays with overhangs
"""
from __future__ import annotations

import os

import numpy as np

from . import spectra as SP
from . import stats as ST
from .measure import ANG_BINS, N_SHAPE

ZETA_REF = (("qEW", 1.25, "--"), ("qKPZ", 0.633, "-."))


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _cm():
    from matplotlib import cm
    return cm


def a_values(W):
    return sorted({float(a) for a in W["aniso"]})


def L_values(W):
    return sorted({int(L) for L in W["L"]})


def color_a(a):
    return _cm().viridis(0.9 * a)


def color_L(L, Ls):
    return _cm().plasma(0.85 * Ls.index(L) / max(1, len(Ls) - 1))


def sel(W, L=None, a=None, mask=None):
    m = np.ones(len(W["L"]), bool)
    if L is not None:
        m &= W["L"] == L
    if a is not None:
        m &= W["aniso"] == a
    if mask is not None:
        m &= mask
    return m


def pool(W, kind, n_min, **k):
    key = "psi" if kind == "psi" else "prof"
    cache = W.setdefault("_pools", {})
    ck = (kind, k.get("L"), k.get("a"), None if k.get("mask") is None else k["mask"].tobytes())
    if ck not in cache:
        cache[ck] = SP.Pooled([W[key][i] for i in np.nonzero(sel(W, **k))[0]], n_min)
    return cache[ck]


def size_name(W):
    return "L" if W["geometry"] == "flat" else "R_eff"


def save(fig, out, name):
    fig.tight_layout()
    fn = os.path.join(out, name)
    fig.savefig(fn, dpi=100)
    _plt().close(fig)
    return fn


# ---------------------------------------------------------------------------
# configurations
# ---------------------------------------------------------------------------
def configs(sets, out, n_show=4):
    plt = _plt()
    rows = []
    for W in sets:
        for a in a_values(W):
            if W["geometry"] == "flat" and a == 0 and any(r[0]["geometry"] == "flat" and r[1] == 0 for r in rows):
                continue                                    # isotropic flat walls: shown once
            rows.append((W, a))
    extra = 2 if any(W["geometry"] == "bubble" for W in sets) else 1
    fig = plt.figure(figsize=(2.6 * n_show + 4.2 * extra, 3.9 * len(rows)))
    gs = fig.add_gridspec(len(rows), n_show + extra, width_ratios=[1] * n_show + [1.7] * extra)
    for r, (W, a) in enumerate(rows):
        Lm = max(L_values(W))
        idx = np.nonzero(sel(W, Lm, a))[0][:n_show]
        b = np.deg2rad(W["burgers_angle"])
        for c in range(n_show):
            ax = fig.add_subplot(gs[r, c])
            ax.set_xticks([]); ax.set_yticks([])
            if c >= len(idx):
                ax.axis("off"); continue
            i = idx[c]
            cv = W["curve"][i]
            if W["geometry"] == "flat":
                Ly = W["Ly"][i]
                ax.fill(np.r_[cv[:, 0], Lm, Lm, 0, 0], np.r_[cv[:, 1], cv[0, 1], 0, 0, cv[0, 1]], color="0.6", lw=0)
                ax.plot(cv[:, 0], cv[:, 1], color=color_a(a), lw=0.6)
                lim = (0, Lm, 0, Ly)
            else:
                ax.fill(cv[:, 0], cv[:, 1], color="0.6", lw=0)
                ax.plot(cv[:, 0], cv[:, 1], color=color_a(a), lw=0.6)
                # the bubble with a margin (small fixed-R0 bubbles would be dots in the box)
                x0, x1, y0, y1 = cv[:, 0].min(), cv[:, 0].max(), cv[:, 1].min(), cv[:, 1].max()
                h = 0.6 * max(x1 - x0, y1 - y0) + 4
                xc, yc = (x0 + x1) / 2, (y0 + y1) / 2
                lim = (max(0, xc - h), min(Lm, xc + h), max(0, yc - h), min(Lm, yc + h))
                ax.text(0.02, 0.02, f"{lim[1] - lim[0]:.0f} of {Lm}", transform=ax.transAxes, fontsize=6)
            ax.set_xlim(lim[0], lim[1]); ax.set_ylim(lim[2], lim[3]); ax.set_aspect("equal")
            ln = 0.12 * (lim[1] - lim[0])
            x0, y0 = lim[0] + 0.15 * (lim[1] - lim[0]), lim[2] + 0.1 * (lim[3] - lim[2])
            ax.annotate("", xy=(x0 + ln * np.cos(b), y0 + ln * np.sin(b)), xytext=(x0, y0),
                        arrowprops=dict(arrowstyle="->", color="red", lw=1.4))
            ax.set_title(f"h_d={W['h_d'][i]:.4f}", fontsize=8)
            if c == 0:
                ax.set_ylabel(f"{W['label']}\na = {a:g}, L = {Lm}", fontsize=8)
        # orientation of the wall line relative to b
        ax = fig.add_subplot(gs[r, n_show])
        Ls = L_values(W)
        cen = 0.5 * (ANG_BINS[1:] + ANG_BINS[:-1])
        for L in Ls:
            m = sel(W, L, a)
            if m.sum():
                ax.plot(cen, W["ang"][m].mean(0) / np.diff(ANG_BINS), "o-", ms=3, color=color_L(L, Ls),
                        label=f"L={L} (n={m.sum()})")
        ax.set_xlim(0, 90); ax.set_xticks([0, 30, 60, 90])
        ax.set_xlabel("angle wall line – b (deg)\n0: screw-like, 90: edge-like", fontsize=8)
        ax.set_ylabel("P (per degree)"); ax.legend(fontsize=6)
        if W["geometry"] == "bubble":
            ax = fig.add_subplot(gs[r, n_show + 1], projection="polar")
            th = 2 * np.pi * (np.arange(N_SHAPE) + 0.5) / N_SHAPE
            for L in Ls:
                m = sel(W, L, a)
                if m.sum():
                    s = W["shape"][m].mean(0)
                    ax.plot(np.r_[th, th[:1]], np.r_[s, s[:1]], color=color_L(L, Ls), lw=1)
            ax.plot([0, 0], [0, 1.3], color="red", lw=1)
            ax.set_ylim(0, 1.35); ax.set_yticks([0.5, 1.0]); ax.tick_params(labelsize=6)
            ax.set_title("<R(θ)>/R_eff, θ from b", fontsize=8)
    fig.suptitle("critical walls (grey: φ > 0 domain; red: Burgers vector b)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fn = os.path.join(out, "configs.png")
    fig.savefig(fn, dpi=100)
    plt.close(fig)
    return fn


# ---------------------------------------------------------------------------
# roughness
# ---------------------------------------------------------------------------
def roughness(sets, out, cfg, txt):
    plt = _plt()
    qlg = cfg.get("qlarge")          # large-scale window q <= qlarge (flat walls: past the qEW -> qKPZ crossover)
    ncol = 5 if qlg else 4
    fig, axs = plt.subplots(len(sets), ncol, figsize=(5.7 * ncol, 4.6 * len(sets)), squeeze=False)
    for row, W in enumerate(sets):
        n_min, qmax = cfg[W["geometry"]]["n_min"], cfg["qmax"]
        As, Ls = a_values(W), L_values(W)
        Lm = max(Ls)
        txt.append(f"== roughness: {W['label']}")
        txt.append(f"   zeta from pooled spectra, modes n >= {n_min} up to q = {qmax:g}; "
                   f"intrinsic S_psi ~ q^(1-2 zeta) | profile ({'u(x)' if W['geometry'] == 'flat' else 'R(theta)'}) "
                   f"S ~ q^-(1+2 zeta) | intrinsic split at the median travel (memory)"
                   + (f" | large scales q <= {qlg:g}: intrinsic, profile" if qlg else ""))
        # S_psi ~ q^(1 - 2 zeta): q^(2 zeta - 1) S_psi is flat for qEW (col 0) and qKPZ (col 1)
        spec = {a: pool(W, "psi", n_min, L=Lm, a=a) for a in As}
        spec = {a: (p, p.mean_err()) for a, p in spec.items() if p.n >= 2}
        for col, (rname, rz, _) in enumerate(ZETA_REF):
            ax, ex = axs[row, col], 2 * rz - 1
            ys = []
            for a, (p, (m, e)) in spec.items():
                ax.errorbar(p.qc, m * p.qc ** ex, e * p.qc ** ex, fmt="o-", ms=2.5, lw=0.8, color=color_a(a),
                            label=f"a={a:g} (n={p.n})")
                y = m * p.qc ** ex
                ys.append(y[(p.qc <= qmax) & np.isfinite(y) & (y > 0)])
            ys = np.concatenate(ys) if ys else np.zeros(0)
            if ys.size:                    # y range from the fit region: the plateau, not the lattice cutoff
                ax.set_ylim(ys.min() / 2, ys.max() * 2)
            ax.axvspan(qmax, 3, color="0.93", zorder=-1)
            if qlg:
                ax.axvline(qlg, color="0.5", ls="--", lw=0.7)
            ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("q (arc length)")
            ax.set_ylabel(rf"$q^{{{ex:.2f}}} S_\psi$  (flat: {rname}, ζ = {rz:g})")
            ax.set_title(f"{W['label']}, L = {Lm}: compensated for {rname}", fontsize=9); ax.legend(fontsize=7)
        ax = axs[row, 2]
        for a, mk in ((As[0], "o"), (As[-1], "s")):
            for L in Ls:
                p = pool(W, "psi", n_min, L=L, a=a)
                if p.n < 4:
                    continue
                z, ze = SP.local_zeta(p, "psi")
                ok = p.qc < 0.3
                ax.errorbar(p.qc[ok], z[ok], ze[ok], fmt=mk + ("-" if a == As[-1] else ":"), color=color_L(L, Ls),
                            ms=2.5, lw=0.7, label=f"a={a:g}, L={L}")
        for name, z, ls in ZETA_REF:
            ax.axhline(z, color="k", ls=ls, lw=0.7, label=name)
        ax.axvspan(qmax, 0.3, color="0.93", zorder=-1)
        if qlg:
            ax.axvline(qlg, color="0.5", ls="--", lw=0.7)
        ax.set_xscale("log"); ax.set_ylim(0, 2); ax.set_xlabel("q"); ax.set_ylabel(r"local $\zeta_{eff}$ (intrinsic)")
        ax.set_title("local exponent (window: factor 2 in q)", fontsize=9); ax.legend(fontsize=5.5, ncol=2)
        ax = axs[row, 3]
        for a in As:
            zi, zp, li, lp = [], [], [], []
            for L in Ls:
                p = pool(W, "psi", n_min, L=L, a=a)
                if p.n < 4:
                    for z_ in (zi, zp, li, lp):
                        z_.append((np.nan, np.nan))
                    continue
                z1 = SP.fit(p, qmax, "psi")
                z2 = SP.fit(pool(W, "prof", n_min, L=L, a=a), qmax, "u")
                tr = W["travel"]
                med = np.nanmedian(tr[sel(W, L, a)])
                zlo = SP.fit(pool(W, "psi", n_min, L=L, a=a, mask=tr < med), qmax, "psi")
                zhi = SP.fit(pool(W, "psi", n_min, L=L, a=a, mask=tr >= med), qmax, "psi")
                zi.append(z1); zp.append(z2)
                line = (f"   a={a:4.2f} L={L:5d} n={p.n:4d}: intrinsic {z1[0]:.3f}±{z1[1]:.3f} | profile "
                        f"{z2[0]:.3f}±{z2[1]:.3f} | travel<{med:.0f}: {zlo[0]:.3f}±{zlo[1]:.3f}, "
                        f">: {zhi[0]:.3f}±{zhi[1]:.3f}")
                if qlg:
                    l1 = SP.fit(p, qlg, "psi")
                    l2 = SP.fit(pool(W, "prof", n_min, L=L, a=a), qlg, "u")
                    li.append(l1); lp.append(l2)
                    line += f" | q<={qlg:g}: {l1[0]:.3f}±{l1[1]:.3f}, {l2[0]:.3f}±{l2[1]:.3f}"
                txt.append(line)
            panels = [(axs[row, 3], zi, zp)] + ([(axs[row, 4], li, lp)] if qlg else [])
            for ax_, a1, a2 in panels:
                z, ze = np.array(a1).T
                ax_.errorbar(Ls, z, ze, fmt="o-", color=color_a(a), label=f"a={a:g}")
                z, ze = np.array(a2).T
                ax_.errorbar(Ls, z, ze, fmt="x:", color=color_a(a), alpha=0.6)
        for c, qq in ((3, qmax), (4, qlg)):
            if c >= ncol:
                continue
            ax = axs[row, c]
            for name, zz, ls in ZETA_REF:
                ax.axhline(zz, color="k", ls=ls, lw=0.7)
            ax.set_xscale("log"); ax.set_xlabel("L"); ax.set_ylabel("ζ  (o intrinsic, x profile)")
            ax.set_title(f"roughness exponent (n ≥ {n_min}, q ≤ {qq:g})" + (", large scales" if c == 4 else ""),
                         fontsize=9)
            ax.legend(fontsize=7)
    return save(fig, out, "roughness.png")


# ---------------------------------------------------------------------------
# threshold
# ---------------------------------------------------------------------------
def threshold(sets, out, cfg, txt):
    plt = _plt()
    fig, axs = plt.subplots(len(sets), 3, figsize=(17, 4.4 * len(sets)), squeeze=False)
    for row, W in enumerate(sets):
        As = a_values(W)
        model = W["geometry"]
        lmin = cfg[W["geometry"]]["lmin"]
        txt.append(f"== threshold: {W['label']}   (<h_d> = h_c + A L^(-1/nu)"
                   f"{' + B/L' if model == 'bubble' else ''}, L >= {lmin})")
        hcs = {}
        for a in As:
            m = sel(W, a=a)
            rows, fits = ST.threshold(W["L"][m], W["h_d"][m], model, lmin)
            Ls = np.array([r["L"] for r in rows])
            mean = np.array([r["mean"] for r in rows]); err = np.array([r["err"] for r in rows])
            var = np.array([r["var"] for r in rows]); verr = np.array([r["var_err"] for r in rows])
            c = color_a(a)
            axs[row, 0].errorbar(Ls, mean, err, fmt="o", color=c, label=f"a={a:g}")
            axs[row, 1].errorbar(Ls, var, verr, fmt="o-", color=c, label=f"a={a:g}")
            txt.append(f"   a={a:4.2f}: " + "  ".join(f"L={r['L']}: {r['mean']:.5f}±{r['err']:.5f}" for r in rows))
            line = []
            if "nu_var" in fits:
                line.append(f"nu(var) = {fits['nu_var'][0]:.2f}±{fits['nu_var'][1]:.2f}")
            for name in ("qEW", "qKPZ", "free"):
                f = fits.get(f"hc_{name}")
                if f is None:
                    continue
                hcs.setdefault(name, []).append((a, f["hc"], f["hc_err"]))
                tag = name if f["nu"] else f"free, 1/nu = {f['inv_nu']:.2f}±{f['inv_nu_err']:.2f}"
                line.append(f"h_c[{tag}] = {f['hc']:.5f}±{f['hc_err']:.5f} (chi2/dof {f['chi2']:.1f})")
                if name == "qKPZ" or (name == "qEW" and "hc_qKPZ" not in fits):
                    Lx = np.geomspace(Ls.min(), Ls.max() * 4, 50)
                    axs[row, 0].plot(Lx, ST.eval_fit(f, Lx), ls="--" if name == "qKPZ" else ":", color=c, lw=0.8)
            txt.append("          " + ";  ".join(line))
        ax = axs[row, 0]
        ax.set_xscale("log"); ax.set_xlabel("L"); ax.set_ylabel("<h_d>")
        ax.set_title(f"{W['label']}: mean threshold (dashed: fit, ν = 1.733)", fontsize=9); ax.legend(fontsize=7)
        ax = axs[row, 1]
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("L"); ax.set_ylabel("var(h_d)")
        ax.set_title("threshold fluctuations ~ L^(-2/ν)  (qEW ν = 4/3, qKPZ 1.733)", fontsize=9); ax.legend(fontsize=7)
        ax = axs[row, 2]
        for (name, v), mk in zip(hcs.items(), "osd"):
            a, h, e = np.array(v).T
            ax.errorbar(a + 0.01 * "osd".index(mk), h, e, fmt=mk + "-", label=f"ν {name}")
        ax.set_xlabel("anisotropy a"); ax.set_ylabel("extrapolated h_c")
        ax.set_title("h_c (L → ∞)", fontsize=9); ax.legend(fontsize=7)
    return save(fig, out, "threshold.png")


# ---------------------------------------------------------------------------
# length and overhangs
# ---------------------------------------------------------------------------
def length(sets, out, txt):
    plt = _plt()
    fig, axs = plt.subplots(len(sets), 3, figsize=(17, 4.4 * len(sets)), squeeze=False)
    for row, W in enumerate(sets):
        flat = W["geometry"] == "flat"
        As, Ls = a_values(W), L_values(W)
        txt.append(f"== length: {W['label']}   P ∝ {size_name(W)}^alpha "
                   f"({'P/L' if flat else 'P / 2 pi R_eff'}; overhangs = fraction of "
                   f"{'columns' if flat else 'rays'} crossing the wall more than once)")
        for a in As:
            m = sel(W, a=a)
            (al, ale), loc, tab = ST.length_scaling(W["size"][m], W["P"][m])
            c = color_a(a)
            norm = tab[:, 0] if flat else 2 * np.pi * tab[:, 0]
            axs[row, 0].errorbar(tab[:, 0], tab[:, 1], tab[:, 2], fmt="o-", ms=3, color=c,
                                 label=f"a={a:g}: α = {al:.3f}±{ale:.3f}")
            axs[row, 1].errorbar(tab[:, 0], tab[:, 1] / norm, tab[:, 2] / norm, fmt="o-", ms=3, color=c, label=f"a={a:g}")
            mv = [(W["multi"][sel(W, L, a)].mean(), ST.sem(W["multi"][sel(W, L, a)])) for L in Ls if sel(W, L, a).sum()]
            Lk = [L for L in Ls if sel(W, L, a).sum()]
            axs[row, 2].errorbar(Lk, *np.array(mv).T, fmt="o-", color=c, label=f"a={a:g}")
            txt.append(f"   a={a:4.2f}: alpha = {al:.3f}±{ale:.3f}; local "
                       + " ".join(f"{x:.2f}±{y:.2f}" for x, y in loc))
            txt.append("          " + "  ".join(f"{size_name(W)}={s:.0f}: {p / n:.3f}" for s, p, n in
                                                  zip(tab[:, 0], tab[:, 1], norm)))
            txt.append("          overhangs: " + "  ".join(f"L={L}: {v:.3f}" for L, (v, _) in zip(Lk, mv)))
        x = "L" if flat else r"$R_{eff} = \sqrt{A/\pi}$"
        for k, yl in enumerate(("P", "P / L" if flat else r"$P / 2\pi R_{eff}$")):
            ax = axs[row, k]
            ax.set_xscale("log"); ax.set_xlabel(x); ax.set_ylabel(yl)
            if k == 0:
                ax.set_yscale("log")
            ax.set_title(f"{W['label']}: wall length", fontsize=9); ax.legend(fontsize=7)
        ax = axs[row, 2]
        ax.set_xscale("log"); ax.set_xlabel("L"); ax.set_ylabel("fraction multivalued")
        ax.set_title("overhangs (KFR)", fontsize=9); ax.legend(fontsize=7)
    return save(fig, out, "length.png")
