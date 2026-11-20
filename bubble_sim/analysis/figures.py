"""Figures and text report shared by analyze_flat.py and analyze_bubbles.py.

A "set" is one run after data.load_* and measure.measure_run: a dict of per-wall arrays with
'label' and 'geometry' ('flat' | 'bubble'). Every figure has one row per set.

    configs.png     example critical walls at the largest L (+ b), the distribution of the angle
                    between the wall line and b for every L, and (bubbles) the mean shape
    roughness.png   pooled S_psi at the largest L, local zeta_eff(q), zeta vs L (intrinsic and
                    profile)
    threshold.png   <h_d>(L) with the fits, var(h_d)(L) -> nu, extrapolated h_c vs a
    hs.png          spectrum of the height along the arc length h(s) (flat: y(s), bubble: |r(s) - centroid|),
                    S_h ~ q^-(1+2 zeta): compensated for qEW / qKPZ, local zeta, zeta vs L
    length.png      wall length P vs size, P / size, fraction of columns / rays with overhangs
    rdensity.png    roughness density (Caballero & Giamarchi): -2 ln chi vs q^2 B, B(r) per a for
                    every L (dashed: exact q -> 0 limit), B at fixed r vs L
    length_dir.png  (bubbles) wall length by direction from b: Pdir / Sdir (contour over smooth
                    shape) vs R_eff at the edge flanks and screw tips, alpha per direction
"""
from __future__ import annotations

import os

import numpy as np

from . import rdensity as RD
from . import spectra as SP
from . import stats as ST
from .measure import ANG_BINS, DIR_BINS, N_SHAPE, N_SMOOTH

ZETA_REF = (("qEW", 1.25, "--"), ("qKPZ", 0.633, "-."))
REF_COLOR = {"qEW": "black", "qKPZ": "red"}     # apart from each other and from the plasma / viridis curves


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


def pool(W, kind, n_min, key=None, **k):
    """Pooled spectra at fixed (L, a): key 'psi' (tangent angle), 'prof' (profile u(x) / R(theta))
    or 'hs' (height along the arc length); default from kind: psi -> 'psi', u -> 'prof'."""
    key = key or ("psi" if kind == "psi" else "prof")
    cache = W.setdefault("_pools", {})
    ck = (key, k.get("L"), k.get("a"), None if k.get("mask") is None else k["mask"].tobytes())
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
            ax.axhline(z, color=REF_COLOR[name], ls=ls, lw=1.0, label=name)
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
                ax.axhline(zz, color=REF_COLOR[name], ls=ls, lw=1.0)
            ax.set_xscale("log"); ax.set_xlim(min(Ls) / 1.25, max(Ls) * 1.25)   # also when all ζ are nan
            ax.set_xlabel("L"); ax.set_ylabel("ζ  (o intrinsic, x profile)")
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
# h(s): height along the arc length
# ---------------------------------------------------------------------------
def hs_spectrum(sets, out, cfg, txt):
    """S_h(q) of h(s), s the arc length (single-valued with overhangs; dh/ds = sin psi, so
    S_h ≈ S_psi / q^2 only where the slopes are small). S_h ~ q^-(1+2 zeta)."""
    plt = _plt()
    qlg = cfg.get("qlarge")
    ncol = 5 if qlg else 4
    fig, axs = plt.subplots(len(sets), ncol, figsize=(5.7 * ncol, 4.6 * len(sets)), squeeze=False)
    for row, W in enumerate(sets):
        n_min, qmax = cfg[W["geometry"]]["n_min"], cfg["qmax"]
        As, Ls = a_values(W), L_values(W)
        Lm = max(Ls)
        hname = "y(s)" if W["geometry"] == "flat" else "|r(s) - centroid|"
        txt.append(f"== h(s) spectrum: {W['label']}   h = {hname} along the arc length s, S_h ~ q^-(1+2 zeta), "
                   f"modes n >= {n_min}: zeta for q <= {qmax:g}" + (f" | q <= {qlg:g}" if qlg else ""))
        spec = {a: pool(W, "u", n_min, key="hs", L=Lm, a=a) for a in As}
        spec = {a: (p, p.mean_err()) for a, p in spec.items() if p.n >= 2}
        for col, (rname, rz, _) in enumerate(ZETA_REF):
            ax, ex = axs[row, col], 1 + 2 * rz
            ys = []
            for a, (p, (m, e)) in spec.items():
                ax.errorbar(p.qc, m * p.qc ** ex, e * p.qc ** ex, fmt="o-", ms=2.5, lw=0.8, color=color_a(a),
                            label=f"a={a:g} (n={p.n})")
                y = m * p.qc ** ex
                ys.append(y[(p.qc <= qmax) & np.isfinite(y) & (y > 0)])
            ys = np.concatenate(ys) if ys else np.zeros(0)
            if ys.size:
                ax.set_ylim(ys.min() / 2, ys.max() * 2)
            ax.axvspan(qmax, 3, color="0.93", zorder=-1)
            if qlg:
                ax.axvline(qlg, color="0.5", ls="--", lw=0.7)
            ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("q (conjugate to the arc length s)")
            ax.set_ylabel(rf"$q^{{{ex:.2f}}} S_h$  (flat: {rname}, ζ = {rz:g})")
            ax.set_title(f"{W['label']}, L = {Lm}: h(s) compensated for {rname}", fontsize=9); ax.legend(fontsize=7)
        ax = axs[row, 2]
        for a, mk in ((As[0], "o"), (As[-1], "s")):
            for L in Ls:
                p = pool(W, "u", n_min, key="hs", L=L, a=a)
                if p.n < 4:
                    continue
                z, ze = SP.local_zeta(p, "u")
                ok = p.qc < 0.3
                ax.errorbar(p.qc[ok], z[ok], ze[ok], fmt=mk + ("-" if a == As[-1] else ":"), color=color_L(L, Ls),
                            ms=2.5, lw=0.7, label=f"a={a:g}, L={L}")
        for name, z, ls in ZETA_REF:
            ax.axhline(z, color=REF_COLOR[name], ls=ls, lw=1.0, label=name)
        ax.axvspan(qmax, 0.3, color="0.93", zorder=-1)
        if qlg:
            ax.axvline(qlg, color="0.5", ls="--", lw=0.7)
        ax.set_xscale("log"); ax.set_ylim(0, 2); ax.set_xlabel("q"); ax.set_ylabel(r"local $\zeta_{eff}$ from $S_h$")
        ax.set_title("local exponent (window: factor 2 in q)", fontsize=9); ax.legend(fontsize=5.5, ncol=2)
        for a in As:
            zs, ls_ = [], []
            for L in Ls:
                p = pool(W, "u", n_min, key="hs", L=L, a=a)
                if p.n < 4:
                    zs.append((np.nan, np.nan)); ls_.append((np.nan, np.nan))
                    continue
                z1 = SP.fit(p, qmax, "u")
                z2 = SP.fit(p, qlg, "u") if qlg else (np.nan, np.nan)
                zs.append(z1); ls_.append(z2)
                txt.append(f"   a={a:4.2f} L={L:5d} n={p.n:4d}: zeta_h {z1[0]:.3f}±{z1[1]:.3f}"
                           + (f" | q<={qlg:g}: {z2[0]:.3f}±{z2[1]:.3f}" if qlg else ""))
            for c, v in ((3, zs), (4, ls_)):
                if c < ncol:
                    z, ze = np.array(v).T
                    axs[row, c].errorbar(Ls, z, ze, fmt="o-", color=color_a(a), label=f"a={a:g}")
        for c, qq in ((3, qmax), (4, qlg)):
            if c >= ncol:
                continue
            ax = axs[row, c]
            for name, zz, ls in ZETA_REF:
                ax.axhline(zz, color=REF_COLOR[name], ls=ls, lw=1.0)
            ax.set_xscale("log"); ax.set_xlim(min(Ls) / 1.25, max(Ls) * 1.25)
            ax.set_xlabel("L"); ax.set_ylabel("ζ from S_h")
            ax.set_title(f"h(s) roughness exponent (n ≥ {n_min}, q ≤ {qq:g})" + (", large scales" if c == 4 else ""),
                         fontsize=9)
            ax.legend(fontsize=7)
    return save(fig, out, "hs.png")


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


# ---------------------------------------------------------------------------
# bubbles: length by direction
# ---------------------------------------------------------------------------
def length_dir(sets, out, txt):
    """P ∝ R^alpha is dominated by the roughest direction; split by the polar direction about
    the centroid (folded from b: 0 = edge flanks, 90 = screw tips). In each direction bin the
    excess Pdir / Sdir over the smooth shape (modes n < N_SMOOTH) goes as R^(alpha - 1):
    constant (alpha = 1) for zeta < 1 (qKPZ), R^(zeta - 1) for zeta > 1 (qEW)."""
    plt = _plt()
    sets = [W for W in sets if W["geometry"] == "bubble"]
    if not sets:
        return None
    fig, axs = plt.subplots(len(sets), 2, figsize=(12, 4.4 * len(sets)), squeeze=False)
    mid = 0.5 * (DIR_BINS[1:] + DIR_BINS[:-1])
    nb = len(mid)
    for row, W in enumerate(sets):
        txt.append(f"== length by direction: {W['label']}   excess Pdir / Sdir over the smooth shape "
                   f"(n < {N_SMOOTH}) ∝ R_eff^(alpha - 1); direction folded from b: "
                   f"0 = edge flanks, 90 = screw tips")
        txt.append("          bins: " + "  ".join(f"[{lo:g},{hi:g})" for lo, hi in zip(DIR_BINS[:-1], DIR_BINS[1:])))
        for a in a_values(W):
            m = sel(W, a=a)
            c = color_a(a)
            al = []
            for j in range(nb):
                r = W["Pdir"][m, j] / W["Sdir"][m, j]
                (sl, se), loc, tab = ST.length_scaling(W["size"][m], r)
                al.append((1 + sl, se, loc, tab))
            for j, ls, nm in ((0, "-", "edge"), (nb - 1, "--", "screw")):
                tab = al[j][3]
                axs[row, 0].errorbar(tab[:, 0], tab[:, 1], tab[:, 2], fmt="o" + ls, ms=3, color=c,
                                     label=f"a={a:g} {nm}: α = {al[j][0]:.3f}±{al[j][1]:.3f}")
            axs[row, 1].errorbar(mid + 0.6 * a, [x[0] for x in al], [x[1] for x in al], fmt="o-", color=c,
                                 label=f"a={a:g}")
            txt.append(f"   a={a:4.2f}: alpha by direction " + "  ".join(f"{x[0]:.3f}±{x[1]:.3f}" for x in al))
            for j, nm in ((0, "edge "), (nb - 1, "screw")):
                loc, tab = al[j][2], al[j][3]
                txt.append(f"          {nm}: " + "  ".join(f"R={t[0]:.0f}: {t[1]:.3f}" for t in tab))
                txt.append(f"          {nm} local alpha: " + " ".join(f"{1 + x:.2f}±{y:.2f}" for x, y in loc))
        ax = axs[row, 0]
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(r"$R_{eff}$")
        ax.set_ylabel("Pdir / Sdir")
        ax.set_title(f"{W['label']}: excess length, edge flanks [0,{DIR_BINS[1]:g})° (solid), "
                     f"screw tips [{DIR_BINS[-2]:g},90]° (dashed)", fontsize=8)
        ax.legend(fontsize=6)
        ax = axs[row, 1]
        for name, z, ls in ZETA_REF:
            ax.axhline(max(1.0, z), color=REF_COLOR[name], ls=ls, lw=1.0, label=f"{name}: α = {max(1.0, z):g}")
        ax.set_xlabel("direction from b (deg): 0 = edge flank, 90 = screw tip"); ax.set_ylabel("α")
        ax.set_title("α per direction (1 + slope of Pdir / Sdir)", fontsize=9); ax.legend(fontsize=7)
    return save(fig, out, "length_dir.png")


# ---------------------------------------------------------------------------
# roughness density (Caballero & Giamarchi)
# ---------------------------------------------------------------------------
R_FIX = 32          # B(R_FIX) - B(0) vs L: ∝ L^(2 zeta - 2) for zeta > 1 (anomalous), flat for zeta < 1
# large-r window of the local slope, as fractions of rmax (flat: L/2; bubble: pi R_eff, where B
# saturates above ~2 pi R_eff / (2 N_SMOOTH) since the smooth shape is subtracted)
R_WIN = dict(flat=(1 / 16, 1 / 4), bubble=(1 / 64, 1 / 16))


def rd_pool(W, L, a, sector=None):
    """chi and Bdir averaged over the walls at (L, a), with r in lattice units (bubbles: arc
    length R_eff dtheta at the mean R_eff); sector: 'edge' / 'screw' (bubbles) for the chi
    averaged over the pairs of rays in that sector only."""
    idx = np.nonzero(sel(W, L, a))[0]
    if not len(idx):
        return None
    sfx = f"_{sector}" if sector else ""
    chi = np.mean([W["chi" + sfx][i] for i in idx], 0)
    bdir = np.mean([W["Bdir" + sfx][i] for i in idx], 0)
    rg = W["rd_r"][idx[0]]
    scale = 1.0 if W["geometry"] == "flat" else 2 * np.pi * W["size"][idx].mean() / W["rd_n"][idx[0]]
    B, y = RD.fit_B(chi)
    return dict(r=rg * scale, B=B, Bdir=bdir, y=y, n=len(idx))


def rdensity(sets, out, txt, sector=None):
    """sector (bubbles): 'edge' or 'screw' -> rdensity_<sector>.png from the sector-restricted chi."""
    plt = _plt()
    if sector:
        sets = [W for W in sets if W["geometry"] == "bubble"]
        if not sets:
            return None
    sname = {"edge": ", edge flanks only (directions [0, 30)° from b)",
             "screw": ", screw tips only (directions [60, 90]° from b)"}.get(sector, "")
    ncol = 2 + max(len(a_values(W)) for W in sets)
    fig, AX = plt.subplots(2 * len(sets), ncol, figsize=(4.6 * ncol, 3.9 * 2 * len(sets)), squeeze=False)
    for row, W in enumerate(sets):
        axs, axl = AX[2 * row:2 * row + 1], AX[2 * row + 1]
        axl[0].axis("off"); axl[-1].axis("off")
        As, Ls = a_values(W), L_values(W)
        flat = W["geometry"] == "flat"
        txt.append(f"== roughness density (Caballero-Giamarchi): {W['label']}{sname}   B(r) = slope of -2 ln chi(q, r) "
                   f"vs q^2 at q -> 0 (= exact <Δu^2> + B(0), u the signed height); "
                   f"local slope d lnB/d ln r of B - B(0) (2 zeta for zeta < 1; 2 for zeta > 1); "
                   f"B({R_FIX}) - B(0) vs L ∝ L^(2 zeta - 2) for zeta > 1; B(0) = overhangs")
        pools = {(L, a): rd_pool(W, L, a, sector) for L in Ls for a in As}
        # -2 ln chi vs q^2 B at the largest L of the largest a
        a0 = As[-1]
        L0 = max(L for L in Ls if pools[(L, a0)] is not None)
        p = pools[(L0, a0)]
        ax = axs[0, 0]
        for j in np.unique(np.linspace(1, len(p["r"]) - 1, 5).astype(int)):
            if np.isfinite(p["B"][j]) and p["B"][j] > 0:
                ax.plot(RD.QS ** 2 * p["B"][j], p["y"][:, j], ".-", ms=3, lw=0.7, label=f"r = {p['r'][j]:.0f}")
        xx = np.geomspace(1e-6, 10, 10)
        ax.plot(xx, xx, "k--", lw=0.8, label="Gaussian: q² B")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1e-6, 1e3); ax.set_ylim(1e-6, 30)
        ax.set_xlabel("q² B(r)"); ax.set_ylabel("-2 ln χ(q, r)")
        ax.set_title(f"{W['label']}{sname}\na = {a0:g}, L = {L0}", fontsize=8); ax.legend(fontsize=6)
        for k, a in enumerate(As):
            ax = axs[0, 1 + k]
            fixed, ref = [], None
            for L in Ls:
                p = pools[(L, a)]
                if p is None:
                    continue
                c = color_L(L, Ls)
                ok = p["r"] > 0
                b0 = p["B"][0]
                d = p["B"] - b0
                rr = p["r"]
                ax.plot(rr[ok], d[ok], "-", color=c, lw=1.3, label=f"L={L} (n={p['n']})")
                ax.plot(rr[ok], p["B"][ok], "-", color=c, lw=0.6, alpha=0.35)
                axl[1 + k].plot(rr[ok], RD.local_slopes(rr, d)[ok], "o-", ms=2.5, lw=1, color=c, label=f"L={L}")
                ok &= np.isfinite(p["B"])
                rmax = rr[ok].max()
                sl_s = RD.local_slope(rr, d, 2, 16)
                sl_m = RD.local_slope(rr, d, 16, 128)
                lo, hi = (f * rmax for f in R_WIN[W["geometry"]])
                sl_l = RD.local_slope(rr, d, lo, hi)
                bf = np.exp(np.interp(np.log(R_FIX), np.log(rr[ok]), np.log(np.maximum(d[ok], 1e-30))))
                fixed.append((L, bf))
                ref = (bf, rr[ok].min(), rmax)
                txt.append(f"   a={a:4.2f} L={L:5d} n={p['n']:4d}: B(0) = {b0:.3g} | slope r 2-16: {sl_s:.2f}  "
                           f"16-128: {sl_m:.2f}  {lo:.0f}-{hi:.0f}: {sl_l:.2f} | B({R_FIX}) - B(0) = {bf:.4g} | "
                           f"max |B - Bdir| / Bdir = {np.nanmax(np.abs(p['B'] - p['Bdir'])[ok] / p['Bdir'][ok]):.1e}")
            if len(fixed) > 1:
                Lf, bf = np.array(fixed).T
                e = ST.loglog_slope(Lf, bf)
                txt.append(f"          B({R_FIX}) - B(0) ∝ L^{e:.2f}  -> zeta = {1 + e / 2:.2f} if anomalous (zeta > 1); "
                           f"local " + " ".join(f"{np.log(bf[i + 1] / bf[i]) / np.log(Lf[i + 1] / Lf[i]):.2f}"
                                                 for i in range(len(Lf) - 1)))
                axs[0, -1].plot(Lf, bf, "o-", color=color_a(a), label=f"a={a:g}: ∝ L^{e:.2f}")
            al = axl[1 + k]
            for name, z, ls in ZETA_REF:
                ex = min(2 * z, 2.0)
                if ref is not None:     # through B - B(0) of the largest L at r = R_FIX
                    r0 = np.geomspace(max(1.0, ref[1]), ref[2] / 2, 20)
                    ax.plot(r0, ref[0] * (r0 / R_FIX) ** ex, color=REF_COLOR[name], ls=ls, lw=1.6,
                            label=f"r^{ex:.3g} ({name}{', anomalous' if z > 1 else ''})")
                al.axhline(ex, color=REF_COLOR[name], ls=ls, lw=1.6, label=f"{ex:.3g} ({name})")
            al.set_xscale("log"); al.set_ylim(0, 2.4); al.set_xlabel("r" if flat else "r = R_eff Δθ")
            al.set_ylabel("d ln[B - B(0)] / d ln r")
            al.set_title(f"a = {a:g}: local slope (window ×2)", fontsize=8); al.legend(fontsize=6)
            al.grid(alpha=0.3)
            ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("r" if flat else "r = R_eff Δθ")
            ax.set_ylabel("B(r) - B(0)")
            ax.set_title(f"a = {a:g}: B(r) - B(0) from χ (faint: B(r), with overhangs)", fontsize=8)
            ax.legend(fontsize=6)
        for k in range(len(As), ncol - 2):
            axs[0, 1 + k].axis("off"); axl[1 + k].axis("off")
        ax = axs[0, -1]
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("L"); ax.set_ylabel(f"B({R_FIX}) - B(0)")
        ax.set_title(f"B({R_FIX}) - B(0): ∝ L^(2ζ-2) if ζ > 1, flat if ζ < 1", fontsize=8); ax.legend(fontsize=7)
    return save(fig, out, f"rdensity_{sector}.png" if sector else "rdensity.png")
