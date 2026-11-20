"""Bubble figures from saved contours, with the flat-interface estimator."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter
import numpy as np

from make_figures import decode, common_bins, spectral_fit, COLORS, MARKERS

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DATA = ROOT / "runs/fss_D0.04/contours.npz"
SIZES = (1024, 2048, 4096, 8192)
RNG = np.random.default_rng(20261123)


def cross(a, b):
    return a[..., 0]*b[..., 1] - a[..., 1]*b[..., 0]


def profiles(curve):
    """Closed contour, made counter-clockwise, with positive phase inside."""
    line = np.asarray(curve, float)
    if np.allclose(line[0], line[-1]):
        line = line[:-1]
    area = .5*np.sum(cross(line, np.roll(line, -1, axis=0)))
    if area < 0:
        line, area = line[::-1], -area
    line = np.vstack((line, line[:1]))
    increments = np.linalg.norm(np.diff(line, axis=0), axis=1)
    keep = np.r_[True, increments > 0]
    line = line[keep]
    c = cross(line[:-1], line[1:])
    centre = np.array([np.sum((line[:-1, 0]+line[1:, 0])*c),
                       np.sum((line[:-1, 1]+line[1:, 1])*c)])/(6*area)
    R = np.sqrt(area/np.pi)
    length = np.r_[0., np.cumsum(np.linalg.norm(np.diff(line, axis=0), axis=1))]
    P = length[-1]
    N = int(np.ceil(P))  # ds <= one lattice spacing
    s = np.arange(N) * P / N
    hs = np.hypot(np.interp(s, length, line[:, 0]) - centre[0],
                  np.interp(s, length, line[:, 1]) - centre[1])
    hs -= hs.mean()
    Ss = P * np.abs(np.fft.rfft(hs) / N)**2
    qs = 2*np.pi*np.fft.rfftfreq(N, d=P/N)
    np.testing.assert_allclose(np.mean(hs**2),
        (2*Ss[1:].sum() - (Ss[-1] if N % 2 == 0 else 0))/P)

    # Rays from the centroid at theta_k = 2 pi (k+1/2)/M, about one per lattice spacing at R.
    M = int(np.ceil(2*np.pi*R))
    d = line - centre
    angle = np.arctan2(d[:, 1], d[:, 0])
    turn = (np.diff(angle) + np.pi) % (2*np.pi) - np.pi
    a0 = angle[:-1]*M/(2*np.pi) - .5
    a1 = a0 + turn*M/(2*np.pi)
    low, high = np.minimum(a0, a1), np.maximum(a0, a1)
    k0 = np.ceil(low).astype(np.int64)
    count = np.maximum(np.ceil(high).astype(np.int64) - k0, 0)
    seg = np.repeat(np.arange(len(a0)), count)
    k = np.repeat(k0, count) + np.arange(count.sum()) - np.repeat(np.cumsum(count)-count, count)
    theta = 2*np.pi*(k+.5)/M
    e = np.stack((np.cos(theta), np.sin(theta)), axis=1)
    step = d[1:][seg] - d[:-1][seg]
    r = cross(d[:-1][seg], step)/cross(e, step)
    sign = np.sign(turn[seg]).astype(int)
    ray = np.mod(k, M)
    h2 = np.bincount(ray, sign*r**2, M)
    winding = np.bincount(ray, sign, M)
    multiplicity = np.bincount(ray, minlength=M)
    np.testing.assert_array_equal(winding, np.ones(M, int))
    assert np.all(multiplicity % 2 == 1)
    np.testing.assert_allclose(np.pi*np.mean(h2), area, rtol=1e-2)
    hx = np.sqrt(h2)
    hx -= hx.mean()
    T = 2*np.pi*R
    Sx = T * np.abs(np.fft.rfft(hx) / M)**2
    qx = 2*np.pi*np.fft.rfftfreq(M, d=T/M)
    np.testing.assert_allclose(np.mean(hx**2),
        (2*Sx[1:].sum() - (Sx[-1] if M % 2 == 0 else 0))/T)
    return dict(line=line, P=P, R=R, centre=centre, x=np.arange(M)*T/M + T/(2*M), hx=hx, s=s, hs=hs,
                qx=qx, Sx=Sx, qs=qs, Ss=Ss, multi=float(np.mean(multiplicity > 1)))


def bubbles(z):
    """Dominant contour (largest enclosed area) of every configuration."""
    best = {}
    for i, (k, A) in enumerate(zip(z["sample"], z["area"])):
        if A > 0 and A > best.get(k, (None, 0))[1]:
            best[k] = (i, A)
    return {k: i for k, (i, _) in best.items()}


def analyse(z, main, a):
    groups, spectra, report = {}, {}, {}
    for L in SIZES:
        ks = np.flatnonzero((z["L"]==L) & (z["aniso"]==a))
        assert len(np.unique(z["seed"][ks])) == len(ks)
        rows = [profiles(decode(z, main[k])) for k in ks]
        for r, k in zip(rows, ks):
            r["seed"] = int(z["seed"][k])
        groups[L] = rows
        Tmin = min(min(r["P"], 2*np.pi*r["R"]) for r in rows)
        spectra[L] = common_bins(rows, Tmin)
        P = np.array([r["P"] for r in rows])
        draws = RNG.integers(len(rows), size=(2000, len(rows)))
        fits = {}
        for kind, (q, values) in spectra[L].items():
            fits[kind] = {f"{lo:g}_{hi:g}": spectral_fit(q, values, lo, hi, draws)
                          for lo, hi in ((0, .03), (0, .1), (.015, .06), (.02, .08))}
        report[L] = dict(n=len(rows), meanP=float(P.mean()), stdP=float(P.std(ddof=1)),
                         meanP_ci95=np.quantile(P[draws].mean(axis=1), [.025, .975]).tolist(),
                         meanR=float(np.mean([r["R"] for r in rows])),
                         multivalued_fraction=float(np.mean([r["multi"] for r in rows])),
                         fits=fits)
        print(a, L, report[L]["n"], report[L]["meanP"], fits["x"]["0_0.03"], fits["s"]["0_0.03"], flush=True)
    Ls = np.array(SIZES)
    Pmean = np.array([report[L]["meanP"] for L in SIZES])
    bootstrap = []
    for _ in range(2000):
        y = [RNG.choice(np.array([r["P"] for r in groups[L]]), len(groups[L]), replace=True).mean()
             for L in SIZES]
        bootstrap.append(np.polyfit(np.log(Ls), np.log(y), 1)[0])
    report["perimeter_fit"] = dict(alpha=float(np.polyfit(np.log(Ls), np.log(Pmean), 1)[0]),
                                   ci95=np.quantile(bootstrap, [.025, .975]).tolist(),
                                   definition="Unweighted fit of log ensemble mean P against log L")
    return groups, spectra, report


def profile_figure(example, name):
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), layout="constrained")
    line = example["line"]
    axes[0].plot(line[:, 0], line[:, 1], color=COLORS[2], lw=1)
    axes[0].plot(*example["centre"], "+", color=".4", ms=6)
    axes[1].plot(example["x"], example["hx"], color=COLORS[0], lw=1)
    axes[2].plot(example["s"], example["hs"], color=COLORS[3], lw=1)
    axes[0].set(xlabel=r"$x$", ylabel=r"$y$", title="(a) Contour", aspect="equal")
    for spine in axes[0].spines.values():
        spine.set_visible(True)
    axes[1].set(xlabel=r"$R\,\theta$", ylabel=r"$h(\theta)-\langle h(\theta)\rangle$",
                title="(b) Area-preserving profile", xlim=(0, 2*np.pi*example["R"]))
    axes[2].set(xlabel=r"$s$", ylabel=r"$h(s)-\langle h(s)\rangle$",
                title="(c) Arc-length profile", xlim=(0, example["P"]))
    limits = [ax.get_ylim() for ax in axes[1:]]
    for ax in axes[1:]:
        ax.set_ylim(min(v[0] for v in limits), max(v[1] for v in limits))
        ax.axhline(0, c=".7", lw=.5, zorder=-1)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=200)
    plt.close(fig)


def scaling_figure(report, spectra, name, second=None):
    """second: optional reference exponent drawn dotted on top of the zeta = 1.25 plateau."""
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), layout="constrained")
    ax = axes[0]
    Ls = np.array(SIZES)
    Pmean = np.array([report[L]["meanP"] for L in SIZES])
    ci = np.array([report[L]["meanP_ci95"] for L in SIZES])
    ax.errorbar(Ls, Pmean, yerr=np.maximum(0, np.stack((Pmean-ci[:, 0], ci[:, 1]-Pmean))),
                fmt="o", c="#275c83", ms=7, capsize=3)
    xx = np.geomspace(880, 10000, 100)
    coefficient = .8*min(Pmean/Ls)
    ax.plot(xx, coefficient*xx, "--", c="black", lw=1.1)
    ax.text(5600, coefficient*5600*.82, r"$L$", ha="center", va="top")
    coefficient = 1.2*max(Pmean/Ls**1.25)
    ax.plot(xx, coefficient*xx**1.25, ":", c="black", lw=1.6)
    ax.text(5600, coefficient*5600**1.25*1.15, r"$L^{1.25}$", ha="center")
    ax.set(xscale="log", yscale="log", xlabel=r"$L$", ylabel=r"$\langle P\rangle$", title="(a) Mean perimeter")
    ax.set_xticks(Ls, labels=[str(L) for L in Ls])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())

    ax = axes[1]
    for L, color, marker in zip(SIZES[2:], COLORS[2:], MARKERS[2:]):
        for kind in ("x", "s"):
            q, values = spectra[L][kind]
            ax.loglog(q, q**3.5*values.mean(axis=0), ls="none", marker=marker, color=color,
                      markerfacecolor=color if kind == "x" else "none", ms=5, markeredgewidth=1)
    q, values = spectra[SIZES[-1]]["s"]
    anchor = np.exp(np.interp(np.log(.03), np.log(q), np.log(values.mean(axis=0))))
    plateau = anchor*.03**3.5
    xx = np.geomspace(.001, .15, 100)
    ax.plot(xx, np.full_like(xx, plateau), "--", color="black", lw=1.3)
    ax.text(.0012, plateau*1.12, r"$\zeta=1.25$", fontsize=10)
    if second is not None:
        ax.plot(xx, plateau*(xx/.03)**(2.5-2*second), ":", c="black", lw=1.6)
        ax.text(.0012, plateau*(.0012/.03)**(2.5-2*second)*.85, fr"$\zeta={second:g}$", fontsize=10, va="top")
    handles = [Line2D([], [], c=c, marker=m, ls="none", label=fr"$L={L}$")
               for L, c, m in zip(SIZES[2:], COLORS[2:], MARKERS[2:])]
    handles.extend([Line2D([], [], c=".25", marker="o", ls="none", label=r"$h(\theta)$"),
                    Line2D([], [], c=".25", marker="o", mfc="none", ls="none", label=r"$h(s)$")])
    ax.legend(handles=handles, frameon=False, fontsize=8, ncol=2, loc="upper left")
    ax.set(xlabel=r"$q$", ylabel=r"$q^{3.5}S(q)$", title="(b) Compensated height spectra")
    ax.set_xticks([.001, .01, .1], labels=["0.001", "0.01", "0.1"])
    ax.xaxis.set_minor_formatter(NullFormatter())

    # As the edge figure: plateau for qKPZ (zeta = 0.63), dotted zeta = 1.25.
    ax = axes[2]
    for L, color, marker in zip(SIZES[2:], COLORS[2:], MARKERS[2:]):
        for kind in ("x", "s"):
            q, values = spectra[L][kind]
            ax.loglog(q, q**2.26*values.mean(axis=0), ls="none", marker=marker, color=color,
                      markerfacecolor=color if kind == "x" else "none", ms=5, markeredgewidth=1)
    q, values = spectra[SIZES[-1]]["s"]
    anchor = np.exp(np.interp(np.log(.02), np.log(q), np.log(values.mean(axis=0))))
    plateau = anchor*.02**2.26
    xx = np.geomspace(.001, .15, 100)
    ax.plot(xx, np.full_like(xx, plateau), "--", color="black", lw=1.3)
    ax.text(.0012, plateau*.85, r"$\zeta=0.63$", fontsize=10, va="top")
    ax.plot(xx, plateau*(xx/.02)**-1.24, ":", c="black", lw=1.6)
    ax.text(.0012, plateau*(.0012/.02)**-1.24*.85, r"$\zeta=1.25$", fontsize=10, va="top")
    ax.legend(handles=handles, frameon=False, fontsize=8, ncol=2, loc="upper right")
    ax.set(xlabel=r"$q$", ylabel=r"$q^{2.26}S(q)$", title="(c) Compensated height spectra")
    ax.set_xticks([.001, .01, .1], labels=["0.001", "0.01", "0.1"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=200)
    plt.close(fig)


def main():
    plt.rcParams.update({"font.size":10, "axes.spines.top":False, "axes.spines.right":False,
                         "font.family":"DejaVu Sans", "pdf.fonttype":42})
    z = np.load(DATA)
    main_contour = bubbles(z)
    report = {}
    for a, names, second in ((0, ("fig7_bubble_profiles", "fig8_bubble_scaling"), None),
                             (1, ("fig9_bubble_anisotropic_profiles", "fig10_bubble_anisotropic_scaling"), .63)):
        groups, spectra, rep = analyse(z, main_contour, a)
        # Prespecified representative: seed zero, L=1024, no aesthetic selection.
        example = next(r for r in groups[1024] if r["seed"] == 0)
        profile_figure(example, names[0])
        scaling_figure(rep, spectra, names[1], second)
        rep["example"] = dict(L=1024, seed=0, P=example["P"], R=example["R"],
                              multivalued_fraction=example["multi"])
        report[f"a={a}"] = rep
        print(a, rep["perimeter_fit"], rep["example"], flush=True)
    report["source"] = str(DATA.relative_to(ROOT))
    report["geometry"] = "bubble, R0 = L/6; h(theta) from rays about the area centroid, period 2 pi R, R = sqrt(A/pi)"
    (OUT/"bubble_statistics.json").write_text(json.dumps(report, indent=2)+"\n")


if __name__ == "__main__":
    main()
