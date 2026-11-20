"""Bubbles without a centre: width of arcs of length ell about their mean line.

Each contour is resampled at uniform arc length (ds <= 1) and cut into arcs of length ell (start
points ell/4 apart). The mean orientation of an arc, the average of its unit tangent, is the
direction of its chord. The reference is the straight line along that direction through the mean
point of the arc; the deviation is the distance from it, w^2(ell) is its mean square, and X is the
chord length. With METHOD = "cubic" the reference is instead a smooth parametric cubic, x(s) and
y(s) fitted separately, the deviation is taken along its normal and X is its length.

Isotropic bubbles (fig11): all arcs; flat isotropic walls (L = 2048) for comparison.
Anisotropic bubbles (fig12, a = 1, b along x): arcs are split by the direction of the normal of
their chord, folded to [0, 90] degrees from b: edge-like arcs (normal within 30 degrees of b) and
screw-like arcs (normal within 30 degrees of the perpendicular), compared with flat edge
(b = 90 degrees) and screw (b = 0) walls at a = 1, L = 2048 (all arcs).
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from make_figures import decode, COLORS, MARKERS
import make_bubble_figures as mb

ROOT = mb.ROOT
OUT = mb.OUT
FLAT = ROOT / "runs/flat_D0.04_b0_A2/walls.npz"
EDGE = ROOT / "runs/flat_D0.04_b90_A2/walls.npz"
DEGREE = 3
METHOD = "line"
ELLS = np.unique(np.rint(np.geomspace(8, 4096, 28)))
SECTORS = dict(edge=(0., 30.), screw=(60., 90.01))
RNG = np.random.default_rng(20261124)


def resample(line):
    """Uniform arc-length points of one period of the contour, ds = P/N <= 1."""
    length = np.r_[0., np.cumsum(np.linalg.norm(np.diff(line, axis=0), axis=1))]
    P = length[-1]
    N = int(np.ceil(P))
    s = np.arange(N)*P/N
    return np.stack((np.interp(s, length, line[:, 0]), np.interp(s, length, line[:, 1])), 1), P/N


def arcs(points, ds, shift, ell, degree=DEGREE, every=1):
    """Arcs of length ell starting every n/4 points (every: keep one arc in that many).
    Returns the arcs (k, 2, n), the reference curves, the normal deviations (k, n), the lengths X
    and the direction of the chord normal folded to [0, 90] degrees from the x axis.
    METHOD "line": the reference is the straight line through the mean point of the arc along its
    mean orientation (the mean unit tangent, i.e. the chord), X is the chord length and the
    deviation is the distance from that line. METHOD "cubic": the smooth parametric cubic."""
    N, n = len(points), int(round(ell/ds))
    extended = np.vstack((points, points[:n] + shift))
    W = sliding_window_view(extended, n, axis=0)[:N:max(1, n//4)*every]
    chord = W[:, :, -1] - W[:, :, 0]
    angle = np.mod(np.degrees(np.arctan2(chord[:, 1], chord[:, 0])) + 90, 180)
    phi = np.minimum(angle, 180 - angle)
    if METHOD == "line":
        X = np.linalg.norm(chord, axis=1)
        tangent = chord/X[:, None]
        normal = np.stack((-tangent[:, 1], tangent[:, 0]), 1)
        centre = W.mean(axis=2)
        relative = W - centre[:, :, None]
        deviation = np.einsum("kin,ki->kn", relative, normal)
        along = np.einsum("kin,ki->kn", relative, tangent)
        ends = np.stack((along.min(axis=1), along.max(axis=1)), 1)
        line = centre[:, :, None] + tangent[:, :, None]*ends[:, None, :]
        return W, line, deviation, X, phi
    t = np.linspace(-1, 1, n)
    V = np.vander(t, degree+1)
    coefficients = np.linalg.lstsq(V, np.moveaxis(W, -1, 0).reshape(n, -1), rcond=None)[0]
    fit = np.moveaxis((V @ coefficients).reshape(n, *W.shape[:2]), 0, -1)
    derivative = np.vander(t, degree)*np.arange(degree, 0, -1)
    tangent = np.moveaxis((derivative @ coefficients[:-1]).reshape(n, *W.shape[:2]), 0, -1)
    speed = np.linalg.norm(tangent, axis=1)
    normal = np.stack((-tangent[:, 1], tangent[:, 0]), 1)/speed[:, None]
    deviation = np.sum((W - fit)*normal, axis=1)
    X = np.trapezoid(speed, t, axis=1)
    return W, fit, deviation, X, phi


def widths(points, ds, shift, P, sector=None):
    """Per contour: mean square deviation and mean fitted length at every ell <= P/4,
    over the arcs whose folded normal direction lies in sector (all arcs if None)."""
    w2, X = np.full(len(ELLS), np.nan), np.full(len(ELLS), np.nan)
    for j, ell in enumerate(ELLS):
        if ell <= P/4:
            _, _, d, x, phi = arcs(points, ds, shift, ell)
            take = np.ones(len(x), bool) if sector is None else (phi >= sector[0]) & (phi < sector[1])
            if take.any():
                w2[j], X[j] = np.mean(d[take]**2), np.mean(x[take])
    return w2, X


def ensemble(w2, X, draws, complete=1.):
    """Mean width, mean fitted length, local exponent d ln w / d ln X and its 95% bootstrap interval,
    at the ell where at least the fraction complete of the contours contribute."""
    keep = np.mean(np.isfinite(w2), axis=0) >= complete
    mean = lambda a: np.nanmean(a[..., keep], axis=-2)
    slope = lambda a, b: np.gradient(.5*np.log(a), np.log(b), axis=-1)
    boot = np.array([slope(mean(w2[i]), mean(X[i])) for i in draws])
    return dict(ell=ELLS[keep], w=np.sqrt(mean(w2)), X=mean(X), zeta=slope(mean(w2), mean(X)),
                zeta_ell=np.gradient(.5*np.log(mean(w2)), np.log(ELLS[keep])),
                ci95=np.quantile(boot, [.025, .975], axis=0), n=len(w2))


def flat_walls(path, a, L=2048):
    z = np.load(path)
    rows = []
    for i in np.flatnonzero((z["L"]==L) & (z["aniso"]==a)):
        curve = decode(z, i)
        line = np.vstack((curve, curve[:1] + (L, 0)))
        line = line[np.r_[True, np.linalg.norm(np.diff(line, axis=0), axis=1) > 0]]
        points, ds = resample(line)
        rows.append(widths(points, ds, np.array([float(L), 0.]), len(points)*ds))
    w2, X = map(np.array, zip(*rows))
    return ensemble(w2, X, RNG.integers(len(rows), size=(2000, len(rows))))


def bubbles(z, main_contour, a, sectors):
    """sectors: {name: sector or None}. Returns {name: {L: ensemble}} and the L = 1024, seed 0 example."""
    result, example = {name: {} for name in sectors}, None
    for L in mb.SIZES:
        ks = np.flatnonzero((z["L"]==L) & (z["aniso"]==a))
        rows = {name: [] for name in sectors}
        for k in ks:
            line = mb.profiles(decode(z, main_contour[k]))["line"]
            points, ds = resample(line[:-1])
            for name, sector in sectors.items():
                rows[name].append(widths(points, ds, np.zeros(2), len(points)*ds, sector))
            if L == 1024 and z["seed"][k] == 0:
                example = (line, points, ds)
        draws = RNG.integers(len(ks), size=(2000, len(ks)))
        for name, sector in sectors.items():
            w2, X = map(np.array, zip(*rows[name]))
            result[name][L] = ensemble(w2, X, draws, 1. if sector is None else .8)
        print(a, L, len(ks), flush=True)
    return result, example


def example_panel(ax, example, colors, count=4):
    """count consecutive arcs of length 256 (all if None); colors(j, phi) gives the color of an arc."""
    line, points, ds = example
    ax.plot(line[:, 0], line[:, 1], color=".75", lw=.8)
    W, fit, _, _, phi = arcs(points, ds, np.zeros(2), 256, every=4)
    count = len(W) if count is None else count
    for j in range(count):
        ax.plot(W[j, 0], W[j, 1], color=colors(j, phi[j]), lw=1)
        ax.plot(fit[j, 0], fit[j, 1], color="black", lw=1.1)
    lo, hi = W[:count].min(axis=(0, 2)), W[:count].max(axis=(0, 2))
    pad = .08*max(hi-lo)
    ax.set(xlim=(lo[0]-pad, hi[0]+pad), ylim=(lo[1]-pad, hi[1]+pad), aspect="equal",
           xlabel=r"$x$", ylabel=r"$y$", title=r"(a) Arcs of length $\ell=256$")
    for spine in ax.spines.values():
        spine.set_visible(True)


def isotropic_figure(result, flat, example):
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), layout="constrained")
    example_panel(axes[0], example, lambda j, phi: (COLORS[2], COLORS[0])[j % 2])

    ax = axes[1]
    for L, color, marker in zip(mb.SIZES, COLORS, MARKERS):
        r = result[L]
        ax.loglog(r["X"], r["w"], ls="none", marker=marker, color=color, ms=5)
    ax.loglog(flat["X"], flat["w"], color=".45", lw=1.1)
    xx = np.geomspace(60, 3000, 100)
    anchor = np.exp(np.interp(np.log(200), np.log(result[8192]["X"]), np.log(result[8192]["w"])))
    ax.plot(xx, 1.6*anchor*(xx/200)**1.25, "--", c="black", lw=1.3)
    ax.text(700, 1.6*anchor*(700/200)**1.25*1.25, r"$\zeta=1.25$", ha="right", fontsize=10)
    ax.plot(xx, .6*anchor*(xx/200)**.63, ":", c="black", lw=1.6)
    ax.text(2500, .6*anchor*(2500/200)**.63*.55, r"$\zeta=0.63$", ha="right", va="top", fontsize=10)
    handles = [Line2D([], [], c=c, marker=m, ls="none", label=fr"$L={L}$")
               for L, c, m in zip(mb.SIZES, COLORS, MARKERS)]
    handles.append(Line2D([], [], c=".45", lw=1.1, label="flat, $L=2048$"))
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="upper left")
    ax.set(xlabel=r"$X$", ylabel=r"$w$", title="(b) Width about the mean line")

    ax = axes[2]
    for L, color, marker in zip(mb.SIZES, COLORS, MARKERS):
        r = result[L]
        ax.errorbar(r["X"], r["zeta"], yerr=np.abs(r["ci95"] - r["zeta"]), ls="none", marker=marker,
                    color=color, ms=5, capsize=2, lw=.8)
    ax.plot(flat["X"], flat["zeta"], color=".45", lw=1.1)
    r = result[8192]
    ax.plot(r["X"], r["zeta_ell"], ls="none", marker=MARKERS[3], color=COLORS[3], mfc="none", ms=5)
    ax.axhline(1.25, ls="--", c="black", lw=1.3)
    ax.axhline(.63, ls=":", c="black", lw=1.6)
    ax.text(6.5, 1.25*1.03, r"$\zeta=1.25$", fontsize=10)
    ax.text(6.5, .63*1.04, r"$\zeta=0.63$", fontsize=10)
    ax.legend(handles=[Line2D([], [], c=".25", marker="o", ls="none", label=r"$d\ln w/d\ln X$"),
                       Line2D([], [], c=COLORS[3], marker=MARKERS[3], mfc="none", ls="none",
                              label=r"$d\ln w/d\ln\ell$, $L=8192$")],
              frameon=False, fontsize=8, loc="upper right")
    ax.set(xscale="log", xlabel=r"$X$", ylabel=r"$d\ln w/d\ln X$", ylim=(.4, 2.6),
           title="(c) Local exponent")
    for ax in axes[1:]:
        ax.xaxis.set_minor_formatter(NullFormatter())
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig11_bubble_arc_width.{ext}", dpi=200)
    plt.close(fig)


def anisotropic_figure(result, flat, example):
    """(a) example arcs, (b) widths, (c) local exponent; filled edge-like, open screw-like arcs."""
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), layout="constrained")
    sector_color = lambda j, phi: ("#2b7bba" if phi < 30 else "#c0392b" if phi >= 60 else ".6")
    example_panel(axes[0], example, sector_color, count=None)
    axes[0].legend(handles=[Line2D([], [], c="#2b7bba", label="edge-like"),
                            Line2D([], [], c="#c0392b", label="screw-like"),
                            Line2D([], [], c=".6", label="other"),
                            Line2D([], [], c="black", label="mean line")],
                   frameon=False, fontsize=7, loc="center left", bbox_to_anchor=(1, .5))
    style = dict(edge=dict(), screw=dict(mfc="none"))
    flat_style = dict(edge="-", screw="--")

    ax = axes[1]
    for name in SECTORS:
        for L, color, marker in zip(mb.SIZES, COLORS, MARKERS):
            r = result[name][L]
            ax.loglog(r["X"], r["w"], ls="none", marker=marker, color=color, ms=5, **style[name])
        ax.loglog(flat[name]["X"], flat[name]["w"], ls=flat_style[name], color=".45", lw=1.1)
    xx = np.geomspace(60, 3000, 100)
    anchor = np.exp(np.interp(np.log(200), np.log(result["edge"][8192]["X"]),
                              np.log(result["edge"][8192]["w"])))
    ax.plot(xx, 1.6*anchor*(xx/200)**1.25, "--", c="black", lw=1.3)
    ax.text(2500, 1.6*anchor*(2500/200)**1.25*.6, r"$\zeta=1.25$", ha="right", va="top", fontsize=10)
    ax.plot(xx, .5*anchor*(xx/200)**.63, ":", c="black", lw=1.6)
    ax.text(2500, .5*anchor*(2500/200)**.63*.55, r"$\zeta=0.63$", ha="right", va="top", fontsize=10)
    handles = [Line2D([], [], c=c, marker=m, ls="none", label=fr"$L={L}$")
               for L, c, m in zip(mb.SIZES, COLORS, MARKERS)]
    handles.extend([Line2D([], [], c=".25", marker="o", ls="none", label="edge-like arcs"),
                    Line2D([], [], c=".25", marker="o", mfc="none", ls="none", label="screw-like arcs"),
                    Line2D([], [], c=".45", lw=1.1, label="flat edge, $L=2048$"),
                    Line2D([], [], c=".45", lw=1.1, ls="--", label="flat screw, $L=2048$")])
    ax.legend(handles=handles, frameon=False, fontsize=7, loc="upper left")
    ax.set(xlabel=r"$X$", ylabel=r"$w$", title="(b) Width about the mean line")

    ax = axes[2]
    for name in SECTORS:
        for L, color, marker in zip(mb.SIZES[2:], COLORS[2:], MARKERS[2:]):
            r = result[name][L]
            ok = np.diff(r["ci95"], axis=0)[0] < .5
            ax.errorbar(r["X"][ok], r["zeta"][ok], yerr=np.abs(r["ci95"] - r["zeta"])[:, ok], ls="none",
                        marker=marker, color=color, ms=5, capsize=2, lw=.8, **style[name])
        ax.plot(flat[name]["X"], flat[name]["zeta"], ls=flat_style[name], color=".45", lw=1.1)
    ax.axhline(1.25, ls="--", c="black", lw=1.3)
    ax.axhline(.63, ls=":", c="black", lw=1.6)
    ax.text(6.5, 1.25*1.03, r"$\zeta=1.25$", fontsize=10)
    ax.text(6.5, .63*1.04, r"$\zeta=0.63$", fontsize=10)
    ax.set(xscale="log", xlabel=r"$X$", ylabel=r"$d\ln w/d\ln X$", ylim=(.4, 2.6),
           title="(c) Local exponent")
    for ax in axes[1:]:
        ax.xaxis.set_minor_formatter(NullFormatter())
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig12_bubble_anisotropic_arc_width.{ext}", dpi=200)
    plt.close(fig)


def table(r):
    return {key: np.asarray(v).tolist() for key, v in r.items()}


def main():
    plt.rcParams.update({"font.size":10, "axes.spines.top":False, "axes.spines.right":False,
                         "font.family":"DejaVu Sans", "pdf.fonttype":42})
    z = np.load(mb.DATA)
    main_contour = mb.bubbles(z)
    report = dict(method=__doc__.strip(), source=[str(p.relative_to(ROOT)) for p in (mb.DATA, FLAT, EDGE)])

    result, example = bubbles(z, main_contour, 0, dict(all=None))
    flat = flat_walls(FLAT, 0)
    isotropic_figure(result["all"], flat, example)
    report["a=0"] = {str(L): table(r) for L, r in result["all"].items()}
    report["a=0"]["flat"] = table(flat)

    result, example = bubbles(z, main_contour, 1, SECTORS)
    flat = dict(edge=flat_walls(EDGE, 1), screw=flat_walls(FLAT, 1))
    anisotropic_figure(result, flat, example)
    report["a=1"] = {name: {str(L): table(r) for L, r in result[name].items()} for name in SECTORS}
    for name in SECTORS:
        report["a=1"][name]["flat"] = table(flat[name])
    (OUT/"bubble_arc_width.json").write_text(json.dumps(report, indent=1)+"\n")


if __name__ == "__main__":
    main()
