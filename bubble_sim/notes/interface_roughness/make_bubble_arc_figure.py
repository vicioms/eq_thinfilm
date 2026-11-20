"""Isotropic bubbles without a centre: width of arcs of length ell about a smooth cubic fit.

Each contour is resampled at uniform arc length (ds <= 1). Every arc of length ell (start points
ell/4 apart) is fitted by a smooth parametric cubic, x(s) and y(s) separately (least squares; a
spline with no interior knot). The deviation is the component of r(s) - r_fit(s) along the unit
normal of the fit, and w^2(ell) is its mean square. A cubic removes the mean, the slope and the
curvature of the arc, so w(X) ~ X^zeta for zeta < 4. X is the length of the fitted cubic.
The flat isotropic walls (L = 2048) are analysed in the same way for comparison.
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
DEGREE = 3
ELLS = np.unique(np.rint(np.geomspace(8, 4096, 28)))
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
    Returns the arcs (k, 2, n), the cubic fits and the normal deviations (k, n)."""
    N, n = len(points), int(round(ell/ds))
    extended = np.vstack((points, points[:n] + shift))
    W = sliding_window_view(extended, n, axis=0)[:N:max(1, n//4)*every]
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
    return W, fit, deviation, X


def widths(points, ds, shift, P):
    """Per contour: mean square deviation and mean fitted length at every ell <= P/4."""
    w2, X = np.full(len(ELLS), np.nan), np.full(len(ELLS), np.nan)
    for j, ell in enumerate(ELLS):
        if ell <= P/4:
            _, _, d, x = arcs(points, ds, shift, ell)
            w2[j], X[j] = np.mean(d**2), np.mean(x)
    return w2, X


def ensemble(w2, X, draws):
    """Mean width, mean fitted length, local exponent d ln w / d ln X and its 95% bootstrap interval."""
    keep = np.all(np.isfinite(w2), axis=0)
    slope = lambda a, b: np.gradient(.5*np.log(a[..., keep]), np.log(b[..., keep]), axis=-1)
    boot = np.array([slope(w2[i].mean(axis=0), X[i].mean(axis=0)) for i in draws])
    return dict(ell=ELLS[keep], w=np.sqrt(w2[:, keep].mean(axis=0)), X=X[:, keep].mean(axis=0),
                zeta=slope(w2.mean(axis=0), X.mean(axis=0)),
                zeta_ell=np.gradient(.5*np.log(w2[:, keep].mean(axis=0)), np.log(ELLS[keep])),
                ci95=np.quantile(boot, [.025, .975], axis=0))


def main():
    plt.rcParams.update({"font.size":10, "axes.spines.top":False, "axes.spines.right":False,
                         "font.family":"DejaVu Sans", "pdf.fonttype":42})
    z = np.load(mb.DATA)
    main_contour = mb.bubbles(z)
    result, example = {}, None
    for L in mb.SIZES:
        ks = np.flatnonzero((z["L"]==L) & (z["aniso"]==0))
        rows = []
        for k in ks:
            line = mb.profiles(decode(z, main_contour[k]))["line"]
            points, ds = resample(line[:-1])
            rows.append(widths(points, ds, np.zeros(2), len(points)*ds))
            if L == 1024 and z["seed"][k] == 0:
                example = (line, points, ds)
        w2, X = map(np.array, zip(*rows))
        result[L] = ensemble(w2, X, RNG.integers(len(ks), size=(2000, len(ks))))
        print(L, len(ks), flush=True)
    zf = np.load(FLAT)
    ids = np.flatnonzero((zf["L"]==2048) & (zf["aniso"]==0))
    rows = []
    for i in ids:
        curve = decode(zf, i)
        line = np.vstack((curve, curve[:1] + (2048, 0)))
        line = line[np.r_[True, np.linalg.norm(np.diff(line, axis=0), axis=1) > 0]]
        points, ds = resample(line)
        rows.append(widths(points, ds, np.array([2048., 0.]), len(points)*ds))
    w2, X = map(np.array, zip(*rows))
    result["flat"] = ensemble(w2, X, RNG.integers(len(ids), size=(2000, len(ids))))

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), layout="constrained")
    ax = axes[0]
    line, points, ds = example
    ax.plot(line[:, 0], line[:, 1], color=".75", lw=.8)
    W, fit, _, _ = arcs(points, ds, np.zeros(2), 256, every=4)
    for j, color in zip(range(4), (COLORS[2], COLORS[0])*2):
        ax.plot(W[j, 0], W[j, 1], color=color, lw=1)
        ax.plot(fit[j, 0], fit[j, 1], color="black", lw=1.1)
    lo, hi = W[:4].min(axis=(0, 2)), W[:4].max(axis=(0, 2))
    pad = .08*max(hi-lo)
    ax.set(xlim=(lo[0]-pad, hi[0]+pad), ylim=(lo[1]-pad, hi[1]+pad), aspect="equal",
           xlabel=r"$x$", ylabel=r"$y$", title=r"(a) Arcs of length $\ell=256$")
    for spine in ax.spines.values():
        spine.set_visible(True)

    ax = axes[1]
    for L, color, marker in zip(mb.SIZES, COLORS, MARKERS):
        r = result[L]
        ax.loglog(r["X"], r["w"], ls="none", marker=marker, color=color, ms=5)
    r = result["flat"]
    ax.loglog(r["X"], r["w"], color=".45", lw=1.1)
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
    ax.set(xlabel=r"$X$", ylabel=r"$w$", title="(b) Width about the cubic")

    ax = axes[2]
    for L, color, marker in zip(mb.SIZES, COLORS, MARKERS):
        r = result[L]
        ax.errorbar(r["X"], r["zeta"], yerr=np.abs(r["ci95"] - r["zeta"]), ls="none", marker=marker,
                    color=color, ms=5, capsize=2, lw=.8)
    r = result["flat"]
    ax.plot(r["X"], r["zeta"], color=".45", lw=1.1)
    r = result[8192]
    ax.plot(r["X"], r["zeta_ell"], ls="none", marker=MARKERS[3], color=COLORS[3], mfc="none", ms=5)
    ax.axhline(1.25, ls="--", c="black", lw=1.3)
    ax.axhline(.63, ls=":", c="black", lw=1.6)
    ax.text(12, 1.25*1.03, r"$\zeta=1.25$", fontsize=10)
    ax.text(12, .63*1.04, r"$\zeta=0.63$", fontsize=10)
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

    report = {str(k): {key: np.asarray(v).tolist() for key, v in r.items()} for k, r in result.items()}
    report["method"] = __doc__.strip()
    report["source"] = [str(mb.DATA.relative_to(ROOT)), str(FLAT.relative_to(ROOT))]
    (OUT/"bubble_arc_width.json").write_text(json.dumps(report, indent=1)+"\n")


if __name__ == "__main__":
    main()
