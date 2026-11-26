"""Figures for the isotropic flat-wall note, using saved contours only."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DATA = ROOT / "runs/flat_D0.04_b0_A2/walls.npz"
SIZES = (256, 512, 1024, 2048)
COLORS = ("#37658b", "#23947d", "#d4932b", "#ac4263")
MARKERS = ("o", "s", "^", "D")
RNG = np.random.default_rng(20261120)


def decode(z, i):
    first, last = z["offsets"][i:i+2]
    steps = z["steps"][first:last].astype(np.int64)
    start = z["start"][i].astype(np.int64)
    return np.vstack((start, start + np.cumsum(steps, axis=0))) / float(z["scale"])


def profiles(curve, L):
    """The winding contour has positive phase below it and winds towards +x."""
    line = np.vstack((curve, curve[:1] + (L, 0)))
    increments = np.linalg.norm(np.diff(line, axis=0), axis=1)
    keep = np.r_[True, increments > 0]
    line = line[keep]
    length = np.r_[0., np.cumsum(np.linalg.norm(np.diff(line, axis=0), axis=1))]
    P = length[-1]
    N = int(np.ceil(P))  # ds <= one lattice spacing
    s = np.arange(N) * P / N
    hs = np.interp(s, length, line[:, 1])
    hs -= hs.mean()
    Ss = P * np.abs(np.fft.rfft(hs) / N)**2
    qs = 2*np.pi*np.fft.rfftfreq(N, d=P/N)
    np.testing.assert_allclose(np.mean(hs**2),
        (2*Ss[1:].sum() - (Ss[-1] if N % 2 == 0 else 0))/P)

    # Periodically tile segments before intersecting x = j+1/2.
    x = np.arange(L) + .5
    hx, winding, multiplicity = np.zeros(L), np.zeros(L, int), np.zeros(L, int)
    for p, q in zip(line[:-1], line[1:]):
        dx = q[0] - p[0]
        if dx == 0:
            continue
        low, high = sorted((p[0], q[0]))
        k0 = int(np.floor((low-.5)/L))
        k1 = int(np.floor((high-.5)/L))
        for k in range(k0, k1+1):
            xx = x + k*L
            take = (xx >= low) & (xx < high)
            y = p[1] + (xx[take]-p[0]) * (q[1]-p[1])/dx
            sign = 1 if dx > 0 else -1
            hx[take] += sign*y
            winding[take] += sign
            multiplicity[take] += 1
    np.testing.assert_array_equal(winding, np.ones(L, int))
    assert np.all(multiplicity % 2 == 1)
    hx -= hx.mean()
    Sx = np.abs(np.fft.rfft(hx))**2/L
    qx = 2*np.pi*np.fft.rfftfreq(L)
    np.testing.assert_allclose(np.mean(hx**2), (2*Sx[1:-1].sum()+Sx[-1])/L)
    return dict(line=line, P=P, x=x, hx=hx, s=s, hs=hs, qx=qx, Sx=Sx, qs=qs, Ss=Ss,
                multi=float(np.mean(multiplicity>1)))


def common_bins(rows, L, upper=.22):
    """Equal weight per wall and at least one mode per wall and representation."""
    lower = 2*2*np.pi/L
    edges = np.geomspace(lower*(1-1e-12), upper, 70)
    bins, start = [], edges[0]
    for stop in edges[1:]:
        complete = all(np.any((r["q"+k] >= start) & (r["q"+k] < stop))
                       for r in rows for k in ("s", "x"))
        if complete:
            bins.append((start, stop))
            start = stop
    output = {}
    for kind in ("x", "s"):
        values, positions = [], []
        for low, high in bins:
            wall_means, wall_logq = [], []
            for r in rows:
                take = (r["q"+kind] >= low) & (r["q"+kind] < high)
                wall_means.append(np.mean(r["S"+kind][take]))
                wall_logq.append(np.mean(np.log(r["q"+kind][take])))
            values.append(wall_means)
            positions.append(np.exp(np.mean(wall_logq)))
        output[kind] = (np.array(positions), np.array(values).T)
    return output


def spectral_fit(q, values, low, high, draws):
    take = (q >= low) & (q <= high)
    if take.sum() < 3:
        return None
    logq = np.log(q[take])
    weights = (logq-logq.mean()) / np.sum((logq-logq.mean())**2)
    zeta = lambda mean: (-np.log(mean[..., take]) @ weights - 1)/2
    point = float(zeta(values.mean(axis=0)))
    boot = zeta(values[draws].mean(axis=1))
    return dict(zeta=point, ci95=np.quantile(boot, [.025, .975]).tolist(),
                bins=int(take.sum()), q_actual=q[take][[0,-1]].tolist())


MODES = 45/2048   # modes 2 <= n <= MODES*L shown in the spectral panels, q_x = 2 pi n / L <= 0.14


def mode_spectra(rows, L):
    """S(n) = <|h_n|^2> at fixed mode index n, equal weight per contour (|h_n|^2 = S_i(q_n)/T)."""
    n = np.arange(2, int(MODES*L)+1)
    return n, {kind: np.mean([r["S"+kind][n]*r["q"+kind][1]/(2*np.pi) for r in rows], axis=0)
               for kind in ("x", "s")}


def mode_panel(ax, modes, sizes, references, kinds=("x", "s")):
    """Raw spectra versus n, filled h(x) and open h(s); references: (zeta, linestyle, factor)
    lines n^-(1+2 zeta), factor > 1 times above (< 1: below) every point shown, labelled at their
    upper left (above) or lower right (below) end."""
    for L, color, marker in zip(sizes, COLORS[-len(sizes):], MARKERS[-len(sizes):]):
        n, S = modes[L]
        for kind in kinds:
            ax.loglog(n, S[kind], ls="none", marker=marker, color=color,
                      markerfacecolor=color if kind == "x" else "none", ms=5, markeredgewidth=1)
    n = modes[sizes[-1]][0]
    nn = np.geomspace(n[0], n[-1], 100)
    for zeta, ls, factor in references:
        ratio = np.concatenate([S[kind]*m**(1+2*zeta) for m, S in (modes[L] for L in sizes) for kind in kinds])
        line = factor*(ratio.max() if factor > 1 else ratio.min())*nn**-(1+2*zeta)
        ax.plot(nn, line, ls, c="black", lw=1.3 if ls == "--" else 1.6)
        low, high = ax.get_ylim()
        if factor > 1:
            ax.text(nn[0]*1.05, line[0]*1.3, fr"$\zeta={zeta:g}$", fontsize=10, ha="left", va="bottom")
            ax.set_ylim(low, max(high, line[0]*25))
        else:
            ax.text(nn[-1], line[-1]*.7, fr"$\zeta={zeta:g}$", fontsize=10, ha="right", va="top")
            ax.set_ylim(min(low, line[-1]*.04), high)
    ax.set(xlabel=r"Fourier mode $n$", ylabel=r"$S(n)$")
    ticks = [t for t in (2, 4, 8, 16, 32, 64, 128) if t <= n[-1]]
    ax.set_xticks(ticks, labels=[str(t) for t in ticks])
    ax.xaxis.set_minor_formatter(NullFormatter())


def main():
    plt.rcParams.update({"font.size":10, "axes.spines.top":False, "axes.spines.right":False,
                         "font.family":"DejaVu Sans", "pdf.fonttype":42})
    z = np.load(DATA)
    groups, spectra, report = {}, {}, {}
    for L in SIZES:
        ids = np.flatnonzero((z["L"]==L) & (z["aniso"]==0))
        assert len(np.unique(z["seed"][ids])) == len(ids)
        rows = [profiles(decode(z, i), L) for i in ids]
        groups[L] = rows
        spectra[L] = common_bins(rows, L)
        P = np.array([r["P"] for r in rows])
        draws = RNG.integers(len(rows), size=(2000, len(rows)))
        means = P[draws].mean(axis=1)
        fits = {}
        for kind, (q, values) in spectra[L].items():
            fits[kind] = {f"{low:g}_{hi:g}": spectral_fit(q, values, low, hi, draws)
                          for low, hi in ((0, .03), (0, .1), (.015, .06), (.02,.08))}
        report[L] = dict(n=len(rows), meanP=float(P.mean()), stdP=float(P.std(ddof=1)),
                         meanP_ci95=np.quantile(means, [.025, .975]).tolist(),
                         multivalued_fraction=float(np.mean([r["multi"] for r in rows])),
                         fits=fits)
        print(L, report[L], flush=True)

    # Prespecified representative: seed zero, L=1024, no aesthetic selection.
    ids = np.flatnonzero((z["L"]==1024)&(z["aniso"]==0))
    index = int(np.flatnonzero(z["seed"][ids]==0)[0])
    r = groups[1024][index]
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), layout="constrained",
                             gridspec_kw={"width_ratios": (1, 1, 1)})
    line = r["line"]
    axes[0].plot(line[:,0], line[:,1], color=COLORS[2], lw=1)
    axes[1].plot(r["x"], r["hx"], color=COLORS[0], lw=1)
    axes[2].plot(r["s"], r["hs"], color=COLORS[3], lw=1)
    margin = .07*np.ptp(line[:, 1])
    axes[0].set(xlabel=r"$x$", ylabel=r"$y$", title="(a) Contour (zoom)",
                xlim=(0, 1024),
                ylim=(max(0, line[:, 1].min()-margin),
                      min(2048, line[:, 1].max()+margin)))
    axes[0].set_xticks([0, 512, 1024])
    for spine in axes[0].spines.values():
        spine.set_visible(True)
    axes[1].set(xlabel=r"$x$", ylabel=r"$h(x)-\langle h(x)\rangle$",
                title="(b) Area-preserving profile", xlim=(0,1024))
    axes[2].set(xlabel=r"$s$", ylabel=r"$h(s)-\langle h(s)\rangle$",
                title="(c) Arc-length profile", xlim=(0,r["P"]))
    limits = [ax.get_ylim() for ax in axes[1:]]
    ylim = min(v[0] for v in limits), max(v[1] for v in limits)
    for ax in axes[1:]:
        ax.set_ylim(ylim)
        ax.axhline(0, c=".7", lw=.5, zorder=-1)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig1_profiles.{ext}", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(1,2,figsize=(10,4.3),layout="constrained")
    ax = axes[0]
    Ls = np.array(SIZES)
    Pmean = np.array([report[L]["meanP"] for L in SIZES])
    ci = np.array([report[L]["meanP_ci95"] for L in SIZES])
    ax.errorbar(Ls,Pmean,yerr=np.maximum(0,np.stack((Pmean-ci[:,0],ci[:,1]-Pmean))),
                fmt="o",c="#275c83",ms=7,capsize=3)
    xx = np.geomspace(220,2500,100)
    ax.plot(xx,.8*xx,"--",c="black",lw=1.1)
    ax.text(1400,.8*1400*.82,r"$L_x$",ha="center",va="top")
    coefficient = 1.2*max(Pmean/Ls**1.25)
    ax.plot(xx,coefficient*xx**1.25,":",c="black",lw=1.6)
    ax.text(1400,coefficient*1400**1.25*1.15,r"$L_x^{1.25}$",ha="center")
    slope = np.polyfit(np.log(Ls),np.log(Pmean),1)[0]
    bootstrap = []
    for _ in range(2000):
        y=[]
        for L in SIZES:
            p=np.array([r["P"] for r in groups[L]])
            y.append(RNG.choice(p,len(p),replace=True).mean())
        bootstrap.append(np.polyfit(np.log(Ls),np.log(y),1)[0])
    report["perimeter_fit"] = dict(alpha=float(slope),ci95=np.quantile(bootstrap,[.025,.975]).tolist(),
                                  definition="Unweighted fit of log ensemble mean P against log Lx")
    ax.set(xscale="log",yscale="log",xlabel=r"$L_x$",ylabel=r"$\langle P\rangle$",title="(a) Mean perimeter")
    ax.set_xticks(Ls,labels=[str(L) for L in Ls]);ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())

    ax = axes[1]
    mode_panel(ax, {L: mode_spectra(groups[L], L) for L in SIZES[2:]}, SIZES[2:], [(1.25, "--", 2)])
    handles=[Line2D([],[],c=c,marker=m,ls="none",label=fr"$L_x={L}$")
             for L,c,m in zip(SIZES[2:],COLORS[2:],MARKERS[2:])]
    handles.extend([Line2D([],[],c=".25",marker="o",ls="none",label=r"$h(x)$"),
                    Line2D([],[],c=".25",marker="o",mfc="none",ls="none",label=r"$h(s)$")])
    ax.legend(handles=handles,frameon=False,fontsize=8,ncol=2,loc="lower left")
    ax.set(title="(b) Height spectra")
    for ext in ("pdf","png"):
        fig.savefig(OUT/f"fig2_scaling.{ext}",dpi=200)
    plt.close(fig)
    report["example"] = dict(Lx=1024,seed=0,P=r["P"],multivalued_fraction=r["multi"])
    report["source"] = str(DATA.relative_to(ROOT))
    (OUT/"statistics.json").write_text(json.dumps(report,indent=2)+"\n")
    print("PERIMETER",report["perimeter_fit"])

if __name__=="__main__":
    main()
