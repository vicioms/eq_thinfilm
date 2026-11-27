"""Anisotropic flat-edge figures from saved contours, with the isotropic estimator."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter
import numpy as np

from make_figures import decode, profiles, common_bins, spectral_fit, mode_spectra, mode_panel, SIZES, COLORS, MARKERS

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
DATA = ROOT / "runs/flat_D0.04_b90_A2/walls.npz"
RNG = np.random.default_rng(20261121)
SIZES = SIZES + (4096,)
COLORS = COLORS + ("#6b4f9e",)
MARKERS = MARKERS + ("v",)


def mode_fits(rows, L):
    """zeta from least-squares slopes of the raw S(n) = <|h_n|^2> in fixed mode ranges."""
    n, S = mode_spectra(rows, L)
    ranges = {"2_8": (2, 8), "first_quarter": (2, n[-1]//4), "rest": (n[-1]//4+1, n[-1])}
    return {name: dict(n=[int(a), int(b)], **{kind: float((-np.polyfit(np.log(n[(n>=a)&(n<=b)]),
                np.log(S[kind][(n>=a)&(n<=b)]), 1)[0]-1)/2) for kind in ("x", "s")})
            for name, (a, b) in ranges.items()}


def main():
    plt.rcParams.update({"font.size":10, "axes.spines.top":False, "axes.spines.right":False,
                         "font.family":"DejaVu Sans", "pdf.fonttype":42})
    z = np.load(DATA)
    groups, spectra, report = {}, {}, {}
    for L in SIZES:
        ids = np.flatnonzero((z["L"]==L)&(z["aniso"]==1))
        assert len(np.unique(z["seed"][ids])) == len(ids)
        rows = [profiles(decode(z,i),L) for i in ids]
        groups[L] = rows
        spectra[L] = common_bins(rows,L)
        P = np.array([r["P"] for r in rows])
        draws = RNG.integers(len(rows),size=(2000,len(rows)))
        fits = {}
        for kind,(q,values) in spectra[L].items():
            fits[kind] = {f"{lo:g}_{hi:g}":spectral_fit(q,values,lo,hi,draws)
                          for lo,hi in ((0,.03),(0,.1),(.015,.06),(.02,.08))}
        report[L] = dict(n=len(rows), meanP=float(P.mean()), stdP=float(P.std(ddof=1)),
                         meanP_ci95=np.quantile(P[draws].mean(axis=1),[.025,.975]).tolist(),
                         multivalued_fraction=float(np.mean([r["multi"] for r in rows])),
                         fits=fits, mode_fits=mode_fits(rows,L) if L>=1024 else None)

    ids = np.flatnonzero((z["L"]==1024)&(z["aniso"]==1))
    index = int(np.flatnonzero(z["seed"][ids]==0)[0])
    example = groups[1024][index]
    line = example["line"]
    fig,axes = plt.subplots(1,3,figsize=(11.5,3.8),layout="constrained")
    axes[0].plot(line[:,0],line[:,1],color=COLORS[2],lw=1)
    axes[1].plot(example["x"],example["hx"],color=COLORS[0],lw=1)
    axes[2].plot(example["s"],example["hs"],color=COLORS[3],lw=1)
    margin = .07*np.ptp(line[:,1])
    axes[0].set(xlabel=r"$x$",ylabel=r"$y$",title="(a) Contour (zoom)",
                xlim=(0,1024),ylim=(max(0,line[:,1].min()-margin),
                                    min(2048,line[:,1].max()+margin)))
    axes[0].set_xticks([0,512,1024])
    for spine in axes[0].spines.values():
        spine.set_visible(True)
    axes[1].set(xlabel=r"$x$",ylabel=r"$h(x)-\langle h(x)\rangle$",
                title="(b) Area-preserving profile",xlim=(0,1024))
    axes[2].set(xlabel=r"$s$",ylabel=r"$h(s)-\langle h(s)\rangle$",
                title="(c) Arc-length profile",xlim=(0,example["P"]))
    limits=[ax.get_ylim() for ax in axes[1:]]
    for ax in axes[1:]:
        ax.set_ylim(min(v[0] for v in limits),max(v[1] for v in limits))
        ax.axhline(0,c=".7",lw=.5,zorder=-1)
    for ext in ("pdf","png"):
        fig.savefig(OUT/f"fig3_edge_profiles.{ext}",dpi=200)
    plt.close(fig)

    fig,axes = plt.subplots(1,2,figsize=(10,4.3),layout="constrained")
    ax=axes[0]
    Ls=np.array(SIZES)
    Pmean=np.array([report[L]["meanP"] for L in SIZES])
    ci=np.array([report[L]["meanP_ci95"] for L in SIZES])
    ax.errorbar(Ls,Pmean,yerr=np.maximum(0,np.stack((Pmean-ci[:,0],ci[:,1]-Pmean))),
                fmt="o",c="#275c83",ms=7,capsize=3)
    xx=np.geomspace(220,4800,100)
    ax.plot(xx,.8*xx,"--",c="black",lw=1.1)
    ax.text(1400,.8*1400*.82,r"$L_x$",ha="center",va="top")
    coefficient=1.2*max(Pmean/Ls**1.25)
    ax.plot(xx,coefficient*xx**1.25,":",c="black",lw=1.6)
    ax.text(1400,coefficient*1400**1.25*1.15,r"$L_x^{1.25}$",ha="center")
    ax.set(xscale="log",yscale="log",xlabel=r"$L_x$",ylabel=r"$\langle P\rangle$",
           title="(a) Mean perimeter")
    ax.set_xticks(Ls,labels=[str(L) for L in Ls])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())

    ax=axes[1]
    mode_panel(ax,{L:mode_spectra(groups[L],L) for L in SIZES[2:]},SIZES[2:],[(.63, "--", .5), (1.25, ":", 2)],
               colors=COLORS,markers=MARKERS)
    handles=[Line2D([],[],c=c,marker=m,ls="none",label=fr"$L_x={L}$")
             for L,c,m in zip(SIZES[2:],COLORS[2:],MARKERS[2:])]
    handles.extend([Line2D([],[],c=".25",marker="o",ls="none",label=r"$h(x)$"),
                    Line2D([],[],c=".25",marker="o",mfc="none",ls="none",label=r"$h(s)$")])
    ax.legend(handles=handles,frameon=False,fontsize=8,ncol=2,loc="upper right")
    ax.set(title="(b) Height spectra")
    for ext in ("pdf","png"):
        fig.savefig(OUT/f"fig4_edge_scaling.{ext}",dpi=200)
    plt.close(fig)
    report["perimeter_fit"]=dict(alpha=float(np.polyfit(np.log(Ls),np.log(Pmean),1)[0]))
    report["example"]=dict(Lx=1024,seed=0,P=example["P"],multivalued_fraction=example["multi"])
    report["source"]=str(DATA.relative_to(ROOT))
    (OUT/"edge_statistics.json").write_text(json.dumps(report,indent=2)+"\n")
    for L in SIZES:
        print(L,report[L]["meanP"],report[L]["fits"]["x"]["0_0.03"],report[L]["fits"]["s"]["0_0.03"])
    print(report["perimeter_fit"],report["example"])


if __name__=="__main__":
    main()
