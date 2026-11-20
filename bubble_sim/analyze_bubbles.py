#!/usr/bin/env python3
"""Critical bubbles (fss_scan.py, bubble geometry): the whole analysis in one place.

    python analyze_bubbles.py --extract runs/fss_D0.04          # where the scan ran -> contours.npz
    python analyze_bubbles.py runs/fss_D0.04 [runs/fss_D0.04_R32 ...] [--out D] [--qmax 0.1] [--qlarge 0.03]

One row per run in every figure, written to --out (default runs/<common name>_analysis) with report.txt.
The same observables as analyze_flat.py, measured on the closed phi = 0 contour of the bubble:
  configs.png    example bubbles, angle between the wall line and b per L, mean shape <R(theta)>
  roughness.png  pooled tangent-angle spectra S_psi(q) ~ q^(1-2 zeta) (modes n >= 12: the lower
                 modes are the bubble's shape), local zeta, zeta vs L (intrinsic and from the
                 polar profile R(theta), KFR Eq. 9 about the centroid); text: memory split
  threshold.png  <h_d>(L) = h_c + A L^(-1/nu) + B/L (B: Laplace pressure, R0 ∝ L), var(h_d)
  length.png     perimeter P ∝ R_eff^alpha, P / 2 pi R_eff (as NB7_perimeter.ipynb), overhangs
  rdensity.png   roughness density chi(q, r) (Caballero & Giamarchi 2022, analysis/rdensity.py): B(r) from
                 -2 ln chi = q^2 B(r) without a single-valued height (overhangs included: B(0) > 0),
                 B(r) per L, B(r = 32) vs L (∝ L^(2 zeta - 2) when zeta > 1, flat for qKPZ)
  length_dir.png perimeter by direction from b (edge flanks vs screw tips): the excess over the
                 smooth shape ∝ R_eff^(alpha - 1), alpha = 1 for qKPZ, zeta for qEW
With a fixed R0 (fss_scan.py --r0) the critical bubble stays near R0 whatever L: there h_d and
the shapes do not scale with L.
Reference values: qEW zeta = 1.25, nu = 4/3; qKPZ / DPD zeta = 0.633, nu = 1.733.
"""
from __future__ import annotations

import argparse
import os

from analysis import data, figures, measure

CFG = dict(bubble=dict(n_min=12, lmin=256), qmax=0.1, qlarge=0.03)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="*", help="bubble runs to analyse")
    ap.add_argument("--extract", metavar="DIR", help="extract the contours of DIR/configs -> DIR/contours.npz")
    ap.add_argument("--out")
    ap.add_argument("--qmax", type=float, default=CFG["qmax"], help="upper end of the zeta fits")
    ap.add_argument("--qlarge", type=float, default=CFG["qlarge"], help="upper end of the large-scale zeta fits")
    ap.add_argument("--nmin", type=int, default=CFG["bubble"]["n_min"], help="lowest mode in the zeta fits")
    ap.add_argument("--lmin", type=int, default=CFG["bubble"]["lmin"], help="smallest L in the threshold fits")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    args = ap.parse_args()
    if args.extract:
        data.extract_bubbles(args.extract, args.workers)
    if not args.dirs:
        return
    CFG["qmax"], CFG["qlarge"] = args.qmax, args.qlarge
    CFG["bubble"].update(n_min=args.nmin, lmin=args.lmin)
    out = args.out or data.default_out(args.dirs)
    os.makedirs(out, exist_ok=True)
    sets = [measure.measure_run(data.load_bubbles(d), args.workers) for d in args.dirs]
    txt = [f"runs: {', '.join(args.dirs)}"] + [f"  {W['label']}: {len(W['L'])} bubbles" for W in sets]
    figs = [figures.configs(sets, out), figures.roughness(sets, out, CFG, txt),
            figures.threshold(sets, out, CFG, txt), figures.length(sets, out, txt),
            figures.rdensity(sets, out, txt),
            figures.length_dir(sets, out, txt)]
    with open(os.path.join(out, "report.txt"), "w") as f:
        f.write("\n".join(txt) + "\n")
    print("\n".join(txt))
    print("saved", ", ".join(figs + [os.path.join(out, "report.txt")]))


if __name__ == "__main__":
    main()
