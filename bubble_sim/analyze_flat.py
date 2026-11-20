#!/usr/bin/env python3
"""Critical flat walls (fss_scan.py --geometry flat): the whole analysis in one place.

    python analyze_flat.py --extract runs/flat_D0.04_b0_A2     # where the scan ran -> walls.npz
    python analyze_flat.py runs/flat_D0.04_b0_A2 runs/flat_D0.04_b90_A2 [--out D] [--qmax 0.1] [--qlarge 0.03]

One row per run (orientation) in every figure, written to --out (default runs/<common name>_analysis,
e.g. runs/flat_D0.04_A2_analysis) with report.txt:
  configs.png    example walls, angle between the wall line and b (faceting) per L
  roughness.png  pooled tangent-angle spectra S_psi(q) ~ q^(1-2 zeta) (modes n >= 2) compensated for qEW
                 (q^1.5 S_psi) and qKPZ (q^0.27 S_psi), local zeta,
                 zeta vs L (intrinsic and from the KFR column profile u(x)) for q <= qmax and for the
                 large scales q <= qlarge only (edge walls cross over from qEW to qKPZ at ~100-150 sites);
                 text: memory split
  threshold.png  <h_d>(L) = h_c + A L^(-1/nu) with nu = 4/3 (qEW), 1.733 (qKPZ) or free, var(h_d)
  hs.png         spectrum of the height along the arc length, h(s) (y(s)), S_h ~ q^-(1+2 zeta):
                 single-valued with overhangs; compensated for qEW / qKPZ, local zeta, zeta vs L
  length.png     wall length P ∝ L^alpha (alpha = 1 for zeta < 1, zeta for zeta > 1), overhangs
  rdensity.png   roughness density chi(q, r) (Caballero & Giamarchi 2022, analysis/rdensity.py): B(r) from
                 -2 ln chi = q^2 B(r) without a single-valued height (overhangs included: B(0) > 0),
                 B(r) per L, B(r = 32) vs L (∝ L^(2 zeta - 2) when zeta > 1, flat for qKPZ)
Reference values: qEW zeta = 1.25, nu = 4/3; qKPZ / DPD zeta = 0.633, nu = 1.733.
"""
from __future__ import annotations

import argparse
import os

from analysis import data, figures, measure

CFG = dict(flat=dict(n_min=2, lmin=0), qmax=0.1, qlarge=0.03)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="*", help="flat-wall runs to analyse (one per orientation)")
    ap.add_argument("--extract", metavar="DIR", help="extract the walls of DIR/configs -> DIR/walls.npz")
    ap.add_argument("--out")
    ap.add_argument("--qmax", type=float, default=CFG["qmax"], help="upper end of the zeta fits")
    ap.add_argument("--qlarge", type=float, default=CFG["qlarge"],
                    help="upper end of the large-scale zeta fits (past the qEW -> qKPZ crossover, ~100-150 sites)")
    ap.add_argument("--nmin", type=int, default=CFG["flat"]["n_min"], help="lowest mode in the zeta fits")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    args = ap.parse_args()
    if args.extract:
        data.extract_flat(args.extract, args.workers)
    if not args.dirs:
        return
    CFG["qmax"], CFG["qlarge"], CFG["flat"]["n_min"] = args.qmax, args.qlarge, args.nmin
    out = args.out or data.default_out(args.dirs)
    os.makedirs(out, exist_ok=True)
    sets = [measure.measure_run(data.load_flat(d), args.workers) for d in args.dirs]
    txt = [f"runs: {', '.join(args.dirs)}"] + [f"  {W['label']}: {len(W['L'])} walls" for W in sets]
    figs = [figures.configs(sets, out), figures.roughness(sets, out, CFG, txt),
            figures.threshold(sets, out, CFG, txt), figures.hs_spectrum(sets, out, CFG, txt),
            figures.length(sets, out, txt),
            figures.rdensity(sets, out, txt)]
    with open(os.path.join(out, "report.txt"), "w") as f:
        f.write("\n".join(txt) + "\n")
    print("\n".join(txt))
    print("saved", ", ".join(figs + [os.path.join(out, "report.txt")]))


if __name__ == "__main__":
    main()
