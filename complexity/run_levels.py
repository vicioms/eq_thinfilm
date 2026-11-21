"""All forces from one f = 0 enumeration (statistical tilt symmetry).

    python3 run_levels.py --A 0.6 --K 20 32 44 --nsamp 2000

For each sample, run.py would need one enumeration per force plus a bisection
for f_c.  Here a single f = 0 enumeration gives, in law and for every f on a
fine grid (see qcount.count_levels):

  N[s, j], Nmin[s, j]  number of equilibria / minima at f_j
  ldet[s, j]           mean ln|det H| / K over them (nan if none)
  Y2[s, j]             participation sum_b (n_b / N)^2 of their positions
                       u mod P over nub bins: 1/nub if spread, 1 if localised
  fc[s]                max_a W_K(a) / K, distributed as the sample threshold

Output: data_levels/A{A}_K{K}.npz (triangle wave, P = 4, kappa = 1).

--model cos: piecewise-linear interpolation (M nodes per period, random grid
offset) of FLRT's Gaussian disorder V_i(u) = xi1 cos u + xi2 sin u,
<xi^2> = sigma^2, period 2 pi; output data_levels/cos_s{sigma}_M{M}_K{K}.npz.
"""
import argparse
import os
import numpy as np
from multiprocessing import Pool
from qcount import count_levels, count_levels_pl, cos_nodes

P = 4.0
NUB = 32


def work(args):
    seed, K, A, fgrid, model, M = args
    rng = np.random.default_rng(seed)
    levels = fgrid * K
    if model == "tri":
        phase = rng.uniform(0.0, P, K)
        N, Nm, sl, ub, Wmax, Wmin, nodes, over = count_levels(phase, A, P, levels, NUB)
    else:
        G, off = cos_nodes(rng, K, M, A)
        N, Nm, sl, ub, Wmax, Wmin, nodes, over = count_levels_pl(G, off, 2 * np.pi, levels, NUB)
    with np.errstate(invalid="ignore", divide="ignore"):
        ld = np.where(N > 0, sl / np.maximum(N, 1), np.nan)
        Y2 = np.where(N > 0, (ub.astype(float) ** 2).sum(1) / np.maximum(N, 1) ** 2, np.nan)
    return N, Nm, ld, Y2, Wmax / K, Wmin / K, nodes, over


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--A", type=float, required=True,
                    help="triangle amplitude A, or sigma for --model cos")
    ap.add_argument("--model", choices=["tri", "cos"], default="tri")
    ap.add_argument("--M", type=int, default=16)
    ap.add_argument("--K", type=int, nargs="+", required=True)
    ap.add_argument("--nsamp", type=int, default=2000)
    ap.add_argument("--nf", type=int, default=121)
    ap.add_argument("--procs", type=int, default=9)
    ap.add_argument("--seed0", type=int, default=500_000_000)
    ap.add_argument("--fmaxfrac", type=float, default=None,
                    help="force grid up to fmaxfrac * A (default 0.6 tri, 1.0 cos)")
    a = ap.parse_args()
    frac = a.fmaxfrac if a.fmaxfrac is not None else (0.6 if a.model == "tri" else 1.0)
    fgrid = np.linspace(0.0, frac * a.A, a.nf)
    os.makedirs("data_levels", exist_ok=True)
    for K in a.K:
        seeds = [a.seed0 + 10_000_000 * K + s for s in range(a.nsamp)]
        with Pool(a.procs) as pool:
            res = pool.map(work, [(sd, K, a.A, fgrid, a.model, a.M) for sd in seeds],
                           chunksize=max(1, a.nsamp // (40 * a.procs)))
        out = {name: np.array([r[i] for r in res]) for i, name in
               enumerate(["N", "Nmin", "ldet", "Y2", "fc", "wmin", "nodes", "overflow"])}
        if a.model == "tri":
            fn = f"data_levels/A{a.A:g}_K{K}.npz"
        else:
            fn = f"data_levels/cos_s{a.A:g}_M{a.M}_K{K}.npz"
        np.savez(fn, f=fgrid, A=a.A, P=P if a.model == "tri" else 2 * np.pi, K=K,
                 M=a.M, model=a.model, seeds=np.array(seeds), **out)
        print(f"{a.model} A={a.A} M={a.M} K={K}: <lnN(0)>/K={np.mean(np.log(np.maximum(out['N'][:,0],1)))/K:.4f} "
              f"<fc>={out['fc'].mean():.4f} nodes={out['nodes'].mean():.3g} "
              f"overflow={out['overflow'].any()}", flush=True)


if __name__ == "__main__":
    main()
