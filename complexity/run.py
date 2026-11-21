"""Per-sample exact counts N(f), N_min(f) and the sample threshold f_c.

    python3 run.py --A 0.6 --K 8 12 16 20 --nsamp 2000 --nf 13

writes data/A{A}_K{K}.npz with, per sample s and force f_j:
  N[s, j]      total number of equilibria (mod the global shift by P)
  Nmin[s, j]   number of local minima (Morse index 0)
  ldet[s, j]   mean over equilibria of ln|det H| / K   (nan if N = 0)
  fc[s]        sample threshold sup{f : N(f) > 0} (bisection, |error| < 1e-6 A)
  nodes[s, j]  leaves visited (cost)
P = 4 throughout (so k = A, L_c = A^{-2/3}), kappa = 1.
"""
import argparse
import os
import numpy as np
from multiprocessing import Pool
from qcount import count

P = 4.0


def sample_fc(phase, A):
    # N(f) > 0 is monotone in f only on average; the threshold of a finite
    # chain is still well defined by no-passing: a pinned configuration exists
    # iff f <= f_c.  Bisection with early exit.
    lo, hi = 0.0, A
    for _ in range(22):
        mid = 0.5 * (lo + hi)
        r = count(phase, A, P, mid, 8, True, 10, 3.0)
        if r[0] > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def work(args):
    seed, K, A, fgrid, do_fc = args
    rng = np.random.default_rng(seed)
    phase = rng.uniform(0.0, P, K)
    nf = len(fgrid)
    N = np.zeros(nf, np.int64)
    Nmin = np.zeros(nf, np.int64)
    ld = np.full(nf, np.nan)
    nodes = np.zeros(nf, np.int64)
    over = False
    for j, f in enumerate(fgrid):
        r = count(phase, A, P, f, 8, False, 10, 3.0)
        N[j] = r[0]
        Nmin[j] = r[2][0]
        nodes[j] = r[1]
        if r[0] > 0:
            ld[j] = r[4] / r[0]
        over |= r[5]
    fc = sample_fc(phase, A) if do_fc else np.nan
    return N, Nmin, ld, nodes, fc, over


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--A", type=float, required=True)
    ap.add_argument("--K", type=int, nargs="+", required=True)
    ap.add_argument("--nsamp", type=int, default=1000)
    ap.add_argument("--fmax", type=float, default=None,
                    help="largest force on the grid (default: A)")
    ap.add_argument("--nf", type=int, default=13)
    ap.add_argument("--nofc", action="store_true")
    ap.add_argument("--procs", type=int, default=9)
    ap.add_argument("--seed0", type=int, default=0)
    a = ap.parse_args()
    fmax = a.A if a.fmax is None else a.fmax
    fgrid = np.linspace(0.0, fmax, a.nf)
    os.makedirs("data", exist_ok=True)
    for K in a.K:
        seeds = [10_000_000 * K + a.seed0 + s for s in range(a.nsamp)]
        jobs = [(sd, K, a.A, fgrid, not a.nofc) for sd in seeds]
        with Pool(a.procs) as pool:
            res = pool.map(work, jobs, chunksize=max(1, a.nsamp // (20 * a.procs)))
        N = np.array([r[0] for r in res])
        Nmin = np.array([r[1] for r in res])
        ld = np.array([r[2] for r in res])
        nodes = np.array([r[3] for r in res])
        fc = np.array([r[4] for r in res])
        over = any(r[5] for r in res)
        fn = f"data/A{a.A:g}_K{K}.npz"
        np.savez(fn, N=N, Nmin=Nmin, ldet=ld, nodes=nodes, fc=fc, f=fgrid,
                 A=a.A, P=P, K=K, seeds=np.array(seeds), overflow=over)
        lnN0 = np.log(np.maximum(N[:, 0], 1))
        print(f"A={a.A} K={K}: <lnN(0)>/K={lnN0.mean()/K:.4f} "
              f"ln<N(0)>/K={np.log(N[:,0].mean())/K:.4f} "
              f"<fc>={np.nanmean(fc):.4f} nodes={nodes.mean():.3g} "
              f"overflow={over}", flush=True)


if __name__ == "__main__":
    main()
