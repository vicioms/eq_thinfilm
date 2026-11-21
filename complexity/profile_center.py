"""Per-sample profile of the levels W = x at which equilibria sit (f = 0
curve, tilt-symmetry frame): centre c_s and width chi_q, and the spread of
the centres between samples, Var(c_s)/K = chi_c.  Kac-Rice: chi_q + chi_c
must equal <V'^2> = A^2/3 asymptotically (Gaussian part of the profile).

    python3 profile_center.py --A 0.6 --K 12 16 20 24 28 32 36 --nsamp 1500
writes profile_center_A{A}.json
"""
import argparse, json
import numpy as np
from multiprocessing import Pool
from qcount import count_levels

P = 4.0


def work(args):
    seed, K, A, levels = args
    rng = np.random.default_rng(seed)
    ph = rng.uniform(0.0, P, K)
    N = count_levels(ph, A, P, levels, 2)[0].astype(float)
    Z = N.sum()
    c = (N * levels).sum() / Z
    v = (N * (levels - c) ** 2).sum() / Z
    return c, v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--A", type=float, default=0.6)
    ap.add_argument("--K", type=int, nargs="+", required=True)
    ap.add_argument("--nsamp", type=int, default=1500)
    ap.add_argument("--procs", type=int, default=3)
    a = ap.parse_args()
    out = {"A": a.A, "chi": a.A ** 2 / 3, "K": [], "chi_c": [], "chi_c_err": [], "chi_q": [], "chi_q_err": []}
    for K in a.K:
        levels = np.linspace(-a.A * K, a.A * K, 801)
        jobs = [(900_000_000 + 10_000 * K + s, K, a.A, levels) for s in range(a.nsamp)]
        with Pool(a.procs) as p:
            r = np.array(p.map(work, jobs, chunksize=20))
        c, v = r[:, 0], r[:, 1]
        n = len(c)
        cc = np.mean(c ** 2) / K  # <c> = 0 by symmetry
        cce = np.std(c ** 2) / np.sqrt(n) / K
        out["K"].append(K); out["chi_c"].append(cc); out["chi_c_err"].append(cce)
        out["chi_q"].append(v.mean() / K); out["chi_q_err"].append(v.std() / np.sqrt(n) / K)
        print(f"K={K}: chi_c={cc:.4f}+-{cce:.4f}  chi_q={v.mean()/K:.4f}  sum={cc+v.mean()/K:.4f}  chi={a.A**2/3:.4f}", flush=True)
    json.dump(out, open(f"profile_center_A{a.A:g}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
