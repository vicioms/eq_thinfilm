"""Stretching rate of the critical (max-W) configuration versus the
near-threshold growth rate c.

    python3 lyap_leader.py   -> lyap_leader.json

For each sample: W_max = max_a W_K; the two crossings of W_max - 1e-9 are the
fold tip, i.e. the critical configuration; ln|psi_{K-1}| there is the log
stretching of the shooting curve along the critical trajectory.  Also the
mean ln|psi| over the equilibria within A/2 of the top, and over all
equilibria at f = 0 for comparison.
"""
import json
import numpy as np
from multiprocessing import Pool
from qcount import count_levels, record_level

P = 4.0


def work(args):
    seed, K, A = args
    rng = np.random.default_rng(seed)
    ph = rng.uniform(0, P, K)
    Wmax = count_levels(ph, A, P, np.array([0.0]), 2)[4]
    lps = record_level(ph, A, P, Wmax - 1e-9, 10)[4]
    near = record_level(ph, A, P, Wmax - 0.5 * A, 2_000_000)[4]
    bulk = record_level(ph, A, P, 0.0, 2_000_000)[4]
    return float(lps.max()), float(near.mean()), float(bulk.mean())


if __name__ == "__main__":
    out = {}
    for A, Ks in [(0.3, [16, 24, 32, 40, 48]), (0.6, [12, 16, 20, 24, 28, 32]), (1.0, [8, 12, 16, 20, 24])]:
        rows = []
        for K in Ks:
            with Pool(9) as p:
                r = np.array(p.map(work, [(77_000_000 + 1000 * K + s, K, A) for s in range(400)]))
            rows.append(dict(K=K, crit=r[:, 0].mean(), crit_err=r[:, 0].std() / 20, near=r[:, 1].mean(),
                             bulk=r[:, 2].mean()))
            print(f"A={A} K={K}: ln|psi| crit={r[:,0].mean():.3f}  near-top={r[:,1].mean():.3f}  bulk={r[:,2].mean():.3f}", flush=True)
        K = np.array([x["K"] for x in rows])
        sl = {key: float(np.polyfit(K[-4:], [x[key] for x in rows][-4:], 1)[0]) for key in ("crit", "near", "bulk")}
        out[str(A)] = dict(rows=rows, slopes=sl)
        print(f"A={A}: slopes  crit={sl['crit']:.4f}  near-top={sl['near']:.4f}  bulk={sl['bulk']:.4f}", flush=True)
    json.dump(out, open("lyap_leader.json", "w"), indent=1)
