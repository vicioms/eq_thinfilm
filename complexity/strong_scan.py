"""Infinite-chain depinning threshold at strong pinning (k > 2) from the exact
lattice minimax (strong_exact.fc_dp, banded kernel r = 2, error ~1e-2).

    python3 strong_scan.py   -> strong_scan.json
"""
import json
import numpy as np
from multiprocessing import Pool
from strong_exact import fc_dp
from annealed import gle, f_ann

P = 4.0
KS = [10, 16, 20, 40, 80, 160, 320]


def work(a):
    seed, K, A = a
    ph = np.random.default_rng(seed).uniform(0, P, K)
    return fc_dp(ph, A, r=2)


if __name__ == "__main__":
    out = {}
    with Pool(9) as pool:
        for A in (3.0, 6.0, 12.0, 24.0, 48.0):
            rows = []
            for K in KS:
                n = 400 if K <= 80 else 200
                v = np.array(pool.map(work, [(91_000_000 + 1000 * K + s, K, A) for s in range(n)]))
                rows.append(dict(K=K, mean=float(v.mean()), err=float(v.std() / np.sqrt(n)), sd=float(v.std())))
            K = np.array(KS, float); m = np.array([x["mean"] for x in rows]); e = np.array([x["err"] for x in rows])
            sel = K >= 40
            c = np.polyfit(1 / K[sel], m[sel], 1, w=1 / e[sel])      # f_c(K) = f_inf + a/K
            k = A * 4 / P
            fa = f_ann(A, P, gle(1.0, k))
            out[str(A)] = dict(A=A, k=k, rows=rows, fc_inf=float(c[1]), slope_1overK=float(c[0]), f_ann=fa,
                               lead=A - P / 2, shift=float(c[1] - (A - P / 2)))
            print(f"k={k:5.1f}: f_c(inf)={c[1]:.4f}  A-P/2={A-P/2:.1f}  shift={c[1]-(A-P/2):.4f}  "
                  f"k*shift/P={k*(c[1]-(A-P/2))/P:.3f}  f_ann={fa:.4f}  f_ann/f_c={fa/c[1]:.3f}  "
                  f"sd*sqrt(K) at K=320: {rows[-1]['sd']*np.sqrt(320):.3f} (P/sqrt12={P/np.sqrt(12):.3f})", flush=True)
            print("      <f_c(K)>:", " ".join(f"{x['K']}:{x['mean']:.4f}" for x in rows), flush=True)
    json.dump(out, open("strong_scan.json", "w"), indent=1)
