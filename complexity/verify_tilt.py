"""Check of the statistical-tilt mapping.

1. Sample by sample: direct count at force f with phases phi_i equals the
   f = 0 level count at W = f K with phases phi_i + f i (i+1) / 2.
2. In law: <ln N>, P(N = 0) and the distribution of f_c from the direct runs
   of run.py (data/, one enumeration per force, f_c by bisection) against
   the level runs of run_levels.py (data_levels/, independent samples).

    python3 verify_tilt.py      -> tilt_check.json
"""
import json
import os
import numpy as np
from qcount import count, count_levels

os.chdir(os.path.dirname(os.path.abspath(__file__)))
P = 4.0
out = {"exact": [], "law": []}

rng = np.random.default_rng(11)
nbad = 0
ntest = 0
for K in [5, 10, 15, 20, 25, 30]:
    for f in [0.0, 0.03, 0.08, 0.12, 0.18]:
        for A in [0.3, 0.6, 1.0]:
            ph = rng.uniform(0, P, K)
            r = count(ph, A, P, f, 8, False, 10, 3.0)
            i = np.arange(K)
            q = count_levels(np.mod(ph + f * i * (i + 1) / 2, P), A, P, np.array([f * K]), 4)
            ntest += 1
            nbad += int(r[0] != q[0][0] or r[2][0] != q[1][0])
out["exact"] = {"tests": ntest, "mismatches": nbad}
print(f"sample-by-sample: {ntest} tests, {nbad} mismatches")

for A in [0.3, 0.6, 1.0]:
    for K in [16, 24, 32, 40]:
        fd = f"data/A{A:g}_K{K}.npz"
        fl = f"data_levels/A{A:g}_K{K}.npz"
        if not (os.path.exists(fd) and os.path.exists(fl)):
            continue
        d, e = np.load(fd), np.load(fl)
        row = {"A": A, "K": K, "fc_bisect": [float(d["fc"].mean()), float(d["fc"].std() / np.sqrt(len(d["fc"])))],
               "fc_levels": [float(e["fc"].mean()), float(e["fc"].std() / np.sqrt(len(e["fc"])))],
               "f": [], "lnN_direct": [], "lnN_levels": []}
        for j, f in enumerate(d["f"][::3]):
            jj = int(np.argmin(abs(e["f"] - f)))
            if abs(e["f"][jj] - f) > 1e-7:
                continue
            a = np.log(np.maximum(d["N"][:, 3 * j], 1))
            b = np.log(np.maximum(e["N"][:, jj], 1))
            row["f"].append(float(f))
            row["lnN_direct"].append([float(a.mean()), float(a.std() / np.sqrt(len(a)))])
            row["lnN_levels"].append([float(b.mean()), float(b.std() / np.sqrt(len(b)))])
        out["law"].append(row)
        z = (row["fc_levels"][0] - row["fc_bisect"][0]) / np.hypot(row["fc_levels"][1], row["fc_bisect"][1])
        zl = [(y[0] - x[0]) / np.hypot(x[1], y[1]) for x, y in zip(row["lnN_direct"], row["lnN_levels"])]
        print(f"A={A} K={K}: <fc> {row['fc_bisect'][0]:.4f} vs {row['fc_levels'][0]:.4f} (z={z:+.1f}); "
              f"<lnN> z-scores at f={['%.3f' % x for x in row['f']]}: {['%+.1f' % x for x in zl]}")
json.dump(out, open("tilt_check.json", "w"), indent=1)
