"""Fractal dimension of the near-critical set in the shooting parameter.

    python3 fractal_top.py   -> fractal_top.json

For each sample, S = {a in [0, P): W_K(a) >= max_a W_K - delta} is a union of
intervals between an up-crossing (W' = det H > 0, even Morse index) and the
next down-crossing of the level.  N(eps) = number of eps-boxes meeting S.
If c = h(mu_max) = d * lambda (Ledrappier-Young), S is a fractal of
dimension d over scales e^{-lambda K} << eps << 1.
"""
import json
import numpy as np
from multiprocessing import Pool
from qcount import count_levels, record_level_a

P = 4.0
EPS = np.logspace(np.log10(0.5), -14, 60)


def boxes(starts, ends, eps):
    out = np.empty(len(eps))
    for q, e in enumerate(eps):
        lo = np.floor(starts / e).astype(np.int64)
        hi = np.floor(ends / e).astype(np.int64)
        n = 0
        last = -1
        for a, b in zip(lo, hi):
            a = max(a, last + 1)
            if b >= a:
                n += b - a + 1
                last = b
        out[q] = n
    return out


def work(args):
    seed, K, A, delta = args
    rng = np.random.default_rng(seed)
    ph = rng.uniform(0, P, K)
    Wmax = count_levels(ph, A, P, np.array([0.0]), 2)[4]
    uz, idx, ld, n, lp, a = record_level_a(ph, A, P, Wmax - delta * A, 5_000_000)
    o = np.argsort(a)
    a = a[o]; up = (idx[o] % 2 == 0)
    # rotate so the list starts with an up-crossing (S is periodic in a)
    k0 = int(np.argmax(up))
    a = np.concatenate([a[k0:], a[:k0] + P]); up = np.concatenate([up[k0:], up[:k0]])
    st = a[0::2]; en = a[1::2]
    if not (up[0::2].all() and (~up[1::2]).all()):
        return None
    return boxes(st, en, EPS), n, float(lp.mean())


if __name__ == "__main__":
    out = {}
    for A, Ks in [(0.3, [32, 48]), (0.6, [24, 32, 40]), (1.0, [20, 28])]:
        for K in Ks:
            for delta in (0.5, 2.0):
                with Pool(9) as p:
                    res = [r for r in p.map(work, [(55_000_000 + 1000 * K + s, K, A, delta) for s in range(300)]) if r is not None]
                logN = np.mean([np.log(r[0]) for r in res], 0)
                lam = np.mean([r[2] for r in res]) / K
                # local slopes d ln N / d ln(1/eps)
                sl = -np.gradient(logN, np.log(EPS))
                # fit in the window e^{-lam K} * 1e3 < eps < 1e-2
                m = (EPS < 1e-2) & (EPS > 1e3 * np.exp(-0.75 * lam * K))
                d = float(-np.polyfit(np.log(EPS[m]), logN[m], 1)[0]) if m.sum() > 3 else float("nan")
                out[f"A{A:g}_K{K}_d{delta:g}"] = dict(A=A, K=K, delta=delta, eps=EPS.tolist(), logN=logN.tolist(),
                                                     slope=sl.tolist(), d=d, n=len(res), lam_near=lam)
                print(f"A={A} K={K} delta={delta}A: d_box={d:.3f}  ({len(res)} samples, fit window {m.sum()} pts)  "
                      f"local slopes at eps=1e-3,1e-5,1e-7: {sl[np.argmin(abs(EPS-1e-3))]:.3f} {sl[np.argmin(abs(EPS-1e-5))]:.3f} {sl[np.argmin(abs(EPS-1e-7))]:.3f}",
                      flush=True)
    json.dump(out, open("fractal_top.json", "w"), indent=1)
