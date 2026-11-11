#!/usr/bin/env python3
"""
Check the CUDA backend against the CPU kernel.

On a machine with an NVIDIA GPU:      python test_gpu_backend.py
Without a GPU (slow, small sizes):    NUMBA_ENABLE_CUDASIM=1 python test_gpu_backend.py --small

1. single sweeps from the same random state must agree to round-off;
2. a full quasi-static run must give the same h_d (to the bisection precision,
   up to the slightly different convergence-check cadence).
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bubble_depinning as bd
from gpu_backend import CUDABackend


def test_sweeps(L, n=6):
    p = bd.Params(L=L, noise="dislocation", coupling="field", Delta=0.1, dJx=0.3, dJy=0.1, seed=3)
    rng = np.random.default_rng(1)
    r, g, Jx, Jy = bd.make_disorder(p, rng)
    r = 0.3 * rng.uniform(-1, 1, (L, L))                  # exercise bond + field + couplings together
    phi0 = rng.uniform(-1.2, 1.2, (L, L))
    be = CUDABackend(r, g, Jx, Jy, p)
    for fo in (False, True):
        a = phi0.copy()
        d = be.put(phi0)
        for _ in range(n):
            for c in (0, 1):
                bd._sweep(a, r, g, Jx, Jy, 0.05, p.eps0, fo, c)
            be._sweep(d, 0.05, fo, False)
        err = np.abs(be.get(d) - a).max()
        print(f"sweeps L={L} forward_only={fo}: max |phi_gpu - phi_cpu| = {err:.2e}")
        assert err < 1e-9, err


def test_run(L, R0, h_start, dh, n_bisect, **noise):
    kw = dict(L=L, R0=R0, h_start=h_start, dh=dh, n_bisect=n_bisect, seed=0, **noise)
    out = {}
    for name in ("cpu", "cuda"):
        t0 = time.time()
        res = bd.run_depinning(bd.Params(**kw), verbose=False, backend=name, keep_snapshots=False)
        out[name] = res
        print(f"{name:5s}: status={res['status_init']}  h_d={res['h_d']:.6f}  steps={len(res['h'])}  "
              f"area_c={bd.bubble_mask(res['phi_c']).sum()}  ({time.time()-t0:.1f}s)")
    tol = dh / 2 ** n_bisect * 2
    dhd = abs(out["cpu"]["h_d"] - out["cuda"]["h_d"])
    print(f"|h_d(cpu) - h_d(cuda)| = {dhd:.2e}  (bisection resolution {dh / 2 ** n_bisect:.1e})")
    assert out["cpu"]["status_init"] == out["cuda"]["status_init"]
    assert dhd <= tol, dhd


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--small", action="store_true", help="tiny sizes for the CUDA simulator")
    a = ap.parse_args()
    if a.small:
        test_sweeps(16)
        test_run(L=32, R0=10, h_start=0.05, dh=0.005, n_bisect=4, noise="white", coupling="bond", Delta=0.6)
    else:
        test_sweeps(256)
        test_run(L=256, R0=40, h_start=0.0, dh=0.002, n_bisect=10, noise="dislocation", coupling="field", Delta=0.08)
        test_run(L=128, R0=20, h_start=0.02, dh=0.002, n_bisect=10, noise="white", coupling="bond", Delta=0.6)
    print("OK")
