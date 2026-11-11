#!/usr/bin/env python3
"""
Batch scan of bubble depinning: anisotropy a x system size L x disorder samples.

Each task runs the quasi-static protocol of bubble_depinning.py with the
dislocation noise (random-field coupling) and stores one JSON line with h_d and
the shape of the critical bubble. The initial radius scales with the box,
R0 = L * r0_frac, so that bubbles of all sizes are geometrically similar.

Backends
  cpu  : tasks run in parallel processes, one numba thread each (--workers).
  cuda : one worker process per GPU (--gpus N, devices 0..N-1), each running
         its tasks one after the other on its own GPU.

Cluster use: split the task list over a job array with --shard i --nshards n.
Each shard writes its own results_shard<i>.jsonl (no concurrent writes);
fss_analyze.py reads all results*.jsonl files. See slurm_fss.sh.

The scan can be interrupted and restarted: finished (L, a, seed) found in any
results*.jsonl of the output folder are skipped.

  python fss_scan.py --aniso 0 0.25 0.5 0.75 1 --L 64 128 256 --samples 64:48 128:32 256:16
  python fss_scan.py --backend cuda --gpus 4 --L 512 1024 --samples 512:32 1024:16
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))


_GPU = None


def _init_worker(backend, gpu_queue):
    """Pin each worker process to one GPU (cuda) or one CPU thread (cpu)."""
    global _GPU
    sys.path.insert(0, HERE)
    if backend == "cuda":
        _GPU = gpu_queue.get()


def run_task(task):
    if task["opts"]["backend"] == "cpu":
        os.environ["NUMBA_NUM_THREADS"] = "1"
    sys.path.insert(0, HERE)
    import numpy as np
    import bubble_depinning as bd

    L, a, seed, opts = task["L"], task["aniso"], task["seed"], task["opts"]
    p = bd.Params(L=L, noise="dislocation", coupling="field", Delta=opts["Delta"], w=opts["w"],
                  burgers_angle=opts["burgers_angle"], aniso=a, R0=L * opts["r0_frac"],
                  h_start=opts["h_start"], dh=opts["dh"], n_bisect=opts["n_bisect"], seed=seed)
    t0 = time.time()
    tries = []
    for _ in range(12):          # adapt h_start if the initial bubble collapses / escapes
        res = bd.run_depinning(p, verbose=False, backend=opts["backend"], keep_snapshots=False,
                               device=_GPU)
        tries.append((p.h_start, res["status_init"]))
        if res["status_init"] == "pinned":
            break
        if res["status_init"] == "collapsed":
            p.h_start += opts["h_step"]
        else:
            p.h_start = max(0.0, p.h_start - opts["h_step"])
            if p.h_start == 0.0 and tries[-1][0] == 0.0:
                break
    out = dict(L=L, aniso=a, seed=seed, Delta=p.Delta, w=p.w, burgers_angle=p.burgers_angle,
               R0=p.R0, h_start=p.h_start, tries=tries, status=res["status_init"],
               time=time.time() - t0, backend=opts["backend"], gpu=_GPU)
    if res["status_init"] == "pinned":
        an = bd.analyse(res["phi_c"])
        out.update(h_d=float(res["h_d"]), n_steps=len(res["h"]))
        if an is not None:
            out.update({k: float(an[k]) for k in
                        ("area", "R_eff", "aspect", "orientation_deg", "a2", "frac_multivalued", "u_o")})
        if opts["save_configs"]:
            d = os.path.join(opts["outdir"], "configs")
            os.makedirs(d, exist_ok=True)
            np.savez_compressed(os.path.join(d, f"L{L}_a{a}_s{seed}.npz"),
                                phi_c=res["phi_c"].astype(np.float16), h_d=res["h_d"])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--aniso", type=float, nargs="+", default=[0, 0.25, 0.5, 0.75, 1.0])
    ap.add_argument("--L", type=int, nargs="+", default=[64, 128, 256])
    ap.add_argument("--samples", nargs="+", default=["64:48", "128:32", "256:16"],
                    help="L:n pairs (number of disorder samples per L)")
    ap.add_argument("--Delta", type=float, default=0.08)
    ap.add_argument("--w", type=float, default=3.0)
    ap.add_argument("--burgers-angle", type=float, default=0.0)
    ap.add_argument("--r0-frac", type=float, default=1 / 6)
    ap.add_argument("--h-start", type=float, default=0.0)
    ap.add_argument("--h-step", type=float, default=0.01, help="h_start correction on collapse")
    ap.add_argument("--dh", type=float, default=0.002)
    ap.add_argument("--n-bisect", type=int, default=10)
    ap.add_argument("--backend", choices=["cpu", "cuda"], default="cpu")
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1),
                    help="CPU worker processes (cpu backend)")
    ap.add_argument("--gpus", type=int, default=1, help="number of GPUs to use (cuda backend)")
    ap.add_argument("--shard", type=int, default=0, help="index of this shard (job-array task id)")
    ap.add_argument("--nshards", type=int, default=1, help="total number of shards")
    ap.add_argument("--outdir", default=os.path.join(HERE, "runs", "fss"))
    ap.add_argument("--no-configs", action="store_true", help="do not store critical configurations")
    args = ap.parse_args()

    nsamp = {int(k): int(v) for k, v in (s.split(":") for s in args.samples)}
    os.makedirs(args.outdir, exist_ok=True)
    out_file = os.path.join(args.outdir, "results.jsonl" if args.nshards == 1
                            else f"results_shard{args.shard}.jsonl")
    done = set()
    for fn in glob.glob(os.path.join(args.outdir, "results*.jsonl")):
        with open(fn) as f:
            for line in f:
                r = json.loads(line)
                done.add((r["L"], r["aniso"], r["seed"]))

    opts = dict(Delta=args.Delta, w=args.w, burgers_angle=args.burgers_angle, r0_frac=args.r0_frac,
                h_start=args.h_start, h_step=args.h_step, dh=args.dh, n_bisect=args.n_bisect,
                save_configs=not args.no_configs, outdir=args.outdir, backend=args.backend)
    # large L first so that the slow tasks do not end up alone at the end
    # the full, deterministic task list is split round-robin over shards
    all_tasks = [dict(L=L, aniso=a, seed=s, opts=opts)
                 for L in sorted(args.L, reverse=True) for s in range(nsamp.get(L, 16)) for a in args.aniso]
    tasks = [t for t in all_tasks[args.shard::args.nshards] if (t["L"], t["aniso"], t["seed"]) not in done]
    nproc = args.gpus if args.backend == "cuda" else args.workers
    print(f"shard {args.shard}/{args.nshards}: {len(tasks)} tasks to do, {args.backend} x {nproc} -> {out_file}",
          flush=True)

    import multiprocessing as mp
    ctx = mp.get_context("spawn")
    gpu_queue = ctx.Manager().Queue()
    for i in range(args.gpus):
        gpu_queue.put(i)
    if args.backend == "cpu":
        os.environ["NUMBA_NUM_THREADS"] = "1"
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=nproc, mp_context=ctx, initializer=_init_worker,
                             initargs=(args.backend, gpu_queue)) as ex, open(out_file, "a") as f:
        futs = [ex.submit(run_task, t) for t in tasks]
        for i, fu in enumerate(as_completed(futs), 1):
            r = fu.result()
            f.write(json.dumps(r) + "\n"); f.flush()
            print(f"[{i}/{len(tasks)} {time.time()-t0:7.0f}s] L={r['L']} a={r['aniso']} seed={r['seed']} "
                  f"{r['status']} h_d={r.get('h_d', float('nan')):.5f} aspect={r.get('aspect', float('nan')):.2f} "
                  f"orient={r.get('orientation_deg', float('nan')):.0f} ({r['time']:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
