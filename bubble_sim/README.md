# Bubble depinning in a disordered Ginzburg-Landau model

`bubble_depinning.py` simulates the 2D Ginzburg-Landau domain-wall model of
Kolton, Ferrero & Rosso (arXiv:2306.13415) starting from a circular bubble,
with optional anisotropic disorder. The model, the algorithm and every option
are described in the module docstring (`python bubble_depinning.py -h`).

Requires numpy, scipy, matplotlib and numba. Run from this folder:

```bash
python bubble_depinning.py --L 256 --R0 30 --out runs/dislocation/my_run
```

Each depinning run writes `<out>.png` (summary figure) and `<out>.npz`
(disorder, h ramp, area, metastable snapshots, critical configuration).
`.npz` files are git-ignored.

## Disorder options

| `--noise`               | what it is                                                                 | default coupling / Δ |
|-------------------------|----------------------------------------------------------------------------|----------------------|
| `dislocation` (default) | resolved shear stress of a random alloy, spectrum (k·b̂)²/\|k\|·e^{-k²w²} (Geslin et al., SM Eq. 21); `--burgers-angle`, `--w` | random field, Δ = 0.08 |
| `gaussian`              | Gaussian correlations with lengths `--xi-x`, `--xi-y`, rotation `--theta` (`inf` = columnar) | random bond, Δ = 0.6 |
| `white`                 | uncorrelated, isotropic (the original paper)                               | random bond, Δ = 0.6 |

`--coupling field|bond` overrides how the noise enters the equation.
`--dJx`, `--dJy` add random elastic couplings on x / y bonds; `--cx`, `--cy`
set the mean stiffness along x / y.

`--check-noise` only generates the dislocation noise and compares its
correlations with SM Eqs. 12-13 (reproduces SM Fig. 1).

`--mode dynamics --h <field>` integrates the equation in time at fixed h
instead of the quasi-static ramp.

## Finite-size scaling scan

`fss_scan.py` runs anisotropy a × size L × disorder samples (dislocation
noise; `--aniso a` interpolates from isotropic, a = 0, to the letter's
spectrum, a = 1, at fixed radial spectrum and variance). `R0 = L/6` (`--r0` sets a fixed radius instead). Results
go to `runs/fss/results*.jsonl` (+ critical configurations in
`runs/fss/configs/`); restarting skips finished samples.
`fss_analyze.py` makes `runs/fss/summary.csv` and `runs/fss/fss_summary.png`
(⟨h_d⟩ vs 1/L, var(h_d) ~ L^(−2/ν), bubble elongation vs a, overhangs).

`fss_shapes.py` looks at the critical bubbles themselves:
`bubbles.png` (snapshots, mean shape ⟨R(θ)⟩/R_eff, overhang probability vs θ)
and `Sq.png` (interface structure factor S(q) ~ q^−(1+2ζ), ζ vs a, and S(q)
of edge vs screw arcs). The number of rays grows with L (≥ 2 per lattice unit
of wall: 1024 up to L = 256, 8192 at L = 2048). The configurations are git-ignored, so on the cluster
run

```bash
python fss_shapes.py --dir runs/fss_cpu --no-plot   # -> runs/fss_cpu/shapes.npz (tracked, few MB)
```

push `shapes.npz`, and plot anywhere with `python fss_shapes.py --dir runs/fss_cpu`
(`fss_analyze.py` also does it when `shapes.npz` is there).

## GPU / cluster

`gpu_backend.py` is a CUDA version of the relaxation (numba.cuda, one thread
per site and checkerboard colour, fp64). Everything else is shared with the
CPU code.

```bash
python bubble_depinning.py --backend cuda --gpu 0 --L 1024 --R0 170   # single run
python fss_scan.py --backend cuda --gpus 4 --L 512 1024 --samples 512:32 1024:16   # 4 GPUs, one node
sbatch slurm_fss.sh                                                    # job array, 1 GPU per task, R0 = L/6
sbatch slurm_fss.sh 32                                                 # same with fixed R0 = 32
python test_gpu_backend.py                                             # GPU vs CPU check (on the GPU node)
```

Environment on the cluster: `numpy scipy matplotlib numba numba-cuda` and a
CUDA toolkit. Use an fp64-capable GPU (A100/H100/V100). Small systems
(L ≲ 256) are faster on CPU cores; the GPU pays off from L ≈ 512 upward.

With `--shard i --nshards n` each job-array task runs a disjoint slice of the
same task list and writes `results_shard<i>.jsonl`; `fss_analyze.py` merges
all shards. Resubmitting skips finished samples.
