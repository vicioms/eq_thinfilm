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
The analysis is in the [Analysis](#analysis) section below.

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
same task list and writes `results_shard<i>.jsonl`; the analysis merges
all shards. Resubmitting skips finished samples.

## Flat walls (wall along or across b)

`--geometry flat` replaces the bubble by a straight wall along x at height
0.1 Ly (+1 below, −1 above), periodic in x and open in y (bottom rows pinned to
+1). The box is L along the wall and Ly = aspect · L along the motion
(`--aspect`, default 2), so the wall travels well beyond its own width (ζ > 1)
before depinning. The box never rotates: `--burgers-angle` rotates b in the
noise spectrum only. It is pushed up and depins when it reaches the top band, as in the
Kolton–Ferrero–Rosso strip. `--burgers-angle` sets the wall direction: 0 is a
screw wall (line ∥ b), 90 an edge wall (line ⊥ b, normal ∥ b, the hard direction).
This is the clean test of qEW (ζ = 1.25, ν = 4/3) vs qKPZ (ζ = 0.63, ν = 1.73)
for each orientation, which the bubbles cannot do.

```bash
python bubble_depinning.py --geometry flat --burgers-angle 90 --L 256 --Delta 0.04 --coupling field
sbatch slurm_fss_flat.sh 0          # screw walls  -> runs/flat_D0.04_b0_A2
sbatch slurm_fss_flat.sh 90         # edge walls   -> runs/flat_D0.04_b90_A2
sbatch slurm_fss_flat.sh 90 4       # same with a L x 4L box -> runs/flat_D0.04_b90_A4
```

## Analysis

One script per geometry, the same observables and figures for both; shared code
in `analysis/` (`data` extraction and loading, `geometry` curves, `measure` one
wall, `spectra` pooled spectra and ζ fits, `stats` thresholds and power laws,
`figures` plots and report). The critical configurations stay on the cluster
(git-ignored); the walls are extracted there into small tracked files:

```bash
python analyze_flat.py --extract runs/flat_D0.04_b0_A2      # -> walls.npz     (one period of each wall)
python analyze_bubbles.py --extract runs/fss_D0.04          # -> contours.npz  (all φ = 0 contours)
```

push them, and analyse anywhere (one row per run in every figure):

```bash
python analyze_flat.py runs/flat_D0.04_b0_A2 runs/flat_D0.04_b90_A2   # -> runs/flat_D0.04_A2_analysis/
python analyze_bubbles.py runs/fss_D0.04                              # -> runs/fss_D0.04_analysis/
```

| output | content |
|---|---|
| `configs.png` | example critical walls at the largest L with b; distribution of the angle between the wall line and b (0 screw-like, 90 edge-like: faceting) for every L; bubbles: mean shape ⟨R(θ)⟩/R_eff |
| `roughness.png` | tangent-angle spectrum S_ψ(q) ~ q^(1−2ζ) along the arc length (valid with overhangs and for ζ > 1), pooled over the walls of each (L, a); local ζ_eff(q); ζ vs L, intrinsic and from the single-valued profile (u(x) per column, R(θ) per ray: KFR Eq. 9) |
| `threshold.png` | ⟨h_d⟩(L) = h_c + A L^(−1/ν) (+ B/L for bubbles, Laplace pressure) with ν = 4/3 (qEW), 1.733 (qKPZ) or free; var(h_d) ~ L^(−2/ν); h_c vs a |
| `length.png` | wall length P ∝ size^α (α = 1 for ζ < 1, α ≈ ζ for ζ > 1); P/L or P/2πR_eff (as `NB7_perimeter.ipynb`); fraction of columns / rays with overhangs |
| `report.txt` | all numbers, with bootstrap errors over disorder samples, and ζ split at the median distance travelled before the last pinning (memory check) |

Spectra: the modes of all walls of one (L, a) go into one array, binned in log q
with bins merged until each holds at least one mode per wall; fits use modes
n ≥ 2 (flat) or n ≥ 12 (bubbles: lower modes are the bubble's shape) up to
`--qmax` (default 0.1, the noise core). With a fixed R0 (`fss_scan.py --r0`)
the critical bubble stays near R0 at every L, so nothing there scales with L.

The previous scripts (`fss_analyze.py`, `fss_shapes.py`, `fss_contours.py`,
`qew_a0.py`, `kpz_edge.py`, `sides.py`, `fss_flat.py`, `perimeter.py`) are in
`legacy/` with a table of what is and is not ported (`legacy/README.md`).
