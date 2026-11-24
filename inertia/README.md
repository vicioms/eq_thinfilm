# Depinning with inertia

Working note on an elastic line with mass, driven through quenched disorder by a constant force or a parabolic well:

m ü_i + η u̇_i = κ Δu_i + F_i(u_i) + f [+ k0 (w − u_i)]

Read [notes/inertial_depinning.pdf](notes/inertial_depinning.pdf). The paper draft is in `paper/` (a separate ShareLaTeX git clone, ignored here).

## Main results

The units are κ = η = 1, with random-field disorder A = 1 (Larkin length L_c = 1.65) unless stated. The inertia parameter is M = κm/(ηL_c)².

- **Exact constraints.** The static threshold f_c does not depend on m, and f_↓ ≤ f_↑ ≤ f_c. In steady sliding, the pinning force equals the velocity variance: f − ηv = η⟨δu̇²⟩/v.
- **Radiation friction.** At M ≫ 1 each site radiates sound into the string, with impedance Z = √(κm). The absorbed power is P = ⟨F²⟩/2Z. The flow curve has a minimum at the scale f_* = (2η⟨F²⟩)^{1/2}(κm)^{−1/4}.
- **Hot-line closure (`src/sc.py`).** Gaussian fluctuations at a temperature fixed by the exact power balance, with no free parameters. It reproduces the sliding branch to three decimals for m ≥ 64. It gives the end of uniform sliding, f_min(M).
- **Closure tested directly.** Temperature and pinning force match to 1–3%. The fluctuations are near-Gaussian but not thermal in shape, and the closure's numbers carry a systematic uncertainty of about 10%.
- **Long-time threshold ≈ 0.975–0.99 f_min** at every mass tested (m = 64 to 16384); 0.986 f_min at m = 64 from runs of 12000m at L = 65536. The edge between hot and pinned line advances above f_0 ≈ 0.98 f_min and retreats below it (m = 256 and 1024), saturating near the sound speed. Below the threshold, sliding lines stop after τ ∝ (f_0 − f)^{−1}. Between the threshold and f_min the line slides in patches. This patchy sliding ends continuously, like an absorbing-state transition: the active fraction decays as a slow power law just above f_∞ and bends down just below it, with no plateau. Lifetimes reach about 10^4 m within 1% of f_∞, so shorter runs place the threshold too low. Inertially stopped lines are more stable than overdamped-relaxed ones (`marginality.py`). Finite ramps give f_↓ ≈ 0.85 f_min and a ramp-up hysteresis that shrinks as about L^{−1/2}; both are finite-time and finite-size effects.
- **Single-site (DMFT) theory of the hot line** (`src/cavity.py`, `src/gdmft2.py`, `runs/dmft_flow.py`, `analysis/single_site.py`, `analysis/dmft_flow.py`, figures/single_site, dmft_flow). The disorder vertex of the MSRJD action is local, so one site is solved exactly in a self-consistent Gaussian bath (the rest of the line as a linear chain). Two solutions: a numerical cavity (site 0 in its real landscape inside a ring of linear sites) and a closed, deterministic Gaussian DMFT: C = disorder correlator smeared by the true spread B; B from S_ξ through the chain propagator; and the one-loop disorder self-energy Σ(ω) = ∫ R⟨Δ″⟩(1 − e^{iωt}) dt in the medium. With Σ(ω) the theory conserves energy (power balance to < 1%; a constant-η_eff version, `gdmft.py`, ran 18–32% too hot). With no input from simulation it gives the flow curve within 1%, T_kin within 3% and S_ξ(0) (1.6–2.1× bare) within 6%, and reproduces C(t) and B(t). Physical picture: S_ξ(0) = Δ̂(0) × the local time of a site on its own trajectory; for B ∝ t^{2H}, S_ξ(0) = Δ̂(0)/[(1−H)v]. Its uniform hot branch stays stable below f_∞ (shallow minimum 0.087 at m = 1024; none down to v = 0.01 at m = 64, 256), so the end of sliding is a nucleation of pinned patches, not a spinodal of uniform sliding.
- **Roughness of the hot line: ζ = 1/2 (Edwards–Wilkinson)**, with a temperature that depends on scale. The large-scale value T_lw = S_ξ(0)/2η_eff is 1.8–4.3 times T_kin and is predicted within about 10%. S_ξ(0) is the measured zero-frequency single-site force noise, and η_eff = η + ∂F_pin/∂v at fixed T. Random bond suppresses large-scale roughness fifty-fold (Δ̂(0) = 0). Bare second order underestimates T_lw by a factor of about 3. For small M, FRG gives qEW exponents (mass irrelevant with eigenvalue −z); at the end of patchy sliding, τ cannot distinguish directed percolation from qEW.
- **Avalanches under slow well driving** (k0 = 3e-3, L = 512): up to M ≈ 0.4 they are the same as at m = 0, event for event. For M ≳ 6 the distribution is bimodal: cold avalanches, plus slips, events larger than any avalanche of the overdamped line under the same driving (s > 8L), that carry 15–88% of the advance (52–91% with the older cut s > 4L). The slips behave like a spring block with kinetic friction 0.90–1.10 f_min, insensitive to the slip definition, and start from an overload (0.99, 1.20 and 1.59 f_min at m = 64, 256 and 1024). τ is not resolved at ξ ≈ 18.
- **No quenched rare regions at the end of sliding** (`runs/rare.py`, `analysis/rare_regions.py`, figures/rare_regions): long-lived activity below f_∞ is a wandering branching cluster that does not come back when the same disorder is rerun from a displaced line, and the decay below f_∞ has no constant-exponent stretch. Infinite randomness and Griffiths behaviour are ruled out; the disorder is effectively annealed for the activity (it is random in u), leaving clean directed percolation as the candidate.
- **Spreading from a hot seed** (`runs/spread.py`, `analysis/spreading.py`, figures/spreading): on the natural absorbing state (lines that died at 0.974 f_min), with 96 runs per force at L = 16384 over 2000m. The seeds die at 0.984 f_min and invade at 0.989 f_min, so f_∞ = (0.987 ± 0.002) f_min. At 0.198: δ ≈ 0.21 ± 0.05 and θ ≈ 0.30–0.38, consistent with DP (0.159, 0.313). The large run (L = 32768, 480 clean runs at each of 0.1981 and 0.1983, 5000m; `analysis/spreading.py`) gives 2/z = 1.50–1.59 ± 0.03 in every window, so z ≈ 1.31 ± 0.03 (statistical): superdiffusive and not DP (1.58). Changing the seed (16 or 256 sites) or the substrate (quenched at 0.95 f_min) moves 2/z between 1.34 and 1.57 (`runs/thresholds.py`), so conserved DP (1.43) is not excluded; δ ≈ 0.20–0.23 and θ ≈ 0.36–0.44 sit 1.5–3σ above DP, and generalized hyperscaling θ + 2δ ≈ 1/z holds. Likely cause: activity carried by hot edges and sound at up to the sound speed (long-range transport).
- **Avalanche exponents** (soft well, k0 = 5e-4, ξ ≈ 45, L = 2048, plus the moments between the two wells): at m = 0 and 1, τ = 1.1 ± 0.05 (three estimators avoiding the cutoff bump: 1.06, 1.09, 1.17), consistent with qEW 1.11, and s ~ T^1.6–1.7. m = 1 reproduces 86% of the m = 0 events, m = 2 only 34%: the overshoot onset is at M = 0.37–0.74. Above it, the inertial force shift grows from 0.6 to 3.3 times the critical-window width at m = 4 (from 3.7 to 20 at m = 16) as ξ goes from 18 to 45, and s ~ T^2. The crossover to stick–slip moves to smaller M as ξ grows.
- **Large-m law depends on the disorder class.** Debye–Waller smearing by the hot line gives f_min ∝ m^{−(2+a)/2(4+a)}, with a = 1 + 2n for Δ̂(k) ∝ k^{2n}. That is m^{−3/10} for random field and m^{−5/14} for random bond. Here m^{−1/4} is only pre-asymptotic. This is checked with the closure up to m = 65536 and with random-bond runs.
- **Collapse and crossover.** Thresholds for A = 1 and A = 0.5 collapse in M. The shift switches on at M ≈ 0.5, where the Larkin-scale mode becomes underdamped, and joins the hot branch smoothly. The onset is earlier for random bond, whose wells are stiffer, so it is not universal.
- **Parabolic-well driving.** At large m this gives stick–slip around f_min.

## Layout

Every script can be launched from any directory: it imports `src/config.py`, which puts `src/` on the path and moves to this folder.

- `src/`: library code.
  - `config.py`: paths and every analysis threshold (active site, rest, uniform vs patchy, avalanche detection, merging, slips, ramp and lifetime stops), each with the measurement that fixes it; `analysis/thresholds.py` makes the supporting figure.
  - `sim.py`: integrator (numba). Hash-generated disorder, A < 0 selects random bond; exact damped step.
  - `protocols.py`: flow curves, force ramps, quench lifetimes, parabolic-well driving.
  - `linresp.py`: second-order friction with the exact lattice susceptibility.
  - `sc.py`: hot-line (Gaussian) closure, random field and random bond.
  - `gdmft2.py`: Gaussian single-site theory with the self-energy Σ(ω); `gdmft.py`: its smearing routines and the constant-η_eff version.
  - `cavity.py`: numerical single-site cavity.
  - `aval_merge.py`, `tau_fit.py`: avalanche merging and maximum-likelihood fits of s^-τ e^(-s/s_c).
- `runs/`: production runs, writing `data/*.npz` (git-ignored).
  - `sliding.py` (modes `flow thresholds well kymo`), `lifetimes_rb.py` (`life heavy rb`), `lifetimes_fmin.py`, `crossover_gap.py` (`cross heavy gap`), `gap64.py`.
  - `closure_flow.py`, `closure_flow_rb.py`, `check_closure.py`: closure flow curves and its direct test.
  - `edges.py` (moving/pinned edge speed), `patchy.py`, `endslide.py`, `endslide_L.py`, `endslide_t.py` (end of sliding at m = 64), `marginality.py`, `rare.py`, `spread.py` (spreading; `SPREAD_L`, `SPREAD_NSUB`, `SPREAD_NK` select the large run).
  - `sq.py`, `sq_near.py`, `sq_diag.py`: structure factor, C(t) and B(t) of the sliding line.
  - `aval.py`: avalanches under slow well driving (`pilot scan dtcheck tau`).
  - `thresholds.py`: short checks behind the analysis thresholds (speed distributions, spreading at several activity thresholds, seeds and substrates, avalanche detection).
  - `cavity_runs.py`, `cavity_polish.py`, `dmft_flow.py`: single-site theory at the measured points and along the flow curves.
- `analysis/`: one script per topic, writing `figures/*.pdf|png` and `data/stats/*.json`.

  | script | figures | stats |
  |---|---|---|
  | `sliding.py` | flow_curves, coexistence, ramp_thresholds, well_driving | sliding.json |
  | `lifetimes_rb.py` | lifetimes, large_mass_rb | lifetimes_rb.json |
  | `closure_crossover.py` | closure_edges, crossover, long_runs | closure_crossover.json |
  | `end_of_sliding.py` | end_of_sliding | |
  | `roughness.py` | roughness | roughness.json |
  | `single_site.py` | single_site | single_site.json |
  | `dmft_flow.py` | dmft_flow | |
  | `avalanches.py` | avalanches | avalanches.json |
  | `tau.py` | tau | tau.json |
  | `rare_regions.py` | rare_regions | |
  | `spreading.py` | spreading | spreading.json |
  | `thresholds.py` | thresholds | thresholds.json |

- `data/`: simulation output (`*.npz`, git-ignored), `data/logs/` (run logs), `data/stats/` (numbers quoted in the note).
- `figures/`: all figures, named by content.
- `notes/`: the working note `main.tex` (built as `inertial_depinning.pdf`) and the `notes_*.md` on exploratory models.
- `archive/`: superseded or exploratory code with negative or mixed results: front models, `sc_hybrid.py`, `sc_noneq.py`, the constant-η_eff flow scan `gdmft_flow.py`, and the first spreading analysis `analyze9.py` (96-run set, superseded by `analysis/spreading.py`).
- `paper/`: ShareLaTeX clone of the paper (ignored by this repository).

## Reproduce

The commands for each section of the note are in its "Reproduce" table. The full set takes several hours; cap the load with `NUMBA_NUM_THREADS`. After the runs, regenerate every figure and statistics file with

```bash
for f in analysis/*.py; do python3 "$f"; done
cd notes && pdflatex main.tex && pdflatex main.tex && cp main.pdf inertial_depinning.pdf
```
