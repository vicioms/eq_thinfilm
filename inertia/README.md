# Depinning with inertia

Working note on an elastic line with mass, driven through quenched disorder by a constant force or a parabolic well:

m ü_i + η u̇_i = κ Δu_i + F_i(u_i) + f [+ k0 (w − u_i)]

Read [inertial_depinning.pdf](inertial_depinning.pdf).

## Main results

The units are κ = η = 1, with random-field disorder A = 1 (Larkin length L_c = 1.65) unless stated. The inertia parameter is M = κm/(ηL_c)².

- **Exact constraints.** The static threshold f_c does not depend on m, and f_↓ ≤ f_↑ ≤ f_c. In steady sliding, the pinning force equals the velocity variance: f − ηv = η⟨δu̇²⟩/v.
- **Radiation friction.** At M ≫ 1 each site radiates sound into the string, with impedance Z = √(κm). The absorbed power is P = ⟨F²⟩/2Z. The flow curve has a minimum at the scale f_* = (2η⟨F²⟩)^{1/2}(κm)^{−1/4}.
- **Hot-line closure (`sc.py`).** Gaussian fluctuations at a temperature fixed by the exact power balance, with no free parameters. It reproduces the sliding branch to three decimals for m ≥ 64. It gives the end of uniform sliding, f_min(M).
- **Closure tested directly.** Temperature and pinning force match to 1–3%. The fluctuations are near-Gaussian but not thermal in shape, and the closure's numbers carry a systematic uncertainty of about 10%.
- **Long-time threshold ≈ 0.975–0.99 f_min** at every mass tested (m = 64 to 16384); 0.986 f_min at m = 64 from runs of 12000m at L = 65536. The edge between hot and pinned line advances above f_0 ≈ 0.98 f_min and retreats below it (m = 256 and 1024), saturating near the sound speed. Below the threshold, sliding lines stop after τ ∝ (f_0 − f)^{−1}. Between the threshold and f_min the line slides in patches. This patchy sliding ends continuously, like an absorbing-state transition: the active fraction decays as a slow power law just above f_∞ and bends down just below it, with no plateau. Lifetimes reach about 10^4 m within 1% of f_∞, so shorter runs place the threshold too low. Inertially stopped lines are more stable than overdamped-relaxed ones (`marginality.py`). Finite ramps give f_↓ ≈ 0.85 f_min and a ramp-up hysteresis that shrinks as about L^{−1/2}; both are finite-time and finite-size effects.
- **Roughness of the hot line: ζ = 1/2 (Edwards–Wilkinson)**, with a temperature that depends on scale. The large-scale value T_lw = S_ξ(0)/2η_eff is 1.8–4.3 times T_kin and is predicted within about 10%. S_ξ(0) is the measured zero-frequency single-site force noise, and η_eff = η + ∂F_pin/∂v at fixed T. Random bond suppresses large-scale roughness fifty-fold (Δ̂(0) = 0). Bare second order underestimates T_lw by a factor of about 3. For small M, FRG gives qEW exponents (mass irrelevant with eigenvalue −z); at the end of patchy sliding, τ cannot distinguish directed percolation from qEW.
- **Avalanches under slow well driving** (k0 = 3e-3, L = 512): up to M ≈ 0.4 they are the same as at m = 0, event for event. For M ≳ 6 the distribution is bimodal: cold avalanches, plus system-spanning slips that carry 50–91% of the advance. The slips behave like a spring block with kinetic friction ≈ f_min, and start from an overload (0.96, 1.17 and 1.56 f_min at m = 64, 256 and 1024). τ is not resolved at ξ ≈ 18.
- **Avalanche exponents** (soft well, k0 = 5e-4, ξ ≈ 45, L = 2048, plus the moments between the two wells): at m = 0 and 1, τ = 1.1 ± 0.05 (three estimators avoiding the cutoff bump: 1.06, 1.09, 1.17), consistent with qEW 1.11, and s ~ T^1.6–1.7. m = 1 reproduces 86% of the m = 0 events, m = 2 only 34%: the overshoot onset is at M = 0.37–0.74. Above it, the inertial force shift grows from 0.6 to 3.3 times the critical-window width at m = 4 (from 3.7 to 20 at m = 16) as ξ goes from 18 to 45, and s ~ T^2. The crossover to stick–slip moves to smaller M as ξ grows.
- **Large-m law depends on the disorder class.** Debye–Waller smearing by the hot line gives f_min ∝ m^{−(2+a)/2(4+a)}, with a = 1 + 2n for Δ̂(k) ∝ k^{2n}. That is m^{−3/10} for random field and m^{−5/14} for random bond. Here m^{−1/4} is only pre-asymptotic. This is checked with the closure up to m = 65536 and with random-bond runs.
- **Collapse and crossover.** Thresholds for A = 1 and A = 0.5 collapse in M. The shift switches on at M ≈ 0.5, where the Larkin-scale mode becomes underdamped, and joins the hot branch smoothly. The onset is earlier for random bond, whose wells are stiffer, so it is not universal.
- **Parabolic-well driving.** At large m this gives stick–slip around f_min.

## Contents

- `sim.py`: integrator (numba). Hash-generated disorder, with A < 0 selecting random bond, and an exact damped step.
- `protocols.py`: flow curves, force ramps, quench lifetimes and parabolic-well driving.
- `linresp.py`: second-order friction with the exact lattice susceptibility.
- `sc.py`, `sc_run.py`, `sc_run2.py`: hot-line closure (random field and random bond).
- `runs.py`, `runs2.py`, `runs3.py`, `runs4.py`: production runs. Output goes to `data/*.npz`, which is git-ignored.
- `edges.py`: hot/pinned edge speed. `patchy.py`: end of sliding at m = 64 (2000m runs). `marginality.py`: stability of inertially stopped versus relaxed lines. `run_gap64.py`: long lifetimes at m = 64, L = 32768. `endslide.py`, `endslide_L.py`, `endslide_t.py`: end of sliding at m = 64 (fine grid, size dependence, activity decay ρ(t) at L = 65536; `lifetimes(..., trace=True)` records ρ per chunk). `check_closure.py`: measured fluctuations of sliding lines against the closure.
- Exploratory, with negative or mixed results: `sc_hybrid.py` (driven shape, power-balance amplitude; fits the patchy low-v branch but gives no sharp threshold), see `notes_hybrid_closure.md`; `sc_noneq.py` (closure with the driven, non-equilibrium noise spectrum; not energy-consistent at low v), see `notes_noneq_closure.md`; `front_model.py` (deterministic front with the closure friction; gets the saturation at c but not the zero below f_min), see `notes_front_model.md`. `front_model_stoch.py` (blocks of size √m with local temperature and noise; reproduces f_min, but gives a jump above f_min instead of the continuous end below it, so quenched heterogeneity is needed), see `notes_front_model_stoch.md`.
- `sq.py`, `sq_near.py`, `sq_diag.py`: structure factor of the hot sliding line, its second-order prediction, and the measured single-site force correlator C(t) and displacement spread B(t).
- `aval.py`, `aval_merge.py`: avalanches under slow well driving for m = 0–1024 (events closer than 4m in time are merged).
- `tau_fit.py`: maximum-likelihood fit of s^-τ e^(-s/s_c) with bootstrap errors.
- `analyze.py` … `analyze7.py`: write `statistics*.json` and `fig1`–`fig13`.
- `main.tex`: the note.

## Reproduce

The full set takes about 3 hours on 10 cores.

```bash
python3 runs.py kymo well flow thresholds
python3 runs2.py rb heavy life
python3 runs3.py
python3 runs4.py cross heavy gap
python3 edges.py
python3 check_closure.py
python3 patchy.py && python3 endslide.py && python3 endslide_L.py
python3 endslide_t.py && python3 endslide_t.py 0.1975 0.1985
python3 sc_run.py
python3 sc_run2.py
python3 analyze.py
python3 analyze2.py
python3 analyze3.py
python3 sq.py run && python3 sq.py theory && python3 sq_near.py 0.10 0.105 0.11
python3 sq_diag.py 1024 0.12 && python3 sq_diag.py 1024 0.10 && python3 sq_diag.py 1024 0.11 && python3 sq_diag.py 256 0.17
python3 analyze4.py
python3 aval.py scan
python3 aval.py tau
python3 analyze5.py
python3 analyze6.py
python3 analyze7.py
pdflatex main.tex && pdflatex main.tex && cp main.pdf inertial_depinning.pdf
```
