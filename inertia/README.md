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
- **Long-time threshold ≈ 0.97–0.98 f_min** at every mass tested (m = 64 to 16384). The edge between hot and pinned line advances above f_0 ≈ 0.98 f_min and retreats below it (m = 256 and 1024), saturating near the sound speed. Below the threshold, sliding lines stop after τ ∝ (f_0 − f)^{−1}. Between the threshold and f_min the line slides in patches, and v falls from ≈ v_min to ≈ 0.2 v_min before sliding dies. Inertially stopped lines are more stable than overdamped-relaxed ones (`marginality.py`). Finite ramps give f_↓ ≈ 0.85 f_min and a ramp-up hysteresis that shrinks as about L^{−1/2}; both are finite-time and finite-size effects.
- **Large-m law depends on the disorder class.** Debye–Waller smearing by the hot line gives f_min ∝ m^{−(2+a)/2(4+a)}, with a = 1 + 2n for Δ̂(k) ∝ k^{2n}. That is m^{−3/10} for random field and m^{−5/14} for random bond. Here m^{−1/4} is only pre-asymptotic. This is checked with the closure up to m = 65536 and with random-bond runs.
- **Collapse and crossover.** Thresholds for A = 1 and A = 0.5 collapse in M. The shift switches on at M ≈ 0.5, where the Larkin-scale mode becomes underdamped, and joins the hot branch smoothly. The onset is earlier for random bond, whose wells are stiffer, so it is not universal.
- **Parabolic-well driving.** At large m this gives stick–slip around f_min.

## Contents

- `sim.py`: integrator (numba). Hash-generated disorder, with A < 0 selecting random bond, and an exact damped step.
- `protocols.py`: flow curves, force ramps, quench lifetimes and parabolic-well driving.
- `linresp.py`: second-order friction with the exact lattice susceptibility.
- `sc.py`, `sc_run.py`, `sc_run2.py`: hot-line closure (random field and random bond).
- `runs.py`, `runs2.py`, `runs3.py`, `runs4.py`: production runs. Output goes to `data/*.npz`, which is git-ignored.
- `edges.py`: hot/pinned edge speed. `patchy.py`: end of sliding at m = 64 (2000m runs). `marginality.py`: stability of inertially stopped versus relaxed lines. `run_gap64.py`: long lifetimes at m = 64, L = 32768. `check_closure.py`: measured fluctuations of sliding lines against the closure.
- Exploratory, with negative or mixed results: `sc_hybrid.py` (driven shape, power-balance amplitude; fits the patchy low-v branch but gives no sharp threshold), see `notes_hybrid_closure.md`; `sc_noneq.py` (closure with the driven, non-equilibrium noise spectrum; not energy-consistent at low v), see `notes_noneq_closure.md`; `front_model.py` (deterministic front with the closure friction; gets the saturation at c but not the zero below f_min), see `notes_front_model.md`.
- `analyze.py`, `analyze2.py`, `analyze3.py`: write `statistics*.json` and `fig1`–`fig9`.
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
python3 sc_run.py
python3 sc_run2.py
python3 analyze.py
python3 analyze2.py
python3 analyze3.py
pdflatex main.tex && pdflatex main.tex && cp main.pdf inertial_depinning.pdf
```
