# Depinning with inertia

Working note on an elastic line with mass, driven through quenched disorder by a constant force or a parabolic well:

m ü_i + η u̇_i = κ Δu_i + F_i(u_i) + f [+ k0 (w − u_i)]

Read [inertial_depinning.pdf](inertial_depinning.pdf).

## Main results

The units are κ = η = 1, with random-field disorder A = 1 (Larkin length L_c = 1.65) unless stated. The inertia parameter is M = κm/(ηL_c)².

- **Exact constraints.** The static threshold f_c does not depend on m, and f_↓ ≤ f_↑ ≤ f_c. In steady sliding, the pinning force equals the velocity variance: f − ηv = η⟨δu̇²⟩/v.
- **Radiation friction.** At M ≫ 1 each site radiates sound into the string, with impedance Z = √(κm). The absorbed power is P = ⟨F²⟩/2Z. The flow curve has a minimum at the scale f_* = (2η⟨F²⟩)^{1/2}(κm)^{−1/4}.
- **Hot-line closure (`sc.py`).** Gaussian fluctuations at a temperature fixed by the exact power balance, with no free parameters. It reproduces the sliding branch to three decimals for m ≥ 64. It gives the end of uniform sliding, f_min(M).
- **Long-time threshold = f_min.** Just below f_min a sliding line breaks into domains and stops after τ ∝ (f_min − f)^{−1}, with only weak dependence on L. Just above it, the line slides indefinitely at the closure's v_min. The transition is discontinuous. Finite ramps give f_↓ ≈ 0.85 f_min and a ramp-up hysteresis that shrinks as about L^{−1/2}; both are finite-time and finite-size effects.
- **Large-m law depends on the disorder class.** Debye–Waller smearing by the hot line gives f_min ∝ m^{−(2+a)/2(4+a)}, with a = 1 + 2n for Δ̂(k) ∝ k^{2n}. That is m^{−3/10} for random field and m^{−5/14} for random bond. Here m^{−1/4} is only pre-asymptotic. This is checked with the closure up to m = 65536 and with random-bond runs.
- **Collapse and crossover.** Thresholds for A = 1 and A = 0.5 collapse in M. Inertia moves the threshold only once f_min(M) < f_c (M ≈ 5).
- **Parabolic-well driving.** At large m this gives stick–slip around f_min.

## Contents

- `sim.py`: integrator (numba). Hash-generated disorder, with A < 0 selecting random bond, and an exact damped step.
- `protocols.py`: flow curves, force ramps, quench lifetimes and parabolic-well driving.
- `linresp.py`: second-order friction with the exact lattice susceptibility.
- `sc.py`, `sc_run.py`, `sc_run2.py`: hot-line closure (random field and random bond).
- `runs.py`, `runs2.py`, `runs3.py`: production runs. Output goes to `data/*.npz`, which is git-ignored.
- `analyze.py`, `analyze2.py`: write `statistics*.json` and `fig1`–`fig6`.
- `main.tex`: the note.

## Reproduce

The full set takes about 1.5 hours on 10 cores.

```bash
python3 runs.py kymo well flow thresholds
python3 runs2.py rb heavy life
python3 runs3.py
python3 sc_run.py
python3 sc_run2.py
python3 analyze.py
python3 analyze2.py
pdflatex main.tex && pdflatex main.tex && cp main.pdf inertial_depinning.pdf
```
