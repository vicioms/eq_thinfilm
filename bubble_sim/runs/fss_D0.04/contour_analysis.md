# Critical bubble walls: contour analysis (series D = 0.04, R0 = L/6)

Data: 7460 critical configurations, L = 128 to 4096, anisotropy a = 0 to 1, b along x.
Code: `fss_contours.py` (validation: `python fss_contours.py --selftest`).
Figures: `contours_overhangs.png`, `contours_roughness.png`, `contours_anisotropy.png`.
All numbers: `contour_report.txt`.

## Literature used

- **Kolton, Ferrero, Rosso**, *Depinning free of the elastic approximation*, PRB 108, 174201 (2023), arXiv:2306.13415.
  - Same GL model with a flat wall in a periodic box and uncorrelated random-bond disorder.
  - h_d is the largest field with a metastable state. The wall is the φ = 0 level set.
  - Roughness uses the mean of the crossings (their Eq. 9) and S(q) (Eq. 10). Overhang size is u_o (Eq. 11).
  - The fraction ρ of walls with overhangs collapses as a function of L/ℓ_o, with ℓ_o ~ Δ^(−2.2).
  - Below ℓ_o they find qEW behaviour (ζ ≈ 1.2–1.25). Above ℓ_o they find a "self-similar, invasion-percolation" regime. Projected effective exponents fall to about 0.9, then 0.5, which they say "must be interpreted carefully".
- **Guyonnet, Agoritsas, Paruch, Bustingorry**, arXiv:1904.11726, and **López, Rodríguez**, cond-mat/9706102.
  - For ζ > 1 (super-rough walls) the local width and B(r) saturate at ζ_loc = 1.
  - Only S(q), the global width, or slope statistics measure ζ.
- **Zhou, Zheng**, PRE 82, 031139 (arXiv:1202.1349): the "envelope" (outermost crossing) and "magnetization" interfaces give different exponents when there are overhangs.
- **Domenichini et al.**, arXiv:2312.16802 (bubble domains): the polar parametrization u(θ) only works while the wall is single valued in θ.
- **Geslin et al. letter (SM Eq. 17)**, used to build the noise: a straight segment parallel to b (screw) feels zero integrated force, an edge segment the maximum. I could not find this letter online, so I take the statement from its use in `bubble_depinning.py`.

## Method choices, and why

1. **Analyse the contour, not R(θ).** Every bubble at L ≥ 2048 has overhangs. R(θ) (mean or outermost crossing) then jumps, which adds a q⁻² tail that grows with L. Instead I use the tangent angle ψ(s) of the φ = 0 line along its arc length s, with δψ = ψ − 2πs/P.
   - For a single-valued wall, δψ = du/ds, so S_ψ(q) = q² S_u(q) ~ q^(1−2ζ).
   - ψ stays defined through overhangs, and S_ψ stays valid for ζ > 1.
2. **Fit window: modes n ≥ 12 and q ≤ 0.1.**
   - At each L the lowest ~5 modes are depleted (the drop at the left of each curve).
   - Above q ≈ 0.12, ζ_eff rises steeply because of the noise cut-off exp(−k²w²) with w = 3 (grey band in the plots).
   - The old window, up to q = 1/w = 0.33, sat partly inside that cut-off.
3. **Orientation-resolved observables** use each wall point's local normal on the wall smoothed over 32 sites. Classes: edge = normal within 30° of b; screw = within 30° of b⊥. They do not use polar quadrants, and there is no mean-shape subtraction to get wrong.
4. **Synthetic validation.**
   - Contours of single-valued walls with known ζ (0.5, 1, 1.25) are recovered within the fit error (±0.03 to ±0.12, depending on the window).
   - Total winding is exactly 1 turn, and E = 0 for isotropic walls.
   - The ruler length of a circle is exact, and the support function of an ellipse is correct.

## Results

### 1. The drop of ζ with L was an artifact of the polar projection

| a | contour S_ψ, L = 1024 | 2048 | 4096 | 4096, q_ov < q < 0.1 | 4096, q < q_ov | polar R_mean, 4096 |
|---|---|---|---|---|---|---|
| 0 | 1.18(3) | 1.18(1) | 1.22(1) | 1.26(1) | 1.10(6) | 0.64(2) |
| 0.5 | 1.11(3) | 1.18(2) | 1.17(1) | 1.20(2) | 1.07(6) | 0.85(2) |
| 1 | 1.13(4) | 1.12(2) | 1.15(1) | 1.11(2) | 1.18(8) | 0.95(2) |

- **The contour spectra at different L overlap** wherever q is above each curve's low-mode drop, and the exponent no longer drifts with L.
- **The polar spectra are dominated by overhang jumps.** At L = 4096 they sit 10²–10⁴ above S_ψ/q² for q > 0.01 (bottom-middle panel of the roughness figure).
- **At a = 0, ζ = 1.22(1).** That is close to the qEW value 1.25 and to KFR's below-ℓ_o value, without the ~0.05 fit-window systematic seen in the synthetic tests.
- **ζ decreases slightly with a**, to 1.12–1.15 at a = 1. That is well below the 1.43 of the letter (2 − 1/ν with ν = 1.76).
- **Caveat at a > 0:** the whole-wall S_ψ mixes edge and screw wall. The non-circular mean shape (pointed tips at a = 1) adds a deterministic ~q⁻² piece; I estimate at most ~10% at the lowest fitted q and under 4% at q ≈ 0.06.
- **KFR crossover at ℓ_o: weak evidence at best.**
  - At a = 0 the large-scale window gives 1.10(6) against 1.26(1) at small scales, a 2.6σ difference.
  - The local slope ζ_eff(q) is flat at 1.1–1.25 between q ≈ 0.008 and 0.1, with no clear step at q_ov.
  - The divider dimension rises slowly, to d_f(ε) ≈ 1.14 at ε = 256, with no plateau. It is nowhere near the percolation hull value of 7/4.
  - So at these sizes there is no sign that the walls have entered a self-similar, percolation-like regime.

### 2. Overhangs: a fixed wall length per overhang

- **Wall length per overhang:** ℓ_ov = total wall length / number of polar overhang runs. It settles at **ℓ_ov ≈ 340** for L ≥ 2048 and is the same for every a.
- **It is larger at small L** (876 at L = 512, 3600 at 256) because those walls have travelled only 7–57 sites.
- **Poisson check:** with that ℓ_ov, ρ = 1 − exp(−P/ℓ_ov) at L = 512 (P ≈ 1075, ℓ_ov = 876 there) gives 0.71 against 0.68 measured. This is the bubble analogue of KFR's ℓ_o, at our Δ.
- **Overhang size u_o (KFR Eq. 11)** grows faster than L, roughly as L^1.7. KFR report ∝ L.
- **Overhangs are small and rare on the intrinsic measure.** The fraction of wall turned back by more than 90° relative to the wall smoothed over 32 sites is about 2·10⁻⁴; over 128 sites it is about 1–1.7%.
- **Topology:** every bubble is one domain. Only 2 of 7460 contain a hole, and none has a detached loop. KFR's pinch-off loops do not occur here.
- **At a > 0, intrinsic overhangs sit on screw wall** (a = 0.75, L = 4096: 1.3·10⁻³ on screw wall against ~0 on edge wall).

### 3. Anisotropy: edge wall stays pinned, screw tips run

- **Wall orientation.** E = ⟨cos 2φ_n⟩ is +1 if all wall has its normal along b (edge) and −1 if perpendicular (screw).
  - It grows with a and keeps growing with L. At a = 1, E goes 0.21 → 0.56 from L = 128 to 4096, with no saturation.
  - At L = 4096 and a = 1, 74% of the wall is edge and 4% is screw.
  - At a = 0 the length distribution is flat within 5%.
- **Advance by direction.** The support function about the nucleation centre gives how far the most advanced wall point got along and across b. At L = 4096, a = 1: **h(⊥b) − R0 = 1080, h(∥b) − R0 = 47**, a ratio of 23 that keeps growing with L (3.5 at L = 128).
- **Interpretation.** Edge wall barely moves, while screw-oriented wall, which feels no integrated force (letter SM Eq. 17), runs. This is a kinetic-Wulff picture: the fast orientations grow out and leave pointed tips, while the slow, strongly pinned edge orientation takes over the wall.
- **Slope fluctuations below 4–16 sites.** var(ψ − ψ_σ) on screw wall over edge wall:

  | a | 0 | 0.25 | 0.5 | 0.75 | 1 |
  |---|---|---|---|---|---|
  | ratio | 1.00 ± 0.03 | 1.15 | 1.40 | 1.66 | 2.02 |

  - The ratio is the same for σ = 4, 8 and 16 and for L = 2048 and 4096.
  - Edge wall gets smoother as a grows (0.086 → 0.059 rad²); screw wall gets rougher (0.088 → 0.120).
  - At σ = 32 the ratio climbs to 2.35 because the tip curvature (radius ~45 sites) leaks in.
  - **This replaces the earlier "screw ~100× rougher" from polar quadrants, which was a shape and overhang artifact.**

### 4. What this means for h_d and the finite-size scaling

- **Memory of the start is limited to the lowest modes.** The wall travels R_eff − R0 ≈ 0.12 L (L ≥ 512). The bubble-scale roughness is a roughly constant 0.3–0.5 of that distance. In S_ψ only the lowest ~5 modes are depleted, so above n ≈ 12 the configurations are L-independent.
- **The wall comes close to the box edge.** At depinning, the farthest wall point is at 0.87–0.94 of the half-box (L ≥ 1024), so the last stretch before the border is only ~6–13% of L/2.
- **At large a, h_d is the threshold of screw-oriented wall.** It is set by tips running along b⊥ while the edge wall stays almost where it started. The decrease of h_c with a in `fss_summary.png` therefore reflects the weakest orientation, not a pinning strength averaged over the bubble.
  - Since E(L) and the advance ratio have not saturated by L = 4096, the large-a threshold may still drift.
  - A square box also lets an elongated bubble reach the edge earlier along its long axis.
- **Open (2σ):** at a = 0, h(0°) and h(90°) differ by 11% at L = 2048–4096 (a 2σ effect, while E = −0.006 ± 0.01). I read it as noise in an extreme-value statistic.

## Not done / limits

- **Fixed-R0 series:** the R0 = 32 series has no contours yet.
- **Error bars at L = 4096:** they come from 32 samples per a; L = 2048 has 64.
- **ℓ_ov at another Δ:** it is measured only at Δ = 0.04, so the Δ^(−2.2) scaling of ℓ_o can't be tested.
- **Edge-only and screw-only exponents:** ζ for each wall class isn't measured. At a = 1, screw wall is only 4% of the length (tips), too little for a spectrum.

## a = 0: is it quenched Edwards-Wilkinson (qEW)?

Script: `qew_a0.py`. Figure: `qew_a0.png`. Numbers: `qew_a0.txt`. Errors come from bootstrapping over disorder samples.

**Targets.** qEW in d = 1 at depinning has ζ = 1.25 and ν = 1/(2 − ζ) = 4/3 (statistical tilt symmetry). KFR's large-scale invasion-percolation-like regime has the same ν = 4/3 but ζ = 1 (self-similar). So ζ, not ν, tells the two apart.

**Noise caveat.** Even at a = 0 the noise is not KFR's white noise. Its spectrum is ∝ |k| e^(−k²w²), which vanishes at k = 0.
- Along the wall, the force correlations are short-range: the integrated force on a straight segment grows like √ℓ.
- Across the wall, the force correlator Δ(u) integrates to zero (random-bond-like) and has a u⁻² tail.
- Random-bond and random-field disorder are expected to share the qEW fixed point at depinning, but the u⁻² tail makes this a check, not a given.

**ζ.**
- **Compensated plot:** q^1.5 S_ψ(q) is flat from q ≈ 0.02 to 0.13 for L = 2048 and 4096, nearly a decade. Its upper end is the noise core w = 3; its lower end is q_ov = 2π/ℓ_ov.
  - L = 4096: ζ = 1.263 ± 0.015 for q_ov < q < 0.1, and 1.260 ± 0.026 for 0.02 < q < 0.07.
  - L = 2048: ζ = 1.217 ± 0.016 for q_ov < q < 0.1.
- **Fit-window scan:** every window with q_lo ≥ 0.015 and q_hi = 0.07–0.14 gives 1.21–1.30.
- **Combined:** ζ = 1.24 ± 0.02 (stat) ± 0.03 (window), which agrees with qEW.
- **Fitting below q_ov lowers the exponent.** The whole n ≥ 12 window gives 1.18–1.22. With an explicit cut-off model, A q^(1−2ζ) exp(−(q/q_c)²) up to q = 0.3, it gives 1.16–1.17, because the window then includes the large-scale regime.

**Above ℓ_ov ≈ 340 there is a crossover.**
- For q < q_ov the compensated spectrum falls. A fit from n ≥ 6 up to q_ov gives ζ_eff = 1.07 at L = 4096 and 1.09 at L = 2048.
- The two sizes agree on the same q range, so this is not a finite-size effect.
- This is KFR's scenario: qEW below the overhang length, and a lower effective exponent above it.
- At L = 4096 the regime above ℓ_ov spans only a factor ~3 in scale. That is too short to say whether ζ → 1 (invasion percolation).
- The divider dimension stays ≈ 1.1, so there is no sign yet of a fractal, percolation-like hull.

**ν.**
- var(h_d) ~ L^(−2/ν) gives ν = 1.31 ± 0.09 for L ≥ 512, consistent with 4/3.
- The estimate depends on the fit range: 1.19 ± 0.05 for L ≥ 256 (small L still remembers the start) and 1.19 ± 0.13 for L ≥ 1024.
- Pairs of successive L give 1.65 ± 0.42, 1.33 ± 0.40 and 1.07 ± 0.31, so the variance alone can't pin ν better than ±0.1.

**Scaling relation.** 1/(2 − ζ) = 1.28–1.36 for ζ = 1.22–1.26, which matches the measured ν = 1.31 ± 0.09.

**Shape of the h_d distribution.**
- Skewness falls from +0.76 ± 0.13 (L = 128) and +0.45 ± 0.16 (L = 256) to values consistent with 0 for L ≥ 512 (+0.23 ± 0.17, −0.14 ± 0.24, −0.14 ± 0.44, −0.22 ± 0.23).
- The rescaled CDFs collapse for L ≥ 512.
- That is consistent with one fixed point once the wall has left the start. With 32–64 samples at the largest L this is only a consistency check.

**Verdict for a = 0.**
- Between the noise core (~10 sites) and the overhang length ℓ_ov ≈ 340, the walls are qEW: ζ = 1.24 ± 0.04, ν = 1.31 ± 0.09, and the two satisfy ν = 1/(2 − ζ).
- At larger scales the effective roughness exponent drops to ≈ 1.08, as KFR describe. The asymptotic large-scale class is not reached at L ≤ 4096.

**Next, for a > 0.** Testing whether wall ∥ b and wall ⊥ b fall in different classes needs walls of each orientation over long lengths. In the bubbles, screw wall is only 4% of the length at a = 1 and sits at the tips. A flat-wall geometry with the wall along or across b (KFR-like strip) is the clean test.
