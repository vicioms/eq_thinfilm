# Bubble theory

Analytic theory of bubble depinning in the disordered 2D φ⁴ (Ginzburg–Landau) model simulated in [`bubble_sim/`](../bubble_sim):

∂ₜφ = c∇²φ + ε₀[(1+r)φ − φ³] + h + g

The programme is to go from the field theory of φ to an interface theory for the bubble wall, and then run the functional RG on that interface.

## Contents

- [`notes_phi4_interface.md`](notes_phi4_interface.md): MSR action for φ⁴ with random-field, random-bond and dislocation disorder; kink and bubble saddle points and their fluctuation spectrum; collective coordinates R, X and u by Faddeev–Popov; tree-level projection onto the interface action, with the disorder correlators Δ(u); collective pinning of a rigid bubble; what the gapped modes add.
- [`notes_frg.md`](notes_frg.md): functional RG for the wall with the dislocation correlator. Exact amplitude scaling of the threshold and why it fails here (the force is isotropic below ~2w, and L_c ≈ 2w); exact a → −a symmetries; the threshold anisotropy as a finite-size effect, C(L) ∝ L^-0.6; orientation smearing excluded by the wall tangents; screw walls pinned by their own roughness; scale-dependent KPZ coefficient and crossover length; predictions to test.
- `frg_anisotropy.py`: every number in notes_frg.md, from the noise spectrum and the flat runs in `bubble_sim/runs` (thresholds from `results*.jsonl`, wall tangents from `walls.npz`). Reports: `frg_anisotropy_A2.txt`, `frg_anisotropy_A4.txt`.
- `check_projection.py`: numerical checks of every closed-form result in the notes (kink constants, Pöschl–Teller spectrum, RF and RB correlators, Fourier transform of φ₀′, Laplace term from the bubble projection). Runs in a few seconds.

## Run

From the repo root:

```bash
python3 bubble_theory/check_projection.py
python3 bubble_theory/frg_anisotropy.py --A A2
```

## Main results so far

- The φ⁴ → interface projection gives σ, η, 2h, the Laplace term and the disorder correlators exactly (checked numerically).
- With the dislocation noise the infinite-system threshold of a flat wall is isotropic to within a few percent: the force is isotropic below about 2w, and the Larkin length is about 2w.
- The anisotropy seen at finite L is a finite-size effect, odd in a by symmetry, with a single amplitude C(L) ∝ L^-0.6 for both edge and screw and both values of a.
- It is not orientation smearing: edge walls stay straight (⟨cos 2θ⟩ ≈ 0.65), yet their threshold anisotropy is C ≈ 0.1 at L = 2048.
