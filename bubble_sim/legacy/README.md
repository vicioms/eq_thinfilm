# Superseded analysis scripts

Replaced by `analyze_flat.py` and `analyze_bubbles.py` (one script per geometry, shared
code in `analysis/`). Kept because some analyses here are not ported (yet) and the result
files in `runs/fss_D0.04/` (`contour_analysis.md`, `qew_a0.*`, `kpz_edge.*`, `sides.*`, ...)
were made with them. Run them from `bubble_sim/`, e.g. `python legacy/qew_a0.py --dir runs/fss_D0.04`.

| script | what it does | in the new scripts |
|---|---|---|
| `fss_analyze.py` | h_d table, var(h_d) -> ν, bubble elongation | threshold.png |
| `fss_shapes.py` | R(θ) rays (`shapes.npz`), contour extraction (`--contours`) | `analyze_bubbles.py --extract`, profile ζ, mean shape |
| `fss_contours.py` | ψ spectra, winding structure function, divider dimension, overhang lengths, normal-angle distribution, support function | ψ spectra, orientation, overhangs (structure function, divider dimension, support function not ported) |
| `qew_a0.py` | a = 0 vs qEW: fit-window scan, ν = 1/(2 − ζ), h_d distribution collapse | partly (ζ, ν) |
| `kpz_edge.py` | qKPZ tests on the edge sides of the bubbles | not ported (the flat edge walls test qKPZ directly) |
| `sides.py` | ζ of the edge / screw sides of the bubbles | not ported |
| `fss_flat.py` | flat walls: extraction and analysis | `analyze_flat.py` |
| `perimeter.py` | wall length vs size | length.png |
