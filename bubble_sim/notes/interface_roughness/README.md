# Interface roughness notes

Working note on contour geometry and roughness in isotropic and anisotropic disorder, for flat interfaces and future bubble analysis.

Read [roughness_note.pdf](roughness_note.pdf). The current results cover flat isotropic, edge, screw, and 45-degree interfaces, and isotropic and anisotropic bubbles.

## Contents

- `main.tex`, `flat_isotropic.tex`, `flat_anisotropic.tex`, `bubbles.tex`: editable LaTeX note.
- `fig1` and `fig2`: isotropic contour, profiles, perimeter, and spectra.
- `fig3` and `fig4`: anisotropic edge results.
- `fig5` and `fig6`: anisotropic screw results.
- `fig7` and `fig8`: isotropic bubble results.
- `fig9` and `fig10`: anisotropic bubble results.
- `fig15` and `fig16`: anisotropic 45-degree results.
- `fig11` to `fig14`: local width of arcs of bubble contours about their mean line (`fig11`, `fig13`) and about a fitted cubic (`fig12`, `fig14`), isotropic and anisotropic.
- `make_figures.py`, `make_edge_figures.py`, `make_screw_figures.py`, `make_diagonal_figures.py`, `make_bubble_figures.py`, `make_bubble_arc_figure.py`: analysis and plotting scripts using saved contours only.
- `statistics.json`, `edge_statistics.json`, `screw_statistics.json`, `diagonal_statistics.json`, `bubble_statistics.json`, `bubble_arc_width.json`: measured quantities and fitting ranges.
- `provenance.json`: input snapshot and SHA-256 checksums.

Each figure is provided as PDF and PNG. All cases use the same contour extraction, area-preserving h(x), arc-length h(s), and equal-realization spectral averaging. Sizes are 256, 512, 1024, and 2048 with 256, 128, 64, and 32 realizations, respectively. Bubbles use 1024, 2048, 4096, and 8192 with 128, 64, 32, and 16 realizations, area-preserving radius h(theta), and arc-length h(s) = |r(s) - c|. Spectral panels show the two largest sizes: the raw spectra S(n) = <|h_n|^2> at fixed mode index 2 <= n <= 45 L / 2048, equal weight per contour, with reference slopes -(1 + 2 zeta). Exponents quoted in the note are fitted to spectra averaged in common q bins (statistics json files).

The screw and 45-degree perimeter panels compare slopes 1.45 and 1.25. Their spectral panels show reference slopes for zeta = 1.45 and zeta = 1.25. Reference exponents are visual guides, not imposed fits.

## Reproduce figures

Requires Python 3, NumPy, and Matplotlib. With the saved `walls.npz` and `contours.npz` files available under `bubble_sim/runs`, run from this folder:

```bash
python3 make_figures.py
python3 make_edge_figures.py
python3 make_screw_figures.py
python3 make_diagonal_figures.py
python3 make_bubble_figures.py
python3 make_bubble_arc_figure.py
```

These commands analyze saved contours and do not run simulations. The figure scripts resolve input paths relative to `bubble_sim`, regardless of the working directory. Statistics are regenerated with the figures. The bundled results use the input snapshot recorded in `provenance.json`; newer input data may give different results.

## Compile the note

Requires pdfLaTeX and the packages listed in `main.tex`.

```bash
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
cp main.pdf roughness_note.pdf
```
