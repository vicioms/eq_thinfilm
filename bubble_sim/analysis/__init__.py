"""Shared analysis code for the critical walls of fss_scan.py runs.

    data      run metadata (results*.jsonl), extraction of the walls from the critical
              configurations (cluster side) and loading of the stored walls
    geometry  curves: arc-length resampling, tangent angle, area, ray crossings, orientation
    measure   one wall -> its observables (length, spectra, orientation, overhangs)
    spectra   pooled mode spectra, zeta fits, local zeta
    stats     bootstrap, power laws, depinning threshold (nu, h_c)
    figures   the figures and the text report shared by analyze_flat.py / analyze_bubbles.py
"""
