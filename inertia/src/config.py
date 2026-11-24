"""Shared paths and analysis thresholds for the inertia project.

Every script in runs/, analysis/ and archive/ starts with

    import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'src'))
    import config

which puts src/ on the import path and makes the inertia folder the working directory, so
that the relative paths data/..., figures/... resolve wherever the script is launched from.

Units: kappa = eta = 1, disorder amplitude A = 1, disorder correlation length r_f = 1.
Each threshold below is either fixed by a scale of the model or read off a gap in a measured
distribution; analysis/thresholds.py makes the figure (figures/thresholds.pdf) and the numbers
(data/stats/thresholds.json) that justify it, including how much results move when it is varied.
"""
import os
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
FIGURES = ROOT / 'figures'
STATS = DATA / 'stats'
os.chdir(ROOT)

# --- activity -------------------------------------------------------------------------------
# A site is active (moving) when |du/dt| > V_ACT.  In partial sliding the moving sites form a
# peak at |du/dt| ~ 0.1, while pinned sites ring down and fill a flat tail in log|du/dt| down
# to 1e-6, so no value is canonical and the active fraction itself depends on V_ACT
# (m = 64, f = 0.198: rho = 0.55, 0.42, 0.30 for V_ACT = 0.0025, 0.02, 0.05).  What we use it
# for does not: the spreading exponents change by less than a tenth of their errors for V_ACT
# between 0.0025 and 0.1 (thresholds.json, 'spreading').
V_ACT = 0.02
# A line (or a spreading run) is at rest when max_i |du_i/dt| < V_REST = V_ACT / 20: far below
# any moving site, so stopping times do not depend on it.
V_REST = 1e-3

# --- uniform versus partial sliding ---------------------------------------------------------
# In uniform sliding every site moves, and the active fraction equals that of Gaussian
# velocities with the measured mean and variance, rho_G = P(|vbar + dv| > V_ACT).  A state is
# partial when rho < (1 - RHO_TOL) rho_G.  For m >= 64 uniform states have rho / rho_G >= 0.987
# and partial ones <= 0.953, so any RHO_TOL between 0.015 and 0.045 gives the same split.
RHO_TOL = 0.03

# --- avalanches under well driving ----------------------------------------------------------
# Start of an avalanche: mean velocity above C_TH * wdot.  With 1.5 or 10 instead of 3 the
# number of small events changes by up to 35%, the slip share by about one percentage point
# (thresholds.json, 'detection').
C_TH = 3.0
# Inertial ringing splits one event into fragments.  Fragments closer in time than
# MERGE_GAP(m) = 4 m / eta (two velocity relaxation times 2m/eta) are merged.  For m >= 16 the
# gap distribution has a deep trough between ringing (~ 6 sqrt(m)) and true waiting times
# (>~ 1e4), 4m lies in it, and every gap inside the trough merges the same fragments.  For
# m <= 8 ringing overlaps the ordinary pauses of qEW avalanches, and the merged fraction depends
# on the gap (m = 4: 0.25, 0.49, 0.56 for 2m, 4m, 16m; thresholds.json, 'merge_gap...').


def merge_gap(m, eta=1.0):
    return 4.0 * max(m, 1.0) / eta


# A slip (king avalanche) is an event larger than any avalanche of the overdamped line (m = 0)
# under the same driving, well and system size: s > max s(m = 0).  Size alone cannot do better:
# P(s) overlaps the qEW tail at intermediate m, and s / L > 1 already includes large qEW
# avalanches.  The kinetic friction of slips, f_k = (f_on + f_stop) / 2, changes by < 3% between
# this definition and s > 4L; the share of the advance carried by slips depends on it.


def slip_cut(s_m0):
    """Slip threshold from the sizes s_m0 of the overdamped avalanches under the same driving."""
    return float(max(s_m0))

# --- finite-time ramps ----------------------------------------------------------------------
# Ramps measure finite-time thresholds by definition; the long-time threshold comes from
# lifetimes.  Ramp-down: the line has stopped when its centre of mass moves less than
# RAMP_STOP r_f during the second half of a dwell.  Ramp-up: it has depinned when it travels
# RAMP_GO r_f without stopping (far beyond the size of any single avalanche at that L).
RAMP_STOP = 0.02
RAMP_GO = 30.0
# Lifetimes: a quenched line has stopped when its centre of mass advances less than LIFE_STOP r_f
# during one observation chunk (max(200, m)); by then its mean speed is < 1e-4, far below V_ACT.
LIFE_STOP = 0.01
