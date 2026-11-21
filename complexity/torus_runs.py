"""Torus transfer-operator runs used in the note.

    python3 torus_runs.py            -> torus.json, torus_decomposition.json

1. tau scan (n = 256) of the quenched pressure, A = 0.6 and 1, with the
   u-participation of the population (1 = uniform).
2. Large-tau grid convergence (A = 0.6, tau A = 6 and 12, n = 128 ... 1024):
   slope -> f_c and intercept -> c of the frozen asymptote.
3. Gap decomposition (n = 256): memory term <Lambda - ln<m>_phi> and one-step
   Jensen term <ln<m>_phi - <ln m>_phi>, with its second-order estimate.
"""
import json
import os
import numpy as np
from multiprocessing import Pool
from torus import pressure, decomposition
from annealed import gle

os.chdir(os.path.dirname(os.path.abspath(__file__)))
SCAN = [0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 6.0, 9.0, 12.0]
DEC = [0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0]


def lam(A, tA):
    return gle(1.0, A) + (np.log(np.sinh(tA) / tA) if tA > 0 else 0.0)


def scan_job(a):
    A, tA = a
    phi, e, loc = pressure(A, tA / A, n=256, steps=2500, burn=300, record_loc=True)
    return dict(A=A, tauA=tA, phi=phi, err=e, loc=loc, Lambda=lam(A, tA))


def conv_job(a):
    n, tA = a
    st = {128: 3000, 256: 2500, 512: 1500, 1024: 800}[n]
    phi, e = pressure(0.6, tA / 0.6, n=n, steps=st, burn=200)
    return dict(n=n, tauA=tA, phi=phi, err=e)


def dec_job(a):
    A, tA = a
    d = decomposition(A, tA / A, n=256, steps=3000, burn=300, every=10)
    L = lam(A, tA)
    return dict(A=A, tauA=tA, Lambda=L, gap=L - d["phi"], memory=L - d["lnEm"],
                jensen=d["jensen"], jensen2=d["jensen2"])


if __name__ == "__main__":
    with Pool(9) as p:
        scan = p.map(scan_job, [(A, t) for A in (0.6, 1.0) for t in SCAN])
        conv = p.map(conv_job, [(n, t) for n in (128, 256, 512, 1024) for t in (6.0, 12.0)])
        dec = p.map(dec_job, [(A, t) for A in (0.3, 0.6, 1.0) for t in DEC])
    cv = []
    for n in (128, 256, 512, 1024):
        a = next(x for x in conv if x["n"] == n and x["tauA"] == 6.0)
        b = next(x for x in conv if x["n"] == n and x["tauA"] == 12.0)
        sl = (b["phi"] - a["phi"]) / (6 / 0.6)
        cv.append(dict(n=n, h=4.0 / n, slope=sl, slope_err=np.hypot(a["err"], b["err"]) / 10,
                       intercept=b["phi"] - 12 / 0.6 * sl))
    json.dump(dict(scan=scan, convergence=cv), open("torus.json", "w"), indent=1)
    json.dump(dec, open("torus_decomposition.json", "w"), indent=1)
    for c in cv:
        print(f"n={c['n']}: slope={c['slope']:.4f}+-{c['slope_err']:.4f} intercept={c['intercept']:.4f}")
