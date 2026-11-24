"""Merge avalanche fragments closer than a time gap (default config.merge_gap(m) = 4 max(m, 1)) in each sample."""
import numpy as np
from config import merge_gap


def merge(r, gap=None):
    m = float(r['m']); wdot = float(r['wdot'])
    gap = merge_gap(m) if gap is None else gap
    s, T, k, w, fo = r['s'], r['T'], r['sample'], r['w_on'], r['f_on']
    S_, T_, F_ = [], [], []
    for kk in np.unique(k):
        sel = np.where(k == kk)[0]; sel = sel[np.argsort(w[sel])]
        ton = w[sel] / wdot; tend = ton + T[sel]
        cs, cst, ce, cf = s[sel[0]], ton[0], tend[0], fo[sel[0]]
        for j in range(1, len(sel)):
            if ton[j] - ce < gap:
                cs += s[sel[j]]; ce = max(ce, tend[j])
            else:
                S_.append(cs); T_.append(ce - cst); F_.append(cf)
                cs, cst, ce, cf = s[sel[j]], ton[j], tend[j], fo[sel[j]]
        S_.append(cs); T_.append(ce - cst); F_.append(cf)
    return np.array(S_), np.array(T_), np.array(F_)
