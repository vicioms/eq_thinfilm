"""Phase-annealed, amplitude-quenched rate for the PL Gaussian model.

    python3 semi_annealed.py      -> semi_annealed_cos.json

Lambda_q(tau) = lim (1/K) < ln || prod_i T^{tau}_{R_i} || > with the site
amplitudes R_i quenched and the phases averaged (annealed.amplitude_quenched_rate),
and Sigma_sa(f) = inf_tau [Lambda_q(tau) - tau f].  For a fixed amplitude
(triangle model) Sigma_sa = Sigma_ann.
"""
import json
import numpy as np
from multiprocessing import Pool
from annealed import amplitude_quenched_rate

SIGMA, M = 0.6, 16
TAUS = np.linspace(0.0, 2.0, 11)


def one(tau):
    lq, e, la = amplitude_quenched_rate(SIGMA, M, tau=tau, K=60_000, seed=int(1000 * tau) + 1)
    return tau, lq, e, la


if __name__ == "__main__":
    with Pool(4) as p:
        res = p.map(one, TAUS)
    out = {"sigma": SIGMA, "M": M, "tau": [r[0] for r in res], "Lambda_q": [r[1] for r in res],
           "Lambda_q_err": [r[2] for r in res], "Lambda_ann": [r[3] for r in res]}
    json.dump(out, open("semi_annealed_cos.json", "w"), indent=1)
    for r in res:
        print("tau=%.2f Lambda_q=%.5f+-%.5f Lambda_ann=%.5f gap=%.5f" % (r[0], r[1], r[2], r[3], r[3] - r[1]))
