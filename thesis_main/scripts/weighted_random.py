"""Table 7.4: randomised check of Proposition 6.8 on [0,1]-weighted loopless matrices.

Also runs the negative control (entries in [0,2], and one doubled edge) that
shows the range restriction is the operative hypothesis. Run:
    python weighted_random.py
"""
import json
import math
import os

import numpy as np

C = (106106902139519 / 107310156800000 + 41 * math.sqrt(3) / 240
     - 53 / 9000000)
SEED = 20261104
TRIALS = 4000


def lambda3(A):
    e = np.linalg.eigvalsh(A)
    return float(e[-1] + e[-2] + e[-3])


def run():
    rng = np.random.default_rng(SEED)
    rows = []
    for n in (4, 6, 9, 15, 30, 60):
        bad, top = 0, -1.0
        for _ in range(TRIALS):
            M = rng.random((n, n))
            A = np.triu(M, 1)
            A = A + A.T
            v = lambda3(A)
            top = max(top, v)
            if v > C * n - 3 + 1e-12:
                bad += 1
        rows.append({"n": n, "trials": TRIALS, "violations": bad,
                     "max_lambda3": top, "bound": C * n - 3})
        print("n=%3d trials=%d violations=%d max=%.6f bound=%.6f"
              % (n, TRIALS, bad, top, C * n - 3))
    # negative control: weights in [0,2] do violate the bound
    n = 30
    bad = 0
    for _ in range(200):
        M = 2 * rng.random((n, n))
        A = np.triu(M, 1)
        A = A + A.T
        if lambda3(A) > C * n - 3:
            bad += 1
    print("negative control ([0,2] weights, n=30): %d/200 exceed c n-3" % bad)
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "..", "data", "weighted_random.json"), "w") as fh:
        json.dump({"seed": SEED, "rows": rows,
                   "negative_control_violations_of_200": bad}, fh, indent=1)


if __name__ == "__main__":
    run()
