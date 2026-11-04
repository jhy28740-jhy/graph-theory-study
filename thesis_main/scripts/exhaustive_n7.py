"""Exhaustive check of Lambda_3(G) <= c*n-3 over all 2097152 labelled graphs on 7 vertices.

Verifies that every maximiser is isomorphic to K2 v C5 (degree sequence
(6,6,4,4,4,4,4) with a 5-cycle on the degree-four vertices). Reproduces the
n=7 row of Table 7.1. Run:
    python exhaustive_n7.py
"""
import json
import math
import os

import numpy as np

C = (106106902139519 / 107310156800000 + 41 * math.sqrt(3) / 240
     - 53 / 9000000)
N, CHUNK = 1 << 21, 1 << 15


def run():
    n = 7
    ii, jj = [], []
    for i in range(n):
        for j in range(i + 1, n):
            ii.append(i)
            jj.append(j)
    ii, jj = np.array(ii), np.array(jj)
    k = np.arange(21, dtype=np.int64)
    best, winners = -1e9, []
    for start in range(0, N, CHUNK):
        stop = min(start + CHUNK, N)
        bits = np.arange(start, stop, dtype=np.int64)
        mask = (((bits[:, None] >> k[None, :]) & 1)).astype(np.float64)
        A = np.zeros((mask.shape[0], n, n))
        A[:, ii, jj] = mask
        A[:, jj, ii] = mask
        ev = np.linalg.eigvalsh(A)
        lam3 = ev[:, -1] + ev[:, -2] + ev[:, -3]
        top = float(lam3.max())
        if top > best + 1e-9:
            best, winners = top, []
        if top > best - 1e-9:
            winners += [int(bits[t]) for t in np.where(lam3 > best - 1e-9)[0]]
    winners = [b for b in winners if True]
    ok = 0
    for bits in winners:
        A = np.zeros((n, n), dtype=int)
        t = 0
        for i in range(n):
            for j in range(i + 1, n):
                A[i, j] = A[j, i] = (bits >> t) & 1
                t += 1
        deg = A.sum(1)
        apex = np.where(deg == 6)[0]
        rim = np.where(deg == 4)[0]
        good = (len(apex) == 2 and len(rim) == 5
                and A[apex[0], apex[1]] == 1
                and sorted(A[np.ix_(rim, rim)].sum(1).tolist()) == [2] * 5)
        ok += int(good)
    res = {"n": n, "graphs": N, "max_lambda3": best, "cn_minus_3": C * n - 3,
           "margin": C * n - 3 - best, "labelled_maximisers": len(winners),
           "all_isomorphic_to_K2_join_C5": ok == len(winners)}
    print(json.dumps(res, indent=1))
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "..", "data", "exhaustive_n7.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    return res


if __name__ == "__main__":
    run()
