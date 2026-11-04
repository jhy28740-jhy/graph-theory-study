"""Exhaustive check of Lambda_3(G) <= c*n-3 over every labelled graph, n=3..6.

Also identifies the extremal graphs up to isomorphism by brute-force canonical
labelling. Reproduces Table 7.1 rows n=3..6. Run:
    python exhaustive_n3_n6.py
"""
import itertools
import json
import math
import os

import numpy as np

C = (106106902139519 / 107310156800000 + 41 * math.sqrt(3) / 240
     - 53 / 9000000)


def lambda3(A):
    e = np.linalg.eigvalsh(A)
    return float(e[-1] + e[-2] + e[-3])


def canonical(A):
    n = A.shape[0]
    best = None
    for perm in itertools.permutations(range(n)):
        key = tuple(A[np.ix_(perm, perm)].flatten().astype(int).tolist())
        if best is None or key < best:
            best = key
    return best


def run():
    out = {}
    for n in range(3, 7):
        m = n * (n - 1) // 2
        best, classes = -1e9, {}
        for bits in range(1 << m):
            A = np.zeros((n, n), dtype=int)
            k = 0
            for i in range(n):
                for j in range(i + 1, n):
                    if (bits >> k) & 1:
                        A[i, j] = A[j, i] = 1
                    k += 1
            v = lambda3(A.astype(float))
            if v > best + 1e-9:
                best, classes = v, {}
            if v > best - 1e-9:
                key = canonical(A)
                classes[key] = classes.get(key, 0) + 1
        out[n] = {"graphs": 1 << m, "max_lambda3": best, "cn_minus_3": C * n - 3,
                  "margin": C * n - 3 - best, "iso_classes": len(classes),
                  "labelled_maximisers": sum(classes.values()),
                  "class_sizes": sorted(classes.values(), reverse=True)}
        print("n=%d graphs=%d max=%.12f cn-3=%.12f margin=%.12f classes=%d"
              % (n, 1 << m, best, C * n - 3, C * n - 3 - best, len(classes)))
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "..", "data", "exhaustive_n3_n6.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    return out


if __name__ == "__main__":
    run()
