"""Tables 7.2 and 7.3: equal and optimally weighted blow-ups of K1 v C5.

Confirms the closed form of Proposition 7.1,
    Lambda_3(H[q]) = ((sqrt5+sqrt6+3)/6) * 6q - 3,
and the convergence of the optimally weighted family to c_w * n - 3. Run:
    python blowup_families.py
"""
import json
import math
import os

import numpy as np

C = (106106902139519 / 107310156800000 + 41 * math.sqrt(3) / 240
     - 53 / 9000000)
CW = (9 + math.sqrt(5) + 2 * math.sqrt(35)) / 18
H = np.array([[0, 1, 1, 1, 1, 1],
              [1, 0, 1, 0, 0, 1],
              [1, 1, 0, 1, 0, 0],
              [1, 0, 1, 0, 1, 0],
              [1, 0, 0, 1, 0, 1],
              [1, 1, 0, 0, 1, 0]], dtype=float)


def lambda3(A):
    e = np.linalg.eigvalsh(A)
    return float(e[-1] + e[-2] + e[-3])


def equal_blowup(q):
    A = (np.kron(H, np.ones((q, q)))
         + np.kron(np.eye(6), np.ones((q, q)) - np.eye(q)))
    np.fill_diagonal(A, 0.0)
    return A


def weighted_blowup(sizes):
    n = sum(sizes)
    A = np.zeros((n, n))
    off, offs = 0, []
    for s in sizes:
        offs.append(off)
        off += s
    for i in range(6):
        for j in range(6):
            si, sj = offs[i], offs[j]
            if i == j:
                A[si:si + sizes[i], sj:sj + sizes[j]] = 1.0 - np.eye(sizes[i])
            elif H[i, j]:
                A[si:si + sizes[i], sj:sj + sizes[j]] = 1.0
    np.fill_diagonal(A, 0.0)
    return A


def run():
    eq = []
    print("equal blow-ups H[q]")
    for q in (1, 2, 3, 5, 10, 20, 50, 100, 250, 500):
        n = 6 * q
        v = lambda3(equal_blowup(q))
        closed = (math.sqrt(5) + math.sqrt(6) + 3) / 6 * n - 3
        eq.append([q, n, v, v / n, v - (C * n - 3), v - (CW * n - 3),
                   abs(v - closed)])
        print("  q=%3d n=%5d L3=%14.6f  closed-form error=%.3e" %
              (q, n, v, abs(v - closed)))
    b = (7 + 5 * math.sqrt(7)) / 126
    op = []
    print("optimally weighted blow-ups")
    for n in (600, 1200, 3000):
        q = max(1, int(round(b * n)))
        sizes = [n - 5 * q] + [q] * 5
        v = lambda3(weighted_blowup(sizes))
        op.append([n, sizes, v, v / n, v - (C * n - 3), v - (CW * n - 3)])
        print("  n=%5d sizes=%s L3=%.6f  L3-(c_w n-3)=%+.6f" %
              (n, sizes, v, v - (CW * n - 3)))
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "..", "data", "blowups.json"), "w") as fh:
        json.dump({"equal": eq, "weighted": op,
                   "b_star": b, "p_star": 1 - 5 * b}, fh, indent=1)


if __name__ == "__main__":
    run()
