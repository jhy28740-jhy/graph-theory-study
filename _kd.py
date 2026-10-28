# -*- coding: utf-8 -*-
"""kappa3(N) = max{ rho((VV^T)_+) : V^T V = I_3 } 直接流形上升(解析梯度).

grad_Q rho(Q_+) = sign_+(Q) odot (x x^T),  x = Perron(Q_+)
grad_V          = 2 (I - V V^T) G V        (右 O(3)-不变函数的 Stiefel 投影梯度)
回缩: QR.  步长回溯."""
import numpy as np, json, sys

def rho_plus(V):
    Q = V @ V.T
    B = np.maximum(Q, 0.0)
    w, U = np.linalg.eigh(B)
    x = U[:, -1]
    if x.sum() < 0: x = -x
    return float(w[-1]), x, Q

def ascend(N, V0, steps=4000, lr=1.0):
    V, _ = np.linalg.qr(V0)
    r, x, Q = rho_plus(V)
    for t in range(steps):
        S = (Q >= 0).astype(float)
        G = S * np.outer(x, x)
        gv = 2.0 * (G @ V - Q @ (G @ V))
        nrm = np.linalg.norm(gv)
        if nrm < 1e-14: break
        Vn, _ = np.linalg.qr(V + (lr / nrm) * gv)
        rn, xn, Qn = rho_plus(Vn)
        if rn > r + 1e-16:
            V, r, x, Q = Vn, rn, xn, Qn
        else:
            lr *= 0.5
            if lr < 1e-14: break
    return r, V

rho0 = (9 + np.sqrt(5) + 2 * np.sqrt(35)) / 18
rng = np.random.default_rng(20261028)
out = {}
for N in range(3, 13):
    bb = (-1.0, None)
    tries = 300 if N <= 7 else 150
    for s in range(tries):
        r, V = ascend(N, rng.normal(size=(N, 3)), steps=1200)
        if r > bb[0]: bb = (r, V)
    r, V = bb
    Q = V @ V.T; B = np.maximum(Q, 0.0); C = np.maximum(-Q, 0.0)
    w, U = np.linalg.eigh(B); x = U[:, -1]
    if x.sum() < 0: x = -x
    x = np.maximum(x, 0); x /= np.linalg.norm(x)
    alpha = float(x @ Q @ x); beta = float(x @ C @ x)
    out[N] = dict(rho=r, gap=float(rho0 - r), diag=np.sort(np.diag(Q))[::-1].tolist(),
                  x=np.sort(x)[::-1].tolist(), rowplus=np.sort(B.sum(1))[::-1].tolist(),
                  rowminus=np.sort(C.sum(1))[::-1].tolist(), alpha=alpha, beta=beta)
    print("N=%2d  kappa3=%.15f   rho0-kappa3=%+.4e   alpha=%.9f beta=%.9f" % (N, r, rho0 - r, alpha, beta))
    print("     diag(Q)   =", np.round(np.sort(np.diag(Q))[::-1], 9))
    print("     rowsum Q+ =", np.round(np.sort(B.sum(1))[::-1], 9))
    print("     rowsum Q- =", np.round(np.sort(C.sum(1))[::-1], 9))
    print("     Perron x  =", np.round(np.sort(x)[::-1], 9)); sys.stdout.flush()
    if N == 6: np.save("_Q6_opt.npy", Q)
json.dump(out, open("_kappa_direct.json", "w"), indent=1)
mx = max(v["rho"] for v in out.values())
print("\nrho0=%.15f  max_{N=3..12} kappa3(N)=%.15f  rho0-max=%+.4e  #violations=%d"
      % (rho0, mx, rho0 - mx, sum(1 for v in out.values() if v["rho"] > rho0 + 1e-9)))
