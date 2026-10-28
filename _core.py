# -*- coding: utf-8 -*-
"""共用核心：F(A)=max_{p in Delta} Lam3(D_sqrt p A D_sqrt p) 的稳健优化 + KKT 诊断.

2026-10-28 修正：softmax 梯度此前写成 g-<p,g>*1（仅在 p 均匀时正确），
正确公式为 g-(sum g)*p = g - Lam3*p，其中 g=diag(Pi3 M Pi3)。有限差分核验 4e-10。
"""
import numpy as np
from scipy.optimize import minimize


def top3(M):
    ev, U = np.linalg.eigh(M)
    o = np.argsort(ev)[::-1]
    return ev[o], U[:, o]


def fg(A, z):
    """z = 未归一化对数权重。返回 (Lam3, d Lam3/d z, eig, p, Pi3)."""
    z = z - z.max()
    p = np.exp(z); p /= p.sum()
    q = np.sqrt(p)
    M = (q[:, None] * A) * q[None, :]
    ev, U = top3(M)
    Q = U[:, :3] @ U[:, :3].T
    g = np.diag(Q @ M @ Q)          # = diag(Pi3 M Pi3); sum(g)=Lam3
    return float(ev[:3].sum()), g - g.sum() * p, ev, p, Q


def Fmax(A, ns=40, seed=13, polish=True):
    n = A.shape[0]
    rng = np.random.default_rng(seed)
    nf = lambda z: -fg(A, z)[0]
    gf = lambda z: -fg(A, z)[1]
    best = (-9.0, None)
    starts = [np.zeros(n)]
    for i in range(n):
        v = np.full(n, -6.0); v[i] = 6.0; starts.append(v)
    starts += [rng.normal(0, 1.0, n) for _ in range(ns)]
    for s0 in starts:
        r = minimize(nf, s0, jac=gf, method="L-BFGS-B",
                     options=dict(maxiter=4000, ftol=1e-16, gtol=1e-13))
        if -r.fun > best[0]:
            best = (-r.fun, r.x)
    if polish:
        r = minimize(nf, best[1], method="Nelder-Mead",
                     options=dict(maxiter=30000, maxfev=60000, xatol=1e-15, fatol=1e-17))
        if -r.fun > best[0]:
            best = (-r.fun, r.x)
        r = minimize(nf, best[1], method="Powell",
                     options=dict(maxiter=30000, xtol=1e-15, ftol=1e-17))
        if -r.fun > best[0]:
            best = (-r.fun, r.x)
    v, gr, ev, p, Q = fg(A, best[1])
    return v, p, ev, Q, best[1]


def kkt_diag(A, p, ev, Q, tol=1e-9):
    """通用 KKT 诊断: (Q odot A) x = Lam3 x  在 supp(p) 上;  自洽时 Q_+ x = Lam3 x."""
    x = np.sqrt(p)
    L = float(ev[:3].sum())
    S = p > tol
    B = Q * A
    r1 = float(np.max(np.abs((B @ x - L * x)[S]))) if S.any() else 0.0
    Qp = np.maximum(Q, 0.0)
    r2 = float(np.max(np.abs((Qp @ x - L * x)[S]))) if S.any() else 0.0
    sc = int(np.array_equal((Q >= -1e-11).astype(float), A))
    alpha = float(x @ Q @ x)
    beta = float(x @ np.maximum(-Q, 0.0) @ x)
    rho = float(np.max(np.linalg.eigvalsh((Qp + Qp.T) / 2)))
    gk = np.diag(Q @ (np.sqrt(p)[:, None] * A * np.sqrt(p)[None, :]) @ Q)
    return dict(L=L, sc=sc, res_univ=r1, res_perron=r2, alpha=alpha, beta=beta,
                rhoQplus=rho, gerr=float(np.max(np.abs(gk - L * p))),
                npos=int((ev > 1e-10).sum()),
                negmass=float(np.clip(-ev, 0, None).sum()),
                act=int(S.sum()))
