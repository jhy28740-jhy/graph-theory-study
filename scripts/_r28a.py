# -*- coding: utf-8 -*-
# _r28a.py  第 28 轮 A 段：P 的上确界的精确数值探索（检验 P1'：sup P =? (3+sqrt5)/5，5-环最优）
#   P(mu) = sum_{i!=j, c_ij<0} rho_i rho_j c_ij^2 = ∫∫ (cos(phi-phi'))_-^2 dmu dmu'
#   约束: rho>=0, sum rho=2, sum rho e^{2i phi}=0  (二维各向同性)
#   用 SLSQP（显式约束）+ 环初始化 + 随机初始化，避免 _r27c 罚函数法的欠收敛。
import numpy as np
from scipy.optimize import minimize
rng = np.random.default_rng(20281)
sqrt5 = np.sqrt(5.0)
P_champ = (3+sqrt5)/5          # 5-环（冠军）的值 = 1.0472135955
P_bound = 10/(3*np.pi)         # Fourier 上界 = 1.0610329539

def Pval(rho, ang):
    c = np.cos(ang[:, None]-ang[None, :])
    return float((rho[:, None]*rho[None, :]*np.where(c < 0, c**2, 0.0)).sum())

def optimize(N, x0, ntries=1):
    # 变量: rho (N), ang (N). 约束: sum rho=2, sum rho cos2ang=0, sum rho sin2ang=0.
    def negf(x):
        rho = x[:N]; ang = x[N:]
        if np.any(rho < -1e-9):
            return 1e3
        return -Pval(np.maximum(rho, 0), ang)
    cons = [
        {'type': 'eq', 'fun': lambda x: x[:N].sum()-2.0},
        {'type': 'eq', 'fun': lambda x: (x[:N]*np.cos(2*x[N:])).sum()},
        {'type': 'eq', 'fun': lambda x: (x[:N]*np.sin(2*x[N:])).sum()},
    ]
    bnds = [(0, 2.0)]*N+[(None, None)]*N
    r = minimize(negf, x0, method='SLSQP', constraints=cons, bounds=bnds,
                 options={'maxiter': 3000, 'ftol': 1e-14})
    x = r.x; rho = np.maximum(x[:N], 0); ang = x[N:]
    viol = abs(rho.sum()-2)+abs((rho*np.cos(2*ang)).sum())+abs((rho*np.sin(2*ang)).sum())
    return Pval(rho, ang), rho, ang, viol

print("="*96)
print("sup P 的数值探索：环初始化 + 随机初始化（SLSQP，显式约束）")
print("  冠军 5-环 P = (3+sqrt5)/5 = %.12f ;  Fourier 上界 10/(3pi) = %.12f" % (P_champ, P_bound))
print("="*96)
best_overall = -1; best_cfg = None
for N in range(3, 17):
    best = -1; bx = None
    # (a) 环初始化: 等权 rho=2/N, 等角 + 随机旋转/扰动
    for rot in np.linspace(0, 2*np.pi, 12, endpoint=False):
        for pert in [0.0, 0.05, 0.15]:
            rho0 = np.full(N, 2.0/N)*(1+pert*(rng.random(N)-0.5))
            ang0 = rot+2*np.pi*np.arange(N)/N+pert*(rng.random(N)-0.5)
            v, rho, ang, viol = optimize(N, np.concatenate([rho0, ang0]))
            if viol < 1e-6 and v > best:
                best = v; bx = (rho.copy(), ang.copy(), viol)
    # (b) 随机初始化
    nrand = 200 if N <= 10 else 80
    for t in range(nrand):
        rho0 = 2*rng.random(N)/N; ang0 = rng.random(N)*2*np.pi
        v, rho, ang, viol = optimize(N, np.concatenate([rho0, ang0]))
        if viol < 1e-6 and v > best:
            best = v; bx = (rho.copy(), ang.copy(), viol)
    if bx is not None and best > best_overall:
        best_overall = best; best_cfg = (N, bx)
    print("  N=%2d  sup P >= %.12f  (违例 %.1e)  %s" %
          (N, best, bx[2] if bx else -1, "<-- 新冠军" if best >= best_overall-1e-12 and N == (best_cfg[0] if best_cfg else -1) else ""))
print("-"*96)
print("  全局最大 sup P >= %.12f  在 N=%d" % (best_overall, best_cfg[0]))
print("  与 5-环冠军 (3+sqrt5)/5 的差 = %.3e" % (best_overall-P_champ))
print("  与 Fourier 上界 10/(3pi) 的差 = %.3e" % (P_bound-best_overall))
if best_cfg is not None:
    N, (rho, ang, viol) = best_cfg
    order = np.argsort(ang)
    print("  最优构型 (N=%d):" % N)
    print("    rho =", np.array2string(rho[order], precision=6))
    print("    ang =", np.array2string(ang[order], precision=6))
    dphi = np.diff(np.concatenate([ang[order], [ang[order][0]+2*np.pi]]))
    print("    相邻夹角 =", np.array2string(dphi, precision=6), " (和=%.6f)" % dphi.sum())
print("ALL R28A DONE")
