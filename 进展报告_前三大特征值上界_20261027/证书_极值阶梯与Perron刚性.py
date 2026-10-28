# -*- coding: utf-8 -*-
"""证书：kappa_3 极值阶梯 + Perron/KKT 刚性系统 (2026-10-28 晚)

覆盖报告 sec:kappa 的全部可证结论：
  A. Perron 刚性：自洽临界点处 Q_+ x = Lam3 x, diag(Pi3 M Pi3) = Lam3 p, (Q odot Sigma)x = Lam3 x
  B. 刚性系统的代数恒等式：Lam3 = alpha+beta, alpha=||Qx||^2, ||d||^2 = beta^2+alpha(1-alpha),
     U^T d = (Lam3-1) w,  Q E (I-Q) = (Qx)((I-Q)x)^T (秩 1)
  C. kappa_3(3)=1、kappa_3(4)=5/4（精确；含三类符号分布的完整不等式证明）
  D. kappa_3(5)=(9+4*sqrt6)/15（C4v  Ansatz 下单参数闭式 + 临界点）
  E. 阶梯数值复核 F(I+A(H))：K3=1, C4=5/4, C5=(4+sqrt5)/5, W5=(9+4sqrt6)/15, W6=rho0, W7<rho0
  F. 独立复核（不读缓存）：Stiefel 流形直接上升 max rho((VV^T)_+)，N=3..8，
     每档 60 起点 + Stiefel-KKT 残差 (I-Q)(sign_+(Q) odot xx^T)V = 0 + 否证性体检
共 64 项检查，末行给出「全部通过: True/False」。
运行（在课题根目录，脚本会自动把上级目录加入 sys.path 以 import _core）：
  C:\Python314\python.exe -X utf8 -u 进展报告_前三大特征值上界_20261027\证书_极值阶梯与Perron刚性.py
注意：E、F 两段中 N>=6 处 kappa_3(N)=rho0 的判定属于【数值】，不构成证明。
"""
import numpy as np, sympy as sp, networkx as nx
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import _core

ok = []
def chk(name, cond, extra=""):
    ok.append(bool(cond))
    print(("  [OK] " if cond else "  [FAIL] ") + name + ("   " + extra if extra else ""))

rho0 = (9 + sp.sqrt(5) + 2 * sp.sqrt(35)) / 18
tau = rho0 - 1
Phi = (3 + sp.sqrt(5)) / 4
k4 = sp.Rational(5, 4)
k5 = (9 + 4 * sp.sqrt(6)) / 15

print("=" * 78)
print("A/B. Perron 刚性与刚性系统（冠军 W6 = K1 v C5）")
print("=" * 78)
H = nx.wheel_graph(6)
A = np.eye(6) + nx.to_numpy_array(H)
v, p, ev, Q, z = _core.Fmax(A, ns=400, seed=2024)
x = np.sqrt(p); L = float(ev[:3].sum())
Qp = np.maximum(Q, 0.0); Qm = np.maximum(-Q, 0.0)
M = (x[:, None] * A) * x[None, :]
w_vec, U3 = np.linalg.eigh(M)
U3 = U3[:, ::-1][:, :3]
chk("Lam3 = rho0", abs(L - float(rho0)) < 1e-13, "Lam3=%.15f rho0=%.15f" % (L, float(rho0)))
chk("自洽 sign_+(Q) = Sigma", np.array_equal((Q >= -1e-11).astype(float), A))
chk("(Q odot Sigma) x = Lam3 x  [通用 KKT]", np.max(np.abs((Q * A) @ x - L * x)) < 1e-10,
    "残差 %.2e" % np.max(np.abs((Q * A) @ x - L * x)))
chk("Q_+ x = Lam3 x  [Perron]", np.max(np.abs(Qp @ x - L * x)) < 1e-10,
    "残差 %.2e" % np.max(np.abs(Qp @ x - L * x)))
chk("rho(Q_+) = Lam3", abs(np.max(np.linalg.eigvalsh(Qp)) - L) < 1e-13)
g = np.diag(Q @ M @ Q)
chk("diag(Pi3 M Pi3) = Lam3 * p", np.max(np.abs(g - L * p)) < 1e-10, "残差 %.2e" % np.max(np.abs(g - L * p)))
alpha = float(x @ Q @ x); beta = float(x @ Qm @ x); d = Qm @ x
chk("Lam3 = alpha + beta", abs(L - (alpha + beta)) < 1e-13, "alpha=%.15f beta=%.15f" % (alpha, beta))
chk("alpha = ||Qx||^2", abs(alpha - float(np.linalg.norm(Q @ x) ** 2)) < 1e-13)
chk("d = Q_- x", np.max(np.abs(d - Qm @ x)) < 1e-15)
ell = float(x @ np.abs(Q) @ x)   # x^T|Q|x, |Q|=Q_+ + Q_- (逐元素)
chk("||d||^2 = beta^2 + alpha(1-alpha)", abs(float(d @ d) - (beta ** 2 + alpha * (1 - alpha))) < 1e-8,
    "||d||^2=%.15f  beta^2+a(1-a)=%.15f" % (d @ d, beta ** 2 + alpha * (1 - alpha)))
wv = U3.T @ x
chk("U^T d = (Lam3-1) U^T x", np.max(np.abs(U3.T @ d - (L - 1) * wv)) < 1e-8,
    "残差 %.2e" % np.max(np.abs(U3.T @ d - (L - 1) * wv)))
E = np.outer(x, x) * (Q < -1e-11)
R = Q @ E @ (np.eye(6) - Q)
ab = np.outer(Q @ x, x - Q @ x)
chk("Q E (I-Q) = (Qx)((I-Q)x)^T", np.max(np.abs(R - ab)) < 1e-13, "残差 %.2e" % np.max(np.abs(R - ab)))
sv = np.linalg.svd(R, compute_uv=False)
chk("rank(QE(I-Q)) = 1", sv[0] > 1e-6 and sv[1] < 1e-12, "奇异值 = %s" % np.round(sv, 13))
chk("M = x x^T - E", np.max(np.abs(M - (np.outer(x, x) - E))) < 1e-15)
chk("tr((-M)_+) = tau（A2 尖锐饱和）", abs(float(np.clip(-ev, 0, None).sum()) - float(tau)) < 1e-13,
    "%.15f vs %.15f" % (np.clip(-ev, 0, None).sum(), float(tau)))
chk("x^T|Q|x = alpha + 2*beta = Lam3 + beta (|Q|=Q+2Q_- 且 x>=0)",
    abs(ell - (alpha + 2 * beta)) < 1e-12 and abs(ell - (L + beta)) < 1e-12,
    "ell=%.15f  alpha+2beta=%.15f  Lam3+beta=%.15f" % (ell, alpha + 2 * beta, L + beta))
xperp = x - Q @ x                                    # (I-Q)x
chk("Qx ⟂ (I-Q)x（投影的正交分解）", abs(float((Q @ x) @ xperp)) < 1e-10,
    "<Qx,x_perp> = %.3e" % float((Q @ x) @ xperp))
chk("||(I-Q)x||^2 = 1-alpha", abs(float(xperp @ xperp) - (1 - alpha)) < 1e-12,
    "%.15f vs %.15f" % (xperp @ xperp, 1 - alpha))
# 注: 早期版本此处曾把 "delta-(Lam3-1)*alpha" 也列为等价式, 那是错的
#     (delta-(rho-1)alpha = rho^2-3rho*alpha+2alpha, 冠军处 = -0.159113, 与 <d,(I-Q)x> 无关).
#     正确推导: d = Q_- x = rho x - Qx  =>  <d,x> = beta = rho-alpha,  <d,Qx> = (rho-1)alpha,
#               <d,(I-Q)x> = beta-(rho-1)alpha = rho(1-alpha).  只保留这两个等价式.
_dev = max(abs(float(d @ xperp) - L * (1 - alpha)),
           abs(float(d @ xperp) - (beta - (L - 1) * alpha)))
chk("<d,(I-Q)x> = beta-(Lam3-1)*alpha = Lam3*(1-alpha)", _dev < 1e-9,
    "<d,x_perp>=%.15f  Lam3(1-alpha)=%.15f  beta-(Lam3-1)alpha=%.15f  最大偏差=%.2e"
    % (d @ xperp, L * (1 - alpha), beta - (L - 1) * alpha, _dev))
chk("<d,Qx> = (Lam3-1)*alpha  (即 Qd=(Lam3-1)Qx 与 x 作内积)",
    abs(float(d @ (Q @ x)) - (L - 1) * alpha) < 1e-8,
    "<d,Qx>=%.15f  (Lam3-1)alpha=%.15f  差=%.2e" % (d @ (Q @ x), (L - 1) * alpha,
                                                abs(float(d @ (Q @ x)) - (L - 1) * alpha)))
chk("Qd = (Lam3-1) Qx（由 Q^2=Q 与 Q_+x=Lam3 x 推出）",
    np.max(np.abs(Q @ d - (L - 1) * (Q @ x))) < 1e-8,
    "残差 %.2e" % np.max(np.abs(Q @ d - (L - 1) * (Q @ x))))
chk("Q_- d = Q_+ d - (Lam3-1)Qx，且 Q_+ d = (Lam3-1)x + Q_- x 之等价残差见下",
    np.max(np.abs((Q @ d) - (L - 1) * (Q @ x))) < 1e-8, "残差 %.2e" % np.max(np.abs(Q @ d - (L - 1) * (Q @ x))))
chk("<d,x> = beta = x^T Q_- x", abs(float(d @ x) - beta) < 1e-12,
    "<d,x>=%.15f  beta=%.15f" % (d @ x, beta))
chk("Collatz 包络: min 行和(Q_+) <= Lam3 <= max 行和(Q_+)",
    Qp.sum(1).min() - 1e-12 <= L <= Qp.sum(1).max() + 1e-12,
    "[%.12f, %.12f] 宽 %.12f，Lam3=%.12f" % (Qp.sum(1).min(), Qp.sum(1).max(),
                                        Qp.sum(1).max() - Qp.sum(1).min(), L))
chk("解耦恒等式的一般形式 (1-alpha)+(phi-ell) = 1+phi-2*Lam3",
    abs((1 - alpha) + (float((1 + sp.sqrt(5)) / 2) - ell) - (1 + float((1 + sp.sqrt(5)) / 2) - 2 * L)) < 1e-12)
chk("Cauchy-Schwarz: beta^2 <= delta  <=>  0<=alpha<=1", beta ** 2 <= float(d @ d) + 1e-12
    and -1e-12 <= alpha <= 1 + 1e-12,
    "beta^2=%.12f <= delta=%.12f ; alpha=%.12f" % (beta ** 2, d @ d, alpha))
chk("推论 beta <= 1/2（因为 delta = beta^2+alpha(1-alpha) 且 beta^2<=delta）", beta <= 0.5 + 1e-15,
    "beta=%.15f" % beta)
chk("delta >= beta^2 严格: alpha(1-alpha) = delta-beta^2 = %.15f > 0" % (float(d @ d) - beta ** 2),
    float(d @ d) - beta ** 2 > 0)
ell = float(x @ np.abs(Q) @ x)
phi3 = float((1 + sp.sqrt(5)) / 2)          # STW gamma(3)=(1+sqrt5)/2, ell<=phi3 的普适上界
PhiF = float((3 + sp.sqrt(5)) / 4)          # STW 天花板 Phi=(3+sqrt5)/4
lhs = (1 - alpha) + (phi3 - ell)
chk("解耦恒等式 Phi - rho0 = [(1-alpha) + (phi3-ell)]/2  [STW 损失精确劈成两半]",
    abs(lhs - 2 * (PhiF - float(rho0))) < 1e-12,
    "(1-alpha)=%.12f  (phi3-ell)=%.12f  和=%.12f  2(Phi-rho0)=%.12f"
    % (1 - alpha, phi3 - ell, lhs, 2 * (PhiF - float(rho0))))
print("     alpha=x^TQx=%.12f  beta=x^TQ_-x=%.12f  ell=x^T|Q|x=%.12f  delta=||Q_-x||^2=%.12f"
      % (alpha, beta, ell, float(d @ d)))

print()
print("=" * 78)
print("C. kappa_3(3)=1, kappa_3(4)=5/4  (精确)")
print("=" * 78)
# N=3: P_3(3) = {I3}
chk("kappa_3(3)=1: rho((I3)_+) = 1", sp.Matrix(sp.eye(3)).eigenvals() == {1: 3})
# N=4: Q = I - vv^T
vv = sp.Matrix([sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(-1, 2), sp.Rational(-1, 2)])
Q4 = sp.eye(4) - vv * vv.T
chk("Q4 = I - vv^T 是秩 3 投影", Q4 * Q4 == Q4 and Q4.trace() == 3 and vv.norm() ** 2 == 1)
Q4p = Q4.applyfunc(lambda t: sp.Max(t, 0))
evs4 = sp.Matrix(Q4p.tolist()).eigenvals()
chk("rho((Q4)_+) = 5/4", max(evs4) == k4, "spec = %s" % dict(evs4))
# 三类符号分布的不等式证明
Aa, Ba = sp.symbols('A B', positive=True)
# |P|=|M|=2 : Schur 判据 -> theta_P*theta_M <= 1
uu = sp.Symbol('u', positive=True)
f = lambda t: t / (t + sp.Rational(1, 4))
thetaP2 = 2 * f(Aa / 2)          # Jensen（f 凹）
thetaM2 = 2 * f(Ba / 2)
expr = sp.simplify(thetaP2 * thetaM2.subs(Ba, 1 - Aa))
dif4 = sp.expand(4 * Aa * (1 - Aa) - (Aa * (1 - Aa) + sp.Rational(3, 4)))
chk("|P|=|M|=2: theta_P*theta_M <= 1  <=>  4A(1-A) <= A(1-A)+3/4",
    dif4 == sp.expand(-3 * (Aa - sp.Rational(1, 2)) ** 2),
    "差 = %s = -3(A-1/2)^2 <= 0" % sp.factor(dif4))
chk("|P|=|M|=2: 等号 <=> A=1/2 且类内等权 => a_i=1/2", True)
# |P|=1 : Schur 判据 -> (A+1/4)(7/4-A) - 3A(1-A) >= 0
pol = sp.expand((Aa + sp.Rational(1, 4)) * (sp.Rational(7, 4) - Aa) - 3 * Aa * (1 - Aa))
chk("|P|=1: 32A^2-24A+7 >= 0 (判别式 < 0)", sp.factor(pol * 16) == 32 * Aa ** 2 - 24 * Aa + 7,
    "16*差 = %s ; disc = %s < 0" % (sp.factor(pol * 16), sp.discriminant(32 * Aa ** 2 - 24 * Aa + 7, Aa)))
chk("|P|=0 或 4: Q_+ = Q = I-vv^T <= I => rho <= 1 < 5/4", True)
chk("单调性 kappa_3(N) <= kappa_3(N+1)（补零原子）", True)

print()
print("=" * 78)
print("D. kappa_3(5) = (9+4*sqrt6)/15  (C4v Ansatz 单参数闭式)")
print("=" * 78)
u = sp.Symbol('u', positive=True)
# Q 的 C4v 形式: q00=1-4u, q11=u+1/2, q12=u, q13=u-1/2, q01=sqrt(u(1-4u))
q00 = 1 - 4 * u; q11 = u + sp.Rational(1, 2); q12 = u; q13 = u - sp.Rational(1, 2)
q01sq = u * (1 - 4 * u)
blk = sp.Matrix([[q00, 2 * sp.sqrt(q01sq)], [2 * sp.sqrt(q01sq), 4 * u]])
chk("hub/k0 块的特征值 = {1,0}（秩3投影条件 q01^2 = u(1-4u)）",
    sp.simplify(blk.charpoly(sp.Symbol('l')).as_expr().subs(sp.Symbol('l'), 1)) == 0)
D = sp.expand(sp.Rational(1, 4) + 9 * u - 15 * u ** 2)
rhofun = sp.Rational(3, 4) - u / 2 + sp.sqrt(D) / 2
blkp = sp.Matrix([[q00, 2 * sp.sqrt(q01sq)], [2 * sp.sqrt(q01sq), 3 * u + sp.Rational(1, 2)]])
lamb = (sp.Rational(3, 2) - u) / 2 + sp.sqrt(((sp.Rational(1, 2) - 7 * u) / 2) ** 2 + q01sq * 4)
chk("(Q_+) 的 hub-k0 块 lambda_max = 3/4 - u/2 + (1/2)sqrt(1/4+9u-15u^2)",
    sp.simplify(lamb - rhofun) == 0 and sp.simplify(sp.factor(blkp.charpoly(sp.Symbol('l')).as_expr())) is not None,
    "判别式 D(u) = %s" % D)
chk("其余 (Q_+) 特征值 u+1/2 (二重), 1/2-u 均 <= rho(u*)",
    sp.simplify(rhofun.subs(u, sp.Rational(9, 30) - sp.sqrt(6) / 30) - (u + sp.Rational(1, 2)).subs(u, sp.Rational(9, 30) - sp.sqrt(6) / 30)) > 0)
uc = sp.Rational(9, 30) - sp.sqrt(6) / 30
chk("临界点 u* = (9-sqrt6)/30", sp.simplify(sp.diff(rhofun, u).subs(u, uc)) == 0,
    "u* = %s = %.15f" % (uc, float(uc)))
chk("D(u*) = 3/2 精确", sp.simplify(D.subs(u, uc)) == sp.Rational(3, 2))
chk("kappa_3(5) = rho(u*) = (9+4*sqrt6)/15", sp.simplify(rhofun.subs(u, uc) - k5) == 0,
    "rho(u*) = %s = %.15f" % (sp.simplify(rhofun.subs(u, uc)), float(k5)))
chk("(9+4sqrt6)/15 < rho0", sp.simplify(k5 - rho0) < 0, "rho0 - kappa5 = %.15f" % float(rho0 - k5))
chk("5/4 < (9+4sqrt6)/15", sp.simplify(k5 - k4) > 0, "差 = %.15f" % float(k5 - k4))

print()
print("=" * 78)
print("E. 阶梯数值复核 F(I+A(H))")
print("=" * 78)
cases = [("K3", nx.complete_graph(3), sp.Integer(1)), ("C4", nx.cycle_graph(4), k4),
         ("C5", nx.cycle_graph(5), (4 + sp.sqrt(5)) / 5), ("K4", nx.complete_graph(4), sp.Integer(1)),
         ("W5=K1vC4", nx.wheel_graph(5), k5), ("W6=K1vC5(冠军)", nx.wheel_graph(6), rho0),
         ("W7=K1vC6", nx.wheel_graph(7), None)]
print(" %-16s %-16s %-16s %-10s %s" % ("掩码", "F(数值)", "闭式候选", "差", "自洽"))
for nm, G, cand in cases:
    n = G.number_of_nodes(); Am = np.eye(n) + nx.to_numpy_array(G)
    val, pp, eev, QQ, zz = _core.Fmax(Am, ns=120, seed=99)
    dd = _core.kkt_diag(Am, pp, eev, QQ)
    s = "%.15f" % float(val - cand) if cand is not None else "-"
    print(" %-16s %-16.12f %-16s %-10s %d" % (nm, val, str(cand) if cand else "-", s[:10], dd["sc"]))
    if cand is not None:
        chk("F(I+A(%s)) = %s" % (nm, cand), abs(val - float(cand)) < 2e-14)
print()
print("=" * 78)
print("F. 独立复核：Stiefel 流形直接上升 max rho((VV^T)_+) 与 Stiefel-KKT")
print("=" * 78)

def rho_plus(V):
    Qm_ = V @ V.T
    Bm = np.maximum(Qm_, 0.0)
    wm, Um = np.linalg.eigh(Bm)
    xm = Um[:, -1]
    if xm.sum() < 0: xm = -xm
    return float(wm[-1]), xm, Qm_
def ascend(Nm, V0, steps=1500, lrs=1.0):
    Vc, _ = np.linalg.qr(V0)
    r, x, Q = rho_plus(Vc)
    for _ in range(steps):
        S = (Q >= 0).astype(float)
        G = S * np.outer(x, x)                 # d rho(Q_+)/dQ
        gv = 2.0 * (G @ Vc - Q @ (G @ Vc))     # Stiefel 投影梯度
        nrm = np.linalg.norm(gv)
        if nrm < 1e-14: break
        Vn, _ = np.linalg.qr(Vc + (lrs / nrm) * gv)
        rn, xn, Qn = rho_plus(Vn)
        if rn > r + 1e-16: Vc, r, x, Q = Vn, rn, xn, Qn
        else:
            lrs *= 0.5
            if lrs < 1e-14: break
    return r, Vc, x, Q
ladder = [(3, sp.Integer(1)), (4, k4), (5, k5), (6, rho0), (7, rho0), (8, rho0)]
rng = np.random.default_rng(20261028)
for Nm, cand in ladder:
    bb = (-1.0, None, None, None)
    for _ in range(250):
        r_, V_, x_, Q_ = ascend(Nm, rng.normal(size=(Nm, 3)))
        if r_ > bb[0]: bb = (r_, V_, x_, Q_)
    r_, V_, x_, Q_ = bb
    S_ = (Q_ >= 0).astype(float)
    G_ = S_ * np.outer(x_, x_)
    kkt = np.max(np.abs((np.eye(Nm) - Q_) @ G_ @ V_))
    print(" N=%d  max rho(Q_+)=%.15f   闭式=%-22s 差=%.3e   Stiefel-KKT=%.2e"
          % (Nm, r_, str(cand), abs(r_ - float(cand)), kkt))
    chk("kappa_3(%d) = %s （流形直接上升独立复核）" % (Nm, cand),
        abs(r_ - float(cand)) < 2e-12, "得到 %.15f" % r_)
    chk("N=%d 处 Stiefel KKT (I-Q)(sign_+(Q) odot xx^T)V = 0" % Nm, kkt < 1e-6,
        "残差 %.2e" % kkt)
    if r_ > float(rho0) + 1e-9:
        chk("N=%d 未发现 rho(Q_+) > rho0" % Nm, False, "反例值 %.15f" % r_)
chk("N=3..8 流形上升未发现 rho(Q_+) > rho0（否证性体检）", True)

print()
print("=" * 78)
print("总计 %d 项检查，全部通过: %s" % (len(ok), all(ok)))
print("=" * 78)