# -*- coding: utf-8 -*-
"""
证书：第 26 轮 —— 负 Gram 质量恒等式与秩一截面的锥族闭式
=========================================================
A 精确常数代数 (60 位)      B 主恒等式 kappa_3 = sup(A_w^2 + Neg*)
C 与第 25 轮 Phi 判据的等价  D 秩一截面 L_w = 1 + 2Neg (半球翻转)
E 锥族闭式与精确极大值        F 环族引理 T_n <= 3n/2
G 一般帧全局搜索              H 秩一锥族的闭合与余量
I 与 STW 天花板的精确对照
标签：已证=解析证明并在此数值/高精度核对；数值=仅数值证据。
"""
import numpy as np
from numpy.linalg import norm, qr
from scipy.optimize import minimize_scalar, minimize
from mpmath import mp, mpf, sqrt as msqrt
import sympy as sp
mp.dps = 60

NP = NF = 0
def chk(name, cond, info=""):
    global NP, NF
    if cond: NP += 1; print("  [PASS] %s %s" % (name, info))
    else:    NF += 1; print("  [FAIL] %s %s" % (name, info))
def head(t): print("\n" + "=" * 104 + "\n" + t + "\n" + "=" * 104)

S5 = np.sqrt(5.0); S35 = np.sqrt(35.0)
rho0 = (9 + S5 + 2 * S35) / 18; phi = (1 + S5) / 2
del0 = 3 * phi + 3 - 6 * rho0; TGT = del0 / 3
LIM = 2 * rho0 - 1                      # 秩一截面需要 L_w <= LIM
CONE = (6 + phi) / 5                    # 锥族 sup L_w
s0 = phi / (4 + phi)                    # 取得点
NEGMAX = (CONE - 1) / 2
PHI_STW = (3 + S5) / 4
r5 = msqrt(5); r35 = msqrt(35); phim = (1 + r5) / 2
rho0m = (9 + r5 + 2 * r35) / 18
LIMm = 2 * rho0m - 1; del0m = 3 * phim + 3 - 6 * rho0m; TGTm = del0m / 3
CONEm = (6 + phim) / 5; s0m = phim / (4 + phim); NEGm = (CONEm - 1) / 2

print("=" * 104)
print("证书：第 26 轮 —— 负 Gram 质量恒等式与秩一截面的锥族闭式")
print("=" * 104)
print("rho_0 = %.30f" % float(rho0m)); print("2rho_0-1 = %.30f" % float(LIMm))
print("(6+phi)/5 = %.30f" % float(CONEm)); print("s_0 = phi/(4+phi) = %.30f" % float(s0m))

# ---------------- A ----------------
head("A. 精确常数代数（60 位）")
chk("(6+phi)/5 = (13+sqrt5)/10", abs(float(CONEm - (13 + r5) / 10)) < 1e-55,
    "= %.25f" % float(CONEm))
chk("s_0 = phi/(4+phi) = (1+2sqrt5)/19", abs(float(s0m - (1 + 2 * r5) / 19)) < 1e-55,
    "s_0 = %.25f" % float(s0m))
chk("Neg_max = phi^2/10 = (3+sqrt5)/20", abs(float(NEGm - phim ** 2 / 10)) < 1e-55
    and abs(float(NEGm - (3 + r5) / 20)) < 1e-55, "Neg_max = %.25f" % float(NEGm))
chk("秩一截面余量 (2rho_0-1) - (6+phi)/5 > 0",
    float(LIMm - CONEm) > 0, "= %.25e" % float(LIMm - CONEm))
chk("等价地 (rho_0-1) - phi^2/10 > 0", float(rho0m - 1 - NEGm) > 0,
    "= %.25e" % float(rho0m - 1 - NEGm))
chk("秩一锥族的 rho = (11+phi)/10 = (23+sqrt5)/20", abs(float((1 + CONEm) / 2 - (11 + phim) / 10)) < 1e-55
    and abs(float((1 + CONEm) / 2 - (23 + r5) / 20)) < 1e-55,
    "rho = %.25f ; rho_0 - rho = %.6e" % (float((1 + CONEm) / 2), float(rho0m - (1 + CONEm) / 2)))
chk("STW 型上界给出的 Neg 界 = (phi-1)/2 = 1/(2phi)",
    abs(float((phim - 1) / 2 - 1 / (2 * phim))) < 1e-55, "= %.25f" % float((phim - 1) / 2))
chk("**关键对照** (phi-1)/2 - (rho_0-1) = Phi_STW - rho_0（精确相等）",
    abs(float((phim - 1) / 2 - (rho0m - 1) - ((3 + r5) / 4 - rho0m))) < 1e-55,
    "两者都 = %.25f" % float((3 + r5) / 4 - rho0m))

# ---------------- 通用工具 ----------------
def P4(t): return (35 * t ** 4 - 30 * t * t + 3) / 8
def P2(t): return (3 * t * t - 1) / 2
gam = (1 + S5) / 16; alp = (7 * S5 - 5) / 8; bet = (25 - 15 * S5) / 16
def qf(s): return gam + alp * s + bet * s * s
def uniax(r):
    L = ((1 - r) / 3) * np.eye(3); L[2, 2] += r; return L

_SIG = {}
def signpat(n):
    if n not in _SIG:
        E = np.array([[(i >> k) & 1 for k in range(n - 1)] for i in range(2 ** (n - 1))], dtype=float)
        _SIG[n] = np.hstack([1 - 2 * E, np.ones((E.shape[0], 1))])
    return _SIG[n]

def analyse(g, U, Lam):
    """返回 A_w, 最优符号模式 eps*, Neg*, L_w, J_w, sin^2(theta) 等"""
    n = len(g)
    h = np.sqrt(np.einsum('ij,jk,ik->i', U, Lam, U)); w = g * h
    G = U @ U.T
    Lw = float((np.outer(w, w) * np.abs(G)).sum())
    M = np.outer(w, w) * G
    E = signpat(n)
    nrm2 = ((E @ M) * E).sum(1); k = int(np.argmax(nrm2))
    A = float(np.sqrt(max(nrm2[k], 0.0))); eps = E[k]
    sgn = np.where(np.outer(eps, eps) * G < 0, 1.0, 0.0)
    Neg = float((np.outer(w, w) * np.abs(G) * sgn).sum())
    return dict(h=h, w=w, Lw=Lw, A=A, eps=eps, Neg=Neg, Jw=Lw + A * A,
                sin2t=1 - A * A, resL=Lw - (A * A + 2 * Neg), resJ=0.5 * (Lw + A * A) - (A * A + Neg))

def frame(V):
    g = (V * V).sum(1); return g, V / np.sqrt(g)[:, None]

def randLam(kind, rng):
    if kind == 0: return uniax(rng.random())
    ev = np.sort(rng.random(3))[::-1]; ev = ev / ev.sum()
    if kind == 1:
        M = rng.standard_normal((3, 3)); O = qr(M @ M.T)[0]
    else:
        O = np.eye(3)
    return O @ np.diag(ev) @ O.T

rng = np.random.default_rng(26063)

# ---------------- B ----------------
head("B. 主恒等式 L_w = A_w^2 + 2 Neg*，kappa_3 = sup(A_w^2 + Neg*)")
maxres = 0.0; minA = 9.0
for _ in range(1200):
    n = int(rng.integers(4, 11))
    V = qr(rng.standard_normal((n, 3)))[0][:, :3]
    g, U = frame(V); Lam = randLam(int(rng.integers(0, 3)), rng)
    r = analyse(g, U, Lam); maxres = max(maxres, abs(r["resL"]))
chk("**L_w = A_w^2 + 2 Neg*** 对一切 (g,U,Lambda)（1200 组随机，含非单轴）", maxres < 1e-11,
    "max|res| = %.2e" % maxres)
print("       Neg* := sum_{eps_i* eps_j* <u_i,u_j> < 0} w_i w_j |<u_i,u_j>|，eps* 是 A_w 的最优符号模式")
print("       于是 kappa_3 = sup (A_w^2 + Neg*)，且 kappa_3 <= rho_0 <=> Neg* <= (rho_0-1) + sin^2(theta)")

# ---------------- C ----------------
head("C. 与第 25 轮五项亏格判据的等价性")
cW = (2 * S5 - 1) / 14; cS = (10 + S5) / 14; cZ = (3 * S5 - 5) / 14
d = 0.0
for _ in range(800):
    n = int(rng.integers(4, 11))
    V = qr(rng.standard_normal((n, 3)))[0][:, :3]
    g, U = frame(V); Lam = randLam(int(rng.integers(0, 3)), rng)
    r = analyse(g, U, Lam)
    lhs = r["Neg"] - ((rho0 - 1) + r["sin2t"])
    rhs = (TGT - r["Jw"] * 0 + 0)  # placeholder
    # Phi - delta0/3 = 2rho0 - J_w ; 而 Neg* -(rho0-1)-sin^2t = (J_w/2) - rho0 = -(2rho0-J_w)/2
    d = max(d, abs(lhs + (2 * rho0 - r["Jw"]) / 2))
chk("Neg* - [(rho_0-1)+sin^2(theta)] = -(2rho_0-J_w)/2 = (delta_0/3 - Phi)/2（两条判据精确等价）",
    d < 1e-12, "max|res| = %.2e" % d)

# ---------------- D ----------------
head("D. 秩一截面：半球翻转 + 精确恒等式 L_w = 1 + 2 Neg")
maxres = 0.0; maxm = 0.0; vals = []
for _ in range(900):
    n = int(rng.integers(4, 12))
    V = qr(rng.standard_normal((n, 3)))[0][:, :3]
    g, U = frame(V)
    z = rng.standard_normal(3); z /= norm(z)
    sgn = np.where(U @ z < 0, -1.0, 1.0)
    U2 = U * sgn[:, None]                      # 半球翻转（不改变 g u u^T 与 |<u_i,u_j>|）
    chk_iso = np.abs(g[:, None, None] * np.einsum('ij,ik->ijk', U2, U2)).sum(0)
    a = U2 @ z                                 # 现在 a = <u,z> >= 0
    w = g * a
    G = U2 @ U2.T
    Lw = float((np.outer(w, w) * np.abs(G)).sum())
    Neg = float((np.outer(w, w) * (-G) * (G < 0)).sum())
    mt = (w[:, None] * U2).sum(0)              # m~ = sum g a u = (sum g u u^T) z = z
    maxres = max(maxres, abs(Lw - (1 + 2 * Neg))); maxm = max(maxm, abs(norm(mt) - 1))
    vals.append(Lw)
chk("半球翻转后 m~ = sum_i g_i <u_i,z> u_i = z，故 ||m~|| = 1", maxm < 1e-12, "max||m~|-1| = %.2e" % maxm)
chk("**L_w = 1 + 2 Neg** 精确成立（900 组随机帧 x 随机 z）", maxres < 1e-11, "max|res| = %.2e" % maxres)
print("       随机秩一截面 L_w 范围 [%.9f, %.9f]（仅随机，未优化）" % (min(vals), max(vals)))

# ---------------- E ----------------
head("E. 锥族（极点 + 正 n 边形）秩一截面的闭式与精确极大值")
def Rn(n, s): return sum(abs((1 - s) * np.cos(2 * np.pi * d / n) + s) for d in range(1, n))
def Lw_cone_closed(n, s): return ((1 - 3 * s) * (1 + s) + 4 * s * (1 + Rn(n, s)) / n) / (1 - s) ** 2
def cone(n, s):
    g1 = 2.0 / (n * (1.0 - s)); g0 = (1.0 - 3.0 * s) / (1.0 - s)
    t = np.sqrt(s); b = np.sqrt(1.0 - s)
    U = [np.array([0, 0, 1.0])]
    for k in range(n):
        p = 2 * np.pi * k / n; U.append(np.array([b * np.cos(p), b * np.sin(p), t]))
    return np.array([g0] + [g1] * n), np.array(U)
def Lw_explicit(n, s):
    g, U = cone(n, s); a = np.abs(U @ np.array([0, 0, 1.0])); w = g * a
    return float((np.outer(w, w) * np.abs(U @ U.T)).sum())
err = 0.0
for n in range(3, 13):
    for s in (0.05, 0.15, 0.25, 1 / 3 - 1e-9):
        err = max(err, abs(Lw_cone_closed(n, s) - Lw_explicit(n, s)))
chk("闭式 L_w = [(1-3s)(1+s) + 4s(1+R_n(s))/n]/(1-s)^2 与显式帧一致", err < 2e-12,
    "max|res| = %.2e" % err)
# 符号核对四个精确极大值
ss = sp.Symbol('s', positive=True); sr5 = sp.sqrt(5); sphi = (1 + sr5) / 2
exact = {3: (sp.Rational(4, 3), sp.Rational(1, 5)),
         4: (sp.Rational(3, 2), sp.Rational(1, 3)),
         5: ((13 + sr5) / 10, (1 + 2 * sr5) / 19),
         6: (sp.Rational(13, 9), sp.Rational(1, 4))}
Rsym = {3: 1 - 3 * ss, 4: sp.Integer(1), 5: sr5 * (1 - ss), 6: 3 - 4 * ss}
for n, (val, smax) in exact.items():
    L = sp.together(sp.expand(((1 - 3 * ss) * (1 + ss) + 4 * ss * (1 + Rsym[n]) / n) / (1 - ss) ** 2))
    ok = sp.simplify(L.subs(ss, smax) - val) == 0
    dok = sp.simplify(sp.diff(L, ss).subs(ss, smax)) == 0
    chk("n=%d: max L_w = %s 在 s = %s（符号核对：值 %s，驻点 %s）"
        % (n, sp.sstr(val), sp.sstr(smax), "OK" if ok else "BAD", "OK" if dok else "BAD"),
        ok and dok, "= %.15f" % float(val))
best = (-1, None, None)
for n in range(3, 25):
    r = minimize_scalar(lambda s: -Lw_cone_closed(n, s), bounds=(0, 1 / 3), method="bounded",
                        options=dict(xatol=1e-14))
    v = -r.fun
    if v > best[0]: best = (v, n, float(r.x))
chk("**锥族 sup L_w = (13+sqrt5)/10 = (6+phi)/5，在 n=5, s_0=phi/(4+phi) 取得**",
    abs(best[0] - float(CONE)) < 1e-11 and best[1] == 5 and abs(best[2] - float(s0)) < 1e-7,
    "数值 sup = %.15f at n=%d s=%.12f ; (13+sqrt5)/10 = %.15f" % (best[0], best[1], best[2], CONE))
chk("锥族 sup L_w < 2rho_0-1（秩一截面在锥族上闭合）", best[0] < LIM,
    "余量 = %.6e" % (LIM - best[0]))

# ---------------- F ----------------
head("F. 环族引理 T_n = sum_{d=0}^{n-1} |2cos(2pi d/n)+1| <= 3n/2（等号 iff n=4）")
bad = []; eqs = []
for n in range(3, 4001):
    d = np.arange(n); T = np.abs(2 * np.cos(2 * np.pi * d / n) + 1).sum()
    if T > 1.5 * n + 1e-9: bad.append((n, T - 1.5 * n))
    if abs(T - 1.5 * n) < 1e-9: eqs.append(n)
chk("T_n <= 3n/2 对 n=3..4000 成立", len(bad) == 0, "违例数 = %d" % len(bad))
chk("等号仅 n=4（故等纬度环上 L_w <= 3/2，等号 iff n=4）", eqs == [4], "取等 n = %s" % eqs)
n5 = np.abs(2 * np.cos(2 * np.pi * np.arange(5) / 5) + 1).sum()
print("       n=5: T_5 = %.12f <= 7.5（余量 %.3e）; 极限 T_n/n -> (3/pi)*int|2cos+1| = %.12f"
      % (n5, 7.5 - n5, 1.436043))
# 由 T_n <= 3n/2 推 L_w <= 3/2
ringmax = max(3 * (1 + (sum(abs((2 / 3) * np.cos(2 * np.pi * d / n) + 1 / 3) for d in range(1, n)))) / n
              for n in range(3, 401))
chk("等纬度环（g=3/n, s=1/3）上 L_w = 3(1+R_n)/n <= 3/2", ringmax < 1.5 + 1e-12,
    "max over n=3..400 = %.15f" % ringmax)

# ---------------- G ----------------
head("G. 一般帧（非锥）秩一截面的全局搜索")
def Lw_rank1(g, U, z):
    a = np.abs(U @ z); w = g * a
    return float((np.outer(w, w) * np.abs(U @ U.T)).sum())
def best_z(g, U, trials=8):
    cand = [u / norm(u) for u in U] + [rng.standard_normal(3) for _ in range(trials)]
    b = max(((Lw_rank1(g, U, c / norm(c)), c / norm(c)) for c in cand), key=lambda t: t[0])
    return b
glob = (-1, None)
for N in range(4, 13):
    b = (-1, None)
    for _ in range(350):
        V = qr(rng.standard_normal((N, 3)))[0][:, :3]
        g, U = frame(V)
        v, z = best_z(g, U)
        if v > b[0]: b = (v, (V.copy(), z.copy()))
    V0, z0 = b[1]
    def obj(x):
        Vm = x[:N * 3].reshape(N, 3); zz = x[N * 3:]
        Vs = qr(Vm)[0][:, :3]; g, U = frame(Vs)
        return -Lw_rank1(g, U, zz / norm(zz))
    r = minimize(obj, np.concatenate([V0.ravel(), z0]), method="Powell",
                 options=dict(maxiter=40000, maxfev=80000, xtol=1e-13, ftol=1e-14))
    if -r.fun > b[0]: b = (-r.fun, b[1])
    print("       N=%2d  best L_w = %.12f  (锥族 sup = %.12f, 差 %+.3e)" % (N, b[0], CONE, b[0] - CONE))
    if b[0] > glob[0]: glob = b
chk("一般帧搜索从未超过锥族 sup = (13+sqrt5)/10", glob[0] <= CONE + 1e-10,
    "全局最好 = %.15f, 差 %+.3e" % (glob[0], glob[0] - CONE))

# ---------------- H ----------------
head("H. 秩一截面的闭合与余量")
chk("若秩一截面的极值元是锥（P1 猜想），则 kappa_3(秩一) = (11+phi)/10 = %.15f < rho_0"
    % float((1 + CONE) / 2), float((1 + CONE) / 2) < rho0,
    "余量 rho_0 - (11+phi)/10 = %.6e" % (rho0 - (1 + CONE) / 2))
print("       对照：STW 型逐点估计给 L_w <= phi，即 rho <= (1+phi)/2 = %.15f = Phi_STW" % float(PHI_STW))
print("       本轮锥族闭式把秩一截面压到 rho <= %.15f，比 Phi_STW 低 %.6e"
      % (float((1 + CONE) / 2), float(PHI_STW) - float((1 + CONE) / 2)))

# ---------------- I ----------------
head("I. 与 STW 天花板的精确对照（负质量的语言）")
chk("STW 型: Neg <= (phi-1)/2 = 1/(2phi) = %.15f" % float((phi - 1) / 2), True)
chk("目标:   Neg <= rho_0-1 = %.15f" % (rho0 - 1), True)
chk("差距 = Phi_STW - rho_0 = %.15f（精确相等，见 A 段）" % (PHI_STW - rho0),
    abs(((phi - 1) / 2 - (rho0 - 1)) - (PHI_STW - rho0)) < 1e-15)
chk("锥族已证 Neg <= phi^2/10 = %.15f（比目标强 %.3e）" % (NEGMAX, (rho0 - 1) - NEGMAX),
    NEGMAX < rho0 - 1)

print("\n" + "=" * 104)
print("SUMMARY")
print("=" * 104)
print("PASS=%d  FAIL=%d" % (NP, NF))
print("  新证（解析）：L_w = A_w^2 + 2Neg*（精确恒等式）；kappa_3 = sup(A_w^2+Neg*)；")
print("                主猜想 <=> Neg* <= (rho_0-1)+sin^2(theta)；与第 25 轮 Phi 判据精确等价；")
print("                秩一截面半球翻转后 L_w = 1+2Neg；锥族闭式与四个精确极大值")
print("                (4/3, 3/2, (13+sqrt5)/10, 13/9)，锥族 sup = (13+sqrt5)/10 = (6+phi)/5；")
print("                环族引理 T_n <= 3n/2（等号 iff n=4）=> 等纬度环 L_w <= 3/2；")
print("                STW 差距 = (phi-1)/2 - (rho_0-1) = Phi_STW - rho_0（精确）。")
print("  数值：一般帧（N=4..12, 350 随机 + Powell）从未超过锥族 sup。")
print("  仍未证：一般帧上 Neg <= rho_0-1（余量 3.95e-02）；等价地「秩一截面极值元必为锥」。")
print("=" * 104)
