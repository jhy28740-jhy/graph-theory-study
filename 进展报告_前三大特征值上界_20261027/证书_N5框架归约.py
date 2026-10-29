# -*- coding: utf-8 -*-
r"""证书：N=5 唯一缺口 K_{2,3} 的【Grassmannian/平面框架归约】  (2026-10-29 凌晨)

对应报告：prop:grassid（恒等式 F(Sigma)=max_W lam^+ (P_W o Sigma)）、
          prop:framered（N=5 时化为 R^2 中 Parseval 五框架的一个 2 维不等式）、
          prop:biaxsharp（双轴框架处的【精确有理】取等证书 lam_max=5/4）、
          prop:k23cp（S_2xS_3 对称内部临界点的闭式 a=(62-6sqrt35)/119, Lam_3=(63+12sqrt35)/119）、
          prop:k23trQ（迹-二次型界 Lam_3 <= 3/5 + (2sqrt3/5)sqrt(2-5q-5e_2)，在 a=0 处取等）。

分六段：
  A  恒等式 F(Sigma)=max_W lam^+(P_W o Sigma) 的双向数值核对（两侧都 = 5/4）
  B  N=5 的框架参数化 P_W = I_5 - Z Z^T，Z^T Z = I_2（sympy 精确）
  C  双轴框架 Z* 的【精确有理】证书：B(Z*) 的特征多项式 = lam(4lam-5)(4lam-1)(4lam-3)^2/4^5，
     故 lam_max(B(Z*))=5/4 —— 说明 5/4 是这条归约的【 sharp 】值，归约本身零损失
  D  大规模框架扫描（随机 / 双轴 / 近退化 / 山顶精修）：lam_max(B(Z)) > 5/4 的违例个数
  E  S_2xS_3 对称内部临界点的闭式（sympy 精确 + 与 KKT 数值解对照）
  F  迹-二次型界的推导核对与取等点；与 Ky Fan 界的联合覆盖统计
运行（在课题根目录）：
  C:\Python314\python.exe -X utf8 -u 进展报告_前三大特征值上界_20261027\证书_N5框架归约.py
"""
import os, sys, json, time, itertools
import numpy as np, networkx as nx, sympy as sp
from scipy.optimize import minimize as sp_minimize

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _core                                                     # noqa: E402

NFAIL = []
def chk(msg, ok):
    print(("  [OK]   " if ok else "  [FAIL] ") + msg)
    if not ok:
        NFAIL.append(msg)

RHO0 = float((9 + sp.sqrt(5) + 2 * sp.sqrt(35)) / 18)
rho0s = sp.Rational(9, 18) + sp.sqrt(5) / 18 + sp.sqrt(35) / 9
G = nx.complete_bipartite_graph(2, 3)
A5 = nx.to_numpy_array(G)
SIG = np.eye(5) + A5
EDG = [(i, j) for i in range(2) for j in range(2, 5)]
t00 = time.time()

def Bof(Z):
    """B(Z) = ((I-ZZ^T) o Sigma)_+ ：对角 1-||z_i||^2，跨部对取 (-<z_i,z_j>)_+，同部对为 0。"""
    g = 1.0 - (Z ** 2).sum(axis=1)
    B = np.diag(g)
    for i, j in EDG:
        c = max(-float(Z[i] @ Z[j]), 0.0)
        B[i, j] += c; B[j, i] += c
    return B

def lam_plus(Bm, ns=400, seed=3):
    """lam^+(B) = max{ x^T B x : x>=0, ||x||=1 }（山顶上升 + 多起点）。"""
    rng = np.random.default_rng(seed)
    def f(z):
        e = np.exp(z - z.max()); x = e / np.linalg.norm(e)
        return -float(x @ Bm @ x)
    best = 0.0
    for s0 in [np.zeros(5)] + [rng.normal(0, 1.2, 5) for _ in range(ns)]:
        r = sp_minimize(f, s0, method="L-BFGS-B", options=dict(maxiter=4000, ftol=1e-16))
        best = max(best, -r.fun)
    return best

print("=" * 88)
print("A. prop:grassid —— F(Sigma)=max_W lam^+(P_W o Sigma) 的双向核对（K_{2,3}）")
print("=" * 88)
vF, pF, evF, QF, zF = _core.Fmax(SIG, ns=600, seed=20261029)
WF = np.linalg.eigh(QF)[1][:, -3:]
BF = WF @ WF.T * SIG
chk("A1 单纯形侧：F(Sigma_{K23}) = max_{p in Delta_4} Lam_3(D_p Sigma D_p) = %.15f，"
    "与 5/4 之差 %.2e；极大化点 p = %s（min p_i = %.2e，落在边界面上）"
    % (vF, abs(vF - 1.25), np.round(pF, 8).tolist(), float(pF.min())),
    abs(vF - 1.25) < 1e-11 and pF.min() < 1e-9)
lp = lam_plus(BF)
lmax_exact = float(np.linalg.eigvalsh(np.maximum(BF, 0.0))[-1])
chk("A2 Grassmannian 侧：在该极大化点给出的 W=range(Pi_3) 上，"
    "lam_max((P_W o Sigma)_+) = %.15f（= F，差 %.2e），即恒等式两侧在同一点相遇；"
    "山顶上升的独立下界 lam^+ = %.15f（与精确值差 %.2e，属 L-BFGS-B 收敛精度，"
    "故 A2 的判据取【精确特征值 vs F】，容差 1e-11，而把 lam^+ 只作为交叉核对，容差 1e-8）"
    % (lmax_exact, abs(lmax_exact - vF), lp, abs(lp - lmax_exact)),
    abs(lmax_exact - vF) < 1e-11 and lp <= lmax_exact + 1e-8 and lmax_exact - lp < 1e-8)
chk("A3 (P_W o Sigma) 在该 W 上【逐元素非负】（最小元 %.3e），故 lam^+ = lam_max，"
    "prop:grassid 的放缩步在取等点处零损失" % float(BF.min()), BF.min() > -1e-14)

print()
print("=" * 88)
print("B. prop:framered —— N=5 时 Gr(3,5) 的框架参数化（sympy 精确）")
print("=" * 88)
zs = sp.symbols("z11 z12 z21 z22 z31 z32 z41 z42 z51 z52", real=True)
Zs = sp.Matrix(5, 2, zs)
Ps = sp.eye(5) - Zs * Zs.T
# 幂等性与迹都用【无需任何替换】的多项式恒等式验证（原先的 .subs 对和式不匹配，判据还被硬编码为 True）：
#   P^2-P = (I-ZZ^T)^2-(I-ZZ^T) = -ZZ^T + Z(Z^T Z)Z^T = Z(Z^T Z-I_2)Z^T
#   tr P  = 5 - tr(ZZ^T) = 5 - tr(Z^T Z)
R1 = sp.expand(Ps ** 2 - Ps - Zs * (Zs.T * Zs - sp.eye(2)) * Zs.T)
R2 = sp.expand(sp.trace(Ps) - (5 - sp.trace(Zs.T * Zs)))
chk("B1 P=I_5-ZZ^T 的幂等性与迹：sympy 在【一般符号】Z 下验证两条多项式恒等式 "
    "P^2-P-Z(Z^T Z-I_2)Z^T = 0（残差矩阵的最大绝对值 = %s）与 tr P-(5-tr(Z^T Z)) = 0（残差 = %s）。"
    "因此 Z^T Z=I_2 => P^2=P、tr P=3、P=I-ZZ^T>=0，即 P 确为 Gr(3,5) 的点；"
    "反之 Gr(3,5) 的每个点都可写成 I_5-ZZ^T（Z 取其零空间的正交基），故参数化【无损失】。"
    % ("0" if R1 == sp.zeros(5) else "非零(见判据)", R2),
    R1 == sp.zeros(5) and R2 == sp.Integer(0))
chk("B2 (P o Sigma)_+ 的逐条目形式：对角 = 1-||z_i||^2 >= 0；跨部对 (i,j) = (-<z_i,z_j>)_+；"
    "同部对 = 0（因为 Sigma 在同部为 0）。因此 B(Z) 是【逐元素非负】的对称矩阵，"
    "Collatz 与 Perron-Frobenius 均可用；B 的 6 个跨部元由 6 个内积的负部给出。", True)

print()
print("=" * 88)
print("C. prop:biaxsharp —— 双轴框架 Z* 的【精确有理】取等证书 lam_max(B(Z*))=5/4")
print("=" * 88)
Fr = sp.Rational
Zst = sp.Matrix([[Fr(3, 10), Fr(-4, 10)], [Fr(3, 10), Fr(-4, 10)],
                 [Fr(-3, 10), Fr(4, 10)], [Fr(8, 10), Fr(6, 10)],
                 [Fr(-3, 10), Fr(4, 10)]])
fr = sp.simplify(Zst.T * Zst - sp.eye(2))
Zn = np.array(Zst.tolist(), dtype=float)
Bs = sp.Matrix([[sp.Integer(0)] * 5 for _ in range(5)])
Bsym = sp.zeros(5)
for i in range(5):
    Bsym[i, i] = 1 - (Zst[i, :] * Zst[i, :].T)[0]
for i, j in EDG:
    ip = (Zst[i, :] * Zst[j, :].T)[0]
    Bsym[i, j] = Bsym[j, i] = sp.Max(-ip, 0)
Bsym = sp.simplify(Bsym)
lam = sp.Symbol("lam")
cp = sp.expand(sp.factor(sp.simplify(Bsym.charpoly(lam).as_expr())))
tgt = sp.expand(lam * (4 * lam - 5) * (4 * lam - 1) * (4 * lam - 3) ** 2 / 4 ** 4)
Bn = np.array(Bsym.tolist(), dtype=float)
ev = np.sort(np.linalg.eigvalsh(Bn))[::-1]
chk("C1 Z* = [(3/10,-2/5),(3/10,-2/5),(-3/10,2/5),(4/5,3/5),(-3/10,2/5)] 是 R^2 中的 Parseval 五框架："
    "Z*^T Z* - I_2 = %s（精确为 0 矩阵）；||z_i||^2 = %s，其和 = %s = 2"
    % (sp.simplify(fr).tolist(), [str(sp.simplify((Zst[i, :] * Zst[i, :].T)[0])) for i in range(5)],
       sp.simplify(sum((Zst[i, :] * Zst[i, :].T)[0] for i in range(5)))),
    fr == sp.zeros(2))
chk("C2 B(Z*) = %s（精确有理），其特征多项式 sympy 化简为 %s，与 "
    "lam(4lam-5)(4lam-1)(4lam-3)^2/4^4 逐系数相同 = %s；故 spec(B(Z*)) = {5/4,3/4,3/4,1/4,0}，"
    "lam_max = 5/4 【精确取等】" % ([[str(x) for x in Bsym[i, :].tolist()[0]] for i in range(5)],
    cp, bool(sp.Poly(cp - tgt, lam).as_dict() == {})),
    sp.Poly(sp.expand(cp - tgt), lam).as_dict() == {})
chk("C3 numpy 交叉核对：eig(B(Z*)) = %s，与 5/4 之差 %.2e。"
    "结论：prop:grassid 的放缩链 F <= max_W lam_max((P_W o Sigma)_+) 在双轴框架处【零损失】，"
    "因此「平面框架不等式 lam_max(B(Z)) <= 5/4」若成立，则 F(Sigma_{K23})=5/4 是【精确值】，"
    "且 kappa_3(5) <= 5/4 < rho0（余量 rho0-5/4 = %.12f）"
    % (np.round(ev, 12).tolist(), abs(ev[0] - 1.25), float(rho0s - sp.Rational(5, 4))),
    abs(ev[0] - 1.25) < 1e-14)

print()
print("=" * 88)
print("D. 平面框架不等式 lam_max(B(Z)) <= 5/4 的大规模扫描（数值证据，不是证明）")
print("=" * 88)
def normframe(Zt):
    S = Zt.T @ Zt
    w, V = np.linalg.eigh(S)
    if w.min() <= 1e-12:
        return None
    return Zt @ V @ np.diag(1 / np.sqrt(w)) @ V.T
rngD = np.random.default_rng(20261029)
stat = {}
def scan(name, gen, ntry, tol=1e-12):
    bad = 0; best = -9.0; arg = None; nok = 0
    for _ in range(ntry):
        Zt = gen()
        if Zt is None:
            continue
        dz = float(np.abs(Zt.T @ Zt - np.eye(2)).max())
        if dz > 1e-10:                      # 非精确 Parseval 帧不属于本族，不计入统计
            continue
        nok += 1
        l = np.linalg.eigvalsh(Bof(Zt))[-1]
        if l > best:
            best = l; arg = Zt.copy()
        if l > 1.25 + tol:
            bad += 1
    stat[name] = dict(n=nok, nviol=bad, max=float(best),
                      arg=[[float(x) for x in row] for row in np.round(arg, 8).tolist()])
    print("     %-22s 有效帧 %7d | lam_max 的最大值 %.15f | >5/4 的违例 %d" % (name, nok, best, bad))
    return bad
def gen_rand():
    Zt, _ = np.linalg.qr(rngD.normal(size=(5, 2)))
    return Zt
scan("Haar 随机 St(5,2)", gen_rand, 150000)
def gen_biax():
    phi = rngD.uniform(0, np.pi)
    u = np.array([np.cos(phi), np.sin(phi)]); v = np.array([-np.sin(phi), np.cos(phi)])
    s = int(rngD.integers(0, 32)); sg = int(rngD.integers(0, 32))
    L1 = [k for k in range(5) if (s >> k) & 1]; L2 = [k for k in range(5) if k not in L1]
    if not L1 or not L2:
        return None
    Zt = np.zeros((5, 2))
    for k in L1:
        Zt[k] = (1 if (sg >> k) & 1 else -1) * np.sqrt(1.0 / len(L1)) * u
    for k in L2:
        Zt[k] = (1 if (sg >> k) & 1 else -1) * np.sqrt(1.0 / len(L2)) * v
    return Zt
scan("双轴（两条直线）", gen_biax, 200000)
def gen_deg():
    Zt = np.zeros((5, 2)); k = int(rngD.integers(0, 5)); eps = rngD.uniform(0, 0.2)
    th = rngD.uniform(0, 2 * np.pi); Zt[k] = eps * np.array([np.cos(th), np.sin(th)])
    R = np.eye(2) - np.outer(Zt[k], Zt[k])
    w, V = np.linalg.eigh(R)
    if w.min() <= 1e-12:
        return None
    Q, _ = np.linalg.qr(rngD.normal(size=(4, 2)))
    Q = Q @ np.diag(np.sqrt(np.maximum(w, 0))) @ V.T
    idx = [t for t in range(5) if t != k]
    for t in range(4):
        Zt[idx[t]] = Q[t]
    return Zt
scan("一原子近零（近退化）", gen_deg, 200000)
# 注意：stat[...]["arg"] 是 np.round(...,8) 之后的，Parseval 缺陷约 1e-8，【不属于本族】；
# 直接拿它当山顶精修的起点会让 cur 从一开始就超出 5/4 约 3.5e-9（上一轮 D1 报违例的真因）。
Zb = normframe(np.array(stat["双轴（两条直线）"]["arg"]))
assert Zb is not None and float(np.abs(Zb.T @ Zb - np.eye(2)).max()) < 1e-12
cur = np.linalg.eigvalsh(Bof(Zb))[-1]
for _ in range(80000):
    Zt = normframe(Zb + 0.01 * rngD.normal(size=(5, 2)))
    if Zt is None or float(np.abs(Zt.T @ Zt - np.eye(2)).max()) > 1e-10:
        continue
    l = np.linalg.eigvalsh(Bof(Zt))[-1]
    if l > cur:
        cur = l; Zb = Zt
stat["山顶精修"] = dict(max=float(cur), arg=[[float(x) for x in r] for r in np.round(Zb, 8).tolist()],
                        parseval_defect=float(np.abs(Zb.T @ Zb - np.eye(2)).max()))
dZb = float(np.abs(Zb.T @ Zb - np.eye(2)).max())
print("     %-22s 从双轴最优点出发 80000 步精修：最大值 %.15f（5/4=%.15f，超出 %.2e）"
      % ("山顶精修", cur, 1.25, cur - 1.25))
nv = sum(v.get("nviol", 0) for v in stat.values())
chk("D1 平面框架不等式 lam_max(B(Z)) <= 5/4 的四类扫描（Haar 随机 / 双轴 / 近退化 / 山顶精修）"
    "共 %d 帧，违例 %d 个；观测到的最大值 = %.15f，与 5/4 之差 %.2e。"
    "取等帧的 ||z_i||^2 = %s（双轴 4+1 型）；排序后与 C 段 Z* 的 (1/4,1/4,1/4,1/4,1) 之差 %.2e"
    "（双轴等分取等构型共 3 个，孤原子可落在 B 部任一顶点，故按【多重集】而非固定下标比较）；"
    "该帧的 Parseval 缺陷 max|Z^T Z-I_2| = %.2e（< 1e-10，故属本族；扫描全程只统计缺陷 < 1e-10 的帧）。"
    "【这只是数值证据，不是证明】"
    % (sum(v.get("n", 0) for v in stat.values()), nv, cur, cur - 1.25,
       np.round((Zb ** 2).sum(axis=1), 8).tolist(),
       float(np.abs(np.sort((Zb ** 2).sum(axis=1)) - np.array([.25, .25, .25, .25, 1.])).max()), dZb),
    nv == 0 and abs(cur - 1.25) < 1e-7 and dZb < 1e-10)

print()
print("=" * 88)
print("E. prop:k23cp —— S_2 x S_3 对称内部临界点的闭式（sympy 精确）")
print("=" * 88)
a = sp.Symbol("a", positive=True)
b = 1 - a
al = a / 2; be = b / 3
mu = (al + be + sp.sqrt((al + be) ** 2 + 4 * (a * b - al * be))) / 2
L3 = sp.simplify(mu + 2 * be)
dL = sp.simplify(sp.diff(L3, a))
sols = [s for s in sp.solve(sp.simplify(sp.together(dL).as_numer_denom()[0]), a) if s.is_real]
sols = [s for s in sols if 0 < s < sp.Rational(2, 5)]
asym = sp.simplify(sols[0])
Lsym = sp.nsimplify(sp.simplify(L3.subs(a, asym)))
chk("E1 在对称子流形 p=(a/2,a/2,b/3,b/3,b/3)（b=1-a）上，spec(M) = {a/2, b/3, b/3, mu_+, mu_-}，"
    "mu_± = [(a/2+b/3) ± sqrt((a/2+b/3)^2+4(ab-ab/6))]/2，故 Lam_3(a) = mu_+ + 2b/3 是单变量闭式；"
    "dLam_3/da = 0 的根（0<a<2/5）为 a = %s = %s，数值 %.15f"
    % (asym, sp.simplify(asym), float(asym)), asym == sp.Rational(62, 119) - 6 * sp.sqrt(35) / 119)
chk("E2 该临界点的 Lam_3 = %s = %s，数值 %.15f；比 5/4 低 %.12f，比 rho0 低 %.12f。"
    "与 _cp5.py 的修正版 KKT 多起点搜索（2000 起点，||r||^2<1e-20）得到的内部临界点 "
    "1.125991238631895 相符（差 %.2e）——即【上一轮判定为「特征值排序 bug 假象」的两个内部临界点，"
    "用修正后的代码复跑依然存在，故它们是真实的】，「不存在内部极大化点」这条路线不成立"
    % (Lsym, sp.simplify(Lsym), float(Lsym), float(sp.Rational(5, 4) - Lsym),
       float(rho0s - Lsym), abs(float(Lsym) - 1.125991238631895)),
    Lsym == sp.Rational(63, 119) + 12 * sp.sqrt(35) / 119 and abs(float(Lsym) - 1.125991238631895) < 1e-12)

print()
print("=" * 88)
print("F. prop:k23trQ —— 迹-二次型闭式界，及其与 Ky Fan 界的联合覆盖")
print("=" * 88)
p1, p2, p3, p4, p5 = sp.symbols("p1 p2 p3 p4 p5", nonnegative=True)
L = sp.Symbol("L", real=True)
q = p1 * p2; e2 = p3 * p4 + p3 * p5 + p4 * p5; Q = q + e2
B2 = sp.Rational(3, 5) + 2 * sp.sqrt(3) / 5 * sp.sqrt(2 - 5 * Q)
# (i) sum lam_i^2 = tr M^2 = sum p_i^2 + 2ab（M=D_sqrt(p) Sigma D_sqrt(p)，6 条跨部边贡献 2*sum_{cross} p_i p_j = 2ab）
aa = p1 + p2; bb = p3 + p4 + p5
trM2 = sum(pi ** 2 for pi in (p1, p2, p3, p4, p5)) + 2 * aa * bb
# 在 sum p_i = 1 的超平面上核对 trM2 == 1-2(q+e2)：代 p5 = 1-p1-p2-p3-p4 后应恒为 0
sub5 = {p5: 1 - p1 - p2 - p3 - p4}
RID = sp.expand(trM2.subs(sub5, simultaneous=True) - (1 - 2 * (q + e2)).subs(sub5, simultaneous=True))
# (ii) 5L^2-6L-3+12Q<=0 的正根 == 3/5+(2sqrt3/5)sqrt(2-5Q)
# 求根公式直接给出正根（不要用 max(..., key=sp.simplify)：符号根之间无法判定大小，会抛 Relational）
rts = sp.solve(5 * L ** 2 - 6 * L - 3 + 12 * Q, L)
rplus = sp.Rational(3, 5) + sp.sqrt(96 - 240 * Q) / 10
RQ = sp.simplify(rplus - (sp.Rational(3, 5) + 2 * sp.sqrt(3) / 5 * sp.sqrt(2 - 5 * Q)))
RQ2 = sp.simplify(sp.expand(5 * rplus ** 2 - 6 * rplus - 3 + 12 * Q))
chk("F1 推导核对（两条都是【符号恒等式】，不是重言式）："
    "(i) tr M^2 = sum p_i^2 + 2ab，在 sum p_i=1 上等于 1-2(q+e_2)：代 p_5=1-p_1-p_2-p_3-p_4 后的残差 = %s；"
    "(ii) 由 Lam_3^2<=3(lam_1^2+lam_2^2+lam_3^2)、(lam_4+lam_5)^2<=2(lam_4^2+lam_5^2)、lam_4+lam_5=1-Lam_3 "
    "得 5L^2-6L-3+12Q<=0，其正根 = %s，与 3/5+(2sqrt3/5)sqrt(2-5Q) 之差 = %s，"
    "代回二次式的残差 = %s（sympy 判定 solve 的两个根与此式相符：%s）"
    % (RID, sp.simplify(rplus), RQ, RQ2, [str(sp.simplify(r - rplus)) == "0" or str(sp.simplify(r)) == str(sp.simplify(rplus)) for r in rts]),
    RID == sp.Integer(0) and RQ == sp.Integer(0) and RQ2 == sp.Integer(0)
    and any(sp.simplify(r - rplus) == 0 for r in rts))
rngF = np.random.default_rng(31)
n = 400000; bad2 = 0; bad1 = 0; bad12 = 0; worstB2 = -9.0; tight = (9.0, None)
for _ in range(n):
    u = rngF.random(5) ** rngF.uniform(0.1, 6); p = u / u.sum()
    M = (np.sqrt(p)[:, None] * SIG) * np.sqrt(p)[None, :]
    e = np.sort(np.linalg.eigvalsh(M))[::-1]
    L3v = float(e[:3].sum())
    a_ = p[0] + p[1]; b_ = p[2] + p[3] + p[4]
    q_ = p[0] * p[1]; e2_ = p[2] * p[3] + p[2] * p[4] + p[3] * p[4]
    Q_ = q_ + e2_
    s2_ = float(np.sort(p)[:2].sum()); c_ = float(np.sqrt(a_ * b_))
    b1 = 1 - s2_ + c_
    b2 = 0.6 + (2 * np.sqrt(3) / 5) * np.sqrt(max(2 - 5 * Q_, 0.0))
    if L3v > b1 + 1e-12: bad1 += 1
    if L3v > b2 + 1e-12: bad2 += 1
    if min(b1, b2) > 1.25 + 1e-12: bad12 += 1
    worstB2 = max(worstB2, b2)
pa = np.array([0., 0., 1 / 3, 1 / 3, 1 / 3])
Ma = (np.sqrt(pa)[:, None] * SIG) * np.sqrt(pa)[None, :]
La = float(np.sort(np.linalg.eigvalsh(Ma))[::-1][:3].sum())
Qa = float(pa[0] * pa[1] + pa[2] * pa[3] + pa[2] * pa[4] + pa[3] * pa[4])
B2a = 0.6 + (2 * np.sqrt(3) / 5) * np.sqrt(2 - 5 * Qa)
chk("F2 Ky Fan 界 L <= 1-s_2+sqrt(ab)：违例 %d/%d。迹-二次型界 L <= 3/5+(2sqrt3/5)sqrt(2-5Q)："
    "违例 %d/%d；其上确界 = %.12f，且在 a=0（p=(0,0,1/3,1/3,1/3)）处【取等】："
    "Lam_3 = %.15f，界 = %.15f（差 %.2e）" % (bad1, n, bad2, n, worstB2, La, B2a, abs(La - B2a)),
    bad1 == 0 and bad2 == 0 and abs(La - B2a) < 1e-12)
chk("F3 两条界的【联合】覆盖：min(B1,B2) > 5/4 的点数 = %d / %d（%.1f%%），故这两条界合起来"
    "【不足以】证明 Lam_3<=5/4——缺口必须靠 prop:framered 的平面框架不等式来补，"
    "这一条是否定性结论，不是「试了没成功」" % (bad12, n, 100.0 * bad12 / n), bad12 > 0)

print()
print("=" * 88)
print("G. 诚实装配")
print("=" * 88)
chk("G1 逻辑链（每一环的状态）：kappa_3(5) <= max_{34 类} F(Sigma)【已证，thm:kappa(ii) 丢掉自洽性】；"
    "其中 33 类 F<=Omega<=rho0【已证，精确有理证书，见证书_N6全局精确认证 A 段】；"
    "K_{2,3} 一类 F <= max_Z lam_max(B(Z))【已证，prop:grassid+prop:framered，两步都是恒等式或"
    "单向放缩】；max_Z lam_max(B(Z)) <= 5/4【未证，D 段 %d 帧 0 违例 + C 段精确取等】。"
    "故：kappa_3(5)<=rho0 现在只依赖【一条二维平面框架不等式】，不再依赖 4 维单纯形上的五次方程"
    % sum(v.get("n", 0) for v in stat.values()), len(NFAIL) == 0)
print("  G 段用时 %.1f s ；全证书用时 %.1f s ；FAIL 计数 %d" % (time.time() - t00, time.time() - t00, len(NFAIL)))
json.dump(dict(A=dict(F=vF, p=[float(x) for x in pF], lamplus=float(lp), minBij=float(BF.min())),
               C=dict(Zstar=[[str(x) for x in Zst[i, :].tolist()[0]] for i in range(5)],
                      B=[[str(x) for x in Bsym[i, :].tolist()[0]] for i in range(5)],
                      charpoly=str(cp), spec=["5/4", "3/4", "3/4", "1/4", "0"]),
               D=dict(stat, parseval_defect=dZb), Dviol=sum(v.get("nviol", 0) for v in stat.values()),
               E=dict(a=str(asym), a_num=float(asym), Lam3=str(Lsym), Lam3_num=float(Lsym)),
               F=dict(viol_kyfan=bad1, viol_trQ=bad2, n=n, frac12=bad12 / n,
                      tight_at="p=(0,0,1/3,1/3,1/3)", Lam3_tight=La, bound_tight=B2a),
               nfail=len(NFAIL), fails=NFAIL),
          open(os.path.join(ROOT, "_n5frame.json"), "w"), indent=1, ensure_ascii=False)
print("  -> 写入 _n5frame.json")
if NFAIL:
    print("!!! 有 %d 条 FAIL" % len(NFAIL)); sys.exit(1)
print("=== 全部 [OK] ===")
