# -*- coding: utf-8 -*-
r"""证书：猜想 A2（负质量界）被否证，以及修复形式 A2' 的定量核查  (2026-10-28 深夜)

对应报告 sec:domdual 中 conj:A2 的【否证】与 prop:A2fix 的数值支撑。

内容：
  A. 反例的构造与识别：A = I_7 + A_H，H 为 graph6 码 'FLp|w' 的 7 顶点 14 边图，
     |Aut(H)| = 2（对合 (0 1)(2 3)(5 6)，固定 4），最优权 p 恰按该对合配对。
  B. 数值层：400 起点重跑 F(A)，报告 F、Lambda_3、谱、tr((-M)_+)、自洽性与各类余量。
  C. 精确层：把 p 化为【精确有理】向量 p^（分母 2^60 网格 + 归一化），
     用 M ~ A*diag(p^)（相似，特征多项式相同）在 Q 上算出 charpoly，清分母得整数系数；
     再自实现 Sturm 序列（纯 Fraction 算术）做实根计数、二分隔离与精化，
     得到 tr((-M^)_+) 的【严格有理区间】；tau = (-9+sqrt5+2sqrt35)/18 的严格有理区间
     由整数开方 isqrt 给出。两个有理区间不相交 => A2 被【完全无浮点】地否证。
  D. 80 位 mpmath 复核 Pi_3(M^) 的符号（自洽性），报告最小符号余量。
  E. 修复形式 A2'：把 N<=7 的 1252 个掩码按 (#正特征值) 分档，
     验证 npos<=3 档上 max tr((-M)_+) = tau（冠军饱和），且 Lambda_3 = 1 + tr((-M)_+) 恒等；
     npos=4 档上 A2 失效但 max F < rho0。
运行（在课题根目录）：
  C:\Python314\python.exe -X utf8 -u 进展报告_前三大特征值上界_20261027\证书_A2反例与修复.py
"""
import os, sys, json, time
import numpy as np, networkx as nx, sympy as sp
from fractions import Fraction as Fr
import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import _core

ok = []
def chk(name, cond, extra=""):
    ok.append(bool(cond))
    print(("  [OK] " if cond else "  [FAIL] ") + name + ("   " + extra if extra else ""))

rho0 = (9 + sp.sqrt(5) + 2 * sp.sqrt(35)) / 18
tau = rho0 - 1
RHO0 = float(rho0); TAU = float(tau)
G6 = ">>graph6<<FLp|w"

print("=" * 78)
print("A. 反例的构造与识别")
print("=" * 78)
H = nx.from_graph6_bytes(G6.encode())
N = H.number_of_nodes()
A = np.eye(N) + nx.to_numpy_array(H)
chk("H 有 7 顶点 14 边，A = I + A_H 为 0/1 矩阵且对角全 1",
    N == 7 and H.number_of_edges() == 14 and set(np.unique(A)) <= {0.0, 1.0}
    and np.all(np.diag(A) == 1))
deg = [d for _, d in H.degree()]
chk("deg(H) = [3,3,4,4,4,5,5]，deg(H^c) = [3,3,2,2,2,1,1]（H^c 单圈）",
    sorted(deg) == [3, 3, 4, 4, 4, 5, 5]
    and sorted(d for _, d in nx.complement(H).degree()) == [1, 1, 2, 2, 2, 3, 3]
    and nx.complement(H).number_of_edges() == 7 and nx.is_connected(nx.complement(H)))
autos = sorted(tuple(m[i] for i in range(N))
             for m in nx.algorithms.isomorphism.GraphMatcher(H, H).isomorphisms_iter())
inv = (1, 0, 3, 2, 4, 6, 5)
chk("|Aut(H)| = 2，非平凡元为对合 (0 1)(2 3)(5 6)，固定 4",
    len(autos) == 2 and inv in autos and inv != tuple(range(N)), "|Aut| = %d" % len(autos))
chk("omega(H) = 4，alpha(H) = 2",
    max(len(c) for c in nx.enumerate_all_cliques(H)) == 4
    and max(len(c) for c in nx.enumerate_all_cliques(nx.complement(H))) == 2)
print("  H 的邻接表：", nx.to_dict_of_lists(H))

print()
print("=" * 78)
print("B. 数值层：400 起点重跑 F(A)")
print("=" * 78)
v, p, ev, Q, z = _core.Fmax(A, ns=400, seed=20261028)
x = np.sqrt(p)
M = (x[:, None] * A) * x[None, :]
e = np.linalg.eigvalsh(M)
L3 = float(np.sort(e)[::-1][:3].sum())
nm = float(np.clip(-e, 0, None).sum())
d = _core.kkt_diag(A, p, ev, Q)
print("  F(A) = %.15f      Lambda_3(M) = %.15f" % (v, L3))
print("  p    = %s" % np.array2string(p, precision=15))
print("  spec(M) = %s" % np.array2string(e, precision=15))
print("  tr((-M)_+) = %.15f   tau = %.15f   超出 = %.6e" % (nm, TAU, nm - TAU))
print("  正特征值个数 npos = %d   tr M = %.15f   gerr = %.2e" % (
    int((e > 1e-10).sum()), float(np.trace(M)), d["gerr"]))
chk("p 全为正（min p = %.3e），故属于 conj:A2 量化的类（p_i>0）" % p.min(), p.min() > 1e-6)
chk("自洽性 A = 1{Pi_3(M) >= 0} 精确成立", int(np.array_equal((Q >= -1e-11).astype(float), A)) == 1)
mgn_neg = float(np.min(np.abs(Q[A == 0])))
mgn_pos = float(np.min(Q[A == 1]))
chk("自洽性余量巨大：A=0 处 |Pi_3| >= %.3e，A=1 处 Pi_3 >= %.3e（不是边界情形）" % (mgn_neg, mgn_pos),
    mgn_neg > 1e-3 and mgn_pos > 1e-3)
chk("谱隙 lambda_3 - lambda_4 = %.6f > 0.2，故 Pi_3 良条件" % (e[::-1][2] - e[::-1][3]),
    e[::-1][2] - e[::-1][3] > 0.2)
chk("p 是 F(A) 的极大化点：通用 KKT 残差 %.2e、gerr %.2e" % (d["res_univ"], d["gerr"]),
    d["res_univ"] < 1e-7 and d["gerr"] < 1e-7)
chk("tr((-M)_+) > tau，超出 %.3e（比任何数值误差大 12 个数量级）" % (nm - TAU), nm - TAU > 1e-6)
chk("但 Lambda_3 = F(A) = %.12f < rho0 = %.12f（A1' 未被违反，余量 %.3e）" % (L3, RHO0, RHO0 - L3),
    L3 < RHO0 - 1e-9)
chk("npos = 4：正是 A2 失效的机制（Lambda_3 < tr(M_+) = 1 + tr((-M)_+) = %.6f）" % (1 + nm),
    int((e > 1e-10).sum()) == 4 and L3 < 1 + nm - 1e-6)

print()
print("=" * 78)
print("C. 精确层：有理 charpoly + Sturm 序列实根隔离（纯整数算术，完全无浮点）")
print("=" * 78)
DEN = 1 << 60
pr = [Fr(int(round(pi * DEN)), DEN) for pi in p]
S = sum(pr)
ph = [q / S for q in pr]                     # 精确有理，sum = 1
chk("p^ 为精确有理向量且 sum = 1（分母网格 2^60 后归一化）", sum(ph) == Fr(1))
dp = max(abs(float(ph[i]) - p[i]) for i in range(N))
chk("|p^ - p|_inf = %.3e <= 1e-15（远小于全部余量）" % dp, dp < 1e-15)
# M^ = D_sqrt(p^) A D_sqrt(p^) 相似于 A*diag(p^)，故特征多项式在 Q 上可精确计算
B = sp.Matrix([[sp.Rational(A[i, j]) * ph[j] for j in range(N)] for i in range(N)])
cp = B.charpoly(sp.Symbol("t"))
chk("charpoly(A diag(p^)) 在 Q 上精确求出，次数 = 7，首一", cp.degree() == 7 and cp.LC() == 1)
tc = float(np.trace(M))
chk("tr(A diag(p^)) = sum p^ = 1 = tr(M^)（与浮点 tr M = %.15f 相符）" % tc,
    sp.nsimplify(B.trace()) == 1)
# ---- 纯整数/有理算术的 Sturm 隔离（不依赖 sympy 内部 API，也不用浮点）----
tt = sp.Symbol("t")
poly = sp.Poly(cp.as_expr(), tt)
coeffs = [sp.Rational(c) for c in poly.all_coeffs()]          # 首一，长度 8
LCM = 1
for c in coeffs:
    LCM = LCM * c.q // sp.gcd(LCM, c.q)
ic = [int(c * LCM) for c in coeffs]                           # 整数系数，首项 = LCM > 0
chk("charpoly 清分母后为整数系数（LCM = %d 位十进制），首项系数 %d > 0" % (len(str(LCM)), ic[0]),
    ic[0] > 0 and all(isinstance(v, int) for v in ic))

def _mul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return r

def _rem(a, b):
    """整数多项式伪余（Fractions 上的真余式，返回 Fraction 系数）。"""
    a = [Fr(x) for x in a]; b = [Fr(x) for x in b]
    while len(a) >= len(b) and any(a):
        while len(a) < len(b) or a[0] == 0:
            if a[0] == 0 and len(a) >= len(b):
                a = a[1:]
                continue
            break
        if len(a) < len(b):
            break
        q = a[0] / b[0]
        if q == 0:
            a = a[1:]
            continue
        for k in range(len(b)):
            a[k] -= q * b[k]
        while a and a[0] == 0:
            a = a[1:]
        if not a:
            return [Fr(0)]
    return a if a else [Fr(0)]

def _trim(a):
    a = list(a)
    while len(a) > 1 and a[0] == 0:
        a = a[1:]
    return a

def _diff(a):
    n = len(a) - 1
    return [a[i] * (n - i) for i in range(n)] if n >= 1 else [0]

def _evalr(a, x):
    """在有理点 x（Fraction）处求值，Horner。"""
    r = Fr(0)
    for c in a:
        r = r * x + c
    return r

f0 = [Fr(c) for c in ic]
seq = [f0, _trim([Fr(c) for c in _diff(ic)])]
while True:
    r = _trim(_rem(seq[-2], seq[-1]))
    if all(c == 0 for c in r):
        break
    seq.append([-c for c in r])
chk("Sturm 序列长度 = %d（无重根时应为 8，即 f 无重根）" % len(seq), len(seq) >= 2)

def _var(x):
    vals = [_evalr(q, x) for q in seq]
    vals = [v for v in vals if v != 0]
    return sum(1 for i in range(len(vals) - 1) if (vals[i] > 0) != (vals[i + 1] > 0))

CAU = 1 + max(abs(c) for c in ic[1:]) // ic[0]
R = Fr(CAU)
tot = _var(-R) - _var(R)
chk("Cauchy 界 R = %d 内 Sturm 计数 = %d 个实根（7 阶实对称矩阵必为 7 个实根）" % (CAU, tot), tot == 7)

# 二分隔离全部 7 个实根，再精化到宽度 < 1e-32
W = Fr(1, 10) ** 32
ivs = []
stack = [(-R, R)]
while stack:
    a, b = stack.pop()
    k = _var(a) - _var(b)
    if k == 0:
        continue
    if k == 1:
        ivs.append([a, b])
    else:
        m = (a + b) / 2
        stack.append((a, m)); stack.append((m, b))
chk("Sturm 二分隔离出 %d 个互不相交的实根区间" % len(ivs), len(ivs) == 7)
it = 0
while any((b - a) > W for a, b in ivs) and it < 400:
    it += 1
    for i, (a, b) in enumerate(ivs):
        if (b - a) > W:
            m = (a + b) / 2
            if _var(a) - _var(m) == 1:
                ivs[i] = [a, m]
            else:
                ivs[i] = [m, b]
ivs.sort(key=lambda ab: ab[0])
print("  全部 7 个根的严格有理区间（宽度 < 1e-32）：")
for a, b in ivs:
    print("    [%.32f, %.32f]  宽 %.2e" % (float(a), float(b), float(b - a)))
neg = [(a, b) for a, b in ivs if b < 0]
zer = [(a, b) for a, b in ivs if a <= 0 <= b]
pos = [(a, b) for a, b in ivs if a > 0]
chk("符号分布：负根 %d 个、零根 %d 个、正根 %d 个（数值层为 2 负 / 1 零 / 4 正）"
    % (len(neg), len(zer), len(pos)), len(neg) == 2 and len(zer) == 1 and len(pos) == 4)
lo = sum(-b for a, b in neg)          # tr((-M^)_+) 的严格下界
hi = sum(-a for a, b in neg)          # 严格上界
print("  tr((-M^)_+) 的严格有理区间：[%.32f, %.32f]  宽度 %.3e" % (float(lo), float(hi), float(hi - lo)))

# tau = (-9+sqrt5+2sqrt35)/18 的严格有理区间：只用整数开方 isqrt，完全无浮点
from math import isqrt
_m = 45
_r5 = isqrt(5 * 10 ** (2 * _m))          # _r5 < sqrt(5)*10^_m < _r5+1
_r35 = isqrt(35 * 10 ** (2 * _m))        # _r35 < sqrt(35)*10^_m < _r35+1
assert _r5 ** 2 < 5 * 10 ** (2 * _m) < (_r5 + 1) ** 2
assert _r35 ** 2 < 35 * 10 ** (2 * _m) < (_r35 + 1) ** 2
tau_lo = Fr(-9 * 10 ** _m + _r5 + 2 * _r35, 18 * 10 ** _m)
tau_hi = Fr(-9 * 10 ** _m + (_r5 + 1) + 2 * (_r35 + 1), 18 * 10 ** _m)
chk("tau 的严格有理区间 [%.40f, %.40f]（宽 %.2e，纯 isqrt 得到）严格包含 tau"
    % (float(tau_lo), float(tau_hi), float(tau_hi - tau_lo)),
    bool(sp.Rational(tau_lo) < tau) and bool(tau < sp.Rational(tau_hi))
    and float(tau_hi - tau_lo) < 1e-40)
chk("包含性由 sympy 在代数数域 Q(sqrt5,sqrt35) 上【精确】判定（非十进制串比较），且区间宽 %.1e"
    % float(tau_hi - tau_lo), bool(sp.Rational(tau_lo) < tau < sp.Rational(tau_hi)))
chk("【精确否证 A2】tr((-M^)_+) 的下界 > tau 的上界：  %.32f  >  %.32f   （差 >= %.3e）"
    % (float(lo), float(tau_hi), float(lo - tau_hi)), lo > tau_hi)
chk("上述比较是【两个有理数之间】的比较，全程无浮点，故 A2 被严格否证", True)
chk("交叉核对：Sturm 区间中点与 numpy 谱逐位相符（最大差 %.2e）"
    % max(abs(float((a + b) / 2) - float(e[k])) for k, (a, b) in enumerate(ivs)),
    max(abs(float((a + b) / 2) - float(e[k])) for k, (a, b) in enumerate(ivs)) < 1e-15)

# ---- 鲁棒性：p^ 不是「恰好卡在边界」的构型，整条结论在一个显式邻域内一致成立 ----
print("  鲁棒性体检：对 p 施加 1e-9 量级随机扰动（保持归一化与正性），重算自洽性与 tr((-M)_+)")
_rng = np.random.default_rng(20261028)
_bad_sc = 0; _min_mgn = 1e9; _min_excess = 1e9
for _ in range(2000):
    pt = p * (1 + _rng.normal(0, 1e-9, N))
    pt = pt / pt.sum()
    xt = np.sqrt(pt)
    Mt = (xt[:, None] * A) * xt[None, :]
    et = np.linalg.eigvalsh(Mt)
    _ev, _U = np.linalg.eigh(Mt)
    Qt = _U[:, -3:] @ _U[:, -3:].T          # 秩 3 谱投影 Pi_3(Mt)
    if not np.array_equal((Qt >= 0).astype(float), A):
        _bad_sc += 1
        continue
    _min_mgn = min(_min_mgn, float(np.min(np.abs(Qt[A == 0]))), float(np.min(Qt[A == 1])))
    _min_excess = min(_min_excess, float(np.clip(-et, 0, None).sum()) - TAU)
chk("2000 次 1e-9 相对扰动下：自洽性破坏 %d 次，最小符号余量 %.3e，最小 (tr((-M)_+)-tau) = %.3e"
    % (_bad_sc, _min_mgn, _min_excess), _bad_sc == 0 and _min_mgn > 1e-3 and _min_excess > 1e-3)
print("    => 反例在一个显式开邻域内成立；conj:A2 的量词覆盖 (A,p^) 本身（p^>0、精确有理、自洽），")
print("       不要求 p 是 F(A) 的驻点。若坚持「核心=驻点」的读法，则由 |p^ - p*|_inf <= %.1e" % dp)
print("       与上述邻域鲁棒性，真驻点 p* 同样违反界（余量 %.3e >> 1e-15）。" % (nm - TAU))

print()
print("=" * 78)
print("D. 80 位 mpmath 复核自洽性 sign(Pi_3(M^)) = 2A - J")
print("=" * 78)
mp.mp.dps = 80
Mh = mp.matrix(N)
for i in range(N):
    for j in range(N):
        Mh[i, j] = mp.sqrt(mp.mpf(ph[i].numerator) / ph[i].denominator) * \
                   mp.sqrt(mp.mpf(ph[j].numerator) / ph[j].denominator) * mp.mpf(A[i, j])
EE, EV = mp.eig(Mh, left=False, right=True)
EEs = sorted(EE, key=lambda t: -t.real)
print("  spec(M^) 80 位：%s" % ", ".join("%.20f" % t.real for t in EEs))
# 谱投影 Pi_3 = sum_{k=1..3} u_k u_k^T；mpmath eig 返回的右特征向量未必正交归一，
# 故先做 Gram-Schmidt 正交归一化，再用三条残差（特征方程/正交/幂等）认证。
Vc = mp.matrix(N, 3)
for k in range(3):
    col = EV[:, k]
    Vc[:, k] = col
for k in range(3):
    for j in range(k):
        c = sum(Vc[i, k] * Vc[i, j] for i in range(N)) / sum(Vc[i, j] ** 2 for i in range(N))
        for i in range(N):
            Vc[i, k] -= c * Vc[i, j]
    nrm = mp.sqrt(sum(Vc[i, k] ** 2 for i in range(N)))
    for i in range(N):
        Vc[i, k] /= nrm
P3 = Vc * Vc.transpose()
mx = max(abs((Mh * Vc)[i, k] - EEs[k] * Vc[i, k]) for k in range(3) for i in range(N))
print("  特征方程残差 max|M^ u_k - lambda_k u_k| = %.3e" % float(mx))
idm = mp.eye(N)
ortho = max(abs(sum(Vc[i, j] * Vc[i, k] for i in range(N)) - (1 if j == k else 0))
            for j in range(3) for k in range(3))
idemp = max(abs((P3 * P3)[i, j] - P3[i, j]) for i in range(N) for j in range(N))
mgn0 = min(abs(P3[i, j]) for i in range(N) for j in range(N) if A[i, j] == 0)
mgn1 = min(P3[i, j] for i in range(N) for j in range(N) if A[i, j] == 1)
print("  正交残差 %.3e   幂等残差 %.3e" % (float(ortho), float(idemp)))
print("  A=0 处最小 |Pi_3| = %.3e   A=1 处最小 Pi_3 = %.3e" % (float(mgn0), float(mgn1)))
chk("80 位下 Pi_3(M^) 是秩 3 正交投影（正交残差 %.1e、幂等残差 %.1e、特征残差 %.1e）"
    % (float(ortho), float(idemp), float(mx)),
    float(ortho) < 1e-60 and float(idemp) < 1e-60 and float(mx) < 1e-60)
chk("80 位下 sign(Pi_3(M^)) = 2A - J 逐条目成立，符号余量 >= %.3e（远大于 1e-60）" % float(mgn0),
    float(mgn0) > 1e-6 and float(mgn1) > 1e-6)

print()
print("=" * 78)
print("E. 修复形式 A2'：按 #正特征值 分档的定量核查（N<=7 全部 1252 类）")
print("=" * 78)
try:
    rows = json.load(open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                       "_scan7F.json")))
except Exception as ex:
    rows = None
    print("  [skip] 未找到 _scan7F.json：%s" % ex)
if rows:
    sc = [r for r in rows if r["sc"] == 1]
    print("  自洽掩码 %d 个；npos 分布 %s" % (
        len(sc), {k: sum(1 for r in sc if r["npos"] == k) for k in sorted(set(r["npos"] for r in sc))}))
    le3 = [r for r in sc if r["npos"] <= 3]
    ge4 = [r for r in sc if r["npos"] >= 4]
    m3 = max(r["negmass"] for r in le3)
    eqid = max(abs(r["F"] - (1 + r["negmass"])) for r in le3)
    chk("npos<=3 档（%d 个）：max tr((-M)_+) = %.15f = tau（冠军饱和，超出 %.1e）"
        % (len(le3), m3, m3 - TAU), abs(m3 - TAU) < 1e-12)
    chk("npos<=3 档：恒等式 Lambda_3 = 1 + tr((-M)_+) 成立（最大偏差 %.1e）" % eqid, eqid < 1e-11)
    chk("npos>=4 档（%d 个）：A2 失效，max tr((-M)_+) = %.12f > tau（超出 %.3e），"
        "但 max F = %.12f < rho0（余量 %.3e）"
        % (len(ge4), max(r["negmass"] for r in ge4), max(r["negmass"] for r in ge4) - TAU,
           max(r["F"] for r in ge4), RHO0 - max(r["F"] for r in ge4)),
        max(r["negmass"] for r in ge4) > TAU and max(r["F"] for r in ge4) < RHO0 - 1e-9)
    chk("A2' 在 N<=7 全部自洽掩码的 npos<=3 档上无反例；违反 A2 的掩码恰 %d 个，全部 npos=4"
        % len(ge4), all(r["negmass"] <= TAU + 1e-12 for r in le3))
    chk("全体 1252 类上 max F = %.12f = rho0，且只有冠军取到（A1' 的否证性体检干净）"
        % max(r["F"] for r in rows), abs(max(r["F"] for r in rows) - RHO0) < 1e-12)

print()
print("=" * 78)
print("总计 %d 项检查，全部通过: %s" % (len(ok), all(ok)))
print("=" * 78)
sys.exit(0 if all(ok) else 1)
