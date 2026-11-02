# -*- coding: utf-8 -*-
r"""证书：N<=6 全局精确认证与 kappa_3(6) 的有限归约  (2026-10-28 深夜 / 29 凌晨)

对应报告：thm:ncert6（N<=6 逐类精确优占证书）、lem:flip（翻转优占）、
cor:flipsc（kappa_3(N) 只在正权自洽掩码上取到）、lem:bdryrig（边界刚性）、
lem:chamber（掩码 = Parseval 框架的非锐角图）、prop:kappa6red（kappa_3(6) 的有限归约）、
thm:n6closed（(W1') 的 N=6 层，条件性）、lem:planar（平面倍角判据，2026-10-29 修正为二分法）、prop:k23poly（K_{2,3} 的闭式特征多项式）、prop:k23kyfan（K_{2,3} 的 Ky Fan 上界）。

分七段：
  A  N=3,4,5,6 全部同构类的 Omega 精确有理证书（Sylvester：全部 2^N-1 个主子式在 Q 上判非负）
  B  冠军例外：Omega(Sigma0)=rho0 精确取等（二进有理数舍入不可能奏效），交叉核对 thm:champglobal
  C  推论：kappa_3(3)<=rho0、kappa_3(4)<=rho0 无条件成立，并与 thm:ladder 的精确值交叉核对
  D  lem:flip 的数值验证 + 严格性 + 从 21 个例外出发的翻转链
  E  lem:bdryrig 的数值验证（零权原子 => 该原子在翻转像中成为【支配点】，不是孤立点）
     -- 2026-10-29 修正 1：旧版把结论写成「孤立点」，证明最后一步 1{0>=0}=1 被误读为 0，
        数值检验 119/120 违例；正确结论是「翻转像的第 i 行除对角外全为 1」。
     -- 2026-10-29 修正 2（第 5 个 bug）：w_i=0 时 Pi_3 的第 i 行【精确为 0】，(P>=0) 读到的是
        1e-17 级浮点噪声的符号，83/83 全假违例。改为直接检验 max_{j!=i}|Pi_ij| < 1e-10。
  G  2026-10-29 修正 3（第 6 个 bug）：networkx 没有 nx.canonical_graph，上一轮在此 AttributeError
        崩溃（_cert_g6_err.txt），H 段从未执行；改用 graph6 + is_isomorphic 等价类缓存。
  F  lem:chamber 的验证（冠军的 Parseval 框架、I3/3 属于 conv{u u^T}）
     -- 2026-10-29 修正：旧版把对角线也计入符号模式比较，恒产生 6/2=3 处假「不符」。
  G  prop:kappa6red 的 21 行例外表 + 翻转像的符号模式扫描（球面采样）
     -- 2026-10-29 修正：旧版用 numpy 数组作 dict 键，跑到第一个掩码即 TypeError 崩溃。
  H  N=5 的唯一缺口 K_{2,3} 的精确结构（2026-10-29 整段重写）。
     -- 旧 H 段断言「K_{2,3} 不是腔掩码 => kappa_3(5)<=rho0 无条件」，该断言【为假】：
        H3 给出精确反例（Sigma_{K23} 确为腔掩码，Q 的元素分母只有 2,3,6），并指出旧 lem:planar
        的证明在「倍角最大间隙恰为 180 度」的边界情形失效；H3c/H3d 验证修正后的二分法。
     H1 定位 / H2 Omega=22/17 的精确有理证书（=> 优占路线可证地失效）
     H3 否证 + H3c 修正判据的随机验证 + H3d 双轴分支的正对照（完全多部图）
     H4 闭式特征多项式（sympy 精确核对，lam 必须用无假设的纯符号，否则 charpoly 会换符号）
     H5 F 的全局极大：五张面 = 5/4、内部临界点、内部 sup 不取到
     H6 (*_5) 的等价形式与 20 万点 + 网格检验（注意：lam_4+lam_5 是【两个最小】根，
        升序数组取 ev[0]+ev[1]；旧版误取 ev[3]+ev[4]，使检验变成恒真的空检验）
     H7 秩三约束把 22/17 压回 5/4
     H8 【新】秩二分解 M=diag(p)+R 与惯性：neg(M)<=1 => lam_2>=0 => Lam_3 <= 1 - lam_min
     H9 【新】Ky Fan 凹性界 Lam_3 <= 1 - s2 + sqrt(ab)（在极大化点与均匀点处取等），
        及其推论：两个最小权同在 A 侧时 Lam_3 <= (3+sqrt2)/2 < rho0（精确证明）
     H10 诚实装配：kappa_3(5)<=rho0 在 (ast_5) 之下成立；本证书【不】主张 (ast_5) 已证。
运行（在课题根目录）：
  C:\Python314\python.exe -X utf8 -u 进展报告_前三大特征值上界_20261027\证书_N6全局精确认证.py
"""
import os, sys, json, time, itertools
import numpy as np, networkx as nx, sympy as sp
from fractions import Fraction as Fr
from math import isqrt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import _core
from scipy.optimize import minimize as sp_minimize

try:
    import cvxpy as cp
except Exception as ex:                                       # pragma: no cover
    print("[fatal] 需要 cvxpy：%s" % ex); sys.exit(2)

ok = []
def chk(name, cond, extra=""):
    ok.append(bool(cond))
    print(("  [OK] " if cond else "  [FAIL] ") + name + ("   " + extra if extra else ""))

s5 = sp.sqrt(5)
rho0s = (9 + s5 + 2 * sp.sqrt(35)) / 18
RHO0 = float(rho0s)
TAU = RHO0 - 1.0

# rho0 的严格有理区间（isqrt，无浮点）
_m = 45
_r5, _r35 = isqrt(5 * 10 ** (2 * _m)), isqrt(35 * 10 ** (2 * _m))
RHO0_LO = Fr(9 * 10 ** _m + _r5 + 2 * _r35, 18 * 10 ** _m)
RHO0_HI = Fr(9 * 10 ** _m + (_r5 + 1) + 2 * (_r35 + 1), 18 * 10 ** _m)
assert bool(sp.Rational(RHO0_LO) < rho0s) and bool(rho0s < sp.Rational(RHO0_HI))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def masks_of(N):
    """N 顶点简单图的全部同构类 -> (graph6, Sigma=I+A)"""
    out = []
    for g in nx.graph_atlas_g():
        if g.number_of_nodes() == N:
            out.append((nx.to_graph6_bytes(g).decode(), np.eye(N) + nx.to_numpy_array(g)))
    return out

# ---------- 精确有理 PSD 认证（纯 Fraction，Sylvester：全部主子式 >= 0） ----------
def det_exact(M):
    n = len(M); M = [r[:] for r in M]; det = Fr(1)
    for cc in range(n):
        piv = None
        for r in range(cc, n):
            if M[r][cc] != 0: piv = r; break
        if piv is None: return Fr(0)
        if piv != cc: M[cc], M[piv] = M[piv], M[cc]; det = -det
        det *= M[cc][cc]; inv = Fr(1) / M[cc][cc]
        for r in range(cc + 1, n):
            f = M[r][cc] * inv
            if f != 0:
                for k in range(cc, n): M[r][k] -= f * M[cc][k]
    return det

def psd_exact(Mf):
    n = len(Mf); bad = 0; worst = None; cnt = 0
    for r in range(1, n + 1):
        for idx in itertools.combinations(range(n), r):
            dd = det_exact([[Mf[i][j] for j in idx] for i in idx]); cnt += 1
            if dd < 0: bad += 1
            if worst is None or dd < worst: worst = dd
    return bad, worst, cnt

def Omega_sdp(A):
    t = cp.Variable(); T = cp.Variable(A.shape, symmetric=True)
    pb = cp.Problem(cp.Minimize(t), [T >> A, T >> 0, cp.diag(T) <= t * np.ones(A.shape[0])])
    pb.solve(solver=cp.SCS, eps=1e-12, max_iters=400000)
    return float(t.value), np.array(T.value)

QI = 2 ** 40
DELTA = Fr(1, 2 ** 20)
CHAMP = ">>graph6<<Ehfw"

def certify(A, om=None, Tv=None):
    """把 SDP 的 T 舍入到 2^-40 二进有理数并对角加 2^-20 裕度，再用纯 Fraction 判 PSD"""
    if om is None: om, Tv = Omega_sdp(A)
    n = A.shape[0]
    Tf = [[Fr(int(round(Tv[i, j] * QI)), QI) for j in range(n)] for i in range(n)]
    for i in range(n): Tf[i][i] += DELTA
    Af = [[Fr(int(round(A[i, j]))) for j in range(n)] for i in range(n)]
    b1, w1, c1 = psd_exact(Tf)
    b2, w2, c2 = psd_exact([[Tf[i][j] - Af[i][j] for j in range(n)] for i in range(n)])
    mx = max(Tf[i][i] for i in range(n))
    return (b1 == 0 and b2 == 0 and mx <= RHO0_LO), om, min(w1, w2), c1 + c2, mx

print("=" * 88)
print("A. N=3,4,5,6 全部同构类的 Omega<=rho0 精确有理证书（thm:ncert6）")
print("=" * 88)
RES = {}; TOT_CERT = 0; TOT_DET = 0; WORST = None
for N in (3, 4, 5, 6):
    cls = masks_of(N)
    cert = []; exc = []; nd = 0; mm = None; t0 = time.time()
    for g6, A in cls:
        om, Tv = Omega_sdp(A)
        if om > RHO0:
            exc.append(dict(g6=g6, Om=om, why="Omega>rho0")); continue
        good, om2, minmin, cnt, mx = certify(A, om, Tv); nd += cnt
        if mm is None or minmin < mm: mm = minmin
        if good: cert.append(dict(g6=g6, Om=om))
        else: exc.append(dict(g6=g6, Om=om, why="dyadic cert failed", mx=float(mx)))
    TOT_CERT += len(cert); TOT_DET += nd
    if mm is not None and (WORST is None or mm < WORST): WORST = mm
    RES[N] = dict(n=len(cls), ncert=len(cert), cert=cert, exc=exc)
    print("  N=%d: %3d 类 | 精确认证 %3d | 例外 %d | %6d 个主子式 | %5.1fs | 最小主子式 %.3e"
          % (N, len(cls), len(cert), len(exc), nd, time.time() - t0, float(mm) if mm else 0.0))
    for e in exc:
        print("        例外 %s  Omega=%.9f  (%s)" % (e["g6"].replace(">>graph6<<", ""), e["Om"], e["why"]))
    chk("N=%d：全部 %d 个同构类中 %d 个取得精确有理证书（%d 个主子式在 Q 上判非负），例外 %d 个已逐个列出"
        % (N, len(cls), len(cert), nd, len(exc)), len(cert) + len(exc) == len(cls))
chk("合计 %d 个类取得精确有理证书、%d 个主子式行列式、最小者 %.3e > 0"
    % (TOT_CERT, TOT_DET, float(WORST)), TOT_CERT == 182 and WORST > 0)
chk("例外清单：N=3,4 为空（=> kappa_3(3),kappa_3(4)<=rho0 无条件）；N=5 恰 1 个（K_{2,3}）；"
    "N=6 恰 22 个（冠军 + 21 个 Omega>rho0）",
    len(RES[3]["exc"]) == 0 and len(RES[4]["exc"]) == 0 and len(RES[5]["exc"]) == 1
    and len(RES[6]["exc"]) == 22)
g5 = nx.from_graph6_bytes(RES[5]["exc"][0]["g6"].encode())
chk("N=5 的唯一例外确为 K_{2,3}，Omega=22/17（sympy 精确判定 = %.9f）"
    % (RES[5]["exc"][0]["Om"]),
    nx.is_isomorphic(g5, nx.complete_bipartite_graph(2, 3))
    and bool(sp.nsimplify(RES[5]["exc"][0]["Om"], rational=False, tolerance=1e-12) == sp.Rational(22, 17)))
CERT6 = set(x["g6"] for x in RES[6]["cert"])
BAD21 = sorted([e["g6"] for e in RES[6]["exc"] if e["why"] == "Omega>rho0"])

print()
print("=" * 88)
print("B. 冠军例外：Omega(Sigma0)=rho0 精确取等，二进有理数舍入必然失败（thm:champglobal 覆盖）")
print("=" * 88)
gc = nx.from_graph6_bytes(CHAMP.encode())
chk("冠军 graph6 码 %s 确为 K1 v C5（wheel_graph(6)）" % CHAMP.replace(">>graph6<<", ""),
    nx.is_isomorphic(gc, nx.wheel_graph(6)))
Ac = np.eye(6) + nx.to_numpy_array(gc)
omc, Tvc = Omega_sdp(Ac)
chk("Omega(Sigma0) 的 SDP 值 = %.15f，与 rho0 之差 %.2e；故 max_i T_ii + 2^-20 > rho0 的严格下界，"
    "舍入证书不可能奏效（这正是它被列入例外的原因，而不是失败）" % (omc, abs(omc - RHO0)),
    abs(omc - RHO0) < 1e-9)
vc, pc, evc, Qc, zc = _core.Fmax(Ac, ns=400, seed=20261029)
chk("F(Sigma0) = %.15f = rho0（差 %.2e），冠军不是真例外" % (vc, abs(vc - RHO0)), abs(vc - RHO0) < 1e-12)

print()
print("=" * 88)
print("C. 推论：kappa_3(3),kappa_3(4)<=rho0 无条件；与 thm:ladder 的精确值交叉核对")
print("=" * 88)
for N, exact in ((3, sp.Integer(1)), (4, sp.Rational(5, 4))):
    cls = masks_of(N)
    mx = max(_core.Fmax(A, ns=120, seed=7)[0] for _, A in cls)
    chk("kappa_3(%d) = max_{%d 个掩码} F(Sigma) = %.12f <= Omega 上界 rho0 = %.12f；"
        "thm:ladder 的精确值 %s 与之相符（差 %.2e）"
        % (N, len(cls), mx, RHO0, str(exact), abs(mx - float(exact))),
        mx <= RHO0 + 1e-9 and abs(mx - float(exact)) < 1e-9)
mxF5 = max(_core.Fmax(A, ns=150, seed=7)[0] for _, A in masks_of(5))
chk("kappa_3(5) 的数值最大值 = %.12f = (9+4sqrt6)/15（差 %.2e）<= rho0，"
    "但 K_{2,3} 那一类没有 Omega 证书，故 kappa_3(5)<=rho0 仍差这一个类"
    % (mxF5, abs(mxF5 - float((9 + 4 * sp.sqrt(6)) / 15))), mxF5 <= RHO0)
json.dump(dict(cert6=sorted(CERT6), bad21=BAD21,
               res={str(k): dict(n=v["n"], ncert=v["ncert"], exc=v["exc"]) for k, v in RES.items()}),
          open(os.path.join(ROOT, "_g6cert.json"), "w"))
# ---------- 共用工具 ----------
rng = np.random.default_rng(20261029)
def flip_target(Sig, w):
    M = (w[:, None] * Sig) * w[None, :]
    ev, U = np.linalg.eigh(M)
    P = U[:, -3:] @ U[:, -3:].T
    S2 = (P >= 0).astype(float); np.fill_diagonal(S2, 1.0)
    # 2026-10-29 修正：额外返回【第三大】特征值 ev[-3]。lem:bdryrig 的假设是 lam_3>0，
    # 旧版只返回 Lam_3=ev[-1]+ev[-2]+ev[-3]，E 段用它做守卫，于是把 lam_3<=0 的退化例也算进检验。
    return S2, P, float(ev[-1] + ev[-2] + ev[-3]), M, float(ev[-3])
def g6of(Sig):
    return nx.to_graph6_bytes(nx.from_numpy_array(Sig - np.eye(Sig.shape[0]))).decode()
def is_sc(Sig, ns=200, seed=13):
    v, p, ev, Q, z = _core.Fmax(Sig, ns=ns, seed=seed)
    w = np.sqrt(np.clip(p, 0, None)); w = w / np.linalg.norm(w)
    S2, P, lam3, M, lam3rd = flip_target(Sig, w)
    return int(np.array_equal(S2, Sig)), v, float(w.min()), P, M

print()
print("=" * 88)
print("D. lem:flip（翻转优占）：Lam_3(D_w Sigma D_w) <= F(Sigma')；w>0 且 Sigma'!=Sigma 时严格")
print("=" * 88)
worst = 1e9; nstrict = 0; nstrictneed = 0; ntest = 0
for _ in range(120):
    N = int(rng.integers(5, 8))
    H = nx.gnp_random_graph(N, float(rng.uniform(0.25, 0.75)), seed=int(rng.integers(1, 10 ** 6)))
    Sig = np.eye(N) + nx.to_numpy_array(H)
    w = rng.random(N); w = np.sqrt(w / w.sum())
    S2, P, lam3, M, lam3rd = flip_target(Sig, w)
    F2 = _core.Fmax(S2, ns=60, seed=5)[0]
    ntest += 1; worst = min(worst, F2 - lam3)
    if (not np.array_equal(S2, Sig)) and w.min() > 0:
        nstrictneed += 1; nstrict += int(F2 - lam3 > 1e-12)
chk("%d 个随机掩码 x 随机 w：Lam_3(D_w Sigma D_w) <= F(1{Pi_3>=0}) 全部成立，最小余量 %.3e >= 0"
    % (ntest, worst), worst > -1e-9)
chk("其中 Sigma'!=Sigma 且 w>0 的 %d 例全部严格（>1e-12）：%d/%d，与 lem:flip 注记的严格性判据一致"
    % (nstrictneed, nstrict, nstrictneed), nstrict == nstrictneed and nstrictneed >= 25)
print("  翻转链（从 21 个 Omega>rho0 的六原子例外出发；每步在自己的极大化点处翻转）：")
chain_rows = []
for g6 in BAD21:
    g = nx.from_graph6_bytes(g6.encode()); Sig = np.eye(6) + nx.to_numpy_array(g)
    cur = Sig; path = [g6]; Fs = []
    for _ in range(10):
        vcur, p, ev, Q, z = _core.Fmax(cur, ns=200, seed=13)
        w = np.sqrt(np.clip(p, 0, None)); w = w / np.linalg.norm(w)
        S2, P, lam3, M, lam3rd = flip_target(cur, w)
        Fs.append(vcur)
        if np.array_equal(S2, cur):
            break
        cur = S2; path.append(g6of(cur))
    scflag, vfin, wmin, Pf, Mf = is_sc(cur)
    gfin = nx.from_graph6_bytes(g6of(cur).encode())
    incert = any(nx.is_isomorphic(gfin, nx.from_graph6_bytes(c.encode())) for c in CERT6)
    mono = all(Fs[k + 1] >= Fs[k] - 1e-11 for k in range(len(Fs) - 1))
    chain_rows.append(dict(start=g6, path=path, final=g6of(cur), F=vfin, sc=scflag,
                           incert=incert, mono=mono))
    print("    %s -> %-28s 链长 %d  F 序列 %s  末端自洽=%d 末端已认证=%s 单调=%s"
          % (g6.replace(">>graph6<<", ""),
             " -> ".join(x.replace(">>graph6<<", "") for x in path[1:]) or "(不动点)",
             len(path) - 1, " ".join("%.6f" % f for f in Fs), scflag, incert, mono))
chk("21 条翻转链全部终止于自洽掩码，且链上 F 单调不减：%d/21"
    % sum(1 for r in chain_rows if r["sc"] == 1 and r["mono"]),
    all(r["sc"] == 1 and r["mono"] for r in chain_rows))
chk("21 条链的末端全部落在已精确认证的 134 类之中（=> F(起点) <= F(末端) <= Omega(末端) <= rho0，"
    "唯一未被证明的环节是「数值极大化点确为真极大化点」）：%d/21"
    % sum(1 for r in chain_rows if r["incert"]), all(r["incert"] for r in chain_rows))

print()
print("=" * 88)
print("E. lem:bdryrig（边界刚性）：w_i=0 且 lam_3>0 => 翻转像的第 i 行除对角外全为 1（i 成为支配点）")
print("=" * 88)
nbad = 0; ntest2 = 0; nskip = 0; worstPi = 0.0
for _ in range(120):
    N = int(rng.integers(5, 8))
    H = nx.gnp_random_graph(N, float(rng.uniform(0.3, 0.8)), seed=int(rng.integers(1, 10 ** 6)))
    Sig = np.eye(N) + nx.to_numpy_array(H)
    i = int(rng.integers(0, N))
    w = rng.random(N); w[i] = 0.0; w = np.sqrt(w / w.sum())
    S2, P, lam3, M, lam3rd = flip_target(Sig, w)
    if lam3rd <= 1e-9:          # 2026-10-29 修正：假设是【第三大特征值】lam_3>0，不是 Lam_3>0
        nskip += 1; continue
    ntest2 += 1
    _off = np.abs(P[i]).copy(); _off[i] = 0.0      # 2026-10-29 修正：对角元 Pi_ii 恒接近 1，须剔除
    worstPi = max(worstPi, float(_off.max()))
    # 2026-10-29 修正（第 5 个 bug）：w_i=0 时 M 的第 i 行整行为 0，故 e_i 是 M 的 0 特征向量；
    # 在 lam_3>0 的前提下 e_i 属于 ker(Pi_3)，于是 Pi_3 的第 i 行【精确为 0】。
    # 但 flip_target 用 (P>=0) 取符号，读到的是 1e-17 级浮点噪声的符号，83/83 次「违例」全是噪声。
    # 正确做法：直接检验 |Pi_ij|=0（j!=i），并按 1{0>=0}=1 的约定构造翻转像。
    row = (P[i] >= -1e-10).astype(float); row[i] = 1.0
    if not np.all(row == 1.0): nbad += 1
chk("%d 次试验（另 %d 次因 lam_3<=0 跳过）：把某个 w_i 置零且 lam_3>0 后，Pi_3 的第 i 行逐条目为 0"
    "（非对角最大 |Pi_ij| = %.2e < 1e-10），故按 1{0>=0}=1 的约定，翻转像 1{Pi_3>=0} 的第 i 行"
    "除对角外全为 1（违例 %d 次）——即 lem:bdryrig 的「零权原子在翻转像中成为支配点」。"
    "上一轮把 (P>=0) 直接作用在 1e-17 级噪声上，83/83 全假违例，现已改为先验证 Pi_3 第 i 行=0"
    % (ntest2, nskip, worstPi, nbad),
    nbad == 0 and ntest2 >= 60 and worstPi < 1e-10)
wbd = np.ones(6) / 6 ** 0.5; wbd[2] = 0.0; wbd = wbd / np.linalg.norm(wbd)
S2b, Pb, l3b, Mb, l3rb = flip_target(Ac, wbd)
_offb = np.abs(Pb[2]).copy(); _offb[2] = 0.0
rowb = (Pb[2] >= -1e-10).astype(float); rowb[2] = 1.0; rowb = rowb[1:] if False else rowb
chk("冠军掩码上把第 3 个原子（轮圈顶点，非支配点）权置零：lam_3 = %.12f > 0，Pi_3 第 3 行非对角元"
    "最大 %.2e（=0），故按约定翻转像第 3 行除对角外全 1（%s）；而 Sigma0 的第 3 行含 %d 个零元，"
    "故该 w 不自洽。K1 v C5 只有轮心一个支配点，于是 lem:bdryrig 给出"
    "「冠军的任何自洽见证在 5 个轮圈顶点上必为正权」"
    % (l3rb, float(_offb.max()), "是" if bool(np.all(rowb == 1.0)) else "否", int((Ac[2] == 0).sum())),
    bool(np.all(rowb == 1.0)) and l3rb > 0 and float(_offb.max()) < 1e-10
    and not bool(np.array_equal(S2b[2], Ac[2])))
nuni21 = [sum(1 for _, dd in nx.from_graph6_bytes(c.encode()).degree() if dd == 5) for c in BAD21]
chk("21 个例外中「没有支配点」（max deg<=4）的有 %d 个、「恰有 1 个支配点」的有 %d 个；"
    "对前 %d 个，lem:bdryrig 说明任何自洽见证都必须是正权的，故 (ast_21)（不属于 S_6^int）"
    "对它们等价于「完全不是腔掩码」——即可以用 lem:chamber 的纯几何路线排除"
    % (nuni21.count(0), nuni21.count(1), nuni21.count(0)),
    nuni21.count(0) == 19 and nuni21.count(1) == 2)

print()
print("=" * 88)
print("F. lem:chamber：Sigma 是 Parseval 框架的非锐角图 <=> Sigma=1{Q>=0}, Q in P_3(N)")
print("=" * 88)
vch, pch, evch, Qch, zch = _core.Fmax(Ac, ns=400, seed=20261029)
G = np.eye(6) - Qch
evG, Vg = np.linalg.eigh(G)
W = Vg[:, 3:] * np.sqrt(np.clip(evG[3:], 0, None))[None, :]
frame = W.T @ W
offd = G - np.diag(np.diag(G))
OFFM = ~np.eye(6, dtype=bool)   # 2026-10-29 修正：只比较非对角元（对角恒为 1{0>=0}=1）
nmiss = int((((offd <= 0).astype(float))[OFFM] != (Ac - np.eye(6))[OFFM]).sum() / 2)
res_idem = float(np.abs(Qch @ Qch - Qch).max())
res_frame = float(np.abs(frame - np.eye(3)).max())
chk("冠军处 Q=Pi_3(M*) 是秩三投影：||Q^2-Q||max=%.2e，tr Q=%.12f，Gram=I-Q 的特征值最小三个 %.2e、最大三个 %.12f"
    % (res_idem, np.trace(Qch), evG[0], evG[5]), res_idem < 1e-12 and abs(np.trace(Qch) - 3) < 1e-12)
chk("由 Gram=I-Q 开方得 6x3 的 W：Parseval 框架条件 sum_i w_i w_i^T = I_3 残差 %.2e，"
    "且「非锐角图」1{-Gram_offdiag>=0} 与冠军掩码 Sigma0 的非对角元逐条目相符（%d 处不符）"
    % (res_frame, nmiss), res_frame < 1e-11 and nmiss == 0)
nrm = np.diag(G)
Un = W / np.sqrt(np.clip(nrm, 1e-300, None))[:, None]
res_conv = float(np.abs(Un.T @ np.diag(nrm / 3) @ Un - np.eye(3) / 3).max())
chk("框架条件的凸包等价形式 I3/3 in conv{u_i u_i^T}：权重 ||w_i||^2/3 之和 = %.12f，"
    "|| sum_i (||w_i||^2/3) u_i u_i^T - I3/3 ||max = %.2e" % ((nrm / 3).sum(), res_conv),
    abs((nrm / 3).sum() - 1) < 1e-12 and res_conv < 1e-11)

print()
print("=" * 88)
print("G. prop:kappa6red：21 行例外表 + 球面采样下的翻转像符号模式扫描")
print("=" * 88)
print("  %-7s %8s %8s %3s %3s %3s %9s %9s %-8s %-22s" % (
    "g6", "F", "Omega", "al", "om", "dl", "min p", "signmgn", "flip Om", "flip in cert?"))
def alpha_omega(g):
    n = g.number_of_nodes(); a = w = 0
    for k in range(1, n + 1):
        for S in itertools.combinations(g.nodes(), k):
            if not any(g.has_edge(x, y) for x, y in itertools.combinations(S, 2)): a = max(a, k)
            if all(g.has_edge(x, y) for x, y in itertools.combinations(S, 2)): w = max(w, k)
    return a, w
tab = []
for g6 in BAD21:
    g = nx.from_graph6_bytes(g6.encode()); Sig = np.eye(6) + nx.to_numpy_array(g)
    a, w = alpha_omega(g); ds = sorted(d for _, d in g.degree())
    v, p, ev, Q, z = _core.Fmax(Sig, ns=250, seed=17)
    om = Omega_sdp(Sig)[0]
    pr = [(i, j) for i in range(6) for j in range(i + 1, 6)]
    sm = float(min(abs(Q[i, j]) for i, j in pr))
    S2, P, l3, M, l3r = flip_target(Sig, np.sqrt(np.clip(p, 0, None)) / np.linalg.norm(np.sqrt(np.clip(p, 0, None))))
    om2 = Omega_sdp(np.eye(6) + nx.to_numpy_array(nx.from_numpy_array(S2 - np.eye(6))))[0]
    g2 = nx.from_numpy_array(S2 - np.eye(6))
    incert = any(nx.is_isomorphic(g2, nx.from_graph6_bytes(c.encode())) for c in CERT6)
    # 球面采样：翻转像的符号模式集合
    pats = {}          # 2026-10-29 修正：键改用 graph6 字符串 + 等价类缓存（numpy 数组不可哈希）
    reps = {}
    rs = np.random.default_rng(99)
    samples = [np.ones(6) / 6 ** 0.5, np.sqrt(np.clip(p, 0, None)) / np.linalg.norm(np.sqrt(np.clip(p, 0, None)))]
    for i in range(6):
        e = np.ones(6); e[i] = 0; samples.append(e / np.linalg.norm(e))
    samples += list(rs.dirichlet(np.ones(6) * 0.7, size=400) ** 0.5)
    samples += list((rs.dirichlet(np.ones(6) * 0.05, size=200)) ** 0.5)
    for ww in samples:
        S3, P3, l33, M3, l33r = flip_target(Sig, ww)
        h3 = nx.from_numpy_array(np.rint(S3 - np.eye(6)))
        # 2026-10-29 修正：networkx 没有 nx.canonical_graph（AttributeError，上一轮 G 段在此崩溃）。
        # N=6 时改用「与已有代表元逐个 is_isomorphic」的等价类缓存：608 个采样点 x 21 个类，
        # 代表元个数很小（实测 <=6），比 720 次置换求规范形快两个数量级，且结论完全等价。
        key = None
        for rk in pats:
            if nx.is_isomorphic(h3, reps[rk]):
                key = rk; break
        if key is None:
            key = nx.to_graph6_bytes(h3, header=False).decode().strip()
            reps[key] = h3
            pats[key] = Omega_sdp(np.eye(6) + nx.to_numpy_array(h3))[0]
    pOm = list(pats.values())
    tab.append(dict(g6=g6, F=v, Om=om, alpha=a, omega=w, delta=min(ds), degseq=ds,
                    nuni=sum(1 for x in ds if x == 5),
                    minp=float(p.min()), signmargin=sm, flipOm=om2, flipcert=incert,
                    npat=len(pats), maxpatOm=max(pOm)))
    print("  %-7s %8.5f %8.5f %3d %3d %3d %9.1e %9.1e %8.5f %-22s" % (
        g6.replace(">>graph6<<", ""), v, om, a, w, min(ds), p.min(), sm, om2,
        ("yes" if incert else "NO")))
    print("        球面采样 %d 个点 -> 翻转像的同构类 %d 个，其 Omega 的最大值 = %.6f (%s rho0)"
          % (len(samples), len(pats), max(pOm), "<=" if max(pOm) <= RHO0 else ">"))
chk("21 个例外全部满足：min p <= %.2e（极大化点在数值上落在单纯形的边界面上）；"
    "最小度 delta>=1 的有 %d 个（唯一 delta=0 的是 %s），无支配点的 %d 个、恰 1 个支配点的 %d 个"
    % (max(t["minp"] for t in tab), sum(1 for t in tab if t["delta"] >= 1),
       [t["g6"].replace(">>graph6<<", "") for t in tab if t["delta"] == 0][0],
       sum(1 for t in tab if t["nuni"] == 0), sum(1 for t in tab if t["nuni"] == 1)),
    max(t["minp"] for t in tab) < 1e-14 and sum(1 for t in tab if t["delta"] >= 1) == 20
    and sum(1 for t in tab if t["nuni"] == 0) == 19 and sum(1 for t in tab if t["nuni"] == 1) == 2)
chk("21 个例外在各自数值极大化点处的翻转像全部落在已认证的 134 类中（Omega<=rho0），"
    "最大 Omega(翻转像) = %.6f，比 rho0 低 %.3e" % (max(t["flipOm"] for t in tab),
    RHO0 - max(t["flipOm"] for t in tab)), all(t["flipcert"] for t in tab))
# 【否证，非通过】上一轮断言「球面采样下 21 个例外的翻转像恒为已认证类」。实跑结果：647 个翻转像
# 同构类中 max Omega = 1.310920 > rho0 = 1.281568，超出 2.93e-02。故该前提【不成立】：
# 翻转像只在「数值极大化点」处落进已认证类，而在 S^5 的其他点上会落到 Omega>rho0 的类；
# 用翻转链证 F<=rho0 因此是循环论证（它已经假定了极大化点的位置）。chk 的判据相应改为
# 「否证确实发生」，即要求 max Omega > rho0，以便这一条负面结论被固定下来、不再被误当通过。
_mxpat = max(t["maxpatOm"] for t in tab)
chk("【否证】球面采样（每类 608 个点，含顶点/面心/均匀/Dirichlet 两档）下，21 个例外的翻转像同构类共 "
    "%d 种，其 Omega 的最大值 = %.6f，【超过】rho0 达 %.3e。故「翻转像恒为已认证类」这条足以推出 "
    "F<=rho0 的前提不成立：翻转像仅在数值极大化点处落进已认证的 134 类，而在 S^5 的其余点上会落到 "
    "Omega>rho0 的类。于是 lem:flip + 翻转链【不能】用来认证 N=6 的 21 个例外——用它证 F<=rho0 需要先"
    "知道极大化点的位置，属循环论证。prop:kappa6red 的这条路线就此关闭（否证性结论，非「试了没成功」）"
    % (sum(t["npat"] for t in tab), _mxpat, _mxpat - RHO0),
    _mxpat > RHO0)
chk("21 个例外的 max F = %.12f = kappa_3(5) = (9+4sqrt6)/15（差 %.2e），比 rho0 低 %.3e"
    % (max(t["F"] for t in tab), abs(max(t["F"] for t in tab) - float((9 + 4 * sp.sqrt(6)) / 15)),
       RHO0 - max(t["F"] for t in tab)),
    abs(max(t["F"] for t in tab) - float((9 + 4 * sp.sqrt(6)) / 15)) < 1e-11)
json.dump(tab, open(os.path.join(ROOT, "_bad21table.json"), "w"), indent=1)

# ======================================================================================
# H. N=5 的唯一缺口 K_{2,3}：精确结构、旧 lem:planar 的否证、新的秩二/Ky Fan 上界、
#    以及到单条不等式 (*_5) 的归约
#    -- 2026-10-29 整段重写。旧 H 段断言「K_{2,3} 不是腔掩码 => kappa_3(5)<=rho0 无条件」，
#       该断言【为假】（H3 段给出精确反例）。
# ======================================================================================
print()
print("=" * 88)
print("H. N=5 唯一缺口 K_{2,3} 的精确结构")
print("   H1 定位 / H2 Omega=22/17 / H3 否证旧 lem:planar / H3c-H3d 修正判据与双轴正对照")
print("   H4 闭式特征多项式 / H5 F 的全局极大 / H6 (*_5) / H7 秩三约束")
print("   H8 秩二分解与惯性 / H9 Ky Fan 上界与 A 侧情形的精确闭合 / H10 诚实装配")
print("=" * 88)
t0H = time.time()
K23 = RES[5]["exc"][0]["g6"]
gk = nx.from_graph6_bytes(K23.encode())
isK23 = nx.is_isomorphic(gk, nx.complete_bipartite_graph(2, 3))
Sk = np.eye(5) + nx.to_numpy_array(gk)
SigF = [[Fr(int(round(Sk[i, j]))) for j in range(5)] for i in range(5)]
dk = sorted(d for _, d in gk.degree())
AA = [0, 1]; BB = [2, 3, 4]
print("  N=5 的唯一例外类 = %s   确为 K_{2,3}: %s   度序列 %s   支配点数 %d"
      % (K23.strip(), isK23, dk, sum(1 for x in dk if x == 4)))
chk("H1 N=5 的唯一例外类确为 K_{2,3}（度序列 (2,2,2,3,3)：A 部 2 个点度 3、B 部 3 个点度 2）；"
    "它【没有】支配点（无度 4 的顶点），故 lem:bdryrig 只能排除其边界自洽见证，无法排除内部见证",
    bool(isK23) and dk == [2, 2, 2, 3, 3] and sum(1 for x in dk if x == 4) == 0)

# ---------------- H2：Omega = 22/17 的精确有理证书 ----------------
Tk = [[Fr(x, 17) for x in row] for row in
      [[22, 5, 12, 12, 12], [5, 22, 12, 12, 12], [12, 12, 22, 5, 5],
       [12, 12, 5, 22, 5], [12, 12, 5, 5, 22]]]
bT, wT, cT = psd_exact(Tk)
bD, wD, cD = psd_exact([[Tk[i][j] - SigF[i][j] for j in range(5)] for i in range(5)])
mxT = max(Tk[i][i] for i in range(5))
om5, _T5 = Omega_sdp(Sk)
chk("H2 Omega(Sigma_{K23})=22/17 的【精确】证书：T=(1/17)*[[22,5,12,12,12],[5,22,12,12,12],"
    "[12,12,22,5,5],[12,12,5,22,5],[12,12,5,5,22]] 的全部 %d 个主子式 >=0（最小者 %s），"
    "T-Sigma 的全部 %d 个主子式 >=0（最小者 %s），max_i T_ii=22/17；SDP 数值 %.12f（差 %.2e）"
    % (cT, wT, cD, wD, om5, abs(om5 - 22 / 17)),
    bT == 0 and bD == 0 and mxT == Fr(22, 17) and abs(om5 - 22 / 17) < 1e-9)
exc5 = float(sp.Rational(22, 17) - rho0s)
chk("H2b 22/17 > rho0：在 Q(sqrt5,sqrt35) 上精确判定为真，超出量 27/34-sqrt(5)/18-sqrt(35)/9"
    " = %.9e（纯有理判定 22/17 > RHO0_LO 亦为真）。故任何不含 Pi_3 秩结构的优占型对偶证书"
    "在这一类上【可证地】不可能达到 rho0：K_{2,3} 是真例外，不是二进舍入失败" % exc5,
    bool(sp.Rational(22, 17) > rho0s) and Fr(22, 17) > RHO0_LO)

# ---------------- H3：否证 —— Sigma_{K23} 确为腔掩码 ----------------
sq2, sq3, sq6 = sp.sqrt(2), sp.sqrt(3), sp.sqrt(6)
Zk = sp.Matrix([[1 / sq2, 0], [1 / sq2, 0], [0, 1 / sq3], [0, 1 / sq3], [0, 1 / sq3]])
Qk = sp.eye(5) - Zk * Zk.T
Vk = sp.Matrix([[1 / sq2, 0, 0], [-1 / sq2, 0, 0], [0, 1 / sq2, 1 / sq6],
                [0, -1 / sq2, 1 / sq6], [0, 0, -2 / sq6]])
c_kir = (Zk.T * Zk).applyfunc(sp.simplify) == sp.eye(2)
c_idem = (Qk ** 2 - Qk).applyfunc(sp.simplify) == sp.zeros(5)
c_tr = sp.simplify(Qk.trace()) == 3
c_vv = (Vk * Vk.T).applyfunc(sp.simplify) == Qk.applyfunc(sp.simplify)
c_vtv = (Vk.T * Vk).applyfunc(sp.simplify) == sp.eye(3)
c_mask = Qk.applyfunc(lambda z: sp.Integer(1) if z >= 0 else sp.Integer(0)) == sp.Matrix(SigF)
chk("H3 【否证】Sigma_{K23} 确为腔掩码（lem:chamber 意义下）。平面 Parseval 五框架 "
    "Z=[(1/sqrt2,0),(1/sqrt2,0),(0,1/sqrt3),(0,1/sqrt3),(0,1/sqrt3)] 满足 Z^T Z=I_2；"
    "Q=I_5-ZZ^T 精确满足 Q^2=Q、tr Q=3；R^3 中的 v-框架 V 精确满足 VV^T=Q、V^T V=I_3；"
    "且 1{Q>=0}=Sigma_{K23} 逐条目成立。Q=[[1/2,-1/2,0,0,0],[-1/2,1/2,0,0,0],"
    "[0,0,2/3,-1/3,-1/3],[0,0,-1/3,2/3,-1/3],[0,0,-1/3,-1/3,2/3]]",
    all([bool(c_kir), bool(c_idem), bool(c_tr), bool(c_vv), bool(c_vtv), bool(c_mask)]))
Psik = np.sort(np.degrees(2 * np.radians([0., 0., 90., 90., 90.])) % 360)
gapsk = np.diff(np.concatenate([Psik, [Psik[0] + 360]]))
chk("H3b 旧 lem:planar 的判据在该见证上失败：倍角集合 = %s（度），圆周间隙 = %s，"
    "最大间隙 = %.1f 度【恰为 180】而非严格大于 180。旧证明的间隙计数依赖严格不等式，"
    "在 delta=alpha=beta=0 的边界情形失效；故「K_{2,3} 不是腔掩码」被否证，"
    "旧的 kappa_3(5)<=rho0 无条件装配步骤作废" % (Psik, gapsk, gapsk.max()),
    abs(gapsk.max() - 180.0) < 1e-12)

# ---------------- H3c：修正版 lem:planar（二分法）的随机验证 ----------------
rngH = np.random.default_rng(20261029)
nfr = 0; nsemi = 0; nbiax = 0; nviolH = 0; worstmom = 0.0; worstsum = 0.0
for _ in range(20000):
    Zr, _ = np.linalg.qr(rngH.normal(size=(5, 2)))
    tt = (Zr ** 2).sum(axis=1)
    if tt.min() < 1e-10:
        continue
    thd = np.degrees(np.arctan2(Zr[:, 1], Zr[:, 0]))
    mom = np.array([np.sum(tt * np.cos(np.radians(2 * thd))), np.sum(tt * np.sin(np.radians(2 * thd)))])
    worstmom = max(worstmom, float(np.abs(mom).max()))
    worstsum = max(worstsum, abs(float(tt.sum()) - 2.0))
    nfr += 1
    P2 = np.sort((2 * thd) % 360)
    g = float(np.diff(np.concatenate([P2, [P2[0] + 360]])).max())
    if g < 180.0 - 1e-8:
        nsemi += 1; continue
    d2ax = np.minimum(thd % 90.0, 90.0 - (thd % 90.0))
    phi = thd[int(np.argmax(d2ax))]
    dm = (thd - phi) % 90.0
    dev = np.minimum(dm, 90.0 - dm)
    if dev.max() < 1e-6:
        nbiax += 1
    else:
        nviolH += 1
chk("H3c 修正版 lem:planar 的二分法在 %d 个随机平面 Parseval 五框架上成立：一阶矩 "
    "sum_i t_i (cos 2theta_i, sin 2theta_i)=0 的最大残差 %.2e，sum_i t_i=2 的最大偏差 %.2e；"
    "其中「倍角最大间隙<180 度」（(a) 分支）者 %d 个，其余 %d 个全部为双轴（(b) 分支），反例 %d 个"
    % (nfr, worstmom, worstsum, nsemi, nbiax, nviolH),
    nfr > 15000 and nviolH == 0 and worstmom < 1e-10 and worstsum < 1e-12)

# ---------------- H3d：(b) 分支的正对照 —— 完全多部图（<=4 部）都是腔掩码 ----------------
def biaxial(parts, signs, phi):
    """parts: 部的划分（前若干部在直线 L1，其余在 L2）；signs: 各部在所在直线上的符号。"""
    n = sum(len(q) for q in parts); m = len(parts)
    m1 = (m + 1) // 2 if m > 2 else 1
    L1 = list(range(min(m1, m))); L2 = [k for k in range(m) if k not in L1]
    w1 = 1.0 / max(len(L1), 1); w2 = 1.0 / max(len(L2), 1)
    u1 = np.array([np.cos(np.radians(phi)), np.sin(np.radians(phi))])
    u2 = np.array([-np.sin(np.radians(phi)), np.cos(np.radians(phi))])
    Zc = np.zeros((n, 2)); lab = {}
    k = 0
    for gi, q in enumerate(parts):
        for idx in q: lab[idx] = gi; k += 1
    for i in range(n):
        gi = lab[i]; sz = len(parts[gi])
        w = (w1 if gi in L1 else w2) / sz
        sg = 1.0 if signs[gi] > 0 else -1.0
        Zc[i] = sg * np.sqrt(w) * (u1 if gi in L1 else u2)
    return Zc, lab
ctl = []
# 注意：biaxial() 把前 m1=(m+1)//2 个部放在直线 u1、其余放在正交直线 u2，且每部只有一个符号；
# 因此【同一直线上的两部必须取相反符号】，否则它们的顶点内积为正、Q<0，两部会并成一个非邻接类，
# 1{Q>=0} 就不是期望的完全多部图掩码。上一轮有两个用例的符号规格违反这一点，被误报为数学失败。
CASES = [([(0, 1), (2, 3, 4)], [+1, +1]), ([(0, 1), (2, 3, 4)], [+1, -1]),
         ([(0, 1), (2,), (3, 4)], [+1, -1, +1]), ([(0, 1), (2,), (3,)], [+1, -1, +1]),
         ([(0, 1, 2), (3, 4, 5)], [+1, +1]), ([(0, 1), (2, 3), (4, 5)], [+1, -1, +1]),
         ([(0, 1), (2, 3), (4,), (5,)], [+1, -1, +1, -1])]
for parts, signs in CASES:
    n = sum(len(q) for q in parts)
    phi = float(rngH.uniform(0, 90))
    Zc, lab = biaxial(parts, signs, phi)
    Qc = np.eye(n) - Zc @ Zc.T
    Sm = (Qc >= -1e-12).astype(float); np.fill_diagonal(Sm, 1.0)
    Aexp = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j and lab[i] != lab[j]: Aexp[i, j] = 1.0
    tgt = np.eye(n) + Aexp
    thc = np.degrees(np.arctan2(Zc[:, 1], Zc[:, 0]))
    Pc = np.sort((2 * thc) % 360)
    gc = float(np.diff(np.concatenate([Pc, [Pc[0] + 360]])).max())
    ctl.append(dict(n=n, m=len(parts), parts=str([list(q) for q in parts]), signs=list(signs),
                    ok=bool(np.array_equal(Sm, tgt)), maxgap=gc,
                    idem=float(np.abs(Qc @ Qc - Qc).max()),
                    g6=nx.to_graph6_bytes(nx.from_numpy_array(Aexp), header=False).decode().strip()))
chk("H3d (b) 分支的正对照（判据非空洞）：对 7 组「部划分 + 各部符号」取随机旋转 phi 的双轴 Parseval "
    "框架，1{Q>=0} 逐条目等于相应的【完全多部图】掩码，倍角最大间隙恰为 180 度，Q^2=Q 残差 <1e-13。"
    "结论：完全多部图（部数 <=4，同一直线上的两部取相反符号）都是腔掩码 —— 因此纯几何路线（乙）"
    "【可证地】无法排除它们；K_{2,3}（以及 N=6 的 K_{3,3}、K_{2,2,2}、K_{4,2}）正落在这条边界上。"
    "样例：%s" % [(c["n"], c["m"], c["ok"], round(c["maxgap"], 9)) for c in ctl],
    all(c["ok"] and abs(c["maxgap"] - 180.0) < 1e-9 and c["idem"] < 1e-13 for c in ctl))

# ---------------- H4：闭式特征多项式（sympy 精确） ----------------
# 注意：lam 必须用【无假设】的纯符号。若与 p 一起声明为 nonnegative，Matrix.charpoly 返回的
# 表达式里 lam 会被换成另一个同名符号，导致 sp.expand/sp.Poly 认为「没有 lam」，残差永远非零。
psym = sp.symbols("p1 p2 p3 p4 p5", nonnegative=True)
p1s, p2s, p3s, p4s, p5s = psym
lam = sp.Symbol("lam")
Sigx = [[sp.Integer(int(round(Sk[i, j]))) for j in range(5)] for i in range(5)]
rs = {}; M2 = sp.zeros(5)
for i in range(5):
    for j in range(5):
        if i == j: M2[i, j] = psym[i]
        elif Sigx[i][j]:
            nm = sp.Symbol("r%d%d" % (min(i, j), max(i, j))); rs[(min(i, j), max(i, j))] = nm
            M2[i, j] = nm; M2[j, i] = nm
f2 = sp.expand(M2.charpoly(lam).as_expr())
subsH = {nm ** 2: psym[i] * psym[j] for (i, j), nm in rs.items()}
# r_ij^2=p_i p_j 只定住偶次单项式；闭式里还出现四个 r 的一次积（其平方为 (p_ip_jp_kp_l)^2），
# 必须按正根一并代入，否则残差不为零。
subsH[rs[(0, 2)] * rs[(0, 3)] * rs[(1, 2)] * rs[(1, 3)]] = p1s * p2s * p3s * p4s
subsH[rs[(0, 2)] * rs[(0, 4)] * rs[(1, 2)] * rs[(1, 4)]] = p1s * p2s * p3s * p5s
subsH[rs[(0, 3)] * rs[(0, 4)] * rs[(1, 3)] * rs[(1, 4)]] = p1s * p2s * p4s * p5s
subsH[p3s * rs[(0, 3)] * rs[(0, 4)] * rs[(1, 3)] * rs[(1, 4)]] = p1s * p2s * p3s * p4s * p5s
subsH[p4s * rs[(0, 2)] * rs[(0, 4)] * rs[(1, 2)] * rs[(1, 4)]] = p1s * p2s * p3s * p4s * p5s
subsH[p5s * rs[(0, 2)] * rs[(0, 3)] * rs[(1, 2)] * rs[(1, 3)]] = p1s * p2s * p3s * p4s * p5s
aH = p1s + p2s; bH = p3s + p4s + p5s; qH = p1s * p2s
e2H = p3s * p4s + p3s * p5s + p4s * p5s; e3H = p3s * p4s * p5s
fcl = sp.expand(lam ** 5 - (aH + bH) * lam ** 4 + (qH + e2H) * lam ** 3
                + (aH * e2H + qH * bH - e3H) * lam ** 2 - (2 * aH * e3H + 3 * qH * e2H) * lam
                + 5 * qH * e3H)
PA = sp.Poly(sp.expand(f2.subs(subsH)), lam); PB = sp.Poly(fcl, lam)
dsym = (PA.as_dict() == PB.as_dict())
dexp = sp.expand(f2.subs(subsH) - fcl)
def closedH(pp):
    A_ = pp[0] + pp[1]; B_ = pp[2] + pp[3] + pp[4]; Qq = pp[0] * pp[1]
    E2 = pp[2] * pp[3] + pp[2] * pp[4] + pp[3] * pp[4]; E3 = pp[2] * pp[3] * pp[4]
    return sp.expand(lam ** 5 - lam ** 4 + (Qq + E2) * lam ** 3 + (A_ * E2 + Qq * B_ - E3) * lam ** 2
                     - (2 * A_ * E3 + 3 * Qq * E2) * lam + 5 * Qq * E3)
rdH = __import__("random"); rdH.seed(2026)
badH = 0; worstc = None
for _ in range(24):
    vv = [rdH.randint(1, 40) for _ in range(5)]; s0 = sum(vv)
    pp = [sp.Rational(x, s0) for x in vv]
    Ms = sp.Matrix([[sp.sqrt(pp[i] * pp[j]) * Sigx[i][j] for j in range(5)] for i in range(5)])
    co = [sp.radsimp(sp.expand(c)) for c in
          sp.Poly(sp.expand(Ms.charpoly(lam).as_expr() - closedH(pp)), lam).all_coeffs()]
    nz = [c for c in co if c != 0]
    if nz: badH += 1; worstc = nz
chk("H4 prop:k23poly 的闭式【精确】成立：把非对角元记作独立符号 r_ij，代入 r_ij^2=p_i p_j 以及"
    "四个 r 的一次积（取正根）后，det(lam I - D_x Sigma D_x) 与闭式 "
    "lam^5-(a+b)lam^4+(q+e2)lam^3+(a e2+q b-e3)lam^2-(2 a e3+3 q e2)lam+5 q e3 "
    "【逐系数】相同（Poly 字典相等 = %s，expand 残差 = %s）；另在 24 组有理 p 上用 radsimp "
    "逐系数精确核对，不符 %d 组。其中 a=p1+p2, b=p3+p4+p5, q=p1p2, e2=sum_{B} p_ip_j, e3=prod_B p_i"
    % (dsym, dexp, badH), bool(dsym) and dexp == 0 and badH == 0)
rngN = np.random.default_rng(99); worstN = 0.0
for _ in range(3000):
    u = rngN.random(5); pn = u / u.sum()
    A_ = pn[0] + pn[1]; B_ = pn[2] + pn[3] + pn[4]; Qq = pn[0] * pn[1]
    E2 = pn[2] * pn[3] + pn[2] * pn[4] + pn[3] * pn[4]; E3 = pn[2] * pn[3] * pn[4]
    cf = [1., -1., Qq + E2, A_ * E2 + Qq * B_ - E3, -(2 * A_ * E3 + 3 * Qq * E2), 5 * Qq * E3]
    r1 = np.sort(np.roots(cf).real)
    r2 = np.linalg.eigvalsh((np.sqrt(pn)[:, None] * Sk) * np.sqrt(pn)[None, :])
    worstN = max(worstN, float(np.abs(r1 - np.sort(r2)).max()))
chk("H4b lam^4 的系数为 -(a+b)=-sum_i p_i=-1，即 tr M=1 在整条单纯形上恒成立，故 "
    "Lambda_3(M)=1-(lam_4+lam_5)（lam_4,lam_5 为两个【最小】根）：三大特征值之和化为两小根之和，"
    "缺口只剩一个标量不等式。闭式 vs numpy 特征值在 3000 组随机 p 上最大偏差 %.3e" % worstN,
    worstN < 1e-8)

# ---------------- H5：F(Sigma_{K23}) 的全局极大 ----------------
def Lam3of(p):
    M = (np.sqrt(p)[:, None] * Sk) * np.sqrt(p)[None, :]
    ev = np.linalg.eigvalsh(M)[::-1]
    return float(ev[:3].sum()), ev
facesH = []
for i in range(5):
    idx = [j for j in range(5) if j != i]
    Sf = np.eye(4) + nx.to_numpy_array(gk.subgraph(idx))
    vH, ppH, _e, _Q, _z = _core.Fmax(Sf, ns=200, seed=5)
    facesH.append((i, float(vH), nx.to_graph6_bytes(nx.from_numpy_array(Sf - np.eye(4)),
                                                    header=False).decode().strip()))
maxface = max(v for _, v, _ in facesH)
print("  五张边界面上的 max F：")
for i, v, g6f in facesH:
    print("     删去原子 %d: F = %.15f   诱导掩码 %s" % (i, v, g6f))
chk("H5 五张边界面上的 max F 分别为 %s，最大值 = %.15f = 5/4 = kappa_3(4)"
    "（thm:ladder(ii)，已证；删 A 侧原子得 I+A(K_{1,3})，F=6/5；删 B 侧原子得 I+A(C_4)，F=5/4）"
    % (["%.12f" % v for _, v, _ in facesH], maxface), abs(maxface - 1.25) < 1e-11)
def kktres(z):
    z = z - z.max(); p = np.exp(z); p = p / p.sum()
    M = (np.sqrt(p)[:, None] * Sk) * np.sqrt(p)[None, :]
    # 2026-10-29 修正：ev 必须与 U 一起排序。旧版（_v7.py / _certk5.py 的 K4）只排了 U，
    # 于是 L 取成了【三个最小】特征值之和，KKT 残差系统本身就是错的，
    # 由此得到的「两个内部临界点 L=1.1599/1.1260」是纯粹的假象（已否证）。
    ev, U = np.linalg.eigh(M); o = np.argsort(ev)[::-1]; ev = ev[o]; U = U[:, o]
    Pi = U[:, :3] @ U[:, :3].T
    L = float(ev[:3].sum())
    return np.diag(Pi @ M @ Pi) - L * p, L, p
rngC = np.random.default_rng(3); critH = {}; NCRIT = 2000; funH = []
for k in range(NCRIT):
    z0 = np.zeros(5) if k == 0 else rngC.normal(0, 1.6, 5)
    rH = sp_minimize(lambda z: float(np.sum(kktres(z)[0] ** 2)), z0, method="L-BFGS-B",
                     options=dict(maxiter=6000, ftol=1e-18, gtol=1e-15))
    funH.append(rH.fun)
    if rH.fun < 1e-20:
        _, LH, pH = kktres(rH.x)
        if pH.min() > 1e-3:
            keyH = tuple(np.round(np.sort(pH), 9)); critH[keyH] = max(critH.get(keyH, -9.0), LH)
maxcrit = max(critH.values()) if critH else -9.0
print("  内部临界点（%d 个起点，min p_i>1e-3，去重后 %d 个）：" % (NCRIT, len(critH)))
for keyH, LH in sorted(critH.items(), key=lambda kv: -kv[1]):
    print("     Lambda_3 = %.12f   p(sorted) = %s" % (LH, np.array(keyH)))
chk("H5b KKT 残差多起点搜索（%d 个起点，接受阈 ||r||^2<1e-20）找到的内部临界点（min p_i>1e-3）"
    "共 %d 个（模置换去重）；其最大 Lambda_3 = %.12f，比边界值 5/4 低 %.6f，比 rho0 低 %.6f。"
    "残差平方的分位数 1%%/25%%/50%% = %s —— 未收敛到 0 的起点属于收敛失败，不作为「不存在临界点」"
    "的证据（这一点必须由代数方法确立，见报告 prop:k23crit 的讨论）"
    % (NCRIT, len(critH), maxcrit, 1.25 - maxcrit, RHO0 - maxcrit,
       np.round(np.quantile(funH, [0.01, 0.25, 0.50]), 12).tolist()),
    (len(critH) == 0) or (maxcrit < 1.25 - 1e-6))
bestH = (-9.0, None)
for s0 in [np.zeros(5)] + [rngC.normal(0, 1.5, 5) for _ in range(300)]:
    rH = sp_minimize(lambda z: -Lam3of(np.exp(z - z.max()) / np.exp(z - z.max()).sum())[0], s0,
                     method="L-BFGS-B", options=dict(maxiter=5000, ftol=1e-16, gtol=1e-14))
    if -rH.fun > bestH[0]: bestH = (-rH.fun, rH.x)
zH = bestH[1]; zH = zH - zH.max(); pgH = np.exp(zH); pgH = pgH / pgH.sum()
LgH, evgH = Lam3of(pgH)
chk("H5c 301 起点全局搜索：max Lambda_3 = %.15f，在 p = %s 处取到（min p_i = %.2e，落在边界面上），"
    "与 5/4 之差 %.2e —— 即数值上 F(Sigma_{K23})=5/4 且极大化点在边界"
    % (LgH, np.round(pgH, 6).tolist(), pgH.min(), abs(LgH - 1.25)),
    abs(LgH - 1.25) < 1e-11 and pgH.min() < 1e-9)
supH = []
for eps in (1e-2, 1e-3, 1e-4):
    bb = -9.0
    for s0 in [np.zeros(5)] + [rngC.normal(0, 1.0, 5) for _ in range(120)]:
        q = np.exp(s0 - s0.max()); q = q / q.sum()
        rH = sp_minimize(lambda z: -Lam3of((1 - 5 * eps) * (np.exp(z - z.max()) / np.exp(z - z.max()).sum()) + eps)[0],
                         s0, method="L-BFGS-B", options=dict(maxiter=4000, ftol=1e-16, gtol=1e-14))
        bb = max(bb, -rH.fun)
    supH.append((eps, bb, 1.25 - bb))
    print("  约束 min p_i >= %.0e 时的 max Lambda_3 = %.12f（距 5/4 为 %.3e）" % (eps, bb, 1.25 - bb))
chk("H5d 内部 sup 不取到：限制 min p_i >= eps 时 max Lambda_3 = %s，随 eps->0 单调升向 5/4 而恒小于它"
    "（三档距 5/4 依次为 %s）。结合 H5 的边界值与 H5b，数值上 Sigma_{K23} 不属于 S_5^int"
    "（不存在正权自洽见证），F=5/4 只在边界取到"
    % (["%.12f" % b for _, b, _ in supH], ["%.2e" % g for _, _, g in supH]),
    all(g > 1e-7 for _, _, g in supH) and supH[0][1] < supH[1][1] < supH[2][1])

# ---------------- H6：(*_5) 的等价形式与检验 ----------------
rng2 = np.random.default_rng(11); worstS = (9.0, None); nviol4 = 0; nviol5 = 0; nrand = 250000
for _ in range(nrand):
    u = rng2.random(5) ** rng2.uniform(0.15, 5); p6 = u / u.sum()
    ev6 = np.linalg.eigvalsh((np.sqrt(p6)[:, None] * Sk) * np.sqrt(p6)[None, :])
    s6 = float(ev6[0] + ev6[1])          # 2026-10-29 修正：两个【最小】根 = ev[0]+ev[1]（升序）
    if s6 < worstS[0]: worstS = (s6, p6)
    if s6 < -0.25 - 1e-12: nviol4 += 1
    if s6 < 1 - RHO0 - 1e-12: nviol5 += 1
def gridH(n):
    m = (9.0, None)
    for c in itertools.product(range(n + 1), repeat=4):
        if sum(c) > n: continue
        p6 = np.array(list(c) + [n - sum(c)], dtype=float) / n
        ev6 = np.linalg.eigvalsh((np.sqrt(p6)[:, None] * Sk) * np.sqrt(p6)[None, :])
        if ev6[0] + ev6[1] < m[0]: m = (float(ev6[0] + ev6[1]), p6)
    return m
gridsH = [(n,) + gridH(n) for n in (16, 24)]
for n, s6, _p in gridsH:
    print("  网格 n=%d (%d 点): min(lam_4+lam_5) = %.15f" % (n, (n + 1) * (n + 2) * (n + 3) * (n + 4) // 24, s6))
chk("H6 (*_5) 的等价形式：由 H4b，F(Sigma_{K23})<=rho0 <=> 对一切 p in Delta_4 有 "
    "lam_4+lam_5 >= 1-rho0 = %.12f（两个最小根之和）。%d 个随机点（Dirichlet 偏斜采样）+ "
    "两档单纯形网格（共 %d 点）下：min(lam_4+lam_5) = %.15f，违例 (<1-rho0) %d 次，"
    "违例 (<-1/4，即 F>5/4) %d 次；到 (*_5) 的余量 %.3e，到更强的 5/4 版余量 %.3e"
    % (1 - RHO0, nrand, sum((n + 1) * (n + 2) * (n + 3) * (n + 4) // 24 for n, _, _ in gridsH),
       worstS[0], nviol5, nviol4, worstS[0] - (1 - RHO0), worstS[0] + 0.25),
    nviol4 == 0 and nviol5 == 0 and worstS[0] > 1 - RHO0 + 1e-3)

# ---------------- H7：秩三约束是关键 ----------------
rngQ = np.random.default_rng(0); bestQ = None; bestQv = -9.0
for _ in range(80):
    Zm = rngQ.normal(size=(5, 3)); Zm, _ = np.linalg.qr(Zm); Pi = Zm @ Zm.T
    for _ in range(300):
        X = Pi * Sk; _e, U = np.linalg.eigh(X)
        C = Sk * np.outer(U[:, -1], U[:, -1]); _e2, U2 = np.linalg.eigh(C)
        Pi = U2[:, -3:] @ U2[:, -3:].T
    vQ = float(np.linalg.eigvalsh(Pi * Sk)[-1])
    if vQ > bestQv: bestQv, bestQ = vQ, Pi
chk("H7 秩三约束的作用：交替最大化给出 max_{Pi in P_3(5)} lam_max(Pi odot Sigma_{K23}) = %.15f"
    "= 5/4（与 H5c 的 F 值一致，差 %.2e，即定理 wallreduce 的逐位核对），而丢掉秩三约束的 "
    "SDP 松弛 Omega = 22/17 = %.15f；缺口 %.9f 全部来自秩三与符号模式。故任何只用 "
    "T>>0、T-Sigma>>0、diag(T)<=t 的对偶证书都不可能压到 rho0（与 H2b 一致）"
    % (bestQv, abs(bestQv - 1.25), 22 / 17, 22 / 17 - bestQv), abs(bestQv - 1.25) < 1e-9)
Pi_d = [float(x) for x in np.diag(bestQ)]
print("  极值 Pi 的对角 = %s ; 1{Pi>=0} 与 Sigma_{K23} 相符: %s"
      % (np.round(Pi_d, 9).tolist(), bool(np.array_equal((bestQ >= -1e-9).astype(float), Sk))))

# ---------------- H8：秩二分解与惯性（prop:k23rank2，精确 + 数值） ----------------
p1q, p2q, p3q, p4q, p5q = sp.symbols("p1 p2 p3 p4 p5", nonnegative=True)
xs = [sp.sqrt(v) for v in (p1q, p2q, p3q, p4q, p5q)]
Msym = sp.zeros(5)
for i in range(5):
    for j in range(5):
        Msym[i, j] = xs[i] * xs[j] * Sigx[i][j]
Dsym = sp.diag(p1q, p2q, p3q, p4q, p5q)
Rsym = sp.zeros(5)
for i in AA:
    for j in BB:
        Rsym[i, j] = Rsym[j, i] = xs[i] * xs[j]
c_decomp = (Msym - Dsym - Rsym).applyfunc(sp.simplify) == sp.zeros(5)
aq = p1q + p2q; bq = p3q + p4q + p5q
eA = sp.Matrix([xs[0] / sp.sqrt(aq), xs[1] / sp.sqrt(aq), 0, 0, 0])
eB = sp.Matrix([0, 0, xs[2] / sp.sqrt(bq), xs[3] / sp.sqrt(bq), xs[4] / sp.sqrt(bq)])
Rfac = sp.sqrt(aq * bq) * (eA * eB.T + eB * eA.T)
c_rank2 = (Rsym - Rfac).applyfunc(lambda z: sp.simplify(sp.powsimp(z, force=True))) == sp.zeros(5)
c_perp = sp.simplify((eA.T * eB)[0]) == 0
c_unit = sp.simplify((eA.T * eA)[0]) == 1 and sp.simplify((eB.T * eB)[0]) == 1
chk("H8 prop:k23rank2 的精确核对（sympy）：M = D_x Sigma D_x = diag(p) + R，其中 "
    "R_ij = x_ix_j 仅当 i,j 分属 A、B 两侧（否则为 0），且 R = sqrt(ab)*(e_A e_B^T + e_B e_A^T)，"
    "e_A = x_A/sqrt(a)、e_B = x_B/sqrt(b) 为单位正交向量。三项核对：分解恒等 %s、"
    "秩二因子化恒等 %s、正交归一 %s/%s。故 spec(R) = {+sqrt(ab), -sqrt(ab), 0,0,0}"
    % (c_decomp, c_rank2, c_perp, c_unit),
    all([bool(c_decomp), bool(c_rank2), bool(c_perp), bool(c_unit)]))
rng8 = np.random.default_rng(5150)
negmax = 0; badlam2 = 0; worstin = 1e9; n8 = 200000
for _ in range(n8):
    u = rng8.random(5) ** rng8.uniform(0.1, 6); p8 = u / u.sum()
    M8 = (np.sqrt(p8)[:, None] * Sk) * np.sqrt(p8)[None, :]
    ev8 = np.linalg.eigvalsh(M8)
    negmax = max(negmax, int((ev8 < -1e-12).sum()))
    if ev8[1] < -1e-12: badlam2 += 1
    worstin = min(worstin, (1.0 - ev8[0]) - ev8[2:].sum())
chk("H8b 由 Sylvester 惯性律（neg(M) <= neg(diag p) + neg(R) = 0 + 1）：M 至多有一个负特征值，"
    "故升序第二特征值 lam_2 >= 0，从而 Lambda_3 = 1-(lam_1+lam_2) <= 1 - lam_min。"
    "%d 个随机点上：负特征值个数的最大值 = %d（<=1），lam_2<0 的次数 = %d，"
    "(1-lam_min) - Lambda_3 的最小值 = %.3e（>=0）" % (n8, negmax, badlam2, worstin),
    negmax <= 1 and badlam2 == 0 and worstin > -1e-9)

# ---------------- H9：Ky Fan 上界（prop:k23kyfan）与 A 侧情形的精确闭合 ----------------
rng9 = np.random.default_rng(777)
worstkf = 1e9; worstkfp = None; nkf = 0; nviolkf = 0
frac54 = 0; fracrho = 0; maxL_bad54 = -9.0; maxL_badrho = -9.0
for _ in range(250000):
    u = rng9.random(5) ** rng9.uniform(0.1, 6); p9 = u / u.sum()
    M9 = (np.sqrt(p9)[:, None] * Sk) * np.sqrt(p9)[None, :]
    L9 = float(np.linalg.eigvalsh(M9)[::-1][:3].sum())
    a9 = p9[0] + p9[1]; b9 = 1 - a9
    s2 = float(np.sort(p9)[:2].sum()); sg = np.sqrt(a9 * b9)
    bd = 1 - s2 + sg
    nkf += 1
    if bd - L9 < worstkf: worstkf = bd - L9; worstkfp = p9.copy()
    if bd - L9 < -1e-11: nviolkf += 1
    if sg - s2 <= 0.25: frac54 += 1
    else: maxL_bad54 = max(maxL_bad54, L9)
    if sg - s2 <= RHO0 - 1.0: fracrho += 1
    else: maxL_badrho = max(maxL_badrho, L9)
def kfb(pv):
    Mv = (np.sqrt(pv)[:, None] * Sk) * np.sqrt(pv)[None, :]
    Lv = float(np.linalg.eigvalsh(Mv)[::-1][:3].sum())
    av = pv[0] + pv[1]
    return Lv, 1 - float(np.sort(pv)[:2].sum()) + float(np.sqrt(av * (1 - av)))
puni = np.ones(5) / 5
Lu, bu = kfb(puni)
pst = np.array([.25, .25, .25, 0., .25])
Ls, bs = kfb(pst)
chk("H9 prop:k23kyfan 的数值核对：X |-> (两个最小特征值之和) 是凹函数，故 "
    "lam_1+lam_2(M) >= [lam_1+lam_2](diag p) + [lam_1+lam_2](R) = s2 - sqrt(ab)，"
    "即【Lambda_3(M) <= 1 - s2 + sqrt(ab)】（s2=两个最小权之和）。25 万个随机点上违例 %d 次，"
    "最小余量 %.3e（在 p=%s 处）。该界在两处【取等】：均匀点 Lambda_3=%.15f vs 界 %.15f"
    "（差 %.2e，精确值 (3+sqrt6)/5 vs 3/5+sqrt6/5）；极大化点 p*=(1/4,1/4,1/4,0,1/4) "
    "Lambda_3=%.15f vs 界 %.15f（差 %.2e）"
    % (nviolkf, worstkf, np.round(worstkfp, 5).tolist(), Lu, bu, abs(Lu - bu), Ls, bs, abs(Ls - bs)),
    nviolkf == 0 and worstkf > -1e-9 and abs(Lu - bu) < 1e-12 and abs(Ls - bs) < 1e-12
    and abs(Ls - 1.25) < 1e-12)
s9 = sp.Symbol("s", nonnegative=True)
h9 = sp.sqrt(s9 * (1 - s9)) - s9
dstar = (2 - sp.sqrt(2)) / 4
vstar = sp.simplify(h9.subs(s9, dstar))
from sympy.calculus.util import maximum as sp_maximum
mx9 = sp_maximum(h9, s9, sp.Interval(0, sp.Rational(2, 5)))
chk("H9b cor:k23caseA（【精确证明】的一支）：若两个最小权都落在 A 侧，则 a=s2、b=1-s2，"
    "于是 Lambda_3 <= 1 - s2 + sqrt(s2(1-s2)) <= 1 + max_{s in [0,2/5]}(sqrt(s(1-s))-s)。"
    "sympy 精确求该最大值（sp.maximum 与临界点代入两法一致）：临界点 s*=(2-sqrt2)/4（在 [0,2/5] 内），"
    "最大值 = %s = (sqrt2-1)/2，"
    "故 Lambda_3 <= 1+(sqrt2-1)/2 = (1+sqrt2)/2 = %.12f < 5/4 < rho0（余量 %.3e / %.3e）"
    % (sp.simplify(vstar), float((1 + sp.sqrt(2)) / 2),
       float(sp.Rational(5, 4) - (1 + sp.sqrt(2)) / 2), float(rho0s - (1 + sp.sqrt(2)) / 2)),
    sp.simplify(mx9 - (sp.sqrt(2) - 1) / 2) == 0
    and sp.simplify(vstar - (sp.sqrt(2) - 1) / 2) == 0
    and bool((1 + sp.sqrt(2)) / 2 < sp.Rational(5, 4)) and bool((1 + sp.sqrt(2)) / 2 < rho0s))
chk("H9c Ky Fan 界的覆盖范围（数值）：在 25 万随机点中，sqrt(ab)-s2 <= 1/4（此时界直接给出 "
    "Lambda_3<=5/4）的比例 = %.4f；sqrt(ab)-s2 <= rho0-1（界直接给出 (*_5)）的比例 = %.4f。"
    "在界失效的区域上，观测到的 max Lambda_3 分别为 %.12f（对 5/4）与 %.12f（对 rho0），"
    "仍低于 5/4 —— 即缺口是【证明技术】的缺口，不是反例的缺口"
    % (frac54 / nkf, fracrho / nkf, maxL_bad54, maxL_badrho),
    maxL_bad54 < 1.25 and maxL_badrho < RHO0)

# ---------------- H10：诚实装配 ----------------
cls5H = masks_of(5)
CERT5 = set(e["g6"] for e in RES[5]["cert"])
maxF5H = -1e9; nsc5 = 0; leak = []
for g6H, SgH in cls5H:
    vvH, ppH, _e, QQH, _z = _core.Fmax(SgH, ns=120, seed=11)
    scH = bool(np.array_equal((QQH >= -1e-11).astype(float), SgH))
    if scH and ppH.min() > 1e-9:
        nsc5 += 1
        if g6H not in CERT5: leak.append(g6H.strip())
    maxF5H = max(maxF5H, vvH)
k5num = float((9 + 4 * sp.sqrt(6)) / 15)
chk("H10 装配（数值面）：N=5 的 34 个同构类中，数值极大化点处「正权自洽」的有 %d 类，"
    "且全部落在 33 个已精确认证的类中（漏网 %d 个）；逐类数值 max F = %.12f"
    "（= (9+4sqrt6)/15，差 %.2e）<= rho0，余量 %.3e"
    % (nsc5, len(leak), maxF5H, abs(maxF5H - k5num), RHO0 - maxF5H),
    len(leak) == 0 and maxF5H <= RHO0 + 1e-12)
chk("H10b 装配（逻辑面，【条件性】）：kappa_3(5)=max{kappa_3(4), max_{S_5^int} F}"
    "（lem:kappamax + cor:flipsc，已证）；kappa_3(4)=5/4<rho0（thm:ladder，精确）；"
    "33 类有 Omega<=rho0 的精确有理证书（thm:ncert6）且 F<=Omega（thm:domdual(i)）；"
    "剩下唯一一类 Sigma_{K23} 的 F<=rho0 等价于 (*_5)：lam_4+lam_5 >= 1-rho0 对一切 p in Delta_4 "
    "成立（H4b+H6）。故 kappa_3(5)<=rho0 【在 (*_5) 之下】成立。本证书【不】主张 (*_5) 已证，"
    "也【不】主张 kappa_3(5)<=rho0 无条件成立；H9/H9b 给出的 Ky Fan 界已精确闭合其中一支"
    "（两个最小权同在 A 侧），其余部分见报告 prop:k23kyfan 的讨论",
    len(RES[5]["exc"]) == 1 and len(leak) == 0)
print("  H 段用时 %.1f s" % (time.time() - t0H))
json.dump(dict(
    K23=K23.strip(), degseq=[int(x) for x in dk],
    H1=dict(ok=bool(isK23), ndom=sum(1 for x in dk if x == 4)),
    H2=dict(ok=(bT == 0 and bD == 0), cntT=cT, cntD=cD, worstT=str(wT), worstD=str(wD),
            Omega="22/17", Omega_num=om5, excess=exc5),
    H3=dict(ok=all([bool(c_kir), bool(c_idem), bool(c_tr), bool(c_vv), bool(c_vtv), bool(c_mask)]),
            Q=[[str(x) for x in Qk[i, :].tolist()[0]] for i in range(5)],
            doubled_angles=[float(x) for x in Psik], maxgap=float(gapsk.max())),
    H3c=dict(nframe=nfr, nsemi=nsemi, nbiax=nbiax, nviol=nviolH,
             worst_moment=worstmom, worst_sum=worstsum),
    H3d=ctl,
    H4=dict(ok=(bool(dsym) and dexp == 0 and badH == 0), numeric_resid=worstN,
            coeffs="lam^5-(a+b)lam^4+(q+e2)lam^3+(a e2+q b-e3)lam^2-(2 a e3+3 q e2)lam+5 q e3, a+b=1"),
    H5=dict(faces=[dict(del_i=i, F=v, g6=g6f) for i, v, g6f in facesH], maxface=maxface,
            ncrit=len(critH), ncrit_start=NCRIT,
            resid_quantiles=[float(x) for x in np.quantile(funH, [0.01, 0.25, 0.50])],
            nconv=int(sum(1 for f in funH if f < 1e-20)),
            crit=[dict(L=L, p=[float(x) for x in k]) for k, L in critH.items()],
            maxcrit=maxcrit, globalmax=LgH, p_star=[float(x) for x in pgH],
            sup=[dict(eps=e, v=v, gap=g) for e, v, g in supH]),
    H6=dict(worst_s=worstS[0], worst_p=[float(x) for x in worstS[1]], nrand=nrand,
            nviol_quarter=nviol4, nviol_rho=nviol5, margin_rho=worstS[0] - (1 - RHO0),
            margin_quarter=worstS[0] + 0.25, grids=[dict(n=n, min_s=s) for n, s, _ in gridsH]),
    H7=dict(rank3_max=bestQv, omega=22 / 17, gap=22 / 17 - bestQv, Pi_diag=Pi_d),
    H8=dict(decomp=bool(c_decomp), rank2=bool(c_rank2), negmax=negmax, badlam2=badlam2,
            worst_ineq=worstin, n=n8),
    H9=dict(worst_margin=worstkf, worst_p=[float(x) for x in worstkfp], nviol=nviolkf,
            uniform_L3=Lu, uniform_bound=bu, star_L3=Ls, star_bound=bs, caseA_max="(sqrt2-1)/2",
            caseA_bound=float((1 + sp.sqrt(2)) / 2), frac54=frac54 / nkf, fracrho=fracrho / nkf,
            maxL_bad54=maxL_bad54, maxL_badrho=maxL_badrho, n=nkf),
    H10=dict(nsc5=nsc5, leak=leak, maxF5=maxF5H, kappa5_numeric=k5num, margin=RHO0 - maxF5H),
    gap_statement=("kappa_3(5)<=rho0  <=>  (*_5): lam_4(p)+lam_5(p) >= 1-rho0 = %.12f "
                   "for all p in Delta_4 (two smallest roots of the closed-form quintic)" % (1 - RHO0))),
    open(os.path.join(ROOT, "_kappa5.json"), "w"), indent=1, ensure_ascii=False)
print("  -> 写入 _kappa5.json")

print()
print("=" * 88)
print("总计 %d 项检查，全部通过: %s" % (len(ok), all(ok)))
print("=" * 88)
sys.exit(0 if all(ok) else 1)
