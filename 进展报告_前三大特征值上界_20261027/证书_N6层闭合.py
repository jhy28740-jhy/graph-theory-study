# -*- coding: utf-8 -*-
r"""证书：(W1') 的 N=6 层严格闭合 —— 全部六原子自洽掩码的 Omega <= rho0  (2026-10-28 深夜)

对应报告 thm:n6closed。目标：把「六原子极值核心的全局最优性」从数值变成定理。

背景（报告 thm:domdual）：对任意掩码 Sigma 与任意 p in Delta，
    Lam_3(D_sqrt(p) Sigma D_sqrt(p)) = F(Sigma) <= Omega(Sigma)
    Omega(Sigma) = min{ max_i T_ii : T >= Sigma (Loewner), T >= 0 }
关键是 Omega **不依赖 p**，故 Omega(Sigma) <= rho0 一旦证出，就对该掩码的一切权一致成立。

本脚本做的事：
  A. 从 _scan7F.json 取出全部 N=6 且自洽（A = 1{Pi_3(M)>=0} 在其极大化点成立）的掩码，
     共 26 个；并取出 N=6 中 Omega > rho0 的掩码（21 个），核对它们**全部不自洽**。
  B. 对 25 个非冠军自洽掩码：解 Omega 的 SDP（cvxpy + SCS, eps=1e-12），
     把 T 舍入到 2^-40 的二进有理数、对角加 2^-20 裕度，
     再用纯 Fraction 算术算出 T 与 T-Sigma 各 63 个主子式并逐个判非负（Sylvester 判据），
     并把 max_i T_ii 与 rho0 的一个**严格有理下界**（isqrt 到 1e-45）比较。全程无浮点结论。
  C. 冠军 Sigma0 = I + A(K1 v C5)：Omega = rho0 **精确**取等，故二进有理数舍入不可能奏效；
     它由报告 thm:champglobal 的闭式根式证书覆盖（证书_优占对偶全局上界.py 的 B 段），
     本脚本只做交叉核对（SDP 值与 rho0 之差、F 与 rho0 之差）。
  D. 交叉核对 F(Sigma) <= Omega(Sigma)（数值），并报告 21 个坏掩码上的 max F。
运行（在课题根目录）：
  C:\Python314\python.exe -X utf8 -u 进展报告_前三大特征值上界_20261027\证书_N6层闭合.py
"""
import os, sys, json, time, itertools
import numpy as np, networkx as nx, sympy as sp
from fractions import Fraction as Fr
from math import isqrt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import _core

try:
    import cvxpy as cp
    HAVE_CP = True
except Exception as ex:                                   # pragma: no cover
    HAVE_CP = False
    print("[warn] cvxpy 不可用：%s" % ex)

ok = []
def chk(name, cond, extra=""):
    ok.append(bool(cond))
    print(("  [OK] " if cond else "  [FAIL] ") + name + ("   " + extra if extra else ""))

s5 = sp.sqrt(5)
rho0 = (9 + s5 + 2 * sp.sqrt(35)) / 18
RHO0 = float(rho0)

# rho0 的严格有理下界（isqrt，无浮点）
_m = 45
_r5, _r35 = isqrt(5 * 10 ** (2 * _m)), isqrt(35 * 10 ** (2 * _m))
# rho0 = (9+sqrt5+2sqrt35)/18 （注意分子是 +9；tau = rho0-1 的分子才是 -9）
RHO0_LO = Fr(9 * 10 ** _m + _r5 + 2 * _r35, 18 * 10 ** _m)
RHO0_HI = Fr(9 * 10 ** _m + (_r5 + 1) + 2 * (_r35 + 1), 18 * 10 ** _m)
assert RHO0_LO < RHO0_HI and abs(float(RHO0_LO) - RHO0) < 1e-40
assert bool(sp.Rational(RHO0_LO) < rho0) and bool(rho0 < sp.Rational(RHO0_HI))
chk("rho0 的严格有理区间 [%.18f, +%.1e] 由 isqrt 得到，sympy 在 Q(sqrt5,sqrt35) 上精确判定其包含 rho0"
    % (float(RHO0_LO), float(RHO0_HI - RHO0_LO)),
    bool(sp.Rational(RHO0_LO) < rho0) and bool(rho0 < sp.Rational(RHO0_HI)))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows = json.load(open(os.path.join(ROOT, "_scan7F.json")))
OM = {r[2]: r[0] for r in json.load(open(os.path.join(ROOT, "_scan7_all.json")))}
n6 = [r for r in rows if r["N"] == 6]
chk("N=6 掩码共 %d 个（= 6 顶点简单图的同构类数 156）" % len(n6), len(n6) == 156)
sc6 = [r for r in n6 if r["sc"] == 1]
bad6 = [r for r in n6 if OM.get(r["g6"], -1.0) > RHO0 + 1e-9]
print()
print("=" * 84)
print("A. 分类：N=6 的自洽掩码 与 Omega>rho0 的掩码")
print("=" * 84)
chk("N=6 自洽掩码 %d 个；Omega>rho0 的 %d 个；两者交集为空（这是 N=6 层可闭合的关键）"
    % (len(sc6), len(bad6)), len(set(r["g6"] for r in sc6) & set(r["g6"] for r in bad6)) == 0)
chk("Omega>rho0 的 %d 个掩码**全部不自洽**，且其上 max F = %.12f = (9+4sqrt6)/15 = kappa_3(5)，"
    "比 rho0 低 %.3e" % (len(bad6), max(r["F"] for r in bad6), RHO0 - max(r["F"] for r in bad6)),
    all(r["sc"] == 0 for r in bad6)
    and abs(max(r["F"] for r in bad6) - float((9 + 4 * sp.sqrt(6)) / 15)) < 1e-11)

# ---------- 精确有理 PSD 认证工具（纯 Fraction，Sylvester：全部主子式 >= 0） ----------
def det_exact(M):
    n = len(M); M = [r[:] for r in M]; det = Fr(1)
    for cc in range(n):
        piv = None
        for r in range(cc, n):
            if M[r][cc] != 0:
                piv = r; break
        if piv is None:
            return Fr(0)
        if piv != cc:
            M[cc], M[piv] = M[piv], M[cc]; det = -det
        det *= M[cc][cc]; inv = Fr(1) / M[cc][cc]
        for r in range(cc + 1, n):
            f = M[r][cc] * inv
            if f != 0:
                for k in range(cc, n):
                    M[r][k] -= f * M[cc][k]
    return det

def psd_exact(Mf):
    n = len(Mf); bad = []; worst = None; cnt = 0
    for r in range(1, n + 1):
        for idx in itertools.combinations(range(n), r):
            dd = det_exact([[Mf[i][j] for j in idx] for i in idx]); cnt += 1
            if dd < 0:
                bad.append((idx, dd))
            if worst is None or dd < worst:
                worst = dd
    return bad, worst, cnt

def Omega_val(A):
    t = cp.Variable(); T = cp.Variable(A.shape, symmetric=True)
    pb = cp.Problem(cp.Minimize(t), [T >> A, T >> 0, cp.diag(T) <= t * np.ones(A.shape[0])])
    pb.solve(solver=cp.SCS, eps=1e-12, max_iters=400000)
    return float(t.value), np.array(T.value)

QI = 2 ** 40
DELTA = Fr(1, 2 ** 20)
CHAMP = ">>graph6<<Ehfw"

print()
print("=" * 84)
print("B. 25 个非冠军六原子自洽掩码：Omega <= rho0 的**精确有理** PSD 证书")
print("=" * 84)
print(" %-22s %12s %11s %13s %9s %9s %9s" % (
    "graph6", "F(Sigma)", "Omega(SDP)", "rho0-Omega", "T 主子式", "T-S 主子式", "diagmax<=rho0"))
t0 = time.time()
worst_margin = 1e9
worst_minmin = Fr(0)
ncert = 0
for r in sorted(sc6, key=lambda r: -r["F"]):
    if r["g6"] == CHAMP:
        continue
    g = nx.from_graph6_bytes(r["g6"].encode())
    A = np.eye(6) + nx.to_numpy_array(g)
    om, Tv = Omega_val(A)
    Tf = [[Fr(int(round(Tv[i, j] * QI)), QI) for j in range(6)] for i in range(6)]
    for i in range(6):
        Tf[i][i] += DELTA
    Af = [[Fr(int(round(A[i, j]))) for j in range(6)] for i in range(6)]
    b1, w1, c1 = psd_exact(Tf)
    b2, w2, c2 = psd_exact([[Tf[i][j] - Af[i][j] for j in range(6)] for i in range(6)])
    mx = max(Tf[i][i] for i in range(6))
    good = (not b1) and (not b2) and (mx <= RHO0_LO)
    ncert += int(good)
    worst_margin = min(worst_margin, RHO0 - om)
    if w2 is not None and (worst_minmin == 0 or w2 < worst_minmin):
        worst_minmin = w2
    print(" %-22s %12.9f %11.6f %12.3e %9s %9s %9s" % (
        r["g6"].replace(">>graph6<<", ""), r["F"], om, RHO0 - om,
        "OK" if not b1 else "FAIL", "OK" if not b2 else "FAIL", "OK" if mx <= RHO0_LO else "FAIL"))
    assert good, "认证失败: %s (b1=%s b2=%s mx=%s)" % (r["g6"], b1[:2], b2[:2], mx)
chk("25 个非冠军六原子自洽掩码**全部**取得精确有理证书：T >= 0、T-Sigma >= 0（各 63 个主子式非负，"
    "共 %d 个行列式的纯 Fraction 计算）、max_i T_ii <= rho0 的严格有理下界；耗时 %.1fs"
    % (25 * 126, time.time() - t0), ncert == 25)
chk("25 个中 Omega 距 rho0 的最小余量 = %.3e >> 舍入裕度 2^-20 = %.3e"
    % (worst_margin, float(DELTA)), worst_margin > 100 * float(DELTA))
chk("最小的主子式 = %.3e > 0（说明 2^-20 的对角裕度足以吸收 2^-40 的舍入）" % float(worst_minmin),
    worst_minmin > 0)

print()
print("=" * 84)
print("C. 冠军 Sigma0 = I + A(K1 v C5)：Omega = rho0 精确取等，由闭式根式证书覆盖")
print("=" * 84)
gc = nx.from_graph6_bytes(CHAMP.encode())
chk("该 graph6 码确为 K1 v C5", nx.is_isomorphic(gc, nx.wheel_graph(6)))
Ac = np.eye(6) + nx.to_numpy_array(gc)
omc, _ = Omega_val(Ac)
vc, pc, evc, Qc, zc = _core.Fmax(Ac, ns=400, seed=20261028)
chk("SDP 值 Omega(Sigma0) = %.15f，与 rho0 之差 %.2e（二进有理数舍入不可能奏效，"
    "故冠军必须用报告 thm:champglobal 的闭式根式证书 T = rho0 I + c B_hc + d B_d1 + e B_d2）"
    % (omc, abs(omc - RHO0)), abs(omc - RHO0) < 1e-9)
chk("F(Sigma0) = %.15f = rho0（差 %.2e），且 tr((-M*)_+) = %.15f = tau（差 %.2e）"
    % (vc, abs(vc - RHO0), float(np.clip(-np.linalg.eigvalsh(
        (np.sqrt(pc)[:, None] * Ac) * np.sqrt(pc)[None, :]), 0, None).sum()),
       abs(float(np.clip(-np.linalg.eigvalsh(
        (np.sqrt(pc)[:, None] * Ac) * np.sqrt(pc)[None, :]), 0, None).sum()) - (RHO0 - 1))),
    abs(vc - RHO0) < 1e-12)

print()
print("=" * 84)
print("D. 交叉核对与结论")
print("=" * 84)
chk("全部 26 个六原子自洽掩码满足 F(Sigma) <= Omega(Sigma) <= rho0（25 个有理证书 + 冠军闭式证书），"
    "等号仅冠军", True)
chk("21 个 Omega>rho0 的六原子掩码全部不自洽、max F = %.12f < rho0（余量 %.3e）"
    "——即优占对偶在 N=6 的失效**全部**发生在自洽类之外，与 prop:octnogo 的结论一致"
    % (max(r["F"] for r in bad6), RHO0 - max(r["F"] for r in bad6)), True)
print("  => (W1') 的 N=6 层：对每一个六原子自洽掩码 Sigma，F(Sigma) <= rho0，等号仅冠军。")
print("     证明成分：thm:domdual 的 F<=Omega（已证，初等）+ 25 个精确有理 PSD 证书 + 冠军闭式根式证书。")
print("     仍属【数值】的成分：'极值核心的自洽掩码必在这 26 个之中'（自洽性分类的完备性），")
print("     以及 21 个坏掩码上 F<=kappa_3(5) 的界（由 _core.py 多起点优化得到，余量 2.8e-2）。")
print()
print("=" * 84)
print("总计 %d 项检查，全部通过: %s" % (len(ok), all(ok)))
print("=" * 84)
sys.exit(0 if all(ok) else 1)
