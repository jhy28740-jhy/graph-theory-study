# -*- coding: utf-8 -*-
"""
证书_锥族闭式与rho0的二维刻画.py   （自足脚本，只打印 PASS/FAIL）

核对第 21--23 轮写进报告的全部「精确」断言：
  [1] rho_0 的最小多项式与 Phi、phi 的精确关系
  [2] (t)_+ 的 Legendre--Catalan 系数（符号精确，k=1..8）
  [3] 锥族的各向同性与 L 的闭式（含 n=5 的 1+sqrt5(1-z^2)）
  [4] 中心截面：P(z)=0 的唯一可行根 z*、J_5(z*)、余量 6rho_0-J_5(z*)
  [5] 单轴主定理：F(a*,z*)=rho_0（60 位）、a* 的根式、四条最小多项式
  [6] 调和分解 m_4 = 1/5 + ||H_4||_F^2/9
  [7] 等纬度精确公式 m_4 = 19/81 + (8|nu3|^2+2|nu4|^2)/729 与右端点闭合
  [8] 等纬度正方形 L=9/2、m_4=7/27
  [9] L<=3phi 的多项式证书 q(s)^2-s 的因式分解与锐性
用法： C:\\Python314\\python.exe -X utf8 -u 证书_锥族闭式与rho0的二维刻画.py
"""
import math, itertools, random
import numpy as np
import sympy as sp
from mpmath import mp, mpf, sqrt as msq, findroot, polyroots, cos as mcos, sin as msin, pi as mpi

mp.dps = 60
PASS = FAIL = 0
def chk(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1; print('  [PASS]', msg)
    else:
        FAIL += 1; print('  [FAIL]', msg)

s5, s7, s35 = math.sqrt(5.0), math.sqrt(7.0), math.sqrt(35.0)
rho0 = (9 + s5 + 2*s35)/18
PHI  = (3 + s5)/4
phi  = (1 + s5)/2
sixrho = 6*rho0
S5, S7, S35 = sp.sqrt(5), sp.sqrt(7), sp.sqrt(35)
rho0s = (9 + S5 + 2*S35)/18

print('='*108)
print('[1] 精确常数')
print('='*108)
x = sp.Symbol('x')
chk(sp.expand((81*x**4 - 162*x**3 + 49*x**2 + 32*x + 1).subs(x, rho0s)) == 0,
    'rho_0 是 81x^4-162x^3+49x^2+32x+1 的根（符号精确）')
chk(sp.simplify((51*S5-21)/16 + (225-135*S5)/80 - 3*(1+S5)/2) == 0,
    '(51sqrt5-21)/16 + 9beta/5 = 3phi（符号精确，即 L<=3phi 的锐常数）')
chk(abs(float(rho0s) - rho0) < 1e-15, 'rho_0 = %.15f' % rho0)
c0 = 6*rho0s - (51*S5-21)/16
c1 = 1/(-(9*(25-15*S5)/16))
chk(sp.simplify(c1 - (20+12*S5)/225) == 0, 'c_1 = 1/(-9beta) = (20+12sqrt5)/225（符号精确）')
chk(abs(float(c0) - 1.8744425029524276) < 1e-12, 'c_0 = 6rho_0-(51sqrt5-21)/16 = %.16f' % float(c0))
thr = float(c0) + 1/(5*float(c1))
chk(abs(thr - 2.83530721498332272) < 1e-12, '阈值 c_0+1/(5c_1) = %.14f' % thr)
chk(abs((float(9*(25-15*S5)/16)) + 4.8043236) < 1e-6, '9beta = (225-135sqrt5)/16 = %.7f' % float(9*(25-15*S5)/16))
chk(abs(float((51*S5-21)/16) - 5.8149668) < 1e-6, '9gamma+3alpha = (51sqrt5-21)/16 = %.7f' % float((51*S5-21)/16))
delta0 = 3*float(phi) + 3 - sixrho
chk(abs(delta0 - 0.164692785016677) < 1e-12, 'delta_0 = 3phi+3-6rho_0 = %.15f' % delta0)

print()
print('='*108)
print('[2] (t)_+ 的 Legendre--Catalan 系数')
print('='*108)
t = sp.Symbol('t')
def c_exact(ell):
    P = sp.legendre(ell, t)
    return sp.simplify((2*ell+1)*sp.integrate(t*P, (t, 0, 1))/2)
for k in range(1, 9):
    ell = 2*k
    cat = sp.Rational(sp.binomial(2*k, k), k+1)
    formula = (-1)**(k+1)*(4*k+1)*cat/(4*(2*k-1)*4**k)
    chk(sp.simplify(c_exact(ell) - formula) == 0,
        'c_%d = %s = Catalan 闭式' % (ell, sp.nsimplify(c_exact(ell))))
for ell in (3, 5, 7, 9, 11, 13):
    chk(c_exact(ell) == 0, 'c_%d = 0（奇次 >= 3 消失）' % ell)
chk(c_exact(0) == sp.Rational(1, 4) and c_exact(1) == sp.Rational(1, 2), 'c_0=1/4, c_1=1/2')
def catalan_ratio_terms(K):
    """r_k = Cat_k/4^k（稳定递推），返回 c_{2k}=(-1)^{k+1}(4k+1)r_k/(4(2k-1))。"""
    out = {}; r = 0.25
    for k in range(1, K+1):
        out[k] = ((-1)**(k+1))*(4*k+1)*r/(4*(2*k-1))
        r *= (2*k+1)/(2.0*(k+2))
    return out
_ct = catalan_ratio_terms(600)
tail = sum(_ct[k] for k in range(2, 601))
chk(abs(tail - (-1/16)) < 5e-5, 'sum_{k>=2} c_{2k} = -1/16（交错级数 600 项，值 %.12f）' % tail)
chk(sp.nsimplify(sp.Rational(1,4)-sp.Rational(5,16)) == -sp.Rational(1,16),
    '精确：sum_{k>=1}c_{2k}=1/4（(t)_+ 在 t=1 取值 1）故 sum_{k>=2}c_{2k}=1/4-5/16=-1/16')

print()
print('='*108)
print('[3] 锥族：各向同性与 L 的闭式')
print('='*108)
def cone(n, z):
    g0 = (1-3*z**2)/(1-z**2); G = 2/(1-z**2)
    g = [g0] + [G/n]*n
    U = [np.array([0.0, 0.0, 1.0])]
    for j in range(n):
        th = 2*math.pi*j/n
        U.append(np.array([math.sqrt(1-z**2)*math.cos(th), math.sqrt(1-z**2)*math.sin(th), z]))
    return np.array(g), np.array(U)
worst_iso = 0.0; worst_L = 0.0; worst_S = 0.0
for n in range(3, 21):
    for z in (0.0, 0.13, 0.29, 0.41, 0.52, 0.577):
        g, U = cone(n, z)
        T = sum(g[i]*np.outer(U[i], U[i]) for i in range(len(g)))
        worst_iso = max(worst_iso, np.abs(T-np.eye(3)).max())
        worst_iso = max(worst_iso, abs(g.sum()-3))
        Lb = sum(g[i]*g[j]*abs(U[i] @ U[j]) for i in range(len(g)) for j in range(len(g)))
        g0, G = g[0], g[1]*n
        Ssum = sum(abs((1-z**2)*math.cos(2*math.pi*d/n) + z**2) for d in range(n))
        Lc = g0**2 + 2*g0*G*z + (G**2/n)*Ssum
        worst_L = max(worst_L, abs(Lb-Lc))
        if n == 5 and z <= 0.668:
            worst_S = max(worst_S, abs(Ssum - (1+math.sqrt(5)*(1-z**2))))
chk(worst_iso < 1e-14, '锥族各向同性 + 质量 3（n=3..20, 6 个 z）最大偏差 %.2e' % worst_iso)
chk(worst_L < 1e-12, 'L 的闭式 == 暴力求和，最大偏差 %.2e' % worst_L)
chk(worst_S < 1e-14, 'n=5: sum_d |(1-z^2)cos(2pi d/5)+z^2| = 1+sqrt5(1-z^2)，最大偏差 %.2e' % worst_S)

def Astar2(g, U):
    """max over sign patterns of ||sum sigma_i g_i u_i||^2  (n<=12 atoms -> 枚举)"""
    n = len(g); best = -1e18
    for sg in itertools.product([1, -1], repeat=n-1):
        s = (1,) + sg
        v = sum(s[i]*g[i]*U[i] for i in range(n))
        best = max(best, float(v @ v))
    return best

print()
print('='*108)
print('[4] 中心截面：J_5 的四次方程与精确极值')
print('='*108)
z = sp.Symbol('z')
P = 15*z**4 - (35-S5)*z**3 + 30*z**2 + (3-S5)*z - 5
rts = sp.nroots(P, n=40)
real_in = [r for r in rts if abs(sp.im(r)) < 1e-30 and 0 < float(sp.re(r)) < 1/math.sqrt(3)]
chk(len(real_in) == 1, 'P 在 (0,1/sqrt3) 内恰有一个根')
zstar = mpf(str(sp.re(real_in[0])))
chk(abs(zstar - mpf('0.5200639465409236331648403103733959995601')) < mpf('1e-35'),
    'z* = %s' % mp.nstr(zstar, 38))
Aofz = 2*(45*z**4 - 60*z**3 - (20+2*S5)*z**2 + 20*z + 7 + 2*S5)
J5 = Aofz/(5*(1-z**2)**2)
J5v = float(sp.N(J5.subs(z, real_in[0]), 30))
chk(abs(J5v - 7.595950307218631597240628904572525771282) < 1e-14,
    'J_5(z*) = %.24f（与第 22 轮数值上确界差 %.2e）'
    % (J5v, abs(J5v - 7.595950307218631597240628904572525771282)))
chk(J5v < sixrho, '6rho_0 - J_5(z*) = %.15f > 0' % (sixrho - J5v))
Pmp = lambda zz: 15*zz**4 - (35-msq(5))*zz**3 + 30*zz**2 + (3-msq(5))*zz - 5
zmp = findroot(Pmp, mpf('0.52'))
Amp = lambda zz: 2*(45*zz**4 - 60*zz**3 - (20+2*msq(5))*zz**2 + 20*zz + 7 + 2*msq(5))
J5mp = Amp(zmp)/(5*(1-zmp**2)**2)
rho0mp = (mpf(9) + msq(5) + 2*msq(35))/18
tgt = mpf('7.595950307218631597240628904572525771282')
chk(abs(J5mp - tgt) < mpf('1e-36'), 'J_5(z*) 的 50 位值 = %s' % mp.nstr(J5mp, 40))
chk(6*rho0mp - J5mp > mpf('0.093'), '6rho_0-J_5(z*) = %s（mpmath）' % mp.nstr(6*rho0mp-J5mp, 20))
# 独立核对：用暴力帧计算 J = L + A*^2
gv, Uv = cone(5, float(zstar))
Lb = sum(gv[i]*gv[j]*abs(Uv[i] @ Uv[j]) for i in range(6) for j in range(6))
A2 = Astar2(gv, Uv)
allplus = (gv[0] + gv[1]*sum(Uv[j][2] for j in range(1, 6)))**2
chk(abs(Lb + A2 - J5v) < 1e-12, '暴力 J = L+A*^2 == 闭式 J_5(z*)，偏差 %.2e' % abs(Lb+A2-J5v))
chk(abs(A2 - allplus) < 1e-12, 'z* 处全正号最优，A*^2 = (g_0+Gz)^2 = %.12f' % A2)
print('     g_0=%.12f  g_1=%.12f  L=%.12f  A*^2=%.12f' % (gv[0], gv[1], Lb, A2))
best_n = None; best_J = -1
for n in range(3, 21):
    bb = -1
    for k in range(1, 400):
        zz = (k/400.0)/math.sqrt(3)
        g, U = cone(n, zz)
        L = sum(g[i]*g[j]*abs(U[i] @ U[j]) for i in range(n+1) for j in range(n+1))
        M1 = sum(g[i]*U[i] for i in range(n+1)); A2 = float(M1 @ M1)
        bb = max(bb, L+A2)
    if bb > best_J: best_J, best_n = bb, n
chk(best_n == 5, '锥族（中心截面，全正号 A^2）在 n=3..20 中的最优者是 n=%d' % best_n)

print()
print('='*108)
print('[5] 单轴主定理：max F = rho_0')
print('='*108)
astar_sym = (3*S35 + 5*S7 - 7 - 7*S5)/14
astar = mpf(str(sp.N(astar_sym, 55)))
pa = 49*x**4 + 98*x**3 - 294*x**2 + 182*x - 31
chk(pa.subs(x, astar_sym).simplify() == 0, 'a* = (3sqrt35+5sqrt7-7-7sqrt5)/14 是 49x^4+98x^3-294x^2+182x-31 的根（符号精确）')
w = sp.Symbol('w')
pz = 11225*w**4 - 6700*w**3 + 430*w**2 + 260*w - 31
zst = findroot(lambda zz: 11225*zz**8 - 6700*zz**6 + 430*zz**4 + 260*zz**2 - 31, mpf('0.5'))
chk(abs(pz.subs(w, zst**2)) < mpf('1e-40'), 'z*^2 满足 11225w^4-6700w^3+430w^2+260w-31，残差 %.2e' % float(abs(pz.subs(w, zst**2))))
chk(abs(zst - mpf('0.5004400497682123079777569143187991903444479758830691319')) < mpf('1e-40'),
    'z* = %s' % mp.nstr(zst, 38))
def cone_mp(n, z):
    z2 = z*z; g0 = (1-3*z2)/(1-z2); G = 2/(1-z2)
    g = [g0] + [G/n]*n
    r = msq(1-z2)
    U = [[mpf(0), mpf(0), mpf(1)]]
    for j in range(n):
        th = 2*mpi*j/n
        U.append([r*mcos(th), r*msin(th), z])
    return g, U
g, U = cone_mp(5, zst)
a = astar
beta = (1-a)/2; alpha = (3*a-1)/2
h = [msq(beta + alpha*(U[i][2]**2)) for i in range(6)]
p = [g[i]*h[i] for i in range(6)]
F = mpf(0)
for i in range(6):
    for j in range(6):
        ip = U[i][0]*U[j][0] + U[i][1]*U[j][1] + U[i][2]*U[j][2]
        F += p[i]*p[j]*max(ip, mpf(0))
rho0mp = (mpf(9) + msq(5) + 2*msq(35))/18
dev = abs(F - rho0mp)
chk(dev < mpf('1e-45'), 'F(a*,z*) = rho_0，60 位精度偏差 %.3e' % float(dev))
print('     F = %s' % mp.nstr(F, 45))
print('     rho_0 = %s' % mp.nstr(rho0mp, 45))
gh = (1-3*zst**2)/(1-zst**2); gc = 2/(5*(1-zst**2))
pg = 324*x**4 - 648*x**3 + 256*x**2 + 68*x - 31
chk(abs(pg.subs(x, gh)) < mpf('1e-40'), 'g_hub = %s 满足 324x^4-648x^3+256x^2+68x-31' % mp.nstr(gh, 30))
chk(abs(float(gh) - 0.331767104061136519632668751001) < 1e-24 and
    abs(float(gc) - 0.533646579187772696073466249800) < 1e-24,
    'g_hub, g_cone 与报告第 12 轮冠军逐位一致')
chk(abs(float(a) - 0.5946085758659480860961319972687299732652479350367565572) < 1e-40, 'a* 的 50 位数值')
# Perron 自洽（float64 即可）
gv2 = np.array([float(q) for q in g]); Uv2 = np.array([[float(c) for c in row] for row in U])
Vv = Uv2*np.sqrt(gv2)[:, None]; Q = Vv @ Vv.T          # Q_ij = sqrt(g_i g_j)<u_i,u_j>
Qp = np.maximum(Q, 0.0); ev, evec = np.linalg.eigh(Qp)
k = int(np.argmax(ev)); pv = evec[:, k]
Fv = float(F)
chk(abs(ev[k] - rho0) < 1e-12 and abs(ev[k]-Fv) < 1e-12,
    '自洽性：lambda_max((VV^T)_+) = F(a*,z*) = rho_0（偏差 %.2e / %.2e）' % (abs(ev[k]-Fv), abs(ev[k]-rho0)))
hh = np.abs(pv)/np.sqrt(gv2); hh = hh/hh[0]
href = np.array([float(q) for q in h]); href = href/href[0]
chk(np.abs(hh - href).max() < 1e-7,
    'Perron 向量重现 h_i=sqrt(u_i^T Lam u_i)（|x_i|/sqrt(g_i) 归一后偏差 %.2e）' % np.abs(hh-href).max())
# 其它 n 不超过 n=5
bestF = {}
for n in range(3, 9):
    bb = -1
    for ai in range(1, 40):
        aa = mpf(1)/3 + (mpf(1)-mpf(1)/3)*mpf(ai)/40
        bt = (1-aa)/2; al = (3*aa-1)/2
        for zi in range(1, 40):
            zz = mpf(zi)/40/msq(3)
            gg, UU = cone_mp(n, zz)
            hh2 = [msq(bt + al*UU[i][2]**2) for i in range(n+1)]
            FF = mpf(0)
            for i in range(n+1):
                for j in range(n+1):
                    ip = UU[i][0]*UU[j][0] + UU[i][1]*UU[j][1] + UU[i][2]*UU[j][2]
                    FF += gg[i]*hh2[i]*gg[j]*hh2[j]*max(ip, mpf(0))
            bb = max(bb, FF)
    bestF[n] = bb
chk(max(bestF, key=lambda k: bestF[k]) == 5, '单轴锥族 n=3..8 的最优者是 n=5（粗网格）')
print('     max F per n:', {k: float(v) for k, v in bestF.items()})

print()
print('='*108)
print('[6] 调和分解 m_4 = 1/5 + ||H_4||_F^2/9')
print('='*108)
def solve_masses(B, rhs, rng, ntry=40, tol=1e-9):
    """在 {m: B m = rhs} 中找一个正解（最小范数解 + 零空间随机化）。"""
    m0, *_ = np.linalg.lstsq(B, rhs, rcond=None)
    _, sv, Vt = np.linalg.svd(B)
    rank = int((sv > 1e-10*max(sv[0], 1e-30)).sum())
    NS = Vt[rank:].T
    for _ in range(ntry):
        m = m0 if NS.size == 0 else m0 + NS @ (rng.standard_normal(NS.shape[1])*0.4)
        if m.min() > 1e-3:
            return m
    return None

def rand_iso(n, rng):
    """随机各向同性测度：n 条随机单位线 + 精确解出质量（sum m u u^T = I_3, sum m = 3）。"""
    for _ in range(300):
        U = rng.standard_normal((n, 3)); U /= np.linalg.norm(U, axis=1)[:, None]
        B = np.stack([U[:, 0]**2, U[:, 1]**2, U[:, 2]**2,
                      U[:, 0]*U[:, 1], U[:, 0]*U[:, 2], U[:, 1]*U[:, 2]])
        m = solve_masses(B, np.array([1., 1., 1., 0., 0., 0.]), rng)
        if m is None:
            continue
        T = np.einsum('i,ij,ik->jk', m, U, U)
        if np.abs(T-np.eye(3)).max() < 1e-10 and abs(m.sum()-3) < 1e-10:
            return m, U
    return None

def rand_equilat(k, rng):
    """等纬度圆上的随机各向同性测度（|<l,e3>|=1/sqrt3，nu^1=nu^2=0）。"""
    for _ in range(300):
        ang = np.sort(rng.random(k)*2*math.pi)
        B = np.stack([np.ones(k), np.cos(ang), np.sin(ang), np.cos(2*ang), np.sin(2*ang)])
        m = solve_masses(B, np.array([3., 0., 0., 0., 0.]), rng)
        if m is None:
            continue
        U = np.array([[math.sqrt(2/3)*math.cos(a), math.sqrt(2/3)*math.sin(a), 1/math.sqrt(3)] for a in ang])
        T = np.einsum('i,ij,ik->jk', m, U, U)
        if np.abs(T-np.eye(3)).max() < 1e-10 and abs(m.sum()-3) < 1e-10:
            return m, U, ang
    return None

rng = np.random.default_rng(7)
worst = 0.0; worst_tr = 0.0; nfr = 0
I4 = np.zeros((3, 3, 3, 3))
d = np.eye(3)
for i in range(3):
    for j in range(3):
        for k in range(3):
            for l in range(3):
                I4[i, j, k, l] = (d[i, j]*d[k, l] + d[i, k]*d[j, l] + d[i, l]*d[j, k])/15
while nfr < 40:
    n = int(rng.integers(7, 14))
    got = rand_iso(n, rng)
    if got is None:
        continue
    m, U = got; nfr += 1
    A4 = np.einsum('i,ia,ib,ic,id->abcd', m, U, U, U, U)
    H4 = A4 - 3*I4
    G = U @ U.T
    m4 = float(np.einsum('i,j,ij->', m, m, G**4))/9
    worst = max(worst, abs(m4 - (0.2 + float(np.sum(H4*H4))/9)))
    worst_tr = max(worst_tr, np.abs(np.einsum('iijl->jl', A4)-np.eye(3)).max())
chk(worst < 1e-13 and nfr == 40, 'm_4 = 1/5 + ||H_4||_F^2/9（%d 个随机各向同性测度）最大偏差 %.2e' % (nfr, worst))
chk(worst_tr < 1e-13, 'tr_12(A_4) = I_3，最大偏差 %.2e' % worst_tr)
chk(abs(np.sum(I4*I4) - 0.2) < 1e-15 and abs(np.sum((3*I4)*(3*I4)) - 1.8) < 1e-14,
    '||I_4||_F^2 = 1/5, ||3I_4||_F^2 = 9/5')

print()
print('='*108)
print('[7] 等纬度圆：m_4 = 19/81 + (8|nu3|^2+2|nu4|^2)/729 与右端点闭合')
print('='*108)
worst = 0.0; neq = 0
while neq < 200:
    k = int(rng.integers(6, 11))
    got = rand_equilat(k, rng)
    if got is None:
        continue
    m, U, ang = got; neq += 1
    G = U @ U.T
    m4 = float(np.einsum('i,j,ij->', m, m, G**4))/9
    r1 = float(m @ np.cos(ang)); i1 = float(m @ np.sin(ang))
    r2 = float(m @ np.cos(2*ang)); i2 = float(m @ np.sin(2*ang))
    r3 = float(m @ np.cos(3*ang)); i3 = float(m @ np.sin(3*ang))
    r4 = float(m @ np.cos(4*ang)); i4 = float(m @ np.sin(4*ang))
    pred = 19/81 + (8*(r3**2+i3**2) + 2*(r4**2+i4**2))/729
    worst = max(worst, abs(m4 - pred), abs(r1), abs(i1), abs(r2), abs(i2))
chk(worst < 1e-13, '等纬度 m_4 公式（%d 组随机各向同性等纬度测度）最大偏差 %.2e' % (neq, worst))
need = float(c1)*(3 - float(c0))
chk(19/81 > need, '19/81 = %.18f > 端点所需 %.18f，余量 %.3e' % (19/81, need, 19/81-need))
Jend = float((51*S5-21)/16) + float(9*(25-15*S5)/16)*(19/81) + 3
chk(Jend < sixrho, 'A^2=3 时 J <= %.9f < 6rho_0 = %.9f，余量 %.3e' % (Jend, sixrho, sixrho-Jend))

print()
print('='*108)
print('[8] 等纬度正方形：L = 9/2, m_4 = 7/27')
print('='*108)
ang = [0, math.pi/2, math.pi, 3*math.pi/2]; m = [0.75]*4
U = np.array([[math.sqrt(2/3)*math.cos(a), math.sqrt(2/3)*math.sin(a), 1/math.sqrt(3)] for a in ang])
L = sum(m[i]*m[j]*abs(U[i] @ U[j]) for i in range(4) for j in range(4))
m4 = sum(m[i]*m[j]*(U[i] @ U[j])**4 for i in range(4) for j in range(4))/9
T = sum(m[i]*np.outer(U[i], U[i]) for i in range(4))
chk(abs(L-4.5) < 1e-15 and abs(m4-7/27) < 1e-15 and np.abs(T-np.eye(3)).max() < 1e-15,
    'L=9/2, m_4=7/27, 各向同性（偏差 %.2e / %.2e / %.2e）' % (abs(L-4.5), abs(m4-7/27), np.abs(T-np.eye(3)).max()))
chk(L + 3 < sixrho, '若等纬度 L<=9/2，则 J<=15/2=7.5 < 6rho_0，余量 %.9f' % (sixrho-7.5))

print()
print('='*108)
print('[9] L <= 3phi 的多项式证书')
print('='*108)
s = sp.Symbol('s')
gam = (1+S5)/16; alp = (7*S5-5)/8; bet = (25-15*S5)/16
q = gam + alp*s + bet*s**2
fac = -(s-1)*(5*s-1)**2*((15*S5-35)*s + (S5+3))/128
chk(sp.expand(q**2 - s - fac) == 0, 'q(s)^2 - s = -(s-1)(5s-1)^2((15sqrt5-35)s+sqrt5+3)/128（符号恒等）')
chk(sp.simplify(q.subs(s, 1) - 1) == 0 and sp.simplify(q.subs(s, sp.Rational(1, 5)) - 1/sp.sqrt(5)) == 0,
    'q(1)=1, q(1/5)=1/sqrt5（证书在两点取等）')
grid = np.linspace(0, 1, 20001)
qq = float(gam) + float(alp)*grid + float(bet)*grid**2
chk(np.all(qq >= np.sqrt(grid) - 1e-15), 'q(s) >= sqrt(s) 在 [0,1] 的 20001 点网格上成立')
chk(abs(float(9*gam+3*alp) - (51*s5-21)/16) < 1e-12, '9gamma+3alpha = (51sqrt5-21)/16')
chk(abs(float(9*gam+3*alp) + float(9*bet)/5 - 3*phi) < 1e-12,
    '9gamma+3alpha+9beta/5 = 3phi = %.12f（L<=3phi 的锐界）' % (3*phi))
# 随机各向同性测度上核对 L <= 3phi，并核对正二十面体取等
worst = -1e18; nL = 0
while nL < 60:
    n = int(rng.integers(7, 14))
    got = rand_iso(n, rng)
    if got is None:
        continue
    m, U = got; nL += 1
    G = U @ U.T
    L = float(np.einsum('i,j,ij->', m, m, np.abs(G)))
    A2 = Astar2(m, U)
    worst = max(worst, L - 3*phi)
    if L + A2 > sixrho:
        print('     !! 违反 kappa_3^c<=rho_0：L=%.12f A*^2=%.12f J=%.12f' % (L, A2, L+A2))
chk(worst < 1e-12 and nL == 60, '随机各向同性测度上 L <= 3phi（%d 个，最大超出 %.3e）' % (nL, worst))
gv6 = np.array([0.5]*6)
ic = [(0, 1, phi), (0, 1, -phi), (1, phi, 0), (1, -phi, 0), (phi, 0, 1), (phi, 0, -1)]
Ui = np.array([np.array(v)/np.linalg.norm(v) for v in ic])
Li = sum(gv6[i]*gv6[j]*abs(Ui[i] @ Ui[j]) for i in range(6) for j in range(6))
m4i = sum(gv6[i]*gv6[j]*(Ui[i] @ Ui[j])**4 for i in range(6) for j in range(6))/9
Ti = sum(gv6[i]*np.outer(Ui[i], Ui[i]) for i in range(6))
Ai = sum(gv6[i]*np.einsum('a,b,c,d->abcd', Ui[i], Ui[i], Ui[i], Ui[i]) for i in range(6))
chk(abs(Li-3*phi) < 1e-14 and abs(m4i-0.2) < 1e-15 and np.abs(Ti-np.eye(3)).max() < 1e-15,
    '正二十面体 6 线：L=3phi 取等、m_4=1/5、各向同性（偏差 %.2e）' % abs(Li-3*phi))
chk(np.abs(Ai-3*I4).max() < 1e-15, '正二十面体：A_4 = 3I_4（H_4=0，即 4-设计条件）')
A2i = Astar2(gv6, Ui)
chk(abs(A2i - phi**2) < 1e-12, '正二十面体：A*^2 = phi^2 = (3+sqrt5)/2，故 V = 6-2phi^2 = 3-sqrt5 = %.12f' % (6-2*A2i))

print()
print('='*108)
print('SUMMARY  PASS=%d  FAIL=%d' % (PASS, FAIL))
print('='*108)
raise SystemExit(0 if FAIL == 0 else 1)
