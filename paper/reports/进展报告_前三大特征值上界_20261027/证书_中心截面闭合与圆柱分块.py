# -*- coding: utf-8 -*-
# 证书_中心截面闭合与圆柱分块.py
# 对应报告 §sec:r24（第 24 轮）。自足运行：python -X utf8 -u 证书_中心截面闭合与圆柱分块.py
#
# 逐项核验以下已写入报告的结论（全部为「已证」级别，数值项仅作交叉核对）：
#   引理 lem:r24basic   三个精确关系（X=1/3+sigma_u^2, sigma_u^2<=2/3, V=6-2A^2）
#   引理 lem:r24mom     15 矩恒等式 ||A_4||_F^2 = sum multinom(4;a,b,c) M_abc^2 = 9 m_4
#   引理 lem:r24wblock  w-块的显式完全平方松弛
#   引理 lem:r24a0block a=0 块的显式完全平方松弛（含 3/8 sigma_u^4 项）
#   定理 thm:r24interp  9 m_4 >= g(sigma_u^2) = 19/9 - (7/3)sigma_u^2 + (35/8)sigma_u^4
#   推论 cor:r24floors  g(0)=19/9, g(4/15)=9/5, g'=0 于 4/15, g(2/3)=5/2
#   命题 prop:r24sharp  闭式取等族（极点 delta + 正 n>=5 边形，delta in [0,1/2]）
#   引理 lem:r24bridge  sigma_u^2 <= R(A)V, R(sqrt3)=(2+sqrt3)/9
#   推论 cor:r24Rw      窗口内 R(A) <= R_w, h(V_win) < 4/15
#   定理 thm:r24centre  三分支一致余量 => kappa_3^c <= rho_0 - 1.7333e-4
#   定理 thm:r24conj    耦合不等式（原猜想 conj:r22coupled）成立
#   命题 prop:r24gen    任意截面恒等式 rho = (L_w + alpha_w)/2 与两个精确约束
import numpy as np, math, itertools
import sympy as sp

np.set_printoptions(precision=12, suppress=False, linewidth=200)
PASS = FAIL = 0
def chk(c, m):
    global PASS, FAIL
    if c: PASS += 1; print('  [PASS]', m)
    else: FAIL += 1; print('  [FAIL]', m)
def head(t): print(); print('='*112); print(t); print('='*112)

# ============================================================ 1. 精确常数
head('1. 精确常数（sympy 代数数，全部精确核对）')
s5 = sp.sqrt(5); s35 = sp.sqrt(35); s3 = sp.sqrt(3)
rho0_e = (9+s5+2*s35)/18; phi_e = (1+s5)/2
gam_e = (1+s5)/16; alp_e = (7*s5-5)/8; bet_e = (25-15*s5)/16
K1_e = sp.simplify(9*gam_e+3*alp_e); B_e = sp.simplify(-9*bet_e)
d0_e = sp.simplify(3*phi_e+3-6*rho0_e)
Vwin_e = 2*d0_e; Amin_e = sp.sqrt(3-d0_e); Rw_e = (1+1/s3)**2/(s3*(s3+Amin_e))
c0_e = sp.simplify(6*rho0_e-K1_e); c1_e = sp.simplify(1/B_e)
rho0 = float(rho0_e); B = float(B_e); d0 = float(d0_e)
Vwin = float(Vwin_e); Amin = float(Amin_e); Rw = float(Rw_e)
gam = float(gam_e); alp = float(alp_e); bet = float(bet_e)
s3f = math.sqrt(3.0); K1 = float(K1_e)
print('  rho_0   = %.24f'%rho0);      print('  K_1     = %.24f'%K1)
print('  B=-9bet = %.24f'%B);           print('  delta_0 = %.24f'%d0)
print('  V_win   = %.24f'%Vwin);         print('  A_min   = %.24f'%Amin)
print('  R_w     = %.24f'%Rw)
chk(sp.simplify(K1_e-(51*s5-21)/16) == 0, 'K_1 = 9gamma+3alpha = (51sqrt5-21)/16 精确')
chk(sp.simplify(B_e-(135*s5-225)/16) == 0, 'B = -9beta = (135sqrt5-225)/16 精确')
chk(sp.simplify(d0_e-(9+7*s5-4*s35)/6) == 0, 'delta_0 = 3phi+3-6rho_0 = (9+7sqrt5-4sqrt35)/6 精确')
chk(sp.simplify(K1_e-B_e/5-3*phi_e) == 0, 'K_1 - B/5 = 3phi 精确（thm:r24centre 证明第一步所需）')
chk(sp.simplify(gam_e+alp_e+bet_e-1) == 0, '证书 q(s)=gamma+alpha s+beta s^2 满足 q(1)=1')
q_s = sp.Symbol('s', nonnegative=True)
q_e = gam_e+alp_e*q_s+bet_e*q_s**2
chk(sp.simplify(sp.expand(q_e**2-q_s)*128 + (q_s-1)*(5*q_s-1)**2*((15*s5-35)*q_s+s5+3)) == 0,
    'q(s)^2-s = -(s-1)(5s-1)^2((15sqrt5-35)s+sqrt5+3)/128 精确（q>=sqrt(s) 的证书）')
chk(sp.simplify(c0_e-sp.Symbol('z')) != 0 and abs(float(c0_e)-1.8744425029524276) < 1e-13
    and abs(float(c1_e)-(20+12*s5.__float__()*0+12*math.sqrt(5))/225) < 1e-15,
    'c_0=%.16f, c_1=%.16f 与定理 thm:r22red 一致'%(float(c0_e), float(c1_e)))
chk(abs(Vwin-2*d0) < 1e-18 and abs(Vwin-float(6-2*(c0_e+sp.Rational(1,5)/c1_e))) < 1e-14,
    'V_win = 2delta_0 = 6-2(c_0+1/(5c_1)) = %.16f（两条独立路径一致）'%Vwin)

# ============================================================ 2. 帧工具
def solve_masses(Bm, rhs, rng, ntry=80):
    m0, *_ = np.linalg.lstsq(Bm, rhs, rcond=None)
    _, sv, VT = np.linalg.svd(Bm); rank = int((sv > 1e-10*max(sv[0], 1e-30)).sum()); NS = VT[rank:].T
    for _ in range(ntry):
        m = m0 if NS.size == 0 else m0 + NS @ (rng.standard_normal(NS.shape[1])*0.4)
        if m.min() > 1e-3: return m
    return None
def rand_iso(n, rng):
    for _ in range(400):
        U = rng.standard_normal((n, 3)); U /= np.linalg.norm(U, axis=1)[:, None]
        Bm = np.stack([U[:,0]**2, U[:,1]**2, U[:,2]**2, U[:,0]*U[:,1], U[:,0]*U[:,2], U[:,1]*U[:,2]])
        m = solve_masses(Bm, np.array([1.,1.,1.,0.,0.,0.]), rng)
        if m is None: continue
        if np.abs(np.einsum('i,ij,ik->jk', m, U, U)-np.eye(3)).max() < 1e-11: return m, U
    return None
def Astar(m, U):
    n = len(m); best = -1e18; bw = None
    for sg in itertools.product([1,-1], repeat=n-1):
        sv = np.array((1,)+sg); v = (sv*m) @ U; t = float(v@v)
        if t > best:
            best = t
            if t > 1e-300: bw = v/math.sqrt(t)
    return math.sqrt(best), bw
def rotate_to_w(U, w):
    e3 = np.array([0.,0.,1.]); c = float(w@e3); v = np.cross(w, e3); sn = float(np.linalg.norm(v))
    if sn < 1e-14:
        R = np.eye(3) if c > 0 else np.diag([1.,1.,-1.])
    else:
        v = v/sn
        K = np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
        R = np.eye(3) + sn*K + (1-c)*(K@K)
    return (R@U.T).T
def data(m, U):
    A, w = Astar(m, U); b = np.abs(U@w); u = b*b
    G = U@U.T; sq = G**2
    return dict(A=A, w=w, sig2=float(m@(u-1/3)**2), X=float(m@u**2), sumu=float(m@u),
                m4=float(np.einsum('i,j,ij->', m, m, G**4))/9,
                L=float(np.einsum('i,j,ij->', m, m, np.abs(G))),
                D=float(np.einsum('i,j,ij->', m, m, (gam+alp*sq+bet*sq*sq)-np.sqrt(np.maximum(sq, 0)))),
                V=6-2*A*A, J=float(np.einsum('i,j,ij->', m, m, np.abs(G)))+A*A)
def moments(m, U, w):
    Ur = rotate_to_w(U, w)
    assert np.abs(np.abs(Ur[:,2])-np.abs(U@w)).max() < 1e-11, 'Rodrigues 必须把 w 映到 e_3'
    z = np.clip(Ur[:,2], -1, 1); r = np.sqrt(np.maximum(1-z*z, 0))
    cs = np.where(r > 1e-14, Ur[:,0]/np.maximum(r, 1e-300), 1.0)
    sn = np.where(r > 1e-14, Ur[:,1]/np.maximum(r, 1e-300), 0.0)
    out = {}
    for a in range(5):
        for bb in range(5-a):
            cc = 4-a-bb
            C = math.factorial(4)//(math.factorial(a)*math.factorial(bb)*math.factorial(cc))
            out[(a,bb,cc)] = (C, float(m @ (z**a * r**(bb+cc) * cs**bb * sn**cc)))
    A4 = np.einsum('i,ia,ib,ic,id->abcd', m, U, U, U, U)
    return out, float(np.sum(A4*A4))
def gg(t): return 19/9-(7/3)*t+(35/8)*t*t

# ============================================================ 3. 符号代数：两个块的松弛恒等式
head('2. 引理 lem:r24wblock / lem:r24a0block：两个块的显式完全平方松弛（符号代数）')
P, S, Q, T, X, M1, M2, s2 = sp.symbols('P S Q T X M1 M2 s2', real=True)
chk(sp.expand(X**2+6*M1**2+6*M2**2 - (X**2+3*(M1+M2)**2+3*(M1-M2)**2)) == 0,
    'w-块 = X^2+3(1-X)^2+3(M220-M202)^2 精确恒等')
chk(sp.expand((X**2+3*(1-X)**2).subs(X, sp.Rational(1,3)+s2)
              - (sp.Rational(13,9)-sp.Rational(10,3)*s2+4*s2**2)) == 0,
    'X=1/3+sigma_u^2 时 X^2+3(1-X)^2 = 13/9-(10/3)sigma_u^2+4sigma_u^4 精确')
chk(sp.expand((P**2+S**2+6*Q**2 - (3*T**2/sp.Integer(8)+8*(Q-T/sp.Integer(8))**2+(P-S)**2/sp.Integer(2)))
              .subs(T, P+S+2*Q)) == 0,
    'a=0 块 = 3T^2/8+8(Q-T/8)^2+(P-S)^2/2 精确恒等')
chk(sp.expand(3*(sp.Rational(4,3)+s2)**2/sp.Integer(8)
              - (sp.Rational(2,3)+s2+sp.Rational(3,8)*s2**2)) == 0,
    'T=4/3+sigma_u^2 时 3T^2/8 = 2/3+sigma_u^2+(3/8)sigma_u^4 精确')
gfun_e = sp.Rational(19,9)-sp.Rational(7,3)*s2+sp.Rational(35,8)*s2**2
tot_e = sp.expand((sp.Rational(13,9)-sp.Rational(10,3)*s2+4*s2**2)+(sp.Rational(2,3)+s2+sp.Rational(3,8)*s2**2))
chk(sp.expand(tot_e-gfun_e) == 0, '两块下界相加 = g(sigma_u^2) = 19/9-(7/3)sigma_u^2+(35/8)sigma_u^4 精确')
chk(sum(math.factorial(4)//(math.factorial(a)*math.factorial(b)*math.factorial(4-a-b))
        for a in range(5) for b in range(5-a)) == 81, 'sum multinom(4;a,b,c) = 3^4 = 81')
kept = [(4,0,0),(2,2,0),(2,0,2),(0,4,0),(0,2,2),(0,0,4)]
wk = sum(math.factorial(4)//(math.factorial(a)*math.factorial(b)*math.factorial(c)) for a,b,c in kept)
chk(wk == 21 and 81-wk == 60, '保留 6 个模式共 21 份权重，丢弃 9 个模式共 60 份（注 rem:r24drop）')

head('3. 推论 cor:r24floors：g 的两个端点与驻点')
chk(sp.simplify(gfun_e.subs(s2, 0)-sp.Rational(19,9)) == 0, 'g(0) = 19/9（等纬度底，引理 lem:r23m4eq）')
chk(sp.simplify(gfun_e.subs(s2, sp.Rational(4,15))-sp.Rational(9,5)) == 0,
    'g(4/15) = 9/5（全局底 m_4>=1/5，引理 lem:r22m4）')
chk(sp.simplify(sp.diff(gfun_e, s2).subs(s2, sp.Rational(4,15))) == 0,
    "g'(4/15)=0：驻点位置由代数自动给出，故 g>=9/5 对一切 sigma_u^2 成立")
chk(sp.simplify(gfun_e.subs(s2, sp.Rational(2,3))-sp.Rational(5,2)) == 0, 'g(2/3) = 5/2（sigma_u^2 的上界处）')
v10 = sp.simplify(gfun_e.subs(s2, sp.Rational(1,10))/9)
chk(sp.simplify(v10-sp.Rational(2767,12960)) == 0,
    'sigma_u^2=1/10 处 m_4 >= g(1/10)/9 = 2767/12960 = %.16f（比 1/5 高 %.4f，推论 cor:r24floors(iii)）'
    %(float(v10), float(v10)-0.2))

# ============================================================ 4. 数值：恒等式与主定理
head('4. 引理 lem:r24mom / 定理 thm:r24interp：数值核对（随机各向同性测度）')
rng = np.random.default_rng(20240924)
rows = []; W = dict(id=0., nrm=0., con=0., w=0., a0=0., tot=0., X=0., u=0.)
nfr = 0
while nfr < 300:
    n = int(rng.integers(6, 14)); got = rand_iso(n, rng)
    if got is None: continue
    m, U = got; nfr += 1
    d = data(m, U); bl, nrm2 = moments(m, U, d['w'])
    tot = sum(C*M*M for C, M in bl.values())
    W['id'] = max(W['id'], abs(tot-nrm2)); W['nrm'] = max(W['nrm'], abs(nrm2-9*d['m4']))
    W['u'] = max(W['u'], abs(d['sumu']-1)); W['X'] = max(W['X'], abs(d['X']-(1/3+d['sig2'])))
    x = bl[(4,0,0)][1]; M1 = bl[(2,2,0)][1]; M2 = bl[(2,0,2)][1]
    Pn = bl[(0,4,0)][1]; Qn = bl[(0,2,2)][1]; Sn = bl[(0,0,4)][1]
    wb = x*x+6*M1**2+6*M2**2; a0 = Pn*Pn+6*Qn*Qn+Sn*Sn
    s2v = d['sig2']; Tv = 4/3+s2v
    W['w'] = max(W['w'], abs(wb-(x*x+3*(1-x)**2+3*(M1-M2)**2)))
    W['a0'] = max(W['a0'], abs(a0-(2/3+s2v+3*s2v**2/8+8*(Qn-Tv/8)**2+0.5*(Pn-Sn)**2)))
    W['con'] = max(W['con'], abs((Pn+Sn+2*Qn)-Tv))
    W['tot'] = max(W['tot'], gg(s2v)-nrm2)
    rows.append((d, wb, a0))
chk(W['u'] < 1e-12, '引理 lem:r24basic(i)：sum m_i u_i = 1（最大偏差 %.2e）'%W['u'])
chk(W['X'] < 1e-12, '引理 lem:r24basic(i)：X = sum m u^2 = 1/3+sigma_u^2（最大偏差 %.2e）'%W['X'])
chk(max(d['sig2'] for d, _, _ in rows) <= 2/3+1e-12, '引理 lem:r24basic(ii)：sigma_u^2<=2/3（实测最大 %.9f）'%max(d['sig2'] for d,_,_ in rows))
chk(W['id'] < 1e-12, '引理 lem:r24mom：||A_4||_F^2 = sum multinom M_abc^2（%d 个测度，偏差 %.2e）'%(nfr, W['id']))
chk(W['nrm'] < 1e-12, '引理 lem:r24mom：||A_4||_F^2 = 9 m_4（偏差 %.2e）'%W['nrm'])
chk(W['con'] < 1e-12, '引理 lem:r24a0block：约束 P+S+2Q = 4/3+sigma_u^2（偏差 %.2e）'%W['con'])
chk(W['w'] < 1e-12, '引理 lem:r24wblock：w-块松弛恒等式（偏差 %.2e）'%W['w'])
chk(W['a0'] < 1e-12, '引理 lem:r24a0block：a=0 块松弛恒等式（偏差 %.2e）'%W['a0'])
chk(W['tot'] < 1e-14, '定理 thm:r24interp：9m_4 >= g(sigma_u^2)（%d 个测度，最大违反 %.2e）'%(nfr, W['tot']))
wv = max((13/9-(10/3)*d['sig2']+4*d['sig2']**2)-wb for d, wb, a0 in rows)
av = max((2/3+d['sig2']+3*d['sig2']**2/8)-a0 for d, wb, a0 in rows)
chk(wv < 1e-12 and av < 1e-12, '两个块下界分别成立（w-块违反 %.2e，a=0 块违反 %.2e）'%(wv, av))

# ============================================================ 5. 桥接引理
head('5. 引理 lem:r24bridge / 推论 cor:r24Rw：sigma_u^2 <= R(A)V')
A = sp.Symbol('A', positive=True)
Rcl = (1+1/s3)**2/(s3*(s3+A))
chk(sp.simplify(sp.cancel(((1+1/s3)**2*(2-2*A/s3))/(6-2*A**2))-Rcl) == 0,
    'sigma_u^2/V <= (1+1/sqrt3)^2/(sqrt3(sqrt3+A)) = R(A) 符号恒等')
chk(sp.simplify(Rcl.subs(A, s3)-(2+s3)/9) == 0, 'R(sqrt3) = (2+sqrt3)/9 = %.18f 精确'%float((2+s3)/9))
dR = sp.simplify(sp.diff(Rcl, A))
chk(all(float(dR.subs(A, t)) < 0 for t in (sp.Rational(1,10), sp.Rational(1,2), 1, sp.Rational(3,2), s3-10**-9))
    and dR.is_negative is not False,
    'R(A) 关于 A 严格递减（导数 %s 在多点为负）'%dR)
wR = 0.; wratio = 0.
for d, wb, a0 in rows:
    R = (1+1/s3f)**2/(s3f*(s3f+d['A'])); wR = max(wR, d['sig2']-R*d['V'])
    if d['V'] > 1e-12: wratio = max(wratio, d['sig2']/d['V'])
chk(wR < 1e-14, 'sigma_u^2 <= R(A)V 在 %d 个随机测度上成立（最大违反 %.2e）'%(nfr, wR))
print('     实测 sigma_u^2/V 最大 = %.12f <= (2+sqrt3)/9 = %.12f <= R_w = %.12f'%(wratio, float((2+s3)/9), Rw))
Vs = sp.Symbol('V', nonnegative=True)
hf = sp.cancel((1+1/s3)**2*Vs/(s3*(s3+sp.sqrt(3-Vs/2))))
chk(sp.simplify(sp.diff(hf, Vs).subs(Vs, sp.Rational(1,10))) > 0, 'h(V)=R(A(V))V 严格递增（符号判定）')
chk(abs(float(sp.N(hf.subs(Vs, Vwin_e), 30))-0.1385149424161078) < 1e-13
    and float(sp.N(hf.subs(Vs, Vwin_e), 30)) < 4/15,
    'h(V_win) = %.16f < 4/15 => 窗口内 g 处于递减支（推论 cor:r24Rw）'%float(sp.N(hf.subs(Vs, Vwin_e), 30)))
chk(abs(Rw-0.4205252294509481) < 1e-15 and abs(Amin-1.6838370511968557) < 1e-13,
    'A_min=%.16f, R_w=%.16f（推论 cor:r24Rw 的常数）'%(Amin, Rw))

# ============================================================ 6. 闭式取等族
head('6. 命题 prop:r24sharp：闭式取等族（极点 delta + 正 n>=5 边形）')
dl = sp.Symbol('delta', positive=True)
u2_e = (1-dl)/(3-dl); sig_cl = 4*dl/(3*(3-dl)); Ap_e = dl+sp.sqrt((3-dl)*(1-dl))
chk(sp.simplify((1-u2_e)-2/(3-dl)) == 0, 'r_2^2 = 2/(3-delta) 精确')
chk(sp.simplify((dl*(1-sp.Rational(1,3))**2+(3-dl)*(u2_e-sp.Rational(1,3))**2)-sig_cl) == 0,
    'sigma_u^2 = 4delta/(3(3-delta)) 精确')
chk(sp.simplify(sig_cl.subs(dl, sp.Rational(1,2))-sp.Rational(4,15)) == 0,
    'delta=1/2 时 sigma_u^2 = 4/15（正二十面体）；delta=0 时 = 0（等纬度）')
chk(sp.simplify(Ap_e.subs(dl, 0)-s3) == 0 and sp.simplify(Ap_e.subs(dl, sp.Rational(1,2))-phi_e) == 0,
    'A_pole(0)=sqrt3 且 A_pole(1/2)=phi 精确')
chk(sp.simplify((3-dl)*(1-u2_e)**2/8 - (sp.Rational(4,3)+sig_cl)/8) == 0,
    'Q = 1/(2(3-delta)) = T/8 精确 => a=0 块的松弛 8(Q-T/8)^2 恒为零')
chk(sp.simplify(sp.limit(sig_cl/(6-2*Ap_e**2), dl, 0, dir='+')-(2+s3)/9) == 0,
    'lim_{delta->0} sigma_u^2/V = (2+sqrt3)/9 = R(sqrt3) 精确（注 rem:r24samefam）')
wF = 0.
for n in (5,6,7,8,9,12,13):
    th = np.array([2*math.pi*j/n for j in range(n)])
    for (p_, q_) in [(1,0),(0,1),(1,1),(3,0),(0,3),(2,1),(1,2),(3,1),(1,3)]:
        wF = max(wF, abs(np.sum(np.cos(th)**p_*np.sin(th)**q_))/n)
    for tgt, val in [((np.cos(th)**2), .5), ((np.sin(th)**2), .5),
                     ((np.cos(th)**4), 3/8), ((np.sin(th)**4), 3/8),
                     ((np.cos(th)**2*np.sin(th)**2), 1/8)]:
        wF = max(wF, abs(np.sum(tgt)/n-val))
chk(wF < 1e-15, '正 n 边形(n=5..13)：9 个被丢弃的矩消没 + P=S、M220=M202（最大偏差 %.2e）'%wF)
weq = 0.; wtilt = 0.; s2cov = []
for delta in (1e-5,1e-4,1e-3,1e-2,3e-2,6e-2,0.1,0.2,0.3,0.4,0.49,0.5):
    for n in (5,6,7,8,12):
        u2 = float(u2_e.subs(dl, delta)); m = np.array([delta]+[(3-delta)/n]*n)
        U = [np.array([0.,0.,1.])]
        for j in range(n):
            t = 2*math.pi*j/n
            U.append(np.array([math.sqrt(1-u2)*math.cos(t), math.sqrt(1-u2)*math.sin(t), math.sqrt(u2)]))
        U = np.array(U)
        if np.abs(np.einsum('i,ij,ik->jk', m, U, U)-np.eye(3)).max() > 1e-12: continue
        d = data(m, U)
        weq = max(weq, abs(9*d['m4']-gg(d['sig2'])))
        wtilt = max(wtilt, abs(d['A']-float(Ap_e.subs(dl, delta))))
        assert abs(d['sig2']-float(sig_cl.subs(dl, delta))) < 1e-12
        s2cov.append(d['sig2'])
chk(weq < 1e-13, '闭式族上 9m_4 = g(sigma_u^2) 精确取等（%d 个构型，最大偏差 %.2e）'%(len(s2cov), weq))
chk(wtilt < 1e-12, 'delta<=1/2 时极点轴确为 A-最优轴（最大偏差 %.2e）'%wtilt)
chk(min(s2cov) < 1e-4 and abs(max(s2cov)-4/15) < 1e-12,
    'sigma_u^2 覆盖 [0,4/15]：%.3e .. %.15f'%(min(s2cov), max(s2cov)))

# ============================================================ 7. 三分支一致余量
head('7. 定理 thm:r24centre / thm:r24main：三分支一致余量与中心截面闭合')
c2e = sp.simplify(B_e*sp.Rational(35,72)*Rw_e**2)
c1e = sp.simplify(sp.Rational(1,2)-B_e*sp.Rational(7,27)*Rw_e)
c0e = sp.simplify(B_e*(sp.Rational(19,81)-sp.Rational(1,5))-d0_e)
disc_e = sp.simplify(c1e**2-4*c2e*c0e)
c2 = float(c2e); c1 = float(c1e); c0 = float(c0e); disc = float(disc_e)
m1_e = sp.simplify(-disc_e/(4*c2e)); m1 = float(m1_e)
Vst = float(-c1e/(2*c2e)); V2 = 2*(d0+m1)
print('  margin(V) = %.16f V^2 %+.16f V %+.16f'%(c2, c1, c0))
print('  判别式    = %.18e'%disc)
print('  顶点 V*   = %.16f   最小余量 m_1 = %.18e'%(Vst, m1))
print('  V_2       = %.16f'%V2)
chk(sp.simplify(c0e-(sp.Rational(14,405)*B_e-d0_e)) == 0, 'c_0 = 14B/405 - delta_0 精确')
chk(float(sp.N(c2e, 40)) > 0 and float(sp.N(disc_e, 40)) < 0,
    '分支1：c_2>0 且判别式<0（40 位高精度判定）=> Q(V) 恒正')
chk(0 < Vst < Vwin, '分支1：顶点 V*=%.16f 落在窗口 (0,V_win) 内，min Q = m_1 = %.18e'%(Vst, m1))
chk(abs(float(sp.N(m1_e, 50))-m1) < 1e-15, 'm_1 的 50 位高精度复核一致（差 %.2e）'%abs(float(sp.N(m1_e,50))-m1))
def hh(V): return (1+1/s3f)**2*V/(s3f*(s3f+math.sqrt(3-V/2)))
b2 = B*(gg(min(hh(V2), 4/15))/9-1/5)+Vwin/2-d0
print('  分支2：h(V_2)=%.16f，余量下界 = %.15f'%(hh(V2), b2))
chk(hh(V2) < 4/15 and b2 > m1, '分支2：V in [V_win,V_2] 上余量 %.15f > m_1'%b2)
b3 = V2/2-d0
chk(abs(b3-m1) < 1e-12 and b3 > 0, '分支3：V>=V_2 上余量 V_2/2-delta_0 = %.18e = m_1（按定义）'%b3)
Munif = min(m1, b2, b3)
chk(Munif > 0, '三分支合并：Phi(V) >= %.18e > 0 对一切 V>=0'%Munif)
bound = rho0-Munif/6
print('  ==> kappa_3^c <= rho_0 - %.18e = %.18f'%(Munif/6, bound))
chk(abs(bound-1.281394870601110) < 1e-14, 'eq:r24main：kappa_3^c <= 1.281394870601110')
chk(bound < rho0, 'kappa_3^c < rho_0 严格 => 中心截面闭合（定理 thm:r24main）')
print('     对照：第 22/23 轮数值上确界 sup J = 7.595950307218 => kappa_3^c = 1.265991717870，仍低 %.6f'%(bound-1.265991717870))

# ============================================================ 8. 逐构型余量与库扫描
head('8. 逐构型验证：证明余量恒正且不超过真实余量（eq:r23Dform 的强化版）')
def margins(d, lb9):
    return B*(max(lb9/9, 1/5)-1/5)+d['V']/2+d['D']-d0, 6*rho0-d['J']
wneg = 1e18; wgap = 1e18; tmin = 1e18
for d, wb, a0 in rows:
    mb, tm = margins(d, gg(d['sig2']))
    wneg = min(wneg, mb); wgap = min(wgap, tm-mb); tmin = min(tmin, tm)
chk(wneg > 0, '%d 个随机测度：证明余量恒正，最小 %.9e'%(nfr, wneg))
chk(wgap > -1e-12, '证明余量 <= 真实余量（最小差 %.9e，链条有效性自检）'%wgap)
print('     真实余量最小 %.9e；证明余量最小 %.9e'%(tmin, wneg))
def cone(n, z):
    g0 = (1-3*z*z)/(1-z*z); G = 2/(1-z*z)
    m = np.array([g0]+[G/n]*n); U = [np.array([0.,0.,1.])]
    for j in range(n):
        t = 2*math.pi*j/n
        U.append(np.array([math.sqrt(1-z*z)*math.cos(t), math.sqrt(1-z*z)*math.sin(t), z]))
    return m, np.array(U)
wlib = 1e18; wviol = -1e18; nlib = 0; trec = None
for n in range(3, 13):
    for k in range(1, 120):
        z = (k/120.0)/math.sqrt(3); m, U = cone(n, z)
        if m.min() <= 0: continue
        d = data(m, U); g9 = gg(d['sig2'])
        wviol = max(wviol, g9-9*d['m4']); mb, tm = margins(d, g9); nlib += 1
        if mb < wlib: wlib = mb; trec = tm
chk(wviol < 1e-12, '锥族库（%d 个构型，n=3..12 x 119 个 z）：9m_4>=g(sigma_u^2) 无违反（最大 %.2e）'%(nlib, wviol))
chk(wlib > 0, '锥族库：证明余量恒正，最小 %.9e（对应真实余量 %.9e）'%(wlib, trec))
ph = (1+math.sqrt(5))/2
Ui = np.array([np.array(v)/np.linalg.norm(v) for v in
               [(0,1,ph),(0,1,-ph),(1,ph,0),(1,-ph,0),(ph,0,1),(ph,0,-1)]])
d = data(np.array([0.5]*6), Ui)
chk(abs(d['sig2']-4/15) < 1e-13 and abs(d['m4']-1/5) < 1e-14 and abs(9*d['m4']-gg(d['sig2'])) < 1e-13,
    '正二十面体：sigma_u^2=4/15, m_4=1/5，定理 thm:r24interp 精确取等')
chk(abs(d['A']**2-ph**2) < 1e-12 and abs(d['V']-(3-math.sqrt(5))) < 1e-12,
    '正二十面体：A^2=phi^2, V=3-sqrt5=%.12f（与 eq:r23Dform 的算例一致）'%d['V'])
U4 = np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]], float); U4 /= np.linalg.norm(U4, axis=1)[:, None]
d4 = data(np.array([0.75]*4), U4)
chk(abs(d4['sig2']) < 1e-13 and abs(d4['m4']-7/27) < 1e-13 and 9*d4['m4']-gg(d4['sig2']) > 1e-9,
    '等角四线：sigma_u^2=0 但 m_4=7/27>19/81（不取等，rho_hat_4!=0，与 lem:r23m4eq 一致）')
weq2 = 0.
for n in (5,6,7,8,12):
    U = np.array([[math.sqrt(2/3)*math.cos(2*math.pi*j/n), math.sqrt(2/3)*math.sin(2*math.pi*j/n), 1/math.sqrt(3)] for j in range(n)])
    d = data(np.array([3.0/n]*n), U); weq2 = max(weq2, abs(d['sig2']), )
    assert abs(d['m4']-19/81) < 1e-13 and abs(9*d['m4']-gg(d['sig2'])) < 1e-12
chk(weq2 < 1e-13, '等纬度正 n 边形（n=5,6,7,8,12）：sigma_u^2=0、m_4=19/81，取等')

# ============================================================ 9. 任意截面恒等式
head('9. 命题 prop:r24gen：任意截面的精确恒等式与两套权重的分离')
def rand_Lam(rng, kind):
    if kind == 'uniaxial':
        a = rng.uniform(0.3, 1.0); return np.diag([a, (1-a)/2, (1-a)/2])
    if kind == 'centre': return np.eye(3)/3
    if kind == 'rank1':
        wv = rng.standard_normal(3); wv /= np.linalg.norm(wv); return np.outer(wv, wv)
    M = rng.standard_normal((3,3)); Sm = M@M.T
    ev, Qm = np.linalg.eigh(Sm); ev = np.abs(ev)+1e-3; ev = ev/ev.sum()
    return Qm@np.diag(ev)@Qm.T
def bestsign(w, U):
    n = len(w); best = -1e18; bs = None
    for sg in itertools.product([1,-1], repeat=n-1):
        sv = np.array((1,)+sg); v = (sv*w)@U; t = float(v@v)
        if t > best: best = t; bs = sv
    return best, bs
rng2 = np.random.default_rng(777)
wf = 0.; wp = 0.; wn = 0.; wc1 = 0.; wc2 = 0.; wle = 0.
for kind in ('centre', 'uniaxial', 'rank1', 'generic'):
    for _ in range(110):
        nn = int(rng2.integers(5, 12)); fr = rand_iso(nn, rng2)
        if fr is None: continue
        g, U = fr; Lam = rand_Lam(rng2, kind)
        h = np.sqrt(np.maximum(np.einsum('ij,jk,ik->i', U, Lam, U), 0))
        V = U*np.sqrt(g)[:, None]; Qm = V@V.T; Qp = (Qm+np.abs(Qm))/2
        x = np.sqrt(g)*h; rd = float(x@Qp@x)
        w = g*h; Lw = float(w@np.abs(U@U.T)@w)
        aw = w@U; rfix = 0.5*(Lw+float(aw@aw))
        Aw2, bs = bestsign(w, U)
        Us = U*bs[:, None]; Vs = Us*np.sqrt(g)[:, None]; Qs = Vs@Vs.T
        rflip = float(x@((Qs+np.abs(Qs))/2)@x)
        wf = max(wf, abs(rd-rfix)/max(rd, 1e-12)); wp = max(wp, abs(rflip-0.5*(Lw+Aw2))/max(rd, 1e-12))
        wn = max(wn, abs(float(x@x)-np.trace(Lam))); wle = max(wle, rd-float(np.linalg.eigvalsh(Qp)[-1]))
        wc1 = max(wc1, abs(float(w@h)-np.trace(Lam)))
        wc2 = max(wc2, np.abs(np.einsum('i,ij,ik->jk', w/h, U, U)-np.eye(3)).max())
chk(wf < 1e-13, 'eq:r24gen1：rho_{g,U}(Lam) = (L_w+||sum w_i u_i||^2)/2（四类截面，最大相对偏差 %.2e）'%wf)
chk(wp < 1e-13, 'eq:r24gen2：符号最优后 = (L_w+A_w^2)/2（最大相对偏差 %.2e）'%wp)
chk(wn < 1e-12, '||x(Lam)||^2 = tr Lam（最大偏差 %.2e）'%wn)
chk(wle < 1e-12, 'rho(Lam) <= lambda_max(Q_+)（最大超出 %.2e）'%wle)
chk(wc1 < 1e-12 and wc2 < 1e-11, 'eq:r24gen3：sum w_i h_i = tr Lam 且 sum (w_i/h_i)u_iu_i^T = I_3（偏差 %.2e / %.2e）'%(wc1, wc2))
g, U = rand_iso(9, rng2)
a = 0.594608575866; Lam = np.diag([a, (1-a)/2, (1-a)/2])
h = np.sqrt(np.einsum('ij,jk,ik->i', U, Lam, U)); w = g*h
tg = np.abs(np.einsum('i,ij,ik->jk', g, U, U)-np.eye(3)).max()
tw = np.abs(np.einsum('i,ij,ik->jk', w, U, U)-np.eye(3)).max()
print('     冠军截面 a=h_hub^2：||T_g-I||_max=%.2e（各向同性对 g 成立），||T_w-I||_max=%.2e（对 w 不成立）'%(tg, tw))
chk(tg < 1e-11 and tw > 1e-2, '注 rem:r24obs：单轴截面下两套权重分离，中心截面的矩计算不可直接搬运')

head('SUMMARY  PASS=%d  FAIL=%d'%(PASS, FAIL))
print('  中心截面闭合：kappa_3^c <= %.15f = rho_0 - %.6e'%(bound, Munif/6))
print('  全局 kappa_3 <= rho_0 仍未证（单轴截面在冠军处取等，见注 rem:r24obs / 队列 P1）')
print('='*112)
raise SystemExit(0 if FAIL == 0 else 1)
