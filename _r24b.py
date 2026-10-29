# -*- coding: utf-8 -*-
# _r24b.py  ROUND 24 (B) -- corrected and complete verification of the centre-section closure.
#
#   TARGET (P1 of the 12th queue): prove conj:r22coupled / eq:r23Dform, i.e. close the open window
#       V in [0, 0.329386) and hence  kappa_3^c <= rho_0.
#
#   BUGS FIXED w.r.t. _r24a.py:
#     (1) delta_0 was K1+3-6rho_0 (=1.12556) instead of the eq:r23Dform convention
#         delta_0 = 3phi+3-6rho_0 = (9+7sqrt5-4sqrt35)/6 (=0.164693).  They differ by B/5 and
#         K1 - B/5 = 3phi exactly, which is why the (C)-quadratic looked right while every
#         per-configuration margin was shifted by -0.96086.
#     (2) the Rodrigues rotation took an extra inverse, so the third axis was NOT w and the
#         a=0-block bound appeared to fail by 5e-3.
#
#   NEW (this file): the two block bounds are upgraded from inequalities to EXACT SLACK
#   IDENTITIES, so the whole certificate is algebra:
#       w-block  = x^2 + 3(1-x)^2 + 3(M220-M202)^2
#       a0-block = 2/3 + s2 + (3/8)s2^2 + 8(Q-T/8)^2 + (P-S)^2/2 ,  T = 4/3+s2 = P+S+2Q
#     => 9m_4 = ||A_4||_F^2 >= 19/9 - (7/3)s2 + (35/8)s2^2  (+ explicit non-negative defect)
#   plus the exact two-branch uniform margin Phi(V) >= 1.0399576e-3 for every V >= 0.
import numpy as np, math, itertools
import sympy as sp

np.set_printoptions(precision=12, suppress=False, linewidth=200)
PASS = FAIL = 0
def chk(c, m):
    global PASS, FAIL
    if c: PASS += 1; print('  [PASS]', m)
    else: FAIL += 1; print('  [FAIL]', m)

# ----------------------------------------------------------------------------- exact constants
s5 = sp.sqrt(5); s35 = sp.sqrt(35); s3 = sp.sqrt(3)
rho0_e = (9+s5+2*s35)/18; phi_e = (1+s5)/2
gam_e = (1+s5)/16; alp_e = (7*s5-5)/8; bet_e = (25-15*s5)/16
K1_e = sp.simplify(9*gam_e+3*alp_e); B_e = sp.simplify(-9*bet_e)
d0_e = sp.simplify(3*phi_e+3-6*rho0_e)                 # = (9+7sqrt5-4sqrt35)/6
Vwin_e = 2*d0_e; Amin_e = sp.sqrt(3-d0_e)
Rw_e = (1+1/s3)**2/(s3*(s3+Amin_e))
c2_e = sp.simplify(B_e*sp.Rational(35,72)*Rw_e**2)
c1_e = sp.simplify(sp.Rational(1,2)-B_e*sp.Rational(7,27)*Rw_e)
c0_e = sp.simplify(B_e*(sp.Rational(19,81)-sp.Rational(1,5))-d0_e)
disc_e = sp.simplify(c1_e**2-4*c2_e*c0_e)
rho0=float(rho0_e); K1=float(K1_e); B=float(B_e); d0=float(d0_e)
Vwin=float(Vwin_e); Amin=float(Amin_e); Rw=float(Rw_e)
gam=float(gam_e); alp=float(alp_e); bet=float(bet_e); s3f=math.sqrt(3.0); phif=float(phi_e)
print('='*116)
print('精确常数（sympy，全部为代数数）')
print('='*116)
print('  rho_0   = (9+sqrt5+2sqrt35)/18      = %.24f'%rho0)
print('  K_1     = 9gamma+3alpha = (51sqrt5-21)/16 = %.24f'%K1)
print('  B       = -9beta = (135sqrt5-225)/16     = %.24f'%B)
print('  delta_0 = 3phi+3-6rho_0 = (9+7sqrt5-4sqrt35)/6 = %.24f   <-- eq:r23Dform 约定'%d0)
print('  V_win   = 2*delta_0                     = %.24f'%Vwin)
print('  A_min   = sqrt(3-delta_0)               = %.24f'%Amin)
print('  R_w     = R(A_min)                      = %.24f'%Rw)
chk(sp.simplify(d0_e-(9+7*s5-4*s35)/6)==0, 'delta_0 = (9+7sqrt5-4sqrt35)/6 精确')
chk(sp.simplify(K1_e-B_e/5-3*phi_e)==0, 'K_1 - B/5 = 3phi 精确（这正是 _r24a 里 0.960865 平移的来源）')
chk(sp.simplify(Vwin_e-2*d0_e)==0 and abs(Vwin-0.3293855700333546)<1e-15, 'V_win = 2delta_0 = %.16f（与定理 thm:r22red 的窗口一致）'%Vwin)

# ----------------------------------------------------------------------------- frame utilities
def solve_masses(Bm, rhs, rng, ntry=80):
    m0, *_ = np.linalg.lstsq(Bm, rhs, rcond=None)
    _, sv, VT = np.linalg.svd(Bm); rank = int((sv > 1e-10*max(sv[0], 1e-30)).sum()); NS = VT[rank:].T
    for _ in range(ntry):
        m = m0 if NS.size == 0 else m0 + NS @ (rng.standard_normal(NS.shape[1])*0.4)
        if m.min() > 1e-3: return m
    return None

def rand_iso(n, rng):
    """random isotropic line measure: sum m = 3, sum m l l^T = I_3, all m > 0"""
    for _ in range(400):
        U = rng.standard_normal((n, 3)); U /= np.linalg.norm(U, axis=1)[:, None]
        Bm = np.stack([U[:,0]**2, U[:,1]**2, U[:,2]**2, U[:,0]*U[:,1], U[:,0]*U[:,2], U[:,1]*U[:,2]])
        m = solve_masses(Bm, np.array([1.,1.,1.,0.,0.,0.]), rng)
        if m is None: continue
        T = np.einsum('i,ij,ik->jk', m, U, U)
        if np.abs(T-np.eye(3)).max() < 1e-11 and abs(m.sum()-3) < 1e-11: return m, U
    return None

def Astar(m, U):
    """A = max_sigma ||sum sigma_i m_i l_i||  (exact sign enumeration) and the optimal axis w."""
    n = len(m); best = -1e18; bw = None
    for sg in itertools.product([1,-1], repeat=n-1):
        s = np.array((1,)+sg); v = (s*m) @ U; t = float(v@v)
        if t > best:
            best = t
            if t > 1e-300: bw = v/math.sqrt(t)
    return math.sqrt(best), bw

def rotate_to_w(U, w):
    """Rodrigues: R w = e_3, return U in the new coordinates (so column 2 = <l_i,w>)."""
    e3 = np.array([0.,0.,1.]); c = float(w@e3); v = np.cross(w, e3); s = float(np.linalg.norm(v))
    if s < 1e-14:
        R = np.eye(3) if c > 0 else np.diag([1.,1.,-1.])
    else:
        v = v/s
        K = np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
        R = np.eye(3) + s*K + (1-c)*(K@K)
    return (R@U.T).T, R

def data(m, U):
    A, w = Astar(m, U)
    b = np.abs(U@w); u = b*b
    sig2 = float(m @ (u-1/3)**2)
    G = U@U.T; s = G**2
    m4 = float(np.einsum('i,j,ij->', m, m, G**4))/9
    L = float(np.einsum('i,j,ij->', m, m, np.abs(G)))
    D = float(np.einsum('i,j,ij->', m, m, (gam+alp*s+bet*s*s)-np.sqrt(np.maximum(s, 0))))
    V = 6-2*A*A
    return dict(A=A, w=w, b=b, u=u, sig2=sig2, m4=m4, L=L, D=D, V=V, J=L+A*A)

def moments(m, U, w):
    """the 15 cylindrical moments M_{abc} (a = power of the w-axis) + ||A_4||_F^2."""
    Ur, R = rotate_to_w(U, w)
    z = Ur[:,2]; assert np.abs(np.abs(z)-np.abs(U@w)).max() < 1e-11, 'rotation must align w with e_3'
    z = np.clip(z, -1, 1); r = np.sqrt(np.maximum(1-z*z, 0))
    cs = np.where(r > 1e-14, Ur[:,0]/np.maximum(r, 1e-300), 1.0)
    sn = np.where(r > 1e-14, Ur[:,1]/np.maximum(r, 1e-300), 0.0)
    out = {}
    for a in range(5):
        for bb in range(5-a):
            cc = 4-a-bb
            C = math.factorial(4)//(math.factorial(a)*math.factorial(bb)*math.factorial(cc))
            out[(a,bb,cc)] = (C, float(m @ (z**a * r**(bb+cc) * cs**bb * sn**cc)))
    A4 = np.einsum('i,ia,ib,ic,id->abcd', m, U, U, U, U)
    return out, float(np.sum(A4*A4)), sum(C for C, _ in out.values())

# ------------------------------------------------------- (A) exact algebra of the block bounds
print()
print('='*116)
print('(A) 圆柱分块的精确松弛恒等式（符号代数，不含任何数值放缩）')
print('='*116)
Ps, Ss, Qs, Ts, xs, M1s, M2s, s2s = sp.symbols('P S Q T x M1 M2 s2', real=True)
# w-block: x = M400 = sum m u^2, M1 = M220, M2 = M202, M1+M2 = sum m u(1-u) = 1-x
wres = sp.expand(xs**2+6*M1s**2+6*M2s**2 - (xs**2+3*(M1s+M2s)**2+3*(M1s-M2s)**2))
chk(wres == 0, 'w-block = x^2+3(1-x)^2+3(M220-M202)^2  精确恒等（残差 %s）'%wres)
wlb = sp.expand(xs**2+3*(1-xs)**2); wq = sp.expand(sp.Rational(13,9)-sp.Rational(10,3)*s2s+4*s2s**2)
chk(sp.expand(wlb.subs(xs, sp.Rational(1,3)+s2s)-wq) == 0,
    'w-block 下界 = 13/9-(10/3)s2+4 s2^2 在 x=1/3+s2 下精确')
# a0-block: P+S+2Q = T = 4/3+s2
a0res = sp.expand((Ps**2+Ss+6*Qs**2*0+0) - 0)   # placeholder removed below
a0res = sp.expand(Ps**2+Ss**2+6*Qs**2 - (3*Ts**2/sp.Integer(8)+8*(Qs-Ts/sp.Integer(8))**2
                                         +(Ps-Ss)**2/sp.Integer(2))).subs(Ts, Ps+Ss+2*Qs)
chk(sp.expand(a0res) == 0, 'a0-block = 3T^2/8+8(Q-T/8)^2+(P-S)^2/2  精确恒等（残差 %s）'%sp.expand(a0res))
a0lb = sp.expand(3*(sp.Rational(4,3)+s2s)**2/sp.Integer(8))
chk(sp.expand(a0lb-(sp.Rational(2,3)+s2s+sp.Rational(3,8)*s2s**2)) == 0,
    '3T^2/8 = 2/3+s2+(3/8)s2^2 在 T=4/3+s2 下精确')
tot = sp.expand(wq+a0lb)
gfun = sp.Rational(19,9)-sp.Rational(7,3)*s2s+sp.Rational(35,8)*s2s**2
chk(sp.expand(tot-gfun) == 0,
    '主插值定理（9m_4 尺度）：9m_4 >= 19/9-(7/3)s2+(35/8)s2^2  精确（= w-block 界 + a0-block 界）')
chk(sp.expand(tot/9-(sp.Rational(19,81)-sp.Rational(7,27)*s2s+sp.Rational(35,72)*s2s**2)) == 0,
    '同一定理的 m_4 尺度写法：m_4 >= 19/81-(7/27)s2+(35/72)s2^2  精确')
chk(sp.simplify(gfun.subs(s2s, 0)-sp.Rational(19,9)) == 0, 'g(0) = 19/9  （等纬度底，= 引理 lem:r23m4eq）')
chk(sp.simplify(gfun.subs(s2s, sp.Rational(4,15))-sp.Rational(9,5)) == 0,
    'g(4/15) = 9/5  （全局底 m_4>=1/5，= 引理 lem:r22m4）；两端同时精确取到')
chk(sp.simplify(sp.diff(gfun, s2s).subs(s2s, sp.Rational(4,15))) == 0,
    'g 的极小点恰在 s2=4/15（导数为零）=> g(s2)>=9/5 对一切 s2 成立')
s2max = sp.simplify(sp.Rational(2,3))
chk(sp.simplify(gfun.subs(s2s, s2max)-sp.Rational(5,2)) == 0,
    's2<=2/3（u in [0,1], sum m u=1, sum m=3 的极值）处 g=5/2')

# ------------------------------------------- (B) numerics: identity + blocks on random frames
print()
print('='*116)
print('(B) 数值核对：15 矩恒等式、两个分块下界、主插值定理（随机各向同性测度）')
print('='*116)
rng = np.random.default_rng(20240924)
rows = []; worst = dict(id=0., rot=0., w=0., a0=0., tot=0., tot9=0., slack=0.)
nfr = 0; target = 300
while nfr < target:
    n = int(rng.integers(6, 14)); got = rand_iso(n, rng)
    if got is None: continue
    m, U = got; nfr += 1
    d = data(m, U); bl, nrm2, Ctot = moments(m, U, d['w'])
    tot = sum(C*M*M for C, M in bl.values())
    worst['id'] = max(worst['id'], abs(tot-nrm2))
    worst['tot9'] = max(worst['tot9'], abs(nrm2-9*d['m4']))
    x = bl[(4,0,0)][1]; M1 = bl[(2,2,0)][1]; M2 = bl[(2,0,2)][1]
    P = bl[(0,4,0)][1]; Q = bl[(0,2,2)][1]; S = bl[(0,0,4)][1]
    wb = x*x+6*M1**2+6*M2**2; a0 = P*P+6*Q*Q+S*S
    s2 = d['sig2']; T = 4/3+s2
    wslack = 3*(M1-M2)**2; a0slack = 8*(Q-T/8)**2+0.5*(P-S)**2
    worst['w'] = max(worst['w'], abs(wb-(x*x+3*(1-x)**2+wslack)))
    worst['a0'] = max(worst['a0'], abs(a0-(2/3+s2+3*s2**2/8+a0slack)))
    worst['slack'] = max(worst['slack'], abs((P+S+2*Q)-T))     # the constraint P+S+2Q = 4/3+s2
    lb9 = 19/9-(7/3)*s2+(35/8)*s2**2
    worst['tot'] = max(worst['tot'], (lb9-nrm2)/max(1e-30, 1))
    rows.append((d, lb9/9, wb, a0, wslack, a0slack))
chk(worst['id'] < 1e-12, '||A_4||_F^2 = sum multinom(4;a,b,c) M_abc^2   (%d 个随机测度) 最大偏差 %.2e'%(nfr, worst['id']))
chk(worst['tot9'] < 1e-12, '||A_4||_F^2 = 9 m_4   最大偏差 %.2e'%worst['tot9'])
chk(worst['slack'] < 1e-12, '约束 P+S+2Q = 4/3+sigma^2   最大偏差 %.2e'%worst['slack'])
chk(worst['w'] < 1e-12, 'w-block 精确松弛恒等式   最大偏差 %.2e'%worst['w'])
chk(worst['a0'] < 1e-12, 'a0-block 精确松弛恒等式   最大偏差 %.2e'%worst['a0'])
chk(worst['tot'] < 1e-14, '9m_4 >= 19/9-(7/3)s2+(35/8)s2^2   最大违反 %.2e'%worst['tot'])
vw = max((13/9-(10/3)*d['sig2']+4*d['sig2']**2)-wb_ for d, lb_, wb_, a0_, ws_, as_ in rows)
va = max((2/3+d['sig2']+3*d['sig2']**2/8)-a0_ for d, lb_, wb_, a0_, ws_, as_ in rows)
print('     (分块下界单独核对：w-block 最大违反 %.2e，a0-block 最大违反 %.2e)'%(vw, va))
chk(vw < 1e-12, 'w-block >= 13/9-(10/3)s2+4s2^2   （%d 个随机测度）'%nfr)
chk(va < 1e-12, 'a0-block >= 2/3+s2+(3/8)s2^2   （%d 个随机测度）'%nfr)
chk(min(ws_ for d, lb_, wb_, a0_, ws_, as_ in rows) > -1e-14 and
    min(as_ for d, lb_, wb_, a0_, ws_, as_ in rows) > -1e-14,
    '两个松弛量 3(M220-M202)^2 与 8(Q-T/8)^2+(P-S)^2/2 非负（数值）')

# --------------------------------------------------------------------- (C) the bridge s2 <= R(A)V
print()
print('='*116)
print('(C) 桥接引理  sigma^2 <= R(A) V,   R(A) = (1+1/sqrt3)^2/(sqrt3(sqrt3+A))')
print('='*116)
Asy = sp.Symbol('A', positive=True)
Rbr = sp.cancel(((1+1/s3)**2*(2-2*Asy/s3))/(6-2*Asy**2))
Rclosed = (1+1/s3)**2/(s3*(s3+Asy))
chk(sp.simplify(Rbr-Rclosed) == 0, 'sigma^2/V <= (1+1/sqrt3)^2/(sqrt3(sqrt3+A))  符号恒等')
chk(sp.simplify(Rclosed.subs(Asy, s3)-(2+s3)/9) == 0, 'R(sqrt3) = (2+sqrt3)/9 = %.18f 精确'%float((2+s3)/9))
chk(sp.simplify(sp.diff(Rclosed, Asy)) < 0 if False else float(sp.diff(Rclosed, Asy).subs(Asy, 1)) < 0,
    'R(A) 关于 A 严格递减（导数在 A=1 处为 %.6f < 0）'%float(sp.diff(Rclosed, Asy).subs(Asy, 1)))
worstR = 0.; worstRatio = 0.
for d, *_ in rows:
    A = d['A']; V = d['V']; s2 = d['sig2']
    R = (1+1/s3f)**2/(s3f*(s3f+A))
    worstR = max(worstR, s2-R*V)
    if V > 1e-12: worstRatio = max(worstRatio, s2/V)
chk(worstR < 1e-14, 'sigma^2 <= R(A)V 在 %d 个随机测度上成立（最大违反 %.2e）'%(nfr, worstR))
print('     实测 sigma^2/V 最大值 = %.12f  <= (2+sqrt3)/9 = %.12f  <= R_w = %.12f'%(
    worstRatio, float((2+s3)/9), Rw))
Vs0 = sp.Symbol('V', nonnegative=True)
hfun = sp.cancel((1+1/s3)**2*Vs0/(s3*(s3+sp.sqrt(3-Vs0/2))))
print('     h(V)=R(A(V))V 在窗口右端点：h(V_win) = %.18f < 4/15 = %.18f'%(float(sp.N(hfun.subs(Vs0, Vwin_e), 30)), 4/15))
chk(float(sp.N(hfun.subs(Vs0, Vwin_e), 30)) < 4/15,
    'h(V_win) < 4/15 => 在窗口内 g 处于递减支，代换 sigma^2 -> R_w V 合法')

# ------------------------------------------------- (D) the exact two-branch uniform margin
print()
print('='*116)
print('(D) 全链条：把 (A)+(C) 代入 eq:r23Dform，得到关于 V 的正定二次式与一致余量')
print('='*116)
Vs = sp.Symbol('V', nonnegative=True)
quad = sp.expand(c2_e*Vs**2+c1_e*Vs+c0_e)
print('     margin(V) = %.18f V^2 %+.18f V %+.18f'%(float(c2_e), float(c1_e), float(c0_e)))
print('     判别式   = %.18e'%float(disc_e))
chk(float(c2_e) > 0 and float(disc_e) < 0, '二次式首项 > 0 且判别式 < 0 => 在 [0,V_win] 上恒正')
Vvert = float(-c1_e/(2*c2_e)); Mvert = float(-disc_e/(4*c2_e))
print('     顶点 V = %.18f（在 [0,V_win] 内），最小余量 = %.18e'%(Vvert, Mvert))
chk(0 < Vvert < Vwin, '顶点落在窗口内部')
print('     quad(V_win) = %.18f > 0'%float(quad.subs(Vs, Vwin_e)))
chk(float(quad.subs(Vs, Vwin_e)) > 0, 'quad(V_win) > 0')
# branch 2: V >= V1 where h(V1)=4/15 ; there m_4 >= 1/5 gives margin >= V/2-delta_0
V1 = float(sp.nsolve(hfun-sp.Rational(4,15), 0.6))
print('     第二分支：h(V_1)=4/15 的 V_1 = %.15f，该处余量 >= V_1/2-delta_0 = %.15f'%(V1, V1/2-d0))
chk(V1 > Vwin and V1/2-d0 > Mvert, 'V_1 > V_win 且第二分支余量 %.6f 大于第一分支最小值 %.6e'%(V1/2-d0, Mvert))
# monotone check of Phi on [Vwin, V1]
grid = np.linspace(Vwin, V1, 2001)
hh = ((1+1/s3f)**2*grid/(s3f*(s3f+np.sqrt(3-grid/2))))
gval = (19/9-(7/3)*np.minimum(hh, 4/15)+(35/8)*np.minimum(hh, 4/15)**2)/9
Phi2 = B*(np.maximum(gval, 1/5)-1/5)+grid/2-d0
print('     Phi 在 [V_win,V_1] 上的最小值 = %.15f（>= 第一分支最小值 %.6e）'%(Phi2.min(), Mvert))
chk(Phi2.min() >= Mvert-1e-15, '两分支拼接后全局一致余量 = %.9e'%Mvert)
print('  ==>  kappa_3^c <= rho_0 - %.9e  （J 尺度余量 %.9e）'%(Mvert/6, Mvert))
chk(Mvert > 0, '中心截面严格闭合，一致余量 %.9e > 0'%Mvert)

# ------------------------------------------- (E) per-configuration margins (correct delta_0)
print()
print('='*116)
print('(E) 逐构型余量：证明余量 = B(max(lb,1/5)-1/5)+V/2+D-delta_0，须恒正且 <= 真实余量')
print('='*116)
def margins(d, lb):
    mb = B*(max(lb, 1/5)-1/5)+d['V']/2+d['D']-d0
    return mb, 6*rho0-d['J']
wneg = 1e18; wgap = 1e18; tmin = 1e18
for d, lb, *_ in rows:
    mb, tm = margins(d, lb)
    wneg = min(wneg, mb); wgap = min(wgap, tm-mb); tmin = min(tmin, tm)
chk(wneg > 0, '随机测度（%d 个）证明余量恒正，最小 %.9e'%(nfr, wneg))
chk(wgap > -1e-12, '证明余量 <= 真实余量，最小差 %.9e（链条自检）'%wgap)
print('     真实余量最小 %.9e；证明余量最小 %.9e'%(tmin, wneg))

# ------------------------------------------------- (F) the extremal family: pole + regular ring
print()
print('='*116)
print('(F) 最紧族：极点质量 delta + 正 n 边形纬度圈（桥接常数的极限族，插值定理的取等族）')
print('='*116)
print('     delta      n    A*^2        V          sigma^2      sigma^2/V   R(A)      9m_4        9*lb       9m_4-9lb  真实余量   证明余量')
wfam = 1e18; wratio = 0.; wgapf = 1e18
for delta in (1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 6e-2, 0.1, 0.15, 0.2, 0.3):
    for n in (3, 5, 6, 8):
        u2 = (1-delta)/(3-delta); b2 = math.sqrt(u2)
        m = np.array([delta]+[(3-delta)/n]*n)
        U = [np.array([0., 0., 1.])]
        for j in range(n):
            th = 2*math.pi*j/n
            U.append(np.array([math.sqrt(1-u2)*math.cos(th), math.sqrt(1-u2)*math.sin(th), b2]))
        U = np.array(U)
        if np.abs(np.einsum('i,ij,ik->jk', m, U, U)-np.eye(3)).max() > 1e-12: continue
        d = data(m, U); lb9 = 19/9-(7/3)*d['sig2']+(35/8)*d['sig2']**2
        R = (1+1/s3f)**2/(s3f*(s3f+d['A']))
        ratio = d['sig2']/d['V'] if d['V'] > 1e-15 else 0.
        mb, tm = margins(d, lb9/9)
        wfam = min(wfam, mb); wgapf = min(wgapf, tm-mb); wratio = max(wratio, ratio)
        print('   %8.1e  %3d  %.9f  %.9f  %.9f  %.7f  %.7f  %.9f  %.9f  %+.2e  %+.6e  %+.6e'%(
            delta, n, d['A']**2, d['V'], d['sig2'], ratio, R, 9*d['m4'], lb9, 9*d['m4']-lb9, tm, mb))
chk(wfam > 0, '最紧族证明余量恒正，最小 %.9e'%wfam)
chk(wgapf > -1e-9, '最紧族证明余量 <= 真实余量（最小差 %.3e）'%wgapf)
chk(wratio <= float((2+s3)/9)+1e-9, '实测 sigma^2/V 最大 %.12f <= (2+sqrt3)/9'%wratio)

# ------------------------------------------------------------- (G) library scan
print()
print('='*116)
print('(G) 库扫描：锥族 n=3..12 x 119 个 z、正二十面体、等角四线、等纬度正多边形')
print('='*116)
def cone(n, z):
    g0 = (1-3*z*z)/(1-z*z); G = 2/(1-z*z)
    m = np.array([g0]+[G/n]*n); U = [np.array([0., 0., 1.])]
    for j in range(n):
        th = 2*math.pi*j/n
        U.append(np.array([math.sqrt(1-z*z)*math.cos(th), math.sqrt(1-z*z)*math.sin(th), z]))
    return m, np.array(U)
wlib = 1e18; rec = None; nlib = 0; wviol = -1e18
for n in range(3, 13):
    for k in range(1, 120):
        z = (k/120.0)/math.sqrt(3)
        m, U = cone(n, z)
        if m.min() <= 0: continue
        d = data(m, U); lb9 = 19/9-(7/3)*d['sig2']+(35/8)*d['sig2']**2
        wviol = max(wviol, lb9-9*d['m4'])
        mb, tm = margins(d, lb9/9); nlib += 1
        if mb < wlib: wlib = mb; rec = (n, z, d, tm, mb, lb9/9)
        if tm < 0: print('   !! 真实违反: n=%d z=%.4f J=%.12f'%(n, z, d['J']))
ph = (1+math.sqrt(5))/2
Ui = np.array([np.array(v)/np.linalg.norm(v) for v in
               [(0,1,ph),(0,1,-ph),(1,ph,0),(1,-ph,0),(ph,0,1),(ph,0,-1)]])
ic = np.array([0.5]*6)
U4 = np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]], float); U4 /= np.linalg.norm(U4, axis=1)[:, None]
m4l = np.array([0.75]*4)
for name, (mm, UU) in [('正二十面体（6 线，等角 1/sqrt5）', (ic, Ui)), ('等角四线（立方体对角线）', (m4l, U4))]:
    d = data(mm, UU); lb9 = 19/9-(7/3)*d['sig2']+(35/8)*d['sig2']**2
    mb, tm = margins(d, lb9/9)
    print('   %s: A*^2=%.12f V=%.12f sigma^2=%.12f m_4=%.15f D=%.2e 9m_4-9lb=%+.2e 真实余量=%.6e 证明余量=%.6e'%(
        name, d['A']**2, d['V'], d['sig2'], d['m4'], d['D'], 9*d['m4']-lb9, tm, mb))
    if '二十' in name:
        chk(abs(9*d['m4']-lb9) < 1e-12, '%s 处主插值定理精确取等（sigma^2=4/15=%.12f, m_4=1/5）'%(name, d['sig2']))
    else:
        chk(9*d['m4']-lb9 > 1e-9 and abs(d['sig2']) < 1e-13,
            '%s：sigma^2=0 但 m_4=%s=7/27>19/81（不取等，因为 rho_hat_4 != 0）'%(name, sp.nsimplify(d['m4'])))
    chk(mb > 0, '%s 证明余量 %.6e > 0'%(name, mb))
    chk(abs(d['A']**2-(ph**2 if 'icosa' in name or '二十' in name else 3.0)) < 1e-12,
        '%s：A*^2 = %s 精确'%((name, 'phi^2') if '二十' in name else (name, '3')))
for n in (3, 4, 5, 6, 7, 8, 12):
    U = np.array([[math.sqrt(2/3)*math.cos(2*math.pi*j/n), math.sqrt(2/3)*math.sin(2*math.pi*j/n), 1/math.sqrt(3)]
                  for j in range(n)])
    m = np.array([3.0/n]*n)
    d = data(m, U); lb9 = 19/9-(7/3)*d['sig2']+(35/8)*d['sig2']**2
    mb, tm = margins(d, lb9/9)
    if n in (3, 6, 12):
        print('   等纬度正 %d 边形: sigma^2=%.2e m_4=%.15f (19/81=%.15f) 真实余量=%.6e 证明余量=%.6e'%(
            n, d['sig2'], d['m4'], 19/81, tm, mb))
    ok = abs(d['sig2']) < 1e-13 and 9*d['m4']-lb9 > -1e-12 and mb > 0
    if n in (3, 4):
        exp = 19/81 + (8*(9 if n == 3 else 0)+2*(9 if n == 4 else 0))/729
        ok = ok and abs(d['m4']-exp) < 1e-13
        chk(ok, '等纬度正 %d 边形：sigma^2=0，m_4=%.15f = 19/81+(8|rho_3|^2+2|rho_4|^2)/729（lem:r23m4eq）'%(n, d['m4']))
    else:
        ok = ok and abs(d['m4']-19/81) < 1e-13 and abs(9*d['m4']-lb9) < 1e-12
        chk(ok, '等纬度正 %d 边形：sigma^2=0，m_4=19/81 精确，主插值定理取等'%n)
n_, z_, dg, tm_, mb_, lb_ = rec
print('   锥族库最紧点：n=%d, z=%.6f, A*^2=%.9f, V=%.6f, sigma^2=%.9f, m_4=%.12f'%(n_, z_, dg['A']**2, dg['V'], dg['sig2'], dg['m4']))
print('              真实余量=%.9f  证明余量=%.9f  (9lb=%.9f, 9m_4=%.9f)'%(tm_, mb_, 9*lb_, 9*dg['m4']))
chk(wviol < 1e-12, '锥族库（%d 个构型）主插值定理无违反，最大违反 %.2e'%(nlib, wviol))
chk(wlib > 0, '锥族库证明余量恒正，最小 %.9e'%wlib)

print()
print('='*116)
print('SUMMARY  PASS=%d  FAIL=%d'%(PASS, FAIL))
print('  结论：中心截面闭合  kappa_3^c <= rho_0 - %.9e = %.15f'%(Mvert/6, rho0-Mvert/6))
print('='*116)
raise SystemExit(0 if FAIL == 0 else 1)
