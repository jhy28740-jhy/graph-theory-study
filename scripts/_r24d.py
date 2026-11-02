# -*- coding: utf-8 -*-
# _r24d.py  ROUND 24 (D): (1) the extremal family in CLOSED FORM (so thm:r24interp is sharp at
#   every sigma_u^2, not just at the two ends); (2) a rigorous three-branch uniform margin.
import numpy as np, math, itertools
import sympy as sp

PASS = FAIL = 0
def chk(c, m):
    global PASS, FAIL
    if c: PASS += 1; print('  [PASS]', m)
    else: FAIL += 1; print('  [FAIL]', m)

s5=sp.sqrt(5); s35=sp.sqrt(35); s3=sp.sqrt(3)
rho0_e=(9+s5+2*s35)/18; phi_e=(1+s5)/2; B_e=(135*s5-225)/16
d0_e=sp.simplify(3*phi_e+3-6*rho0_e); Vwin_e=2*d0_e; Amin_e=sp.sqrt(3-d0_e)
Rw_e=(1+1/s3)**2/(s3*(s3+Amin_e))
rho0=float(rho0_e); B=float(B_e); d0=float(d0_e); Vwin=float(Vwin_e)
Rw=float(Rw_e); s3f=math.sqrt(3.0)
def gg(t): return 19/9-(7/3)*t+(35/8)*t*t

# ===================================================================== (F2) closed-form extremal family
print('='*116)
print('(F2) 取等族的闭式：极点质量 delta + 正 n 边形（n>=5）纬度圈，u_2=(1-delta)/(3-delta)')
print('='*116)
dl = sp.Symbol('delta', positive=True)
u2_e = (1-dl)/(3-dl); r2_e = sp.simplify(1-u2_e)
chk(sp.simplify(r2_e-2/(3-dl)) == 0, 'r_2^2 = 1-u_2 = 2/(3-delta) 精确')
sig_e = sp.simplify(dl*(1-sp.Rational(1,3))**2 + (3-dl)*(u2_e-sp.Rational(1,3))**2)
sig_closed = 4*dl/(3*(3-dl))
chk(sp.simplify(sig_e-sig_closed) == 0, 'sigma_u^2 = 4delta/(3(3-delta)) 精确（delta in [0,1] 时扫遍 [0,2/3]）')
chk(sp.simplify(sig_closed.subs(dl,1)-sp.Rational(2,3)) == 0, 'delta=1 时 sigma_u^2 = 2/3（达到上界）')
T_e = sp.simplify(sp.Rational(4,3)+sig_closed); Q_e = sp.simplify((3-dl)*r2_e**2*sp.Rational(1,8)/(3-dl)*1)
Qring = sp.simplify((3-dl)*r2_e**2/8)
chk(sp.simplify(Qring-1/(2*(3-dl))) == 0, 'Q = sum m c s = 1/(2(3-delta)) 精确')
chk(sp.simplify(Qring-T_e/8) == 0, 'Q = T/8 精确 => a0-块的松弛 8(Q-T/8)^2 恒为零')
print('     T = 4/3+sigma_u^2 = %s ; T/8 = %s'%(sp.simplify(T_e), sp.simplify(T_e/8)))
# the w-block slack: M220 = M202 for a regular n-gon (n>=3) since sum cos^2 = sum sin^2 = n/2
print('     正 n 边形 (n>=3): sum_j cos^2 = sum_j sin^2 = n/2  => M_220 = M_202，w-块松弛 3(M220-M202)^2 = 0')
print('     正 n 边形 (n>=5): sum_j cos^p sin^q = 0 对一切 1<=p+q<=4, (p,q)!=(2,2)  => 9 个被丢弃的矩全为零')
# Fourier verification for n>=5
worstF = 0.
for n in (5,6,7,8,9,12,13):
    th = np.array([2*math.pi*j/n for j in range(n)])
    for (p,q) in [(1,0),(0,1),(1,1),(3,0),(0,3),(2,1),(1,2),(3,1),(1,3)]:
        worstF = max(worstF, abs(np.sum(np.cos(th)**p*np.sin(th)**q))/n)
    worstF = max(worstF, abs(np.sum(np.cos(th)**2)/n-0.5)); worstF = max(worstF, abs(np.sum(np.sin(th)**2)/n-0.5))
    worstF = max(worstF, abs(np.sum(np.cos(th)**4)/n-3/8)); worstF = max(worstF, abs(np.sum(np.sin(th)**4)/n-3/8))
    worstF = max(worstF, abs(np.sum(np.cos(th)**2*np.sin(th)**2)/n-1/8))
chk(worstF < 1e-15, '正 n 边形(n=5..13) 的 Fourier 消没与均值恒等式，最大偏差 %.2e'%worstF)
print('     极点贡献：r_pole=0，故任何 (b,c)!=(0,0) 的模式在极点上为零')
# numeric end-to-end: 9m_4 == g(sigma_u^2) for the family
def Astar(m,U):
    n=len(m); best=-1e18; bw=None
    for sg in itertools.product([1,-1],repeat=n-1):
        s=np.array((1,)+sg); v=(s*m)@U; t=float(v@v)
        if t>best:
            best=t
            if t>1e-300: bw=v/math.sqrt(t)
    return math.sqrt(best),bw
# closed forms for the family (verified symbolically below)
Ap_e = sp.simplify(dl+sp.sqrt((3-dl)*(1-dl)))
print('     A_pole(delta) = delta+sqrt((3-delta)(1-delta))；A_pole(0)=sqrt3, A_pole(1/2)=phi')
chk(sp.simplify(Ap_e.subs(dl,0)-s3)==0 and sp.simplify(Ap_e.subs(dl,sp.Rational(1,2))-phi_e)==0,
    'A_pole(0)=sqrt(3)（等纬度）且 A_pole(1/2)=phi（delta=1/2 恰是正二十面体）')
print('     V(delta) = 6-2A_pole^2；delta->0 时 sigma_u^2/V -> (2+sqrt3)/9')
lim_ratio = sp.limit(sig_closed/(6-2*Ap_e**2), dl, 0, dir='+')
chk(sp.simplify(lim_ratio-(2+s3)/9)==0, 'lim_{delta->0} sigma_u^2/V = (2+sqrt3)/9 精确（sympy 极限）')
worst_eq = 0.; rows=[]; worst_tilt = 0.
for delta in (1e-5,1e-4,1e-3,1e-2,3e-2,6e-2,0.1,0.15,0.2,0.3,0.4,0.45,0.49,0.5):
    for n in (5,6,7,8,12):
        u2=(1-delta)/(3-delta); b2=math.sqrt(u2)
        m=np.array([delta]+[(3-delta)/n]*n)
        U=[np.array([0.,0.,1.])]
        for j in range(n):
            t=2*math.pi*j/n
            U.append(np.array([math.sqrt(1-u2)*math.cos(t),math.sqrt(1-u2)*math.sin(t),b2]))
        U=np.array(U)
        if np.abs(np.einsum('i,ij,ik->jk',m,U,U)-np.eye(3)).max()>1e-12: continue
        A,w=Astar(m,U); b=np.abs(U@w); u=b*b
        s2=float(m@(u-1/3)**2); G=U@U.T
        m4=float(np.einsum('i,j,ij->',m,m,G**4))/9
        V=6-2*A*A
        worst_eq=max(worst_eq,abs(9*m4-gg(s2)))
        Ap=float(Ap_e.subs(dl,delta)); worst_tilt=max(worst_tilt,abs(A-Ap))
        assert abs(s2-float(sig_closed.subs(dl,delta)))<1e-12, (delta,n,s2)
        rows.append((delta,n,A*A,V,s2,m4))
chk(worst_eq<1e-13,'delta in [0,1/2] 的闭式族上 9m_4 = g(sigma_u^2) 精确取等（%d 个构型，最大偏差 %.2e）'%(len(rows),worst_eq))
chk(worst_tilt<1e-12,'delta<=1/2 时极轴确为 A*-最优轴，A*=A_pole（最大偏差 %.2e）'%worst_tilt)
print('     注意：delta>1/2 时 A*-最优轴会倾斜离开极轴（数值：delta=0.7 处 w_z=0.529，A*=1.6362>A_pole=1.5307），')
print('           故取等族恰好在 delta=1/2 处终止；而该族扫过的 sigma_u^2 in [0,4/15] 正是链条需要的全部范围。')
print('     sigma_u^2 覆盖范围：%.9f .. %.9f（= [0,4/15]，4/15=%.9f）'%(min(r[4] for r in rows),max(r[4] for r in rows),4/15))
print('     sigma_u^2/V 的最大值：%.12f  vs (2+sqrt3)/9 = %.12f（桥接常数也被同一族取到）'%(
    max((r[4]/r[3] for r in rows if r[3]>1e-15)), float((2+s3)/9)))
chk(max((r[4]/r[3] for r in rows if r[3]>1e-15)) <= float((2+s3)/9)+1e-9,
    '桥接引理在同一族上取到极限 (2+sqrt3)/9 => 整条链条在闭式族上逐点取等')

# ===================================================================== (D2) rigorous 3-branch margin
print()
print('='*116)
print('(D2) 严格三分支的一致余量（不依赖网格）')
print('='*116)
V2 = 2*(d0 + 1.039957626350013513e-03)
print('     分支 1：V in [0,V_win]，V_win=%.15f'%Vwin)
c2=float(B_e*sp.Rational(35,72)*Rw_e**2); c1=0.5-float(B_e*sp.Rational(7,27)*Rw_e)
c0=float(B_e*(sp.Rational(19,81)-sp.Rational(1,5))-d0_e)
disc=c1*c1-4*c2*c0; Vvert=-c1/(2*c2); m1=-disc/(4*c2)
print('       quad(V)=%.18f V^2 %+.18f V %+.18f, 判别式=%.18e, 顶点 V=%.18f, 最小值=%.18e'%(c2,c1,c0,disc,Vvert,m1))
chk(c2>0 and disc<0 and 0<Vvert<Vwin, '分支 1：quad 首项>0、判别式<0、顶点在窗口内 => min = %.18e'%m1)
print('     分支 2：V in [V_win,V_2]，V_2=2(delta_0+min_1)=%.15f'%V2)
def h(V): return (1+1/s3f)**2*V/(s3f*(s3f+math.sqrt(3-V/2)))
m2 = B*(gg(min(h(V2),4/15))/9-1/5)+Vwin/2-d0
print('       h(V_2)=%.15f < 4/15=%.15f，g 递减 => Phi >= B(g(h(V_2))/9-1/5)+V_win/2-delta_0 = %.15f'%(
    h(V2),4/15,m2))
chk(h(V2)<4/15 and m2>m1, '分支 2 余量 %.15f > 分支 1 最小值 %.15e'%(m2,m1))
print('     分支 3：V >= V_2，m_4>=1/5 => Phi >= V/2-delta_0 >= V_2/2-delta_0 = %.15e'%(V2/2-d0))
chk(abs(V2/2-d0-m1)<1e-9 and V2/2-d0>0, '分支 3 余量 = %.15e（按 V_2 的定义恰为分支 1 的最小值）'%(V2/2-d0))
Munif=min(m1,m2,V2/2-d0)
print('  ==> 全局一致余量 = %.18e（J 尺度），即 kappa_3^c <= rho_0 - %.18e = %.18f'%(Munif,Munif/6,rho0-Munif/6))
chk(Munif>0,'中心截面闭合，一致余量 %.18e'%Munif)
# symbolic high-precision confirmation of the branch-1 minimum
m1_e = sp.simplify(-(( sp.Rational(1,2)-B_e*sp.Rational(7,27)*Rw_e)**2
                     -4*(B_e*sp.Rational(35,72)*Rw_e**2)*(B_e*(sp.Rational(19,81)-sp.Rational(1,5))-d0_e))
                   /(4*B_e*sp.Rational(35,72)*Rw_e**2))
print('     高精度（50 位）复核分支 1 最小值 = %s'%sp.N(m1_e,50))
chk(abs(float(sp.N(m1_e,50))-m1)<1e-15,'50 位高精度与双精度一致（差 %.2e）'%abs(float(sp.N(m1_e,50))-m1))
print('     kappa_3^c <= %.18f   （数值上确界 1.265991717870，仍低 %.6f）'%(rho0-Munif/6, rho0-Munif/6-1.265991717870))

print()
print('='*116); print('SUMMARY  PASS=%d  FAIL=%d'%(PASS,FAIL)); print('='*116)
raise SystemExit(0 if FAIL==0 else 1)
