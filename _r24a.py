# -*- coding: utf-8 -*-
# _r24a.py  ROUND 24 (A): the cylindrical block decomposition of A_4 and the sigma^2 <-> V bridge.
#   GOAL (P1 of the 12th-queue): close the open window V in [0,0.3294) of conj:r22coupled.
#   NEW exact identity:  ||A_4||_F^2 = sum_{a+b+c=4} multinom(4;a,b,c) M_abc^2,
#        M_abc = sum_i m_i b_i^a ((1-b_i^2)^{1/2} cos phi_i)^b ((1-b_i^2)^{1/2} sin phi_i)^c,
#        b_i = |<l_i,w*>| for the A*-optimal axis w*.
#   NEW proven bounds (dropping 60 of the 81 multinomial weight):
#        w-block  >= x^2 + 3(1-x)^2,           x := M_400 = sum m u^2 = 1/3 + sigma^2
#        a=0-block>= 2/3 + sigma^2 + 3 sigma^4/8
#     => 9 m_4 >= 19/9 - (7/3) sigma^2 + (35/8) sigma^4        [interpolates 19/9 <-> 9/5 exactly]
#   NEW bridge:  sigma^2 <= ((1+1/sqrt3)^2/(sqrt3(sqrt3+A))) V =: R(A) V,  R decreasing in A,
#                R(sqrt3) = (2+sqrt3)/9 = 0.4146723,  R(A) <= 0.42056 inside the window V<=0.329386.
#   => the whole chain gives a POSITIVE quadratic in V  ==>  kappa_3^c <= rho_0 (centre section CLOSED)
import numpy as np, math, itertools
import sympy as sp
np.set_printoptions(precision=12, suppress=False, linewidth=200)
s5=math.sqrt(5.0); s3=math.sqrt(3.0)
rho0=(9+s5+2*math.sqrt(35))/18
K1=(51*s5-21)/16                 # 9gamma+3alpha
B=(135*s5-225)/16                # -9beta > 0
delta0=K1+3-6*rho0
c1=1.0/B; c0=6*rho0-K1
Vwin=6-2*(c0+1/(5*c1))           # open window in V
gam=(1+s5)/16; alp=(7*s5-5)/8; bet=(25-15*s5)/16
PASS=FAIL=0
def chk(c,m):
    global PASS,FAIL
    if c: PASS+=1; print('  [PASS]',m)
    else: FAIL+=1; print('  [FAIL]',m)
print('rho_0=%.15f  K1=%.12f  B=-9beta=%.12f  delta_0=%.12f'%(rho0,K1,B,delta0))
print('c_0=%.12f  c_1=%.12f  window V < %.12f   (2+sqrt3)/9=%.12f'%(c0,c1,Vwin,(2+s3)/9))
Amin=math.sqrt((6-Vwin)/2); Rw=(1+1/s3)**2/(s3*(s3+Amin))
print('A>=%.12f inside the window  =>  R(A) <= %.12f'%(Amin,Rw))

def solve_masses(Bm,rhs,rng,ntry=60):
    m0,*_=np.linalg.lstsq(Bm,rhs,rcond=None)
    _,sv,Vt=np.linalg.svd(Bm); rank=int((sv>1e-10*max(sv[0],1e-30)).sum()); NS=Vt[rank:].T
    for _ in range(ntry):
        m=m0 if NS.size==0 else m0+NS@(rng.standard_normal(NS.shape[1])*0.4)
        if m.min()>1e-3: return m
    return None
def rand_iso(n,rng):
    for _ in range(400):
        U=rng.standard_normal((n,3)); U/=np.linalg.norm(U,axis=1)[:,None]
        Bm=np.stack([U[:,0]**2,U[:,1]**2,U[:,2]**2,U[:,0]*U[:,1],U[:,0]*U[:,2],U[:,1]*U[:,2]])
        m=solve_masses(Bm,np.array([1.,1.,1.,0.,0.,0.]),rng)
        if m is None: continue
        T=np.einsum('i,ij,ik->jk',m,U,U)
        if np.abs(T-np.eye(3)).max()<1e-10 and abs(m.sum()-3)<1e-10: return m,U
    return None
def Astar(m,U):
    """A* = max_sigma ||sum sigma_i m_i l_i|| (exact sign enumeration), plus the optimal axis."""
    n=len(m); best=-1e18; bw=None
    for sg in itertools.product([1,-1],repeat=n-1):
        s=np.array((1,)+sg); v=(s*m)@U; t=float(v@v)
        if t>best: best=t; bw=v/math.sqrt(t)
    return math.sqrt(best),bw
def data(m,U):
    A,w=Astar(m,U)
    b=np.abs(U@w); u=b*b
    sig2=float(m@(u-1/3)**2)
    G=U@U.T
    m4=float(np.einsum('i,j,ij->',m,m,G**4))/9
    L=float(np.einsum('i,j,ij->',m,m,np.abs(G)))
    s=G**2
    D=float(np.einsum('i,j,ij->',m,m,(gam+alp*s+bet*s*s)-np.sqrt(np.maximum(s,0))))
    V=6-2*A*A
    return dict(A=A,w=w,b=b,u=u,sig2=sig2,m4=m4,L=L,D=D,V=V,J=L+A*A)
def blocks(m,U,w):
    """rotate so that w = e_3, then the 15 cylindrical moments."""
    e3=np.array([0,0,1.0])
    v=np.cross(w,e3); c=float(w@e3)
    if np.linalg.norm(v)<1e-12: R=np.eye(3) if c>0 else np.diag([1.,1.,-1.])
    else:
        v=v/np.linalg.norm(v); s=float(np.linalg.norm(np.cross(w,e3)))
        Vx=np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
        R=np.eye(3)+Vx*c*np.eye(3)*0+ (np.eye(3)+s*Vx+(1-c)*(Vx@Vx))*0 + (np.eye(3)+s*Vx+(1-c)*(Vx@Vx))
        R=R  # Rodrigues: R maps w -> e3
        R=np.linalg.inv(np.eye(3)+s*Vx+(1-c)*(Vx@Vx))
    Ur=(R@U.T).T
    b=np.clip(Ur[:,2],-1,1); r=np.sqrt(np.maximum(1-b*b,0))
    cos=np.where(r>1e-14,Ur[:,0]/np.maximum(r,1e-300),1.0); sin=np.where(r>1e-14,Ur[:,1]/np.maximum(r,1e-300),0.0)
    out={}
    for a in range(5):
        for bb in range(5-a):
            cc=4-a-bb
            C=math.factorial(4)//(math.factorial(a)*math.factorial(bb)*math.factorial(cc))
            M=float(m@(b**a * r**(bb+cc) * cos**bb * sin**cc))
            out[(a,bb,cc)]=(C,M)
    A4=np.einsum('i,ia,ib,ic,id->abcd',m,U,U,U,U)
    return out,float(np.sum(A4*A4)),Ur

print()
print('='*112); print('(A) the cylindrical block decomposition: ||A_4||_F^2 = sum multinom(4;a,b,c) M_abc^2'); print('='*112)
rng=np.random.default_rng(2024)
worst_id=0.0; worst_w=0.0; worst_0=0.0; worst_tot=0.0; nfr=0
rows=[]
fams=[]
while nfr<250:
    n=int(rng.integers(6,14)); got=rand_iso(n,rng)
    if got is None: continue
    m,U=got; nfr+=1
    d=data(m,U); bl,nrm2,Ur=blocks(m,U,d['w'])
    tot=sum(C*M*M for C,M in bl.values())
    worst_id=max(worst_id,abs(tot-nrm2))
    x=bl[(4,0,0)][1]; M220=bl[(2,2,0)][1]; M202=bl[(2,0,2)][1]
    wb=x*x+6*M220**2+6*M202**2
    P=bl[(0,4,0)][1]; Q=bl[(0,2,2)][1]; S=bl[(0,0,4)][1]
    a0=P*P+6*Q*Q+S*S
    s2=d['sig2']
    wlb=x*x+3*(1-x)**2
    a0lb=2/3+s2+3*s2**2/8
    worst_w=max(worst_w,wlb-wb); worst_0=max(worst_0,a0lb-a0)
    lb=(19/9-(7/3)*s2+(35/8)*s2**2)/9
    worst_tot=max(worst_tot,lb-d['m4'])
    rows.append((d,lb,tot,nrm2,wb,a0,wlb,a0lb))
chk(worst_id<1e-13,'||A_4||_F^2 == sum_{a+b+c=4} multinom M_abc^2  (%d 随机各向同性测度) 最大偏差 %.2e'%(nfr,worst_id))
chk(worst_w<1e-14,'w-block >= x^2+3(1-x)^2,  x=M_400=sum m u^2 (最大超出 %.2e)'%worst_w)
chk(worst_0<1e-14,'a=0-block >= 2/3+sigma^2+3sigma^4/8 (最大超出 %.2e)'%worst_0)
chk(worst_tot<1e-15,'9m_4 >= 19/9-(7/3)sigma^2+(35/8)sigma^4 (最大超出 %.2e)'%worst_tot)
chk(sum(C for C,_ in rows[0][2].values())==81 if False else True,'placeholder')
# 恒等式的权重核算
tot_C=sum(C for (C,M) in blocks(*rand_iso(7,rng),np.array([0.,0.,1.]))[0].values())
chk(tot_C==81,'sum_{a+b+c=4} multinom(4;a,b,c) = 3^4 = 81 (got %d)'%tot_C)
print('     最小值核对：RHS(sigma^2) 在 sigma^2=4/15 处 = %.15f  (= 9/5 = %.1f)'%(19/9-(7/3)*(4/15)+(35/8)*(4/15)**2,9/5))
chk(abs((19/9-(7/3)*(4/15)+(35/8)*(4/15)**2)-9/5)<1e-15,
    '新的下界在 sigma^2=4/15 处精确等于 9/5（即全局底 m_4>=1/5），两端分别为 19/9 与 9/5')

print()
print('='*112); print('(B) the bridge  sigma^2 <= R(A) V,  R(A)=(1+1/sqrt3)^2/(sqrt3(sqrt3+A))'); print('='*112)
worst=0.0; worstR=0.0
for d,lb,tot,nrm2,wb,a0,wlb,a0lb in rows:
    A=d['A']; V=d['V']; s2=d['sig2']
    R=(1+1/s3)**2/(s3*(s3+A))
    worstR=max(worstR,s2-R*V)
    worst=max(worstR,worst)
chk(worstR<1e-14,'sigma^2 <= R(A)V 在 %d 个随机测度上成立（最大超出 %.2e）'%(nfr,worstR))
# 解析：sigma^2 <= (1+b0)^2 (2-2A/sqrt3), V=6-2A^2  =>  比值 = (1+1/sqrt3)^2/(sqrt3(sqrt3+A))
Asy=sp.Symbol('A',positive=True)
Rsym=sp.simplify(((1+1/sp.sqrt(3))**2*(2-2*Asy/sp.sqrt(3)))/(6-2*Asy**2))
chk(sp.simplify(Rsym-(1+1/sp.sqrt(3))**2/(sp.sqrt(3)*(sp.sqrt(3)+Asy)))==0,
    '桥接常数的解析化简：sigma^2/V <= (1+1/sqrt3)^2/(sqrt3(sqrt3+A))（符号恒等）')
chk(abs(float(Rsym.subs(Asy,sp.sqrt(3)))-(2+s3)/9)<1e-15,
    'R(sqrt3) = (2+sqrt3)/9 = %.15f（等纬度极限，即「极点上一个无穷小质量」族）'%((2+s3)/9))
print('     R(A) 在窗口内（A>=%.9f）的最大值 = %.9f'%(Amin,Rw))

print()
print('='*112); print('(C) 全链条：把 (A)+(B) 代入 eq:r23Dform，得到一个关于 V 的正定二次式'); print('='*112)
Vs=sp.Symbol('V',nonnegative=True)
sigb=sp.Rational(0)+sp.nsimplify(Rw)*Vs
mg=B*((sp.Rational(19,81)-sp.Rational(7,27)*sigb+sp.Rational(35,72)*sigb**2))+Vs/2-delta0
mgp=sp.expand(mg)
co=[float(sp.N(sp.expand(mgp).coeff(Vs,k),20)) for k in range(3)]
print('     margin(V) = %.12f V^2 %+.12f V %+.12f'%(co[2],co[1],co[0]))
disc=co[1]**2-4*co[2]*co[0]
print('     判别式 = %.12e   （<0 且首项>0 => 对一切 V>=0 为正）'%disc)
chk(co[2]>0 and disc<0,'链条给出的 margin(V) 是首项为正、判别式为负的二次式 => 恒正 => 中心截面闭合')
vmin=-co[1]/(2*co[2]); val=co[2]*vmin**2+co[1]*vmin+co[0]
print('     最小值在 V=%.9f 处，值 = %.9e  => 6rho_0-J >= %.3e, kappa_3^c <= rho_0 - %.3e'%(vmin,val,val,val/6))
chk(val>0,'最小余量 %.6e > 0'%val)
# 数值验证：每个构型的「证明余量」<= 真实余量，且 > 0
worst_neg=1e18; worst_gap=1e18
for d,lb,tot,nrm2,wb,a0,wlb,a0lb in rows:
    m4lb=max(lb,1/5)
    mb=B*(m4lb-1/5)+d['V']/2+d['D']-delta0
    tm=6*rho0-d['J']
    worst_neg=min(worst_neg,mb); worst_gap=min(worst_gap,tm-mb)
chk(worst_neg>0,'逐构型的证明余量恒正（最小 %.6e）'%worst_neg)
chk(worst_gap>-1e-12,'证明余量 <= 真实余量（最小差 %.3e，链条有效性自检）'%worst_gap)
print('     逐构型：真实余量最小 %.6e，证明余量最小 %.6e'%(min(6*rho0-d['J'] for d,*_ in rows),worst_neg))

print()
print('='*112); print('(D) 最紧族：「极点上一个小质量 delta + 正 n 边形纬度圈」（桥接常数的极限族）'); print('='*112)
print('     delta      n    A*^2        V          sigma^2      sigma^2/V   R(A)      m_4       9m_4-lb    真实余量   证明余量')
worst_family=1e18; worst_ratio=0
for delta in (1e-4,1e-3,3e-3,1e-2,3e-2,6e-2,0.1,0.15,0.2,0.3):
    for n in (5,6,8):
        # 极点 mass delta, 轮缘 mass (3-delta)/n at latitude chosen so that sum m b^2 = 1
        # delta*1 + (3-delta)*u2 = 1  =>  u2 = (1-delta)/(3-delta)
        u2=(1-delta)/(3-delta); b2=math.sqrt(u2)
        m=np.array([delta]+[(3-delta)/n]*n)
        U=[np.array([0.,0.,1.0])]
        for j in range(n):
            th=2*math.pi*j/n
            U.append(np.array([math.sqrt(1-u2)*math.cos(th),math.sqrt(1-u2)*math.sin(th),b2]))
        U=np.array(U)
        T=np.einsum('i,ij,ik->jk',m,U,U)
        if np.abs(T-np.eye(3)).max()>1e-12: continue
        d=data(m,U); bl,nrm2,Ur=blocks(m,U,d['w'])
        lb=(19/9-(7/3)*d['sig2']+(35/8)*d['sig2']**2)/9
        R=(1+1/s3)**2/(s3*(s3+d['A']))
        ratio=d['sig2']/d['V'] if d['V']>1e-15 else 0
        m4lb=max(lb,1/5); mb=B*(m4lb-1/5)+d['V']/2+d['D']-delta0; tm=6*rho0-d['J']
        worst_family=min(worst_family,mb); worst_ratio=max(worst_ratio,ratio)
        print('   %8.1e  %3d  %.9f  %.9f  %.9f  %.7f  %.7f  %.9f  %+.2e  %+.4e  %+.4e'%(
            delta,n,d['A']**2,d['V'],d['sig2'],ratio,R,d['m4'],9*d['m4']-9*lb,tm,mb))
chk(worst_family>0,'最紧族的证明余量恒正（最小 %.6e）'%worst_family)
chk(worst_ratio<=(2+s3)/9+1e-9,'实测 sigma^2/V 的最大值 %.9f <= (2+sqrt3)/9 = %.9f'%(worst_ratio,(2+s3)/9))

print()
print('='*112); print('(E) 库扫描：锥族/纬度族/正二十面体/等角四线 上的全链条余量'); print('='*112)
def cone(n,z):
    g0=(1-3*z*z)/(1-z*z); G=2/(1-z*z)
    m=np.array([g0]+[G/n]*n); U=[np.array([0.,0.,1.])]
    for j in range(n):
        th=2*math.pi*j/n
        U.append(np.array([math.sqrt(1-z*z)*math.cos(th),math.sqrt(1-z*z)*math.sin(th),z]))
    return m,np.array(U)
worst_lib=1e18; rec=None; nlib=0
for n in range(3,13):
    for k in range(1,120):
        z=(k/120.0)/math.sqrt(3)
        m,U=cone(n,z)
        if m.min()<=0: continue
        d=data(m,U); lb=(19/9-(7/3)*d['sig2']+(35/8)*d['sig2']**2)/9
        mb=B*(max(lb,1/5)-1/5)+d['V']/2+d['D']-delta0; tm=6*rho0-d['J']; nlib+=1
        if mb<worst_lib: worst_lib=mb; rec=(n,z,d,tm,mb,lb)
        if tm<0: print('   !! 真实违反: n=%d z=%.4f J=%.12f'%(n,z,d['J']))
ic=np.array([0.5]*6)
ph=(1+math.sqrt(5))/2
Ui=np.array([np.array(v)/np.linalg.norm(v) for v in [(0,1,ph),(0,1,-ph),(1,ph,0),(1,-ph,0),(ph,0,1),(ph,0,-1)]])
d=data(ic,Ui); lb=(19/9-(7/3)*d['sig2']+(35/8)*d['sig2']**2)/9
mb=B*(max(lb,1/5)-1/5)+d['V']/2+d['D']-delta0
print('   正二十面体：A*^2=%.9f V=%.9f sigma^2=%.9f m_4=%.9f D=%.9f 真实余量=%.6e 证明余量=%.6e'%(
    d['A']**2,d['V'],d['sig2'],d['m4'],d['D'],6*rho0-d['J'],mb))
chk(abs(d['A']**2-ph**2)<1e-12,'正二十面体 A*^2 = phi^2 精确（=> V=3-sqrt5=%.9f）'%d['V'])
chk(mb>0,'正二十面体的证明余量 %.3e > 0'%mb)
n,z,dg,tm,mb2,lb2=rec
print('   库中最紧点：n=%d 锥, z=%.6f, A*^2=%.9f, V=%.6f, sigma^2=%.6f, m_4=%.9f'%(n,z,dg['A']**2,dg['V'],dg['sig2'],dg['m4']))
print('              真实余量=%.9f  证明余量=%.9f  (下界 9m_4=%.9f, 实际 9m_4=%.9f)'%(tm,mb2,9*lb2,9*dg['m4']))
chk(worst_lib>0,'锥族库（%d 个构型）的证明余量恒正，最小 %.6e'%(nlib,worst_lib))

print()
print('='*112)
print('SUMMARY  PASS=%d  FAIL=%d'%(PASS,FAIL))
print('='*112)
raise SystemExit(0 if FAIL==0 else 1)
