# -*- coding: utf-8 -*-
# _r4.py ROUND 18 part G/H/I: Stiefel maximizers, explicit u-KKT formula for Lam~, exact dilution identity D(h)=beta+E_-/2
import numpy as np, math
np.set_printoptions(precision=12, suppress=False, linewidth=200)
exec(open('_o1.py').read().split('print("="*104)')[0])
s5=math.sqrt(5.0); rho0=(9+s5+2*math.sqrt(35))/18; tau0=rho0-1
PASS=FAIL=0
def chk(cond,msg):
    global PASS,FAIL
    if cond: PASS+=1; print('  [PASS]',msg)
    else:    FAIL+=1; print('  [FAIL]',msg)
def core(U,g):
    N=len(g); Kf=U@U.T; P=np.maximum(-Kf,0.0); V=U*np.sqrt(g)[:,None]; Q=V@V.T
    Qp=np.maximum(Q,0.0); w,ev=np.linalg.eigh(Qp); idx=np.argsort(-w); w=w[idx]; ev=ev[:,idx]
    rho=w[0]; x=ev[:,0].copy()
    if x.sum()<0: x=-x
    h=x/np.sqrt(g); c=g*h; m=c@U; d=P@c
    return dict(N=N,rho=rho,x=x,h=h,c=c,m=m,alpha=float(m@m),tau=rho-1.0,D=P@g,d=d,
                beta=float(c@d),Kf=Kf,P=P,g=g,U=U,Q=Q,w=w,ev=ev)
def dil(C): return float((C['g']*C['h']**2*C['D']).sum())
def Eminus(C):
    W=np.where(C['Kf']<0, np.abs(C['Kf'])*np.outer(C['g']*C['h']**2,np.ones(C['N'])),0.0)
    H=C['h']; D2=(H[:,None]-H[None,:])**2
    return float(np.sum(np.where(C['Kf']<0, np.outer(C['g'],C['g'])*np.abs(C['Kf'])*D2,0.0)))
def stiefel_max(N,seed,iters=3000,step=0.05):
    r=np.random.default_rng(seed); Qr,_=np.linalg.qr(r.normal(size=(N,3))); V=Qr[:,:3]
    for it in range(iters):
        Kf=V@V.T; Qp=np.maximum(Kf,0.0); w,ev=np.linalg.eigh(Qp); i0=int(np.argmax(w)); x=ev[:,i0]
        if x.sum()<0: x=-x
        K=(Kf>0).astype(float); G=2.0*(np.outer(x,x)*K)@V
        G=G-V@((V.T@G+G.T@V)/2.0); V=V+step*G
        Qr,_=np.linalg.qr(V); V=Qr[:,:3]
    return V
print('='*104); print('(I) 精确稀释恒等式  D(h) = beta + (1/2) E_-(h),  E_- = sum_{i~j} g_i g_j |<u_i,u_j>| (h_i-h_j)^2'); print('='*104)
Cc=core(Uc,gc)
chk(abs(dil(Cc)-(Cc['beta']+0.5*Eminus(Cc)))<1e-14,'冠军处恒等式成立  D(h)=%.15f beta=%.15f E_-=%.3e'%(dil(Cc),Cc['beta'],Eminus(Cc)))
rng=np.random.default_rng(9090); worst=0.0
for t in range(2000):
    N=int(rng.integers(4,14)); Z=rng.normal(size=(N,3)); Uu=Z/np.linalg.norm(Z,axis=1)[:,None]
    gg=rng.dirichlet(np.ones(N))*3.0; M=(Uu*gg[:,None]).T@Uu
    ev,EV=np.linalg.eigh(M); S=EV@np.diag(ev**-0.5)@EV.T; W2=Uu@S; nr=np.linalg.norm(W2,axis=1)
    Uu=W2/nr[:,None]; gg=gg*nr**2
    C=core(Uu,gg); e=abs(dil(C)-(C['beta']+0.5*Eminus(C)))
    if e>worst: worst=e
chk(worst<1e-12,'稀释恒等式在 2000 个随机构型(N=4..13)上成立 max|残差|=%.3e'%worst)
print('  => tau = tr(Lam~ Xi) - (1-alpha) - (1/2) E_-(h)  是【精确】重构（非松弛），冠军处 E_-=0');
# E_- via Perron equation:  E_- = (1/rho^2) sum_{i~j} g_i g_j |<u_i,u_j>| zeta_ij^2, zeta=<u_i-u_j,m>+d_i-d_j
def Eminus2(C):
    U=C['U']; g=C['g']; m=C['m']; d=C['d']; Kf=C['Kf']; rho=C['rho']
    z=(U@m)[:,None]-(U@m)[None,:]+d[:,None]-d[None,:]
    return float(np.sum(np.where(Kf<0,np.outer(g,g)*np.abs(Kf)*z**2,0.0)))/rho**2
chk(abs(Eminus(Cc)-Eminus2(Cc))<1e-13,'E_- 的第二显式（只含 m,d,几何）成立：%.15f vs %.15f'%(Eminus(Cc),Eminus2(Cc)))
w2=0.0
for t in range(1200):
    N=int(rng.integers(4,14)); Z=rng.normal(size=(N,3)); Uu=Z/np.linalg.norm(Z,axis=1)[:,None]
    gg=rng.dirichlet(np.ones(N))*3.0; M=(Uu*gg[:,None]).T@Uu
    ev,EV=np.linalg.eigh(M); S=EV@np.diag(ev**-0.5)@EV.T; W2=Uu@S; nr=np.linalg.norm(W2,axis=1)
    Uu=W2/nr[:,None]; gg=gg*nr**2
    C=core(Uu,gg); e=abs(Eminus(C)-Eminus2(C))
    if e>w2: w2=e
chk(w2<1e-11,'E_- 第二显式在 1200 个随机构型上成立 max|残差|=%.3e'%w2)
print(); print('='*104); print('(H) u-KKT 显式公式  nu_i=-h_i d_i,  Lam~ u_i = (h_i/rho)(W_i + d_i u_i)'); print('='*104)
def lam_ukkt(C):
    U=C['U']; g=C['g']; h=C['h']; rho=C['rho']; c=g*h; d=C['d']; N=C['N']
    ac=C['Kf']>1e-13
    W=np.array([(c[ac[i]][:,None]*U[ac[i]]).sum(0) if ac[i].any() else np.zeros(3) for i in range(N)])
    Lam=(1.0/rho)*((U*(g*h)[:,None]).T@(W+d[:,None]*U))
    return Lam,W
def ukkt_res(C):
    U=C['U']; g=C['g']; h=C['h']; rho=C['rho']; d=C['d']
    Lam,W=lam_ukkt(C); Lam=0.5*(Lam+Lam.T)
    r1=np.abs((Lam@U.T).T-(h/rho)[:,None]*(W+d[:,None]*U)).max()
    r2=np.abs((U@Lam*U).sum(1)-h**2).max()
    r3=abs(np.trace(Lam)-1.0)
    r4=np.abs(h[:,None]*W-rho*(Lam@U.T).T+ (h*d)[:,None]*U).max()
    r5=np.linalg.eigvalsh(Lam).min()
    return Lam,r1,r2,r3,r4,r5
Lam,r1,r2,r3,r4,r5=ukkt_res(Cc); Lt=np.load('_r1_champ.npz')['Lt']
print('   冠军: |Lam u_i-(h_i/rho)(W_i+d_i u_i)|=%.2e  |u\'Lam u-h^2|=%.2e  |tr Lam-1|=%.2e  |u-KKT|=%.2e  lam_min(Lam)=%.6f'%(r1,r2,r3,r4,r5))
print('   eig(Lam)=',np.array2string(np.linalg.eigvalsh(Lam)[::-1],precision=10),'  |Lam-Lam~(lstsq fit)|=%.2e'%np.abs(Lam-Lt).max())
chk(max(r1,r2,r3,r4)<1e-11,'冠军: u-KKT 显式公式给出同一个 Lam~，且 tr=1、h^2=u\'Lam u 全部自洽 (max=%.2e)'%max(r1,r2,r3,r4))
chk(r5>-1e-10,'冠军: Lam~ >= 0 (lam_min=%.3e)'%r5)
print(); print('='*104); print('(G) Stiefel 梯度上升极大元 (N=6..12)：稀释恒等式 / u-KKT 自洽 / alpha+D(h)'); print('='*104)
for N in [6,7,8,9,10,12]:
    best=None
    for sd in range(4):
        V=stiefel_max(N,4242+N*37+sd)
        g=(V**2).sum(1); U=V/np.sqrt(g)[:,None]
        C=core(U,g)
        if best is None or C['rho']>best[0]['rho']: best=(C,U,g)
    C=best[0]; Lam,a1,a2,a3,a4,a5=ukkt_res(C)
    val=dil(C)+C['alpha']
    print('   N=%2d rho-rho0=%+.3e min g=%.2e alpha+D(h)-rho0=%+.3e  E_-/2=%.3e  uKKT残差=%.2e  tr(Lam)-1=%+.2e  lam_min(Lam)=%.3e'%(
        N,C['rho']-rho0,C['g'].min(),val-rho0,0.5*Eminus(C),max(a1,a2,a4),a3,a5))
    chk(abs(dil(C)-(C['beta']+0.5*Eminus(C)))<1e-10 and max(a1,a2,a4)<1e-8 and abs(a3)<1e-8,'N=%d: 恒等式+u-KKT 自洽'%N)
    if N>=7: chk(a5>-1e-7,'N=%d: Lam~ >= 0（分支 II 所需的半定性成立）'%N)
print(); print('第四段汇总：PASS=%d  FAIL=%d'%(PASS,FAIL))
