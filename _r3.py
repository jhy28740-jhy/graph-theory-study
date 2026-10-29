# -*- coding: utf-8 -*-
# _r3.py ROUND 18 part F/G/H: strong dilution inequality alpha + D(h) <= rho0 ; rank-1 maximizer ; u-KKT symmetry system
import numpy as np, math
from scipy.optimize import minimize
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
    h=x/np.sqrt(g); c=g*h; m=c@U; alpha=float(m@m); d=P@c
    return dict(N=N,rho=rho,x=x,h=h,c=c,m=m,alpha=alpha,tau=rho-1.0,D=P@g,d=d,
                beta=float(c@d),E=float((np.maximum(-Q,0.0)**2).sum()),Kf=Kf,P=P,g=g,U=U,Q=Q,w=w,ev=ev)
def whiten(U,g):
    M=(U*g[:,None]).T@U; ev,EV=np.linalg.eigh(M); ev=np.maximum(ev,1e-12)
    S=EV@np.diag(ev**-0.5)@EV.T; W=U@S; nrm=np.linalg.norm(W,axis=1)
    return W/nrm[:,None], g*nrm**2
def dil(C):  # D(h) = sum g_i h_i^2 D_i
    return float((C['g']*C['h']**2*C['D']).sum())
print('='*104); print('(F) 强稀释不等式  alpha + D(h_Perron) <= rho0  （比主猜想强，冠军处取等）'); print('='*104)
Cc=core(Uc,gc); S=dil(Cc)+Cc['alpha']
print('  冠军: alpha=%.15f  D(h)=%.15f  和=%.15f  rho0=%.15f  余量=%.3e'%(Cc['alpha'],dil(Cc),S,rho0,S-rho0))
chk(abs(S-rho0)<1e-12,'冠军处 alpha+D(h) = rho0【精确取等】')
chk(abs(dil(Cc)-Cc['beta'])<1e-12,'且 D(h)=beta（稀释引理在冠军处取等）')
rng=np.random.default_rng(31337); worst=-1e9; wN=0; nviol=0; ntested=0; worst_rho=0
for t in range(1500):
    N=int(rng.integers(4,15)); Z=rng.normal(size=(N,3)); Uu=Z/np.linalg.norm(Z,axis=1)[:,None]
    gg=rng.dirichlet(np.ones(N)*float(rng.choice([0.4,1.0,3.0])))*3.0
    Uu,gg=whiten(Uu,gg)
    if gg.min()<1e-8: continue
    C=core(Uu,gg); ntested+=1
    val=dil(C)+C['alpha']
    if val-worst>-1e-15: worst=val; wN=N; worst_rho=C['rho']
    if val>rho0+1e-10: nviol+=1
chk(nviol==0,'强稀释不等式在 %d 个随机帧构型(N=4..14)上成立，max(alpha+D(h))=%.12f (N=%d) vs rho0=%.12f'%(ntested,worst,wN,rho0))
print('  该最坏构型的真实 rho = %.12f'%worst_rho)
# champion perturbations (most adversarial: near the maximizer)
print('  --- 冠军的随机微扰（最靠近极大元的方向）---')
wp=-1e9; wn=0
for t in range(1200):
    eps=float(rng.choice([1e-3,1e-2,3e-2,1e-1]))
    Up=Uc+eps*rng.normal(size=Uc.shape); Up=Up/np.linalg.norm(Up,axis=1)[:,None]
    gp=np.maximum(gc*(1+eps*rng.normal(size=gc.shape)),1e-6)
    Up,gp=whiten(Up,gp)
    C=core(Up,gp); val=dil(C)+C['alpha']
    if val>wp: wp=val; wn=C['N']
chk(wp<=rho0+1e-10,'冠军微扰 1200 次：max(alpha+D(h))=%.12f <= rho0=%.12f (余量 %.3e)'%(wp,rho0,wp-rho0))
print(); print('='*104); print('(G) Stiefel 梯度上升极大元上的强稀释不等式 (N=6..12)'); print('='*104)
def stiefel_max(N,seed,iters=4000):
    r=np.random.default_rng(seed); V0=r.normal(size=(3,N)); Qr,_=np.linalg.qr(V0.T); V=Qr[:,:3].T
    for it in range(iters):
        Kf=V.T@V; Qp=np.maximum(Kf,0.0); w,ev=np.linalg.eigh(Qp); i0=np.argmax(w); x=ev[:,i0]
        if x.sum()<0: x=-x
        K=(Kf>0).astype(float); G=2.0*(np.outer(x,x)*K)@V.T
        G=G-V@(V.T@G+G.T@V)/2.0
        V=V+0.05*G
        Qr,Rr=np.linalg.qr(V.T); V=Qr[:,:3].T
    return V.T
for N in [6,7,8,9,10,12]:
    best=None
    for sd in range(6):
        V=stiefel_max(N,1000+N*10+sd)
        g=(V**2).sum(1); U=V/np.sqrt(g)[:,None]
        C=core(U,g)
        if best is None or C['rho']>best['rho']: best=C
    val=dil(best)+best['alpha']
    print('   N=%2d rho-rho0=%+.3e min g=%.2e alpha=%.9f D(h)=%.9f alpha+D(h)=%.12f  余量=%+.3e  D(h)-beta=%.2e'%(
        N,best['rho']-rho0,best['g'].min(),best['alpha'],dil(best),val,val-rho0,dil(best)-best['beta']))
    chk(val<=rho0+1e-9,'N=%d 极大元满足强稀释不等式'%N)
print(); print('='*104); print('(H) u-KKT 显式公式 Lam~ u_i = (h_i/rho)(W_i + d_i u_i) 与对称性系统残差'); print('='*104)
def ukkt_resid(C):
    U=C['U']; g=C['g']; h=C['h']; rho=C['rho']; c=g*h; d=C['d']
    ac=U@U.T>1e-13
    W=np.array([(c[ac[i]][:,None]*U[ac[i]]).sum(0) if ac[i].any() else np.zeros(3) for i in range(C['N'])])
    Lam=(1.0/rho)*((U*(g*h)[:,None]).T@ (W+(d[:,None]*U)))
    Lam=0.5*(Lam+Lam.T)
    res1=np.abs(Lam@U.T-U.T@(W+d[:,None]*U)*0.0)  # placeholder
    # exact: compare Lam u_i  vs  (h_i/rho)(W_i+d_i u_i)
    A=(Lam@U.T).T
    B=(h/rho)[:,None]*(W+d[:,None]*U)
    r1=np.abs(A-B).max()
    r2=np.abs((U@Lam*U).sum(1)-h**2).max()
    r3=abs(np.trace(Lam)-1)
    nu=-h*d
    r4=np.abs(h[:,None]*W-2*(rho/2)*(Lam@U.T).T-nu[:,None]*U).max()
    return Lam,r1,r2,r3,r4
for nm,C in [('冠军',Cc)]:
    Lam,r1,r2,r3,r4=ukkt_resid(C)
    Lt=np.load('_r1_champ.npz')['Lt']
    print('   %s: |Lam u_i-(h_i/rho)(W_i+d_i u_i)|=%.3e  |u\'Lam u-h^2|=%.3e  |tr Lam-1|=%.3e  |u-KKT|=%.3e'%(nm,r1,r2,r3,r4))
    print('   eig(Lam_uKKT)=',np.array2string(np.linalg.eigvalsh(Lam)[::-1],precision=10),' |Lam-Lam~(fit)|=%.3e'%np.abs(Lam-Lt).max())
    chk(max(r1,r2,r3,r4)<1e-11,'%s: u-KKT 显式公式给出同一个 Lam~（四重残差 %.2e）'%(nm,max(r1,r2,r3,r4)))
# symmetry system  h_i R_ij - h_j R_ji = h_i d_j - h_j d_i + (h_j d_j - h_i d_i)<u_i,u_j>
def symres(C):
    U=C['U']; g=C['g']; h=C['h']; d=C['d']; Kf=C['Kf']; N=C['N']
    Wt=np.maximum(-Kf,0.0)*np.outer(g*h,np.ones(N))   # Wt[i,k]= g_k h_k (<u_i,u_k>)_-  for k~i
    R=np.einsum('ik,jk->ij',Wt,U)                     # R_ij = sum_{k~i} g_k h_k |<u_i,u_k>| <u_j,u_k>
    L=h[:,None]*R-h[None,:]*R.T
    Rt=h[:,None]*d[None,:]-h[None,:]*d[:,None]+(np.outer(d,h**0)*0+ (np.outer(np.ones(N),h*d)-np.outer(h*d,np.ones(N))))*Kf
    return np.abs(L-Rt).max()
print('   对称性系统残差(冠军) = %.3e'%symres(Cc))
chk(symres(Cc)<1e-11,'u-KKT 的对称性等价形式在冠军处成立（残差 %.2e）'%symres(Cc))
print(); print('第三段汇总：PASS=%d  FAIL=%d'%(PASS,FAIL))
