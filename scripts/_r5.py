# -*- coding: utf-8 -*-
# _r5.py ROUND 18 part H(fixed)/J: correct u-KKT (nu_i=0), and refined (STAR') with the Dirichlet term kept.
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
    h=x/np.sqrt(g); c=g*h; m=c@U; d=P@c
    return dict(N=N,rho=rho,x=x,h=h,c=c,m=m,alpha=float(m@m),tau=rho-1.0,D=P@g,d=d,
                beta=float(c@d),Kf=Kf,P=P,g=g,U=U,Q=Q,w=w,ev=ev)
def dil(C): return float((C['g']*C['h']**2*C['D']).sum())
def Em(h,U,g,Kf):
    return float(np.sum(np.where(Kf<0,np.outer(g,g)*np.abs(Kf)*(h[:,None]-h[None,:])**2,0.0)))
def stiefel_max(N,seed,iters=3000,step=0.05):
    r=np.random.default_rng(seed); Qr,_=np.linalg.qr(r.normal(size=(N,3))); V=Qr[:,:3]
    for it in range(iters):
        Kf=V@V.T; Qp=np.maximum(Kf,0.0); w,ev=np.linalg.eigh(Qp); i0=int(np.argmax(w)); x=ev[:,i0]
        if x.sum()<0: x=-x
        K=(Kf>0).astype(float); G=2.0*(np.outer(x,x)*K)@V
        G=G-V@((V.T@G+G.T@V)/2.0); V=V+step*G
        Qr,_=np.linalg.qr(V); V=Qr[:,:3]
    return V
print('='*104); print('(H-fix) u-KKT 的正确形式：nu_i=0，Lam~ u_i = (h_i/rho) W_i，Lam~ = (1/rho) sum_{acute} g_i g_j h_i h_j u_i u_j^T'); print('='*104)
def lam_ukkt(C):
    U=C['U']; g=C['g']; h=C['h']; rho=C['rho']; c=g*h; N=C['N']; ac=C['Kf']>1e-13
    W=np.array([(c[ac[i]][:,None]*U[ac[i]]).sum(0) if ac[i].any() else np.zeros(3) for i in range(N)])
    Lam=(1.0/rho)*np.einsum('i,ij,ik->jk',g*h,W,U)     # (1/rho) sum_i g_i h_i W_i u_i^T
    return Lam,W
def ukkt_res(C):
    U=C['U']; h=C['h']; rho=C['rho']; Lam,W=lam_ukkt(C)
    r1=np.abs((Lam@U.T).T-(h/rho)[:,None]*W).max()          # Lam~ u_i = (h_i/rho) W_i
    r2=np.abs((U@Lam*U).sum(1)-h**2).max()                  # h_i^2 = u_i' Lam~ u_i
    r3=abs(np.trace(Lam)-1.0)
    r4=np.abs(h[:,None]*W-rho*(Lam@U.T).T).max()            # R15-3
    r5=np.linalg.eigvalsh(Lam).min()
    r6=np.abs(Lam-Lam.T).max()
    return Lam,r1,r2,r3,r4,r5,r6
Cc=core(Uc,gc); Lt=np.load('_r1_champ.npz')['Lt']
Lam,r1,r2,r3,r4,r5,r6=ukkt_res(Cc)
print('   冠军: |Lam u_i-(h_i/rho)W_i|=%.2e |u\'Lam u-h^2|=%.2e |tr Lam-1|=%.2e |R15-3|=%.2e 反对称=%.2e lam_min=%.9f'%(r1,r2,r3,r4,r6,r5))
print('   eig(Lam~)=',np.array2string(np.linalg.eigvalsh(Lam)[::-1],precision=10),'  |Lam~-lstsq fit|=%.2e'%np.abs(Lam-Lt).max())
chk(max(r1,r2,r3,r4,r6)<1e-11,'冠军: 显式公式 Lam~=(1/rho)sum_{acute} g_ig_jh_ih_j u_iu_j^T 满足全部四条 KKT 关系 (max=%.2e)'%max(r1,r2,r3,r4,r6))
chk(r5>-1e-10,'冠军: Lam~ >= 0, lam_min=%.9f'%r5)
chk(abs(np.trace(Lam)-1)<1e-12,'tr(Lam~)=1 是【恒等式】(=sum g_i h_i^2)，不依赖临界性')
print(); print('='*104); print('(J) 精化 (STAR\')：Psi(Lam) = tr(Lam Xi)+alpha(Lam)-1-tau0-(1/2)E_-(h(Lam)) <= 0 ?'); print('='*104)
def makePsi(U,g,D,Kf):
    Xi=(U*(g*D)[:,None]).T@U
    def Psi(A):
        Lm=A.reshape(3,3); Lm=Lm@Lm.T; t=np.trace(Lm)
        if t<=1e-14: return -1e9
        Lm=Lm/t
        q=np.maximum(np.einsum('ij,jk,ik->i',U,Lm,U),0.0); h=np.sqrt(q)
        m=(U*(g*h)[:,None]).sum(0)
        return float(np.trace(Lm@Xi))+float(m@m)-1.0-tau0-0.5*Em(h,U,g,Kf)
    return Psi
def maxPsi(U,g,D,Kf,ntry=14,seed=0):
    Psi=makePsi(U,g,D,Kf); rng=np.random.default_rng(seed); best=-1e9; bA=None
    starts=[np.eye(3)/math.sqrt(3)]+[np.eye(3)[k].reshape(3,1)@np.eye(3)[k].reshape(1,3) for k in range(3)]
    for _ in range(ntry): starts.append(rng.normal(size=(3,3)))
    for A0 in starts:
        for meth in ('L-BFGS-B','Nelder-Mead'):
            try:
                r=minimize(lambda A:-Psi(A),A0.ravel(),method=meth,
                           options=(dict(maxiter=8000,maxfev=30000) if meth=='Nelder-Mead' else dict(maxiter=8000)))
                if -r.fun>best: best=-r.fun; bA=r.x.copy()
            except Exception: pass
    Lm=bA.reshape(3,3); Lm=Lm@Lm.T; Lm/=np.trace(Lm)
    return best,Lm
print('  --- 冠军几何 ---')
b,Lm=maxPsi(Cc['U'],Cc['g'],Cc['D'],Cc['Kf'],ntry=24,seed=3)
Lam_c,_,_,_,_,_,_=ukkt_res(Cc)
print('   max Psi = %+.12f  (需 <=0)   argmax eig(Lam)=%s'%(b,np.array2string(np.linalg.eigvalsh(Lm)[::-1],precision=8)))
print('   在 KKT 点 Lam~ 处 Psi = %+.3e'%(makePsi(Cc['U'],Cc['g'],Cc['D'],Cc['Kf'])(np.linalg.cholesky(Lam_c+1e-18*np.eye(3)).ravel())))
chk(b<=1e-9,'(STAR\') 在冠军几何上成立：max_Lam Psi = %+.3e <= 0'%b)
print('  --- Stiefel 极大元几何 (N=7..12) 与随机构型 ---')
rng=np.random.default_rng(5150); nbad=0; nb=0; worst=-1e9
tests=[]
for N in [6,7,8,9,10,12]:
    V=stiefel_max(N,777+N); g=(V**2).sum(1); U=V/np.sqrt(g)[:,None]; tests.append(('max N=%d'%N,U,g))
for t in range(18):
    N=int(rng.integers(5,13)); Z=rng.normal(size=(N,3)); Uu=Z/np.linalg.norm(Z,axis=1)[:,None]
    gg=rng.dirichlet(np.ones(N))*3.0; M=(Uu*gg[:,None]).T@Uu
    ev,EV=np.linalg.eigh(M); S=EV@np.diag(ev**-0.5)@EV.T; W2=Uu@S; nr=np.linalg.norm(W2,axis=1)
    tests.append(('rand N=%d t=%d'%(N,t),W2/nr[:,None],gg*nr**2))
for nm,U,g in tests:
    C=core(U,g); b,Lm=maxPsi(U,g,C['D'],C['Kf'],ntry=8,seed=hash(nm)%9999)
    nb+=1; worst=max(worst,b)
    if b>1e-8: nbad+=1
    Ck,_,_,_,_,_,_=ukkt_res(C)
    psi_kkt=makePsi(U,g,C['D'],C['Kf'])(np.linalg.cholesky(Ck+1e-15*np.eye(3)).ravel()) if np.linalg.eigvalsh(Ck).min()>1e-13 else float('nan')
    print('   %-14s rho-rho0=%+.2e  max Psi=%+.6e  Psi(Lam~_KKT)=%+.3e  E_-(h_true)/2=%.3e'%(nm,C['rho']-rho0,b,psi_kkt,0.5*Em(C['h'],U,g,C['Kf'])))
chk(nbad==0,'(STAR\') 在 %d 个几何(极大元+随机)上全部成立，max Psi = %+.3e'%(nb,worst))
print(); print('第五段汇总：PASS=%d  FAIL=%d'%(PASS,FAIL))
