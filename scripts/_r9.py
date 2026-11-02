# -*- coding: utf-8 -*-
# _r9.py  R18+: analytic concavity of rho(Lam) (geometric-mean representation), Lam~>=0 at critical points,
#         Lam~ is the GLOBAL argmax, tangent bound rho(Lam) <= rho*tr(Lam).
import numpy as np, math
from scipy.optimize import minimize
np.set_printoptions(precision=12, suppress=False, linewidth=200)
exec(open('_o1.py').read().split('print("="*104)')[0])
s5=math.sqrt(5.0); rho0=(9+s5+2*math.sqrt(35.0))/18
PASS=FAIL=0
def chk(c,m):
    global PASS,FAIL
    if c: PASS+=1; print('  [PASS]',m)
    else: FAIL+=1; print('  [FAIL]',m)
def core(U,g):
    N=len(g); Kf=U@U.T; V=U*np.sqrt(g)[:,None]; Q=V@V.T; Qp=np.maximum(Q,0.0)
    w,ev=np.linalg.eigh(Qp); idx=np.argsort(-w); w=w[idx]; ev=ev[:,idx]
    rho=w[0]; x=ev[:,0].copy()
    if x.sum()<0: x=-x
    h=x/np.sqrt(g); c=g*h; m=c@U
    return dict(N=N,rho=rho,x=x,h=h,c=c,m=m,alpha=float(m@m),tau=rho-1.0,Kf=Kf,g=g,U=U,Q=Q,Qp=Qp)
def whiten(U,g):
    M=(U*g[:,None]).T@U; ev,EV=np.linalg.eigh(M)
    if ev.min()<=1e-10: return None,None
    S=EV@np.diag(ev**-0.5)@EV.T; W=U@S; nr=np.linalg.norm(W,axis=1)
    return W/nr[:,None], g*nr**2
def lam_ukkt(C):
    U=C['U']; g=C['g']; h=C['h']; rho=C['rho']; c=g*h; N=C['N']; ac=C['Kf']>1e-13
    W=np.array([(c[ac[i]][:,None]*U[ac[i]]).sum(0) if ac[i].any() else np.zeros(3) for i in range(N)])
    return (1.0/rho)*np.einsum('i,ij,ik->jk',g*h,W,U)
def rhoLam(U,g,L):
    h=np.sqrt(np.maximum(np.einsum('ij,jk,ik->i',U,L,U),0.0))
    m=(U*(g*h)[:,None]).sum(0); return float(m@m)+float((np.outer(g*h,g*h)*np.maximum(-(U@U.T),0.0)).sum())
def rhoLam_geo(U,g,L):   # geometric-mean representation  sum_ij (Q+)_ij sqrt(L_i L_j),  L_i = g_i u_i'Lam u_i
    Qp=np.maximum(np.sqrt(np.outer(g,g))*(U@U.T),0.0)
    Lin=g*np.einsum('ij,jk,ik->i',U,L,U)
    xx=np.sqrt(np.maximum(Lin,0.0))
    return float((Qp*np.outer(xx,xx)).sum())
def stiefel_max(N,seed,iters=3000,step=0.05):
    r=np.random.default_rng(seed); Qr,_=np.linalg.qr(r.normal(size=(N,3))); V=Qr[:,:3]
    for it in range(iters):
        Kf=V@V.T; Qp=np.maximum(Kf,0.0); w,ev=np.linalg.eigh(Qp); x=ev[:,int(np.argmax(w))]
        if x.sum()<0: x=-x
        G=2.0*(np.outer(x,x)*(Kf>0).astype(float))@V; G=G-V@((V.T@G+G.T@V)/2.0); V=V+step*G
        Qr,_=np.linalg.qr(V); V=Qr[:,:3]
    return V
print('='*104); print('(L1) rho(Lam) 的几何平均表示  rho(Lam)=sum_ij (Q+)_ij sqrt(L_i L_j),  L_i=g_i u_i\'Lam u_i  => 凹（解析）'); print('='*104)
Cc=core(Uc,gc); rng=np.random.default_rng(1717); w=0.0
for t in range(600):
    A=rng.normal(size=(3,3)); L=A@A.T; L/=np.trace(L)
    w=max(w,abs(rhoLam(Uc,gc,L)-rhoLam_geo(Uc,gc,L)))
print('  冠军几何上两种写法最大差 = %.3e'%w); chk(w<1e-12,'几何平均表示与原定义一致 (max %.2e)'%w)
def conc(U,g,n,seed):
    r=np.random.default_rng(seed); best=1e9
    for t in range(n):
        A1=r.normal(size=(3,3)); A2=r.normal(size=(3,3)); L1=A1@A1.T; L1/=np.trace(L1); L2=A2@A2.T; L2/=np.trace(L2)
        s=float(r.random()); best=min(best,rhoLam(U,g,(1-s)*L1+s*L2)-((1-s)*rhoLam(U,g,L1)+s*rhoLam(U,g,L2)))
    return best
wc=conc(Uc,gc,800,3); print('  冠军凹性余量 %+.4e'%wc); chk(wc>-1e-11,'冠军几何上 rho(Lam) 凹（解析已证，数值余量 %+.2e）'%wc)
print(); print('='*104); print('(L2) 临界点上 Lam~ >= 0（可容许性）与 Lam~ 是 {tr=1} 上的全局极大点'); print('='*104)
def maxrho(U,g,ntry=18,seed=0):
    f=lambda A:-rhoLam(U,g,(lambda M:M@M.T/np.trace(M@M.T))(A.reshape(3,3)))
    r=np.random.default_rng(seed); best=-1e9; bA=None
    for A0 in [np.eye(3)/math.sqrt(3)]+[r.normal(size=(3,3)) for _ in range(ntry)]:
        for meth in ('L-BFGS-B','Nelder-Mead'):
            try:
                rr=minimize(f,A0.ravel(),method=meth,options=(dict(maxiter=9000,maxfev=30000) if meth=='Nelder-Mead' else dict(maxiter=9000)))
                if -rr.fun>best: best=-rr.fun; bA=rr.x.copy()
            except Exception: pass
    M=bA.reshape(3,3); return best,M@M.T/np.trace(M@M.T)
rows=[]
for N in [6,7,8,9,10,11,12]:
    V=stiefel_max(N,555+N*13); g=(V**2).sum(1); U=V/np.sqrt(g)[:,None]; C=core(U,g)
    Lt=lam_ukkt(C); lmin=np.linalg.eigvalsh(Lt).min()
    res=max(np.abs((Lt@U.T).T-(C['h']/C['rho'])[:,None]*np.array([( (C['g']*C['h'])[(C['Kf'][i]>1e-13)][:,None]*U[(C['Kf'][i]>1e-13)]).sum(0) for i in range(N)])).max(),
            np.abs((U@Lt*U).sum(1)-C['h']**2).max(), abs(np.trace(Lt)-1))
    b,Lb=maxrho(U,g,16,N)
    tang=0.0
    r2=np.random.default_rng(9000+N)
    for t in range(400):
        A=r2.normal(size=(3,3)); L=A@A.T; tang=max(tang,rhoLam(U,g,L)-C['rho']*np.trace(L))
    print('   N=%2d rho-rho0=%+.2e lam_min(Lam~)=%+.9f KKT残差=%.1e max_Lam rho(Lam)=%.12f (差 %+.2e) |Lam_arg-Lam~|=%.1e 切向界违背=%.1e'%(
        N,C['rho']-rho0,lmin,res,b,b-C['rho'],np.abs(Lb-Lt).max(),tang))
    rows.append((N,lmin,res,b-C['rho'],np.abs(Lb-Lt).max(),tang))
    chk(lmin>-1e-9 and res<1e-8,'N=%d: Lam~>=0 (lam_min=%+.3e) 且满足 KKT'%(N,lmin))
    chk(abs(b-C['rho'])<1e-7,'N=%d: max_Lam rho(Lam) = rho（Lam~ 是全局极大点），差 %.2e'%(N,b-C['rho']))
    chk(tang<1e-9,'N=%d: 切向界 rho(Lam) <= rho*tr(Lam) 在 400 个随机 Lam 上成立（最大违背 %.2e）'%(N,tang))
print(); print('第八段汇总：PASS=%d FAIL=%d'%(PASS,FAIL))
