# -*- coding: utf-8 -*-
# _r2.py  ROUND 18 part D:  sup_{Lam>=0,tr=1} Phi(Lam) = sum_i g_i D_i (u_i'Lam u_i) + ||sum_i g_i h_i u_i||^2 - 1
#                          h_i = sqrt(u_i'Lam u_i).   Conjecture (STAR): sup over geometries = tau0.
import numpy as np, math
from scipy.optimize import minimize
np.set_printoptions(precision=12, suppress=False, linewidth=190)
exec(open('_o1.py').read().split('print("="*104)')[0])
s5=math.sqrt(5.0); rho0=(9+s5+2*math.sqrt(35))/18; tau0=rho0-1
def core(U,g):
    Kf=U@U.T; P=np.maximum(-Kf,0.0); V=U*np.sqrt(g)[:,None]; Q=V@V.T; Qp=np.maximum(Q,0.0)
    w,ev=np.linalg.eigh(Qp); idx=np.argsort(-w); w=w[idx]; ev=ev[:,idx]
    rho=w[0]; x=ev[:,0].copy()
    if x.sum()<0: x=-x
    h=x/np.sqrt(g); c=g*h; m=c@U
    return dict(rho=rho,x=x,h=h,c=c,m=m,alpha=float(m@m),tau=rho-1.0,D=P@g,beta=float(c@(P@c)),
                E=float((np.maximum(-Q,0.0)**2).sum()),Kf=Kf,P=P,g=g,U=U,Q=Q,
                varpi=float(np.sum(np.where(Kf<0,np.outer(g*h**2,g*h**2),0.0))))
def whiten(U,g):
    M=(U*g[:,None]).T@U; ev,EV=np.linalg.eigh(M); ev=np.maximum(ev,1e-12)
    S=EV@np.diag(ev**-0.5)@EV.T; W=U@S; nrm=np.linalg.norm(W,axis=1)
    return W/nrm[:,None], g*nrm**2
def makePhi(U,g,D):
    Gu=(U* (g*D)[:,None]).T@U                 # Xi
    def Phi(A):
        Lam=A.reshape(3,3); Lam=Lam@Lam.T; t=np.trace(Lam)
        if t<=1e-14: return -1e9
        Lam=Lam/t
        q=np.einsum('ij,jk,ik->i',U,Lam,U); q=np.maximum(q,0.0)
        h=np.sqrt(q); m=(U*(g*h)[:,None]).sum(0)
        return float(np.trace(Lam@Gu))+float(m@m)-1.0
    return Phi,Gu
def maxPhi(U,g,D,ntry=14,seed=0):
    Phi,_=makePhi(U,g,D); rng=np.random.default_rng(seed); best=-1e9; bA=None
    starts=[np.eye(3)/math.sqrt(3)]
    for k in range(3):
        E=np.zeros((3,3)); E[k,k]=1.0; starts.append(E)
    for _ in range(ntry): starts.append(rng.normal(size=(3,3)))
    for A0 in starts:
        for meth in ('L-BFGS-B','Nelder-Mead'):
            try:
                r=minimize(lambda A:-Phi(A),A0.ravel(),method=meth,
                           options=dict(maxiter=6000,maxfev=20000) if meth=='Nelder-Mead' else dict(maxiter=6000))
                if -r.fun>best: best=-r.fun; bA=r.x.copy()
            except Exception: pass
    Lam=(bA.reshape(3,3)); Lam=Lam@Lam.T; Lam/=np.trace(Lam)
    return best,Lam
PASS=FAIL=0
def chk(cond,msg):
    global PASS,FAIL
    if cond: PASS+=1; print('  [PASS]',msg)
    else:    FAIL+=1; print('  [FAIL]',msg)
print('='*104); print('(D) 冠军几何上 max_{Lam>=0,tr=1} Phi(Lam)  ==?  tau0 （检验 (STAR) 的尖锐性与刚性）'); print('='*104)
Cc=core(Uc,gc); Xi=(Cc['U']*(Cc['g']*Cc['D'])[:,None]).T@Cc['U']
print('  冠军真实 tau = %.15f   tau0 = %.15f'%(Cc['tau'],tau0))
b,Lam=maxPhi(Cc['U'],Cc['g'],Cc['D'],ntry=22,seed=1)
print('  max Phi = %.15f   (tau0 = %.15f, 超出 %+.3e)'%(b,tau0,b-tau0))
print('  argmax eig(Lam) =',np.array2string(np.linalg.eigvalsh(Lam)[::-1],precision=10))
Lt=np.load('_r1_champ.npz')['Lt']
print('  冠军 KKT Lam~ eig =',np.array2string(np.linalg.eigvalsh(Lt)[::-1],precision=10))
print('  |Lam_argmax - Lam~|_max = %.3e'%np.abs(Lam-Lt).max())
chk(b<=tau0+1e-9,'(STAR) 在冠军几何上成立且【取等】 max Phi = %.13f <= tau0 = %.13f'%(b,tau0))
chk(np.abs(Lam-Lt).max()<3e-4,'唯一极大元就是 KKT 矩阵 Lam~（刚性） max|.|=%.3e'%np.abs(Lam-Lt).max())
print(); print('='*104); print('(E) 随机几何上的 (STAR) 检验：max_Lam Phi(Lam) <= tau0 ?'); print('='*104)
rng=np.random.default_rng(777); rows=[]; viol=[]
for t in range(60):
    N=int(rng.integers(5,13)); Z=rng.normal(size=(N,3)); Uu=Z/np.linalg.norm(Z,axis=1)[:,None]
    gg=rng.dirichlet(np.ones(N)*1.0)*3.0
    Uu,gg=whiten(Uu,gg)
    if gg.min()<1e-7: continue
    C=core(Uu,gg); b,Lam=maxPhi(Uu,gg,C['D'],ntry=8,seed=t)
    rows.append((N,C['rho'],b,b-tau0,C['tau']))
    if b>tau0+1e-8: viol.append((t,N,C['rho'],b,b-tau0))
    print('   N=%2d  rho=%.9f (rho-rho0=%+.2e)  tau=%.9f  maxPhi=%.9f  maxPhi-tau0=%+.3e'%(N,C['rho'],C['rho']-rho0,C['tau'],b,b-tau0))
print('  违背 (STAR) 的几何数 = %d / %d'%(len(viol),len(rows)))
chk(len(viol)==0,'(STAR) 在 %d 个随机几何上全部成立'%len(rows))
np.save('_r2_rows.npy',np.array(rows))
print(); print('第二段汇总：PASS=%d  FAIL=%d'%(PASS,FAIL))
