# -*- coding: utf-8 -*-
# _r6.py ROUND 18 part K: algebraic rho(Lambda); grad = rho*I; homogeneity; rho(Lam)<=rho; CONCAVITY; broad scan
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
    N=len(g); Kf=U@U.T; V=U*np.sqrt(g)[:,None]; Q=V@V.T; Qp=np.maximum(Q,0.0)
    w,ev=np.linalg.eigh(Qp); idx=np.argsort(-w); w=w[idx]; ev=ev[:,idx]
    rho=w[0]; x=ev[:,0].copy()
    if x.sum()<0: x=-x
    h=x/np.sqrt(g); c=g*h; m=c@U; P=np.maximum(-Kf,0.0)
    return dict(N=N,rho=rho,x=x,h=h,c=c,m=m,alpha=float(m@m),tau=rho-1.0,D=P@g,d=P@c,
                beta=float(c@(P@c)),Kf=Kf,P=P,g=g,U=U,Q=Q)
def rhoLam(U,g,Lam):
    q=np.maximum(np.einsum('ij,jk,ik->i',U,Lam,U),0.0); h=np.sqrt(q)
    Kf=U@U.T; P=np.maximum(-Kf,0.0); m=(U*(g*h)[:,None]).sum(0)
    return float(m@m)+float((np.outer(g*h,g*h)*P).sum())
def norm1(Lam):
    t=np.trace(Lam); return Lam/t if t>0 else Lam
print('='*104); print('(K1) 代数泛函 rho(Lambda)=||sum g_i h_i u_i||^2 + sum g_i g_j h_i h_j (<u_i,u_j>)_- ,  h_i=sqrt(u_i\'Lam u_i)'); print('='*104)
Cc=core(Uc,gc); Lt=np.load('_r1_champ.npz')['Lt']
print('  冠军 rho = %.15f   rho(Lam~) = %.15f   差=%.3e'%(Cc['rho'],rhoLam(Uc,gc,Lt),Cc['rho']-rhoLam(Uc,gc,Lt)))
chk(abs(Cc['rho']-rhoLam(Uc,gc,Lt))<1e-12,'rho(Lam~) = rho（KKT 矩阵把代数泛函送回真值）')
# gradient identity  d rho / d Lam = rho * I_3
eps=1e-6; G=np.zeros((3,3))
for a in range(3):
    for b in range(3):
        E=np.zeros((3,3)); E[a,b]=1.0; E[b,a]=1.0 if a!=b else 1.0
        sc=0.5 if a!=b else 1.0
        G[a,b]=(rhoLam(Uc,gc,Lt+eps*E)-rhoLam(Uc,gc,Lt-eps*E))/(2*eps)*sc
print('  FD 梯度 d rho/d Lam =',np.array2string(G,precision=9))
print('  rho*I_3            =',np.array2string(Cc['rho']*np.eye(3),precision=9))
chk(np.abs(G-Cc['rho']*np.eye(3)).max()<1e-6,'【新恒等式】nabla_Lam rho(Lam) = rho(Lam)*I_3  max|FD-rho I|=%.3e'%np.abs(G-Cc['rho']*np.eye(3)).max())
# homogeneity
for t in [0.3,0.7,1.0,1.6,4.0]:
    pass
hh=[(t,rhoLam(Uc,gc,t*Lt)-t*rhoLam(Uc,gc,Lt)) for t in [0.25,0.5,2.0,7.0]]
chk(max(abs(v) for _,v in hh)<1e-11,'rho(Lam) 关于 Lam 是【一次正齐次】的 rho(tLam)=t rho(Lam)  max|.|=%.3e'%max(abs(v) for _,v in hh))
# rho(Lam) <= rho for tr Lam = 1
rng=np.random.default_rng(4711); worst=-1e9
for t in range(3000):
    A=rng.normal(size=(3,3)); Lm=A@A.T; Lm=norm1(Lm)
    v=rhoLam(Uc,gc,Lm)-Cc['rho']
    if v>worst: worst=v
    if t<2000: continue
chk(worst<1e-9,'rho(Lam) <= rho 对一切 Lam>=0,tr=1（冠军几何，3000 次）max(rho(Lam)-rho)=%.3e'%worst)
print(); print('='*104); print('(K2) rho(Lam) 的【凹性】检验（决定该路线能否凸化/给对偶证书）'); print('='*104)
def concavity_scan(U,g,ncase=400,seed=0):
    r=np.random.default_rng(seed); worst=1e9
    for t in range(ncase):
        A1=r.normal(size=(3,3)); A2=r.normal(size=(3,3))
        L1=norm1(A1@A1.T); L2=norm1(A2@A2.T); s=float(r.random())
        Lm=(1-s)*L1+s*L2
        marg=rhoLam(U,g,Lm)-((1-s)*rhoLam(U,g,L1)+s*rhoLam(U,g,L2))
        if marg<worst: worst=marg
    return worst
wc=concavity_scan(Uc,gc,600,1)
print('  冠军几何: min [rho((1-s)L1+sL2) - (1-s)rho(L1)-s rho(L2)] = %+.6e'%wc)
chk(wc>-1e-10,'rho(Lam) 在冠军几何上【凹】（600 组随机对，最小凹性余量 %+.3e）'%wc)
allw=1e9
for t in range(12):
    N=int(rng.integers(5,12)); Z=rng.normal(size=(N,3)); Uu=Z/np.linalg.norm(Z,axis=1)[:,None]
    gg=rng.dirichlet(np.ones(N))*3.0; M=(Uu*gg[:,None]).T@Uu
    ev,EV=np.linalg.eigh(M); S=EV@np.diag(ev**-0.5)@EV.T; W2=Uu@S; nr=np.linalg.norm(W2,axis=1)
    Uu=W2/nr[:,None]; gg=gg*nr**2
    w2=concavity_scan(Uu,gg,150,100+t); allw=min(allw,w2)
    print('   rand N=%2d  凹性余量 = %+.4e'%(N,w2))
chk(allw>-1e-10,'rho(Lam) 在 12 个随机几何上也凹（全局最小余量 %+.3e）=> 每个几何上 max_Lam rho(Lam) 是【凸优化】'%allw)
print(); print('='*104); print('(K3) 广谱扫描：max_{Lam>=0,tr=1} rho(Lam) <= rho0 ?  (含冠军微扰=最不利方向)'); print('='*104)
def maxrhoLam(U,g,ntry=16,seed=0):
    Kf=U@U.T; f=lambda A:-rhoLam(U,g,norm1((lambda M:M@M.T)(A.reshape(3,3))))
    r=np.random.default_rng(seed); best=-1e9; bA=None
    st=[np.eye(3)/math.sqrt(3),Lt]+[r.normal(size=(3,3)) for _ in range(ntry)]
    for A0 in st:
        for meth in ('L-BFGS-B','Nelder-Mead'):
            try:
                rr=minimize(f,A0.ravel(),method=meth,options=(dict(maxiter=9000,maxfev=30000) if meth=='Nelder-Mead' else dict(maxiter=9000)))
                if -rr.fun>best: best=-rr.fun; bA=rr.x.copy()
            except Exception: pass
    M=bA.reshape(3,3); return best,norm1(M@M.T)
b,Lm=maxrhoLam(Uc,gc,24,5)
print('  冠军几何: max rho(Lam) = %.15f  (rho0=%.15f, 差 %+.3e)  argmax eig=%s'%(b,rho0,b-rho0,np.array2string(np.linalg.eigvalsh(Lm)[::-1],precision=8)))
print('  |argmax - Lam~| = %.3e'%np.abs(Lm-Lt).max())
chk(b<=rho0+1e-10,'冠军几何上 max_Lam rho(Lam) = rho0【精确取到】')
chk(np.abs(Lm-Lt).max()<2e-5,'唯一极大点是 KKT 矩阵 Lam~（刚性）')
worst=-1e9; wn=None; nbad=0; nt=0
for t in range(40):
    if t%2==0:
        eps=float(rng.choice([1e-3,1e-2,4e-2])); Up=Uc+eps*rng.normal(size=Uc.shape); Up=Up/np.linalg.norm(Up,axis=1)[:,None]
        gp=np.maximum(gc*(1+eps*rng.normal(size=gc.shape)),1e-6); nm='champ-pert eps=%.0e'%eps
    else:
        N=int(rng.integers(5,13)); Z=rng.normal(size=(N,3)); Up=Z/np.linalg.norm(Z,axis=1)[:,None]
        gp=rng.dirichlet(np.ones(N))*3.0; nm='rand N=%d'%N
    M=(Up*gp[:,None]).T@Up; ev,EV=np.linalg.eigh(M); S=EV@np.diag(np.maximum(ev,1e-12)**-0.5)@EV.T
    W2=Up@S; nr=np.linalg.norm(W2,axis=1); Up=W2/nr[:,None]; gp=gp*nr**2
    if gp.min()<1e-8: continue
    bb,_=maxrhoLam(Up,gp,10,t); nt+=1
    if bb>worst: worst=bb; wn=nm
    if bb>rho0+1e-9: nbad+=1
    print('   %-22s max rho(Lam)=%.12f   max-rho0=%+.3e'%(nm,bb,bb-rho0))
chk(nbad==0,'广谱扫描 %d 个几何（含冠军微扰）全部满足 max_Lam rho(Lam) <= rho0；全局最大 = %.12f (%s)'%(nt,worst,wn))
print(); print('第六段汇总：PASS=%d  FAIL=%d'%(PASS,FAIL))
