# -*- coding: utf-8 -*-
# _r1.py  ROUND 18: obtuse-dilution majorant + (STAR) reformulation tr(Lam Xi)-(1-alpha) <= tau0
import numpy as np, math
np.set_printoptions(precision=12, suppress=False, linewidth=190)
exec(open('_o1.py').read().split('print("="*104)')[0])
def fit_Lt(U,h2):
    rows=[[u[0]*u[0],u[1]*u[1],u[2]*u[2],2*u[0]*u[1],2*u[0]*u[2],2*u[1]*u[2]] for u in U]
    A=np.array(rows); sol,_,rank,_=np.linalg.lstsq(A,np.array(h2),rcond=None)
    Lt=np.array([[sol[0],sol[3],sol[4]],[sol[3],sol[1],sol[5]],[sol[4],sol[5],sol[2]]])
    return Lt,np.abs(A@sol-np.array(h2)).max(),rank
PASS=FAIL=0
def chk(cond,msg):
    global PASS,FAIL
    if cond: PASS+=1; print('  [PASS]',msg)
    else:    FAIL+=1; print('  [FAIL]',msg)
def core(U,g):
    N=len(g); Kf=U@U.T; P=np.maximum(-Kf,0.0)
    V=U*np.sqrt(g)[:,None]; Q=V@V.T; Qp=np.maximum(Q,0.0); Qm=np.maximum(-Q,0.0)
    w,ev=np.linalg.eigh(Qp); idx=np.argsort(-w); w=w[idx]; ev=ev[:,idx]
    rho=w[0]; x=ev[:,0].copy()
    if x.sum()<0: x=-x
    h=x/np.sqrt(g); c=g*h; m=c@U; alpha=float(m@m); d=P@c; beta=float(c@d)
    D=P@g; E=float((Qm**2).sum()); delta=(g*d)@U
    I=np.eye(N); wx=(I-Q)@x; Qmx=Qm@x
    varpi=float(np.sum(np.where(Kf<0, np.outer(c**2,c**2),0.0)))
    return dict(N=N,rho=rho,x=x,h=h,c=c,m=m,alpha=alpha,tau=rho-1.0,d=d,D=D,beta=beta,delta=delta,
                E=E,Q=Q,Qp=Qp,Qm=Qm,Kf=Kf,P=P,g=g,U=U,w=w,ev=ev,varpi=varpi,wx=wx,Qmx=Qmx)
def whiten(U,g,rng=None):
    M=(U*g[:,None]).T@U; ev,EV=np.linalg.eigh(M); ev=np.maximum(ev,1e-12)
    S=EV@np.diag(ev**-0.5)@EV.T; W=U@S; nrm=np.linalg.norm(W,axis=1)
    return W/nrm[:,None], g*nrm**2
print('='*104); print('(A) 恒等式族：delta = tau*m,  <m,delta>=alpha*tau,  (I-Q)Q_-x = rho(I-Q)x,  beta=sum g_i h_i d_i'); print('='*104)
C=core(Uc,gc); I=np.eye(C['N'])
chk(abs(C['delta']-C['tau']*C['m']).max()<1e-12,'delta = tau*m  (冠军)  max|.|=%.3e'%abs(C['delta']-C['tau']*C['m']).max())
chk(abs(float(C['m']@C['delta'])-C['alpha']*C['tau'])<1e-12,'<m,delta> = alpha*tau   %.15f vs %.15f'%(float(C['m']@C['delta']),C['alpha']*C['tau']))
chk(abs((I-C['Q'])@C['Qmx']-C['rho']*C['wx']).max()<1e-12,'(I-Q)Q_-x = rho*(I-Q)x   max|.|=%.3e'%abs((I-C['Q'])@C['Qmx']-C['rho']*C['wx']).max())
chk(abs(float(C['Qmx']@C['Qmx'])-((1-C['alpha'])*(1+2*C['tau'])+C['tau']**2))<1e-12,'||Q_-x||^2 = (1-alpha)(1+2tau)+tau^2   %.15f vs %.15f'%(float(C['Qmx']@C['Qmx']),(1-C['alpha'])*(1+2*C['tau'])+C['tau']**2))
chk(abs(C['Q']@C['Qmx']-C['tau']*(C['Q']@C['x'])).max()<1e-12,'Q Q_-x = tau*Qx   max|.|=%.3e'%abs((C['Q']@C['Qmx']-C['tau']*(C['Q']@C['x']))).max())
chk(abs(C['beta']-float((C['g']*C['h']*C['d']).sum()))<1e-13,'beta = sum_i g_i h_i d_i   %.15f'%C['beta'])
chk(abs(C['beta']-(C['tau']+1-C['alpha']))<1e-12,'tau-形式 tau = beta-(1-alpha)   %.15f'%C['tau'])
chk(abs(C['beta']**2-C['varpi']*C['E'])<1e-12,'thm:betaAA 在冠军处【取等】 beta^2 = varpi*e   %.15f vs %.15f'%(C['beta']**2,C['varpi']*C['E']))
print('  冠军: rho=%.12f alpha=%.12f beta=%.12f tau=%.12f e=%.12f varpi=%.12f'%(C['rho'],C['alpha'],C['beta'],C['tau'],C['E'],C['varpi']))
print('  D_i =',np.array2string(C['D'],precision=12)); print('  d_i =',np.array2string(C['d'],precision=12))
print(); print('='*104); print('(B) 钝性稀释引理  beta <= Dh := sum_i g_i h_i^2 D_i，等号 <=> h 在 G_- 每连通支上为常值'); print('='*104)
Dh=float((C['g']*C['h']**2*C['D']).sum())
chk(C['beta']<=Dh+1e-12,'冠军处取等：beta=%.15f  Dh=%.15f  差=%.3e'%(C['beta'],Dh,C['beta']-Dh))
adj=C['Kf']<-1e-12; comp=[-1]*C['N']; nc=0
for s in range(C['N']):
    if comp[s]>=0: continue
    st=[s]; comp[s]=nc
    while st:
        v=st.pop()
        for u2 in range(C['N']):
            if adj[v,u2] and comp[u2]<0: comp[u2]=nc; st.append(u2)
    nc+=1
groups=[sorted([k for k in range(C['N']) if comp[k]==cc]) for cc in range(nc)]
constcomp=all(C['h'][gp].max()-C['h'][gp].min()<1e-10 for gp in groups)
print('  G_- 连通支 = %d : %s'%(nc,groups))
chk(constcomp,'等号条件成立：h 在每个钝连通支上为常值 (枢纽孤立 + C5 支)')
rng=np.random.default_rng(20261029); worst=-1e9; ntested=0
for t in range(3000):
    N=int(rng.integers(4,14)); Z=rng.normal(size=(N,3)); Uu=Z/np.linalg.norm(Z,axis=1)[:,None]
    gg=rng.dirichlet(np.ones(N))*3.0
    Uu2,gg2=whiten(Uu,gg)
    if gg2.min()<=1e-9: continue
    Cc=core(Uu2,gg2); ntested+=1
    diff=float((Cc['g']*Cc['h']**2*Cc['D']).sum())-Cc['beta']
    if -diff>worst: worst=-diff
chk(worst<1e-11,'钝性稀释引理在 %d 个随机帧构型(N=4..13)上成立，最大违背 %.3e'%(ntested,worst))
print(); print('='*104); print('(C) (STAR) 精确重构：Phi(Lam)=tr(Lam Xi)+alpha-1，冠军处 Phi(Lam~)=tau0【精确取等】'); print('='*104)
Xi=(C['U']*(C['g']*C['D'])[:,None]).T@C['U']
print('  eig(Xi) =',np.array2string(np.linalg.eigvalsh(Xi)[::-1],precision=12),'  tr(Xi)=%.12f'%np.trace(Xi))
Lt,resid,rank=fit_Lt(C['U'],C['h']**2)
print('  冠军 KKT 矩阵 Lam~ eig=',np.array2string(np.linalg.eigvalsh(Lt)[::-1],precision=12),' tr=%.15f resid=%.2e'%(np.trace(Lt),resid))
chk(abs(np.trace(Lt)-1)<1e-10,'tr(Lam~)=1')
chk(np.abs((C['U']@Lt*C['U']).sum(1)-C['h']**2).max()<1e-10,'h_i^2 = u_i^T Lam~ u_i  max|.|=%.3e'%np.abs((C['U']@Lt*C['U']).sum(1)-C['h']**2).max())
chk(np.linalg.eigvalsh(Lt).min()>-1e-12,'Lam~ >= 0   lam_min=%.3e'%np.linalg.eigvalsh(Lt).min())
PhiC=float(np.trace(Lt@Xi)+C['alpha']-1)
chk(abs(PhiC-C['tau'])<1e-10,'【核心】Phi(Lam~) = tr(Lam~ Xi)+alpha-1 = tau = tau0   %.15f vs %.15f'%(PhiC,C['tau']))
chk(abs(float(np.trace(Lt@Xi))-C['beta'])<1e-10,'tr(Lam~ Xi) = Dh = beta  %.15f'%float(np.trace(Lt@Xi)))
print('  => (STAR) 是【尖锐】路线：冠军处精确取等，且 Phi 只依赖 (g,u,Lam)，不含特征向量方程')
np.savez('_r1_champ.npz',Lt=Lt,Xi=Xi,U=Uc,g=gc)
print(); print('第一段汇总：PASS=%d  FAIL=%d'%(PASS,FAIL))
