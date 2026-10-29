# -*- coding: utf-8 -*-
# _o1.py  ROUND 15 (part 1): trace-1 SCF rigidity, acute-moment identity (**), sigma-identities,
#                            spectral floor n_+, nu_- control, coincident-atom merging
import numpy as np, math, sympy as sp
np.set_printoptions(precision=12, suppress=False, linewidth=190)
s5=math.sqrt(5.0); s7=math.sqrt(7.0); s35=math.sqrt(35.0)
rho0=(9+s5+2*s35)/18; tau=rho0-1; PHI=(3+s5)/4; MOHAR=(1+math.sqrt(3))/2
exec(open('_n1.py').read().split('print("="*100)')[0])

def champion():
    Vhh=14.0/(-7*s5-7+5*s7+3*s35); Vhc=math.sqrt(7.0/5.0); Vcc=-3.5-7*s5/10+3*s35/10+5*s7/2
    Ah=math.sqrt(Vhh); Ac=Vhc/Ah; Rc=math.sqrt(Vcc-Ac**2)
    uh=(91-25*s7)/126; uc=(7+5*s7)/126; wv=np.array([uh]+[uc]*5)
    cols=[np.array([Ah,0.0,0.0])]
    for j in range(5):
        th=2*math.pi*j/5; cols.append(np.array([Ac,Rc*math.cos(th),Rc*math.sin(th)]))
    A0=np.array([math.sqrt(wv[k])*cols[k][rr] for rr in range(3) for k in range(6)]).reshape(3,6)
    P=A0.T@np.linalg.inv(A0@A0.T)@A0
    w,Vs=np.linalg.eigh(P); idx=np.argsort(-w); Vs=Vs[:,idx]; w=w[idx]
    V=Vs[:,:3]*np.sqrt(np.maximum(w[:3],0))
    g=np.linalg.norm(V,axis=1)**2; U=V/np.sqrt(g)[:,None]
    return U,g,P
Uc,gc,Pc=champion()

def unpack(U,g,tol=1e-12):
    N=U.shape[0]; V=U*np.sqrt(g)[:,None]; Q=V@V.T; Qp=np.maximum(Q,0.0)
    w,ev=np.linalg.eigh(Qp); idx=np.argsort(-w); w=w[idx]; ev=ev[:,idx]
    rho=w[0]; x=ev[:,0].copy()
    if x.sum()<0: x=-x
    h=x/np.sqrt(g); c=g*h
    Kf=U@U.T; K=np.maximum(Kf,0.0)
    p=K@c; wp=c@U
    alpha=float(wp@wp); ell=float(c@(np.abs(Kf)@c)); beta=float(c@(np.maximum(-Kf,0.0)@c))
    sig=np.maximum(-Kf,0.0)@c
    npos=int((w>tol).sum()); n0=int((np.abs(w)<=tol).sum()); nneg=int((w<-tol).sum())
    return dict(N=N,Q=Q,Qp=Qp,w=w,ev=ev,rho=rho,x=x,h=h,c=c,K=K,Kf=Kf,p=p,wp=wp,alpha=alpha,
                ell=ell,beta=beta,sig=sig,g=g,U=U,Psi=float(c.sum()),H=float((c/p).sum()),
                npos=npos,n0=n0,nneg=nneg,E=np.maximum(-Q,0.0))

print("="*104); print("(1) R15-2  trace-1 rigidity:  h_i^2 = u_i^T Lt u_i ,  tr(Lt)=1 ,  Lt = (2/rho)*Lambda_SCF"); print("="*104)
def fit_Lt(U,h2):
    rows=[]
    for i in range(U.shape[0]):
        u=U[i]; rows.append([u[0]*u[0],u[1]*u[1],u[2]*u[2],2*u[0]*u[1],2*u[0]*u[2],2*u[1]*u[2]])
    A=np.array(rows); sol,_,rank,_=np.linalg.lstsq(A,np.array(h2),rcond=None)
    Lt=np.array([[sol[0],sol[3],sol[4]],[sol[3],sol[1],sol[5]],[sol[4],sol[5],sol[2]]])
    return Lt,np.abs(A@sol-np.array(h2)).max(),rank
d=unpack(Uc,gc); Lt,resid,rank=fit_Lt(Uc,d['h']**2)
evLt=np.linalg.eigvalsh(Lt)[::-1]
print("  champion: fit residual=%.2e  rank=%d  tr(Lt)=%.15f"%(resid,rank,np.trace(Lt)))
print("  eig(Lt)=",np.array2string(evLt,precision=12))
print("  2*lam_min+lam_max-1 = %.3e    (=> lam_min=(1-lam_max)/2)"%(2*evLt[2]+evLt[0]-1))
print("  lam_max(Lt) - h_hub^2 = %.3e"%(evLt[0]-d['h'][0]**2))
rr=np.random.default_rng(5); Z=rr.normal(size=(60000,3)); Z/=np.linalg.norm(Z,axis=1)[:,None]
thg=np.einsum('ij,jk,ik->i',Z,Lt,Z)
print("  min_{S^2} u^T Lt u = %.12f > 0  => the relaxation 'theta>=0 on ALL of S^2' is TIGHT at the champion"%thg.min())
# exact closed forms
S5,S7,S35=sp.sqrt(5),sp.sqrt(7),sp.sqrt(35)
g0S=(9-4*S5+S35)/18; gS=(45+4*S5-S35)/90; rho0S=(9+S5+2*S35)/18
a2=sp.Rational(91,126)-25*S7/126
hhub2=sp.simplify(a2/g0S)
print("  EXACT  h_hub^2 = lam_max(Lt) = %s = %.15f"%(sp.simplify(sp.rationalize(hhub2) if hasattr(sp,'rationalize') else hhub2),float(hhub2)))
print("         numerical lam_max = %.15f  diff=%.2e"%(evLt[0],abs(float(hhub2)-evLt[0])))
lmin=sp.simplify((1-hhub2)/2)
print("  EXACT  lam_min(Lt) = (1-h_hub^2)/2 = %s = %.15f  (numeric %.15f, diff %.2e)"%(lmin,float(lmin),evLt[2],abs(float(lmin)-evLt[2])))
print("  EXACT  B = lam_max-lam_min = (3*h_hub^2-1)/2 = %s = %.15f"%(sp.simplify(hhub2-lmin),float(hhub2-lmin)))
print("  check tr = 2*lam_min+lam_max -1 = %s"%sp.simplify(2*lmin+hhub2-1))

print(); print("="*104); print("(2) R15-3  acute-moment identity   h_i*W_i = rho*Lt*u_i ,   W_i = sum_{j:<u_i,u_j>>0} c_j u_j"); print("="*104)
def star_residual(U,g,Lt):
    dd=unpack(U,g); N=dd['N']; c=dd['c']; rho=dd['rho']; h=dd['h']
    acute=(U@U.T)>0.0
    W=np.array([(c[acute[i]][:,None]*U[acute[i]]).sum(axis=0) for i in range(N)])
    return np.abs(h[:,None]*W-rho*(U@Lt.T)).max(), dd, W
err,dd,W=star_residual(Uc,gc,Lt)
print("  champion  max|h_i W_i - rho Lt u_i| = %.3e   (rho=%.12f)"%(err,dd['rho']))
# corollary (b): |Lt^{-1/2}(w-S_i)| = rho  and  <zeta,S_i> = |S_i|^2_{Lt^-1}/2
Li=np.linalg.inv(Lt); Lh=np.linalg.cholesky(Li); zeta=Li@dd['wp']
Kf=Uc@Uc.T; obt=Kf<0
Svec=np.array([(dd['c'][obt[i]][:,None]*Uc[obt[i]]).sum(axis=0) if obt[i].any() else np.zeros(3) for i in range(6)])
nrm=np.array([math.sqrt(float((dd['wp']-Svec[i])@Li@(dd['wp']-Svec[i]))) for i in range(6)])
print("  |Lt^-1/2 (w-S_i)| =",np.array2string(nrm,precision=12),"  rho=%.12f  maxdev=%.2e"%(dd['rho'],np.abs(nrm-dd['rho']).max()))
print("  <zeta,S_i> - 0.5*S_i^T Lt^-1 S_i =",np.array2string(np.array([float(zeta@Svec[i])-0.5*float(Svec[i]@Li@Svec[i]) for i in range(6)]),precision=3))
fa=[i for i in range(6) if not obt[i].any()]
print("  fully-acute atoms = %s ;  u_hub=%s ;  zeta/|zeta|=%s"%(fa,np.array2string(Uc[fa[0]],precision=9),np.array2string(zeta/np.linalg.norm(zeta),precision=9)))
print("  rho/|zeta| = %.12f  vs  h_hub = %.12f   (=> all fully-acute atoms share ONE direction & ONE h)"%(dd['rho']/np.linalg.norm(zeta),dd['h'][fa[0]]))
print("  sigma_i =",np.array2string(dd['sig'],precision=12)); print("  S_i norms =",np.array2string(np.linalg.norm(Svec,axis=1),precision=12))

print(); print("="*104); print("(3) R15-4  sigma-identities (valid for EVERY configuration, no criticality needed)"); print("="*104)
rr2=np.random.default_rng(2024)
tests=[('champion',Uc,gc)]
for N in (4,5,6,7,8,9,10,12,14):
    for s in range(3):
        fm=random_frame_measure(N,seed=9100+17*N+s)
        if fm is not None: tests.append(('rand N=%d#%d'%(N,s),fm[0],fm[1]))
w1=w2=w3=w4=0.0
print("   %-15s %10s %12s %12s %12s"%("config","rho","|sc1|","|vec|","|beta|"))
for name,U,g in tests:
    dd=unpack(U,g); sig=dd['sig']
    a1=float(dd['g']@sig**2); b1=dd['beta']**2+dd['alpha']*(1-dd['alpha'])
    a2=(dd['g']*sig)@U; b2=(dd['rho']-1)*dd['wp']
    a3=float((dd['g']*dd['h']*sig).sum()); b3=dd['beta']
    a4=dd['rho']; b4=dd['alpha']+dd['beta']
    w1=max(w1,abs(a1-b1)); w2=max(w2,float(np.abs(a2-b2).max())); w3=max(w3,abs(a3-b3)); w4=max(w4,abs(a4-b4))
    if name=='champion' or name in ('rand N=7#0','rand N=9#0','rand N=14#0'):
        print("   %-15s %10.6f %12.2e %12.2e %12.2e"%(name,dd['rho'],abs(a1-b1),np.abs(a2-b2).max(),abs(a3-b3)))
print("  WORST over %d configs:  sum g sig^2 id %.2e | vector id %.2e | beta=sum g h sig %.2e | rho=alpha+beta %.2e"%(len(tests),w1,w2,w3,w4))

print(); print("="*104); print("(4) R15-5  spectral floor  n_+(Q_+) >= 1+(3-rho)^2/(3-rho^2)"); print("="*104)
print("   %-15s %10s %5s %4s %5s %9s %12s"%("config","rho","n_+","n_0","n_-","floor","margin"))
viol=0; wm=9e9
for name,U,g in tests:
    dd=unpack(U,g); rho=dd['rho']; fl=1+(3-rho)**2/(3-rho**2); m=dd['npos']-fl
    wm=min(wm,m)
    if m<-1e-9: viol+=1
    if name in ('champion','rand N=7#0','rand N=9#0','rand N=14#0'):
        print("   %-15s %10.6f %5d %4d %5d %9.4f %+12.4f"%(name,rho,dd['npos'],dd['n0'],dd['nneg'],fl,m))
print("  violations=%d  worst margin=%+.6f"%(viol,wm))
for lbl,r in (("rho0",rho0),("Phi",(3+s5)/4),("Mohar",(1+math.sqrt(3))/2),("5/4",1.25)):
    print("   floor at %-6s rho=%.9f : 1+(3-rho)^2/(3-rho^2) = %.6f  => n_+ >= %d"%(lbl,r,1+(3-r)**2/(3-r**2),math.ceil(1+(3-r)**2/(3-r**2))-1e-12))
dd=unpack(Uc,gc)
print("  champion: spec(Q_+)=",np.array2string(dd['w'],precision=9)," n_+=%d n_0=%d n_-=%d (floor=4 -> TIGHT)"%(dd['npos'],dd['n0'],dd['nneg']))

print(); print("="*104); print("(5) R15-6  Weyl: Q_+=Q+E, Q>=0, E=(-Q)_+ entrywise, trE=0  =>  nu_-(Q_+) <= nu_-(E)"); print("="*104)
bad=0
for name,U,g in tests:
    dd=unpack(U,g); nE=int((np.linalg.eigvalsh(dd['E'])<-1e-12).sum())
    if dd['nneg']>nE: bad+=1; print("   VIOLATION",name,dd['nneg'],nE)
print("  configs=%d violations=%d"%(len(tests),bad))
dd=unpack(Uc,gc); eE=np.linalg.eigvalsh(dd['E'])
print("  champion: nu_-(Q_+)=%d  nu_-(E)=%d  eig(E)=%s"%(dd['nneg'],int((eE<-1e-12).sum()),np.array2string(eE,precision=9)))
print("  tr(E)=%.2e ; E support == obtuse graph ; eig(E) matches g_p*|<pent,pentagram>|*eig(A(C5)):"%np.trace(dd['E']))
val=-float((Uc[1]@Uc[3]))*gc[1]
print("     predicted %s"%np.array2string(np.sort(np.r_[val*2,val*2*math.cos(2*math.pi/5),val*2*math.cos(2*math.pi/5),val*2*math.cos(4*math.pi/5),val*2*math.cos(4*math.pi/5),0.0]),precision=9))
print("  => if Q>=0 entrywise then E=0, nu_-(Q_+)=0, and cor:nuneg forces N<=6 at any interior local max")

print(); print("="*104); print("(6) R15-8  coincident atoms: merging preserves rho EXACTLY, splitting adds a zero eigenvalue"); print("="*104)
md=0.0
for trial in range(6):
    fm=random_frame_measure(7,seed=400+trial)
    if fm is None: continue
    U,g=fm; U2=U.copy(); U2[1]=U2[0].copy()
    d1=unpack(U2,g)
    gm=np.delete(g,1); gm[0]=g[0]+g[1]; Um=np.delete(U2,1,axis=0)
    d2=unpack(Um,gm)
    md=max(md,abs(d1['rho']-d2['rho']))
    print("   trial%d rho_split=%.13f rho_merged=%.13f diff=%.2e ; smallest |eig(Q+_split)|=%.2e"
          %(trial,d1['rho'],d2['rho'],abs(d1['rho']-d2['rho']),np.abs(d1['w']).min()))
print("  WORST |rho_split-rho_merged| = %.3e  => rho(mu) depends only on the MEASURE (thm:nfree is safe as stated)"%md)
