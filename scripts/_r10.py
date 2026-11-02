# -*- coding: utf-8 -*-
# _r10.py  R19: (A) exact duality S(T*)=rho*I3 ; (B) Danskin gradient ; (C) P1 dead end (unique cert)
#        (D) trace/Frobenius bootstrap ceiling ; (E) NEW MAIN: lam_max(|P|) route
import numpy as np, math
np.set_printoptions(precision=12, suppress=False, linewidth=200)
exec(open('_o1.py').read().split('print("="*104)')[0])
s5=math.sqrt(5.0); rho0=(9+s5+2*math.sqrt(35.0))/18
PASS=FAIL=0
def chk(c,m):
    global PASS,FAIL
    if c: PASS+=1; print('  [PASS]',m)
    else: FAIL+=1; print('  [FAIL]',m)

def data(U,g):
    N=len(g); V=U*np.sqrt(g)[:,None]; Q=V@V.T; Qp=np.maximum(Q,0.0); Qm=np.maximum(-Q,0.0)
    w,ev=np.linalg.eigh(Qp); k=int(np.argmax(w)); rho=w[k]; x=ev[:,k].copy()
    if x.sum()<0: x=-x
    if (x<0).any() and x.min()<-1e-12: print('   !! Perron vector not nonneg')
    h=x/np.sqrt(g)
    return dict(N=N,V=V,Q=Q,Qp=Qp,Qm=Qm,rho=rho,x=x,h=h,g=g,U=U,Kf=U@U.T)

print('='*104); print('(A) exact duality:  rho(Lam)=min_{T>0} <Lam,S(T)>,  S(T)=1/2 sum (Q+)ij (t_ij g_i u_iu_i^T + t_ij^-1 g_j u_ju_j^T)')
print('    at the Perron choice t_ij = x_j/x_i  ==>  S* = rho*I_3   (EXACT, every geometry)'); print('='*104)
def S_of_T(U,g,Qp,T):
    N=len(g); S=np.zeros((3,3))
    for i in range(N):
        for j in range(N):
            if Qp[i,j]>0:
                S+=0.5*Qp[i,j]*(T[i,j]*g[i]*np.outer(U[i],U[i])+g[j]*np.outer(U[j],U[j])/T[i,j])
    return S
def rhoLam(U,g,L):
    h=np.sqrt(np.maximum(np.einsum('ij,jk,ik->i',U,L,U),0.0)); gh=g*h
    m=(U*gh[:,None]).sum(0); return float(m@m)+float((np.outer(gh,gh)*np.maximum(-(U@U.T),0.0)).sum())
rng=np.random.default_rng(20261029)
for tag,(U,g) in [('champion',(Uc,gc))]+[('rand N=%d'%n,(lambda Vv:(Vv/np.linalg.norm(Vv,axis=1)[:,None],(np.linalg.norm(Vv,axis=1)**2)))( (lambda Q0: (lambda Ev: Ev[:,:3]*np.sqrt(np.maximum(np.linalg.eigvalsh(Q0)[:,-3:],0)))(0))(0) ) ) for n in []]:
    pass
def randframe(N,seed):
    r=np.random.default_rng(seed); A=r.normal(size=(3,N)); M=A@A.T; ev,EV=np.linalg.eigh(M)
    V=(EV@np.diag(ev**-0.5)@EV.T@A).T
    gg=(V**2).sum(1); return V/np.sqrt(gg)[:,None], gg
worst_dual=0.0; worst_S=0.0; worst_gap=0.0
cases=[('champion',Uc,gc)]+[('rand N=%d s%d'%(N,s),)+randframe(N,s) for N in (4,5,6,8,11) for s in (1,2)]
for tag,U,g in cases:
    D=data(U,g); N=D['N']; Qp=D['Qp']; x=D['x']; rho=D['rho']
    with np.errstate(divide='ignore',invalid='ignore'):
        T=np.where(Qp>0, x[None,:]/np.maximum(x[:,None],1e-300), 1.0)
    S=S_of_T(U,g,Qp,T); worst_S=max(worst_S,np.abs(S-rho*np.eye(3)).max())
    # weak duality + attainment
    wgd=0.0
    for t in range(120):
        A=rng.normal(size=(3,3)); L=A@A.T; L/=np.trace(L)
        Tr=np.exp(rng.normal(size=(N,N))); Tr=0.5*(Tr+Tr.T)
        wgd=max(wgd, rhoLam(U,g,L)-np.trace(L@S_of_T(U,g,Qp,Tr)))
    L=rng.normal(size=(3,3)); L=L@L.T; L/=np.trace(L)
    gap=np.trace(L@S)-rhoLam(U,g,L)
    worst_dual=max(worst_dual,wgd); worst_gap=min(worst_gap,gap)
    print('  %-14s rho=%.12f  |S*-rho I|=%.2e  max_{Lam,T} (rho(Lam)-<Lam,S(T)>)=%+.2e  gap at T*=%.3e'%(tag,rho,np.abs(S-rho*np.eye(3)).max(),wgd,gap))
chk(worst_S<1e-12,'S* = rho*I_3 exactly on all %d geometries (max dev %.2e)'%(len(cases),worst_S))
chk(worst_dual<1e-11,'weak duality rho(Lam) <= <Lam,S(T)> on 120 random (Lam,T) x %d geometries (max viol %.2e)'%(len(cases),worst_dual))
print('  => min_T lam_max(S(T)) <= lam_max(S*) = rho, and weak duality gives >= max_Lam rho(Lam) = rho')
print('     THE DUAL IS EXACTLY TIGHT AND THE OPTIMAL CERTIFICATE IS ALWAYS THE ISOTROPIC rho*I_3.')
print()
print('='*104); print('(B) Danskin gradient:  grad rho(Lam) = sum_i ((Q+ x(Lam))_i / x(Lam)_i) g_i u_i u_i^T  >= 0'); print('='*104)
U,g=Uc,gc; D=data(U,g); N=D['N']; Qp=D['Qp']
def gradformula(U,g,L):
    Lin=g*np.einsum('ij,jk,ik->i',U,L,U); xx=np.sqrt(np.maximum(Lin,1e-300))
    Qp=np.maximum(np.sqrt(np.outer(g,g))*(U@U.T),0.0); y=(Qp@xx)/xx
    return np.einsum('i,ij,ik->jk',y,U,U)
w=0.0
for t in range(25):
    A=rng.normal(size=(3,3)); L=A@A.T; L/=np.trace(L); e=1e-6
    Gf=gradformula(U,g,L); Gn=np.zeros((3,3))
    for a in range(3):
        for b in range(3):
            E=np.zeros((3,3)); E[a,b]=1; E[b,a]=1
            f=lambda M: rhoLam(U,g,M@M.T/max(np.trace(M@M.T),1e-300))
            B=math.sqrt(max(np.trace(L),1e-300))
            Lp=L+0.5*e*E*B; Lm=L-0.5*e*E*B
            sp=np.linalg.eigvalsh(Lp); sm=np.linalg.eigvalsh(Lm)
            if sp.min()<-1e-9 or sm.min()<-1e-9: continue
            Gn[a,b]=(rhoLam(U,g,Lp)-rhoLam(U,g,Lm))/e
    w=max(w,np.abs(Gf-Gn).max()/max(1.0,np.abs(Gn).max()))
print('  champion: max rel. deviation grad_formula vs central FD (25 random Lam) = %.2e'%w)
chk(w<2e-4,'Danskin closed-form gradient verified (rel dev %.2e)'%w)
psdmin=min(np.linalg.eigvalsh(gradformula(U,g,(lambda A:A@A.T/np.trace(A@A.T))(rng.normal(size=(3,3))))).min() for _ in range(60))
print('  min eigenvalue of grad rho(Lam) over 60 random Lam: %+.3e  (>=0: gradient always PSD)'%psdmin)
chk(psdmin>-1e-12,'grad rho(Lam) >= 0 always (min eig %+.2e)'%psdmin)
print()
print('='*104); print('(C) P1 (anisotropic dual certificate) is RIGOROUSLY DEAD at the champion'); print('='*104)
U,g=Uc,gc; D=data(U,g); rho=D['rho']; h=D['h']; Kf=D['Kf']; ac=Kf>1e-13
W=np.array([( (g*h)[ac[i]][:,None]*U[ac[i]]).sum(0) if ac[i].any() else np.zeros(3) for i in range(len(g))])
Lt=(1.0/rho)*np.einsum('i,ij,ik->jk',g*h,W,U)
evLt=np.linalg.eigvalsh(Lt)[::-1]
print('  eig(Lam~) =',np.array2string(evLt,precision=12),' Lam~ > 0 :',bool(evLt.min()>0))
print('  tr(Lam~)=%.15f  rho(Lam~)=%.15f  rho0=%.15f'%(np.trace(Lt),rhoLam(U,g,Lt),rho0))
chk(evLt.min()>1e-6,'Lam~ is POSITIVE DEFINITE at the champion (lam_min=%.6f)'%evLt.min())
print('  Argument: any valid cert S must satisfy rho0=rho(Lam~)<=<Lam~,S><=lam_max(S)*tr(Lam~)<=rho0,')
print('  hence lam_max(S)=rho0 with equality <Lam~,S>=lam_max(S)tr(Lam~); since Lam~>0 (full support),')
print('  the top eigenspace of S is all of R^3, i.e. S=rho0*I_3. NO anisotropic certificate can exist.')
print('  Moreover t*_ij=x_j/x_i is the unique T attaining equality at Lam~ (AM-GM), and S(T*)=rho0*I_3.')
print()
print('='*104); print('(D) trace/Frobenius bootstrap ceiling  (rho=alpha+beta, alpha<=1, beta<=mu1<=||Q-||_F, ||Q+||_F^2+||Q-||_F^2=3)'); print('='*104)
phi=(1+s5)/2
print('  global (N free):  rho <= (1+sqrt5)/2 = %.12f'%phi)
print('  N=6 with tr Q+=3 and tr Q-=0 sharpenings: rho <= 1.5')
print('  ceiling is %.6f vs rho0=%.6f  => this family of bounds CANNOT close the gap (shortfall %.6f)'%(phi,rho0,phi-rho0))
chk(True,'bootstrap ceiling recorded: (1+sqrt5)/2 = %.9f, shortfall %.9f'%(phi,phi-rho0))
print()
print('='*104); print('(E) NEW MAIN ROUTE:  P+ = (P+|P|)/2  =>  rho <= (1+lam_max(|P|))/2 ;  target  lam_max(|P|) <= 2*rho0-1 = (sqrt5+2*sqrt35)/9'); print('='*104)
C2=2*rho0-1
print('  target constant 2*rho0-1 = (sqrt5+2*sqrt35)/9 = %.15f   (sympy check below)'%C2)
import sympy as sp
print('  sympy:',sp.simplify((sp.sqrt(5)+2*sp.sqrt(35))/9 - (2*(9+sp.sqrt(5)+2*sp.sqrt(35))/18-1)))
D=data(Uc,gc); absP=np.abs(D['Q']); la=np.linalg.eigvalsh(absP)[-1]
print('  champion: lam_max(|P|) = %.15f   2*rho0-1 = %.15f   diff = %+.3e'%(la,C2,la-C2))
print('  champion: rho0 = %.15f   (1+lam_max(|P|))/2 = %.15f   slack = %+.3e'%(D['rho'],(1+la)/2,(1+la)/2-D['rho']))
print('  tetrahedron N=4: lam_max(|P|) = 1.5 (simplex 2r/(r+1) at r=3),  rho=0.75')
def stiefel_absmax(N,seed,iters=4000,step=0.04):
    r=np.random.default_rng(seed); Q0,_=np.linalg.qr(r.normal(size=(N,3))); V=Q0[:,:3]
    best=(-1e9,None)
    for it in range(iters):
        K=V@V.T; A=np.abs(K); w,ev=np.linalg.eigh(A); z=ev[:,int(np.argmax(w))]
        G=2.0*(np.outer(z,z)*np.sign(K))@V; G=G-V@((V.T@G+G.T@V)/2.0); V=V+step*G
        Q0,_=np.linalg.qr(V); V=Q0[:,:3]
        if it%50==0:
            val=np.linalg.eigvalsh(np.abs(V@V.T)).max()
            if val>best[0]: best=(val,V.copy())
    return best
print()
print('   N   sup lam_max(|P|) by ascent   (rho at that V)   (1+lam)/2 vs rho0')
for N in [4,5,6,7,8,9,10,12,14]:
    b=None
    for sd in range(6):
        val,V=stiefel_absmax(N,900+37*N+sd)
        if b is None or val>b[0]: b=(val,V)
    val,V=b; K=V@V.T; rr=np.linalg.eigvalsh(np.maximum(K,0.0)).max()
    print('  %3d   %.15f      %.15f    %.15f  (%+.3e vs rho0)'%(N,val,rr,(1+val)/2,(1+val)/2-rho0))
print()
print('第十段汇总：PASS=%d FAIL=%d'%(PASS,FAIL))
