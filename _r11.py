# -*- coding: utf-8 -*-
# _r11.py  R19: spherical-harmonic (Legendre) decomposition of P+ ; exact c_ell ; champion budget split
import numpy as np, math, sympy as sp
np.set_printoptions(precision=12, suppress=False, linewidth=200)
exec(open('_o1.py').read().split('print("="*104)')[0])
s5=math.sqrt(5.0); rho0=(9+s5+2*math.sqrt(35.0))/18; phi=(1+s5)/2
t=sp.symbols('t')
def cl(ell):
    P=sp.legendre(ell,t); return sp.simplify((2*ell+1)/2*sp.integrate(t*P,(t,0,1)))
print('='*104); print('(1) exact Legendre coefficients  (t)_+ = sum_ell c_ell P_ell(t) on [-1,1]'); print('='*104)
C={}
for ell in range(0,19):
    C[ell]=cl(ell)
    if ell<=14: print('   c_%-2d = %-22s = %+.15f'%(ell,sp.sstr(C[ell]),float(C[ell])))
print('   odd ell>=3 vanish:',all(C[l]==0 for l in [3,5,7,9,11,13,15,17]))
print('   check partial sums at t=0.37,-0.62:',end=' ')
for L in (2,4,6,10,18):
    f=sum(float(C[l])*sp.legendre(l,t) for l in range(0,L+1))
    print('L=%d: %+.5f/%+.5f'%(L,float(f.subs(t,sp.Rational(37,100))),float(f.subs(t,sp.Rational(-62,100)))),end='  ')
print()
print('   EXACT: c0=1/4, c1=1/2, c2=5/16, c4=-3/32, c6=13/256 ->',sp.simplify(C[0]),sp.simplify(C[1]),sp.simplify(C[2]),sp.simplify(C[4]),sp.simplify(C[6]))

def data(U,g):
    V=U*np.sqrt(g)[:,None]; Q=V@V.T; Qp=np.maximum(Q,0.0)
    w,ev=np.linalg.eigh(Qp); k=int(np.argmax(w)); rho=w[k]; x=ev[:,k].copy()
    if x.sum()<0: x=-x
    return V,Q,Qp,rho,x
print(); print('='*104); print('(2) champion budget: rho = sum_ell c_ell * x^T A_ell x,  (A_ell)_ij = sqrt(g_i g_j) P_ell(<u_i,u_j>)  >= 0 (PSD kernel)'); print('='*104)
U,g=Uc,gc; V,Q,Qp,rho,x=data(U,g); Kf=U@U.T; N=len(g); sg=np.sqrt(g)
print('   rho = %.15f   rho0 = %.15f   alpha = x^T P x = %.15f'%(rho,rho0,float(x@(Q@x))))
print('   (sum_i sqrt(g_i) x_i)^2 = %.15f   [lam_max(A_0)=sum g_i=3]'%float((sg@x)**2))
rows=[]; tot=0.0
for ell in range(0,25):
    Pl=np.array([[float(sp.legendre(ell,t).subs(t,sp.nsimplify(Kf[i,j],rational=True,maxdigits=14))) if abs(Kf[i,j])<1-1e-12 else float(sp.legendre(ell,t).subs(t,sp.Rational(1))) for j in range(N)] for i in range(N)])
    A=np.outer(sg,sg)*Pl
    q=float(x@A@x); lm=float(np.linalg.eigvalsh(0.5*(A+A.T)).max()); cf=float(C[ell]) if ell in C else float(cl(ell))
    tot+=cf*q; rows.append((ell,cf,q,lm,cf*q))
    if ell<=14: print('   ell=%-2d c=%+.10f  x^T A x=%+.12f  lam_max(A)=%+.9f  contribution=%+.12f'%(ell,cf,q,lm,cf*q))
print('   sum of contributions (ell<=24) = %.15f  vs rho = %.15f  (diff %.2e)'%(tot,rho,tot-rho))
s01=rows[0][4]+rows[1][4]; srest=tot-s01
print()
print('   ell<=1 part  c0*x^T A0 x + c1*x^T A1 x = %.15f'%(s01))
print('   ell<=1 CEILING (c0*lam_max(A0)+c1*lam_max(A1) = 3/4+1/2) = %.15f  = 5/4'%(0.75+0.5))
print('   ell>=2 remainder at the champion = %.15f'%(srest))
print('   rho0 - 5/4 = %.15f      rho0 - (ell<=1 actual) = %.15f'%(rho0-1.25,rho0-s01))
print('   ell=2 term = %+.12f ; ell=4 term = %+.12f ; ell=6 = %+.12f ; ell=8 = %+.12f'%(rows[2][4],rows[4][4],rows[6][4],rows[8][4]))
print()
print('='*104); print('(3) how much of the ell>=2 remainder can be captured by lam_max bounds alone?'); print('='*104)
pos=sum(r[4] for r in rows if r[4]>0); neg=sum(r[4] for r in rows if r[4]<0)
print('   sum of POSITIVE contributions = %+.12f ; sum of NEGATIVE = %+.12f'%(pos,neg))
print('   termwise ceiling sum_ell c_ell lam_max(A_ell) over ell<=24 =')
tc=sum(r[1]*r[3] for r in rows if r[1]>0)
print('      %.12f   (>> rho0 = %.12f) => termwise lam_max bounds are USELESS'%(tc,rho0))
print()
print('='*104); print('(4) same split for other geometries (is the ell<=1 part always <= 5/4? and where is rho large?)'); print('='*104)
def randframe(Nn,seed):
    r=np.random.default_rng(seed); A=r.normal(size=(3,Nn)); M=A@A.T; ev,EV=np.linalg.eigh(M)
    V=(EV@np.diag(ev**-0.5)@EV.T@A).T; gg=(V**2).sum(1); return V/np.sqrt(gg)[:,None],gg
def icosahedral():
    gr=(1+math.sqrt(5))/2; pts=[]
    for s1 in (1,-1):
        for s2 in (1,-1):
            pts.append(np.array([0.0,s1*1.0,s2*gr])); pts.append(np.array([s1*1.0,s2*gr,0.0])); pts.append(np.array([s2*gr,0.0,s1*1.0]))
    P=np.array(pts); P=P/np.linalg.norm(P,axis=1)[:,None]
    M=P.T@P; ev,EV=np.linalg.eigh(M); V=(EV@np.diag(ev**-0.5)@EV.T@P.T).T
    gg=(V**2).sum(1); return V/np.sqrt(gg)[:,None],gg
Ui,gi=icosahedral(); Ki=Ui@Ui.T
print('   icosahedral ETF: |<u_i,u_j>| off-diag = %.15f  (1/sqrt5=%.15f)'%(np.abs(Ki[np.triu_indices(6,1)]).mean(),1/math.sqrt(5)))
Vd,Qd,Qpd,rd,xd=data(Ui,gi); print('   icosahedral ETF: lam_max(|P|) = %.15f  (=phi? diff %.2e)   rho = %.15f'%(np.linalg.eigvalsh(np.abs(Qd)).max(),np.linalg.eigvalsh(np.abs(Qd)).max()-phi,rd))
for tag,(Ux,gx) in [('icosahedron',(Ui,gi)),('tetrahedron',randframe(4,3)),('rand N=6',(lambda f:f)(randframe(6,11))),('rand N=9',randframe(9,12))]:
    Vx,Qx,Qpx,rxx,xx=data(Ux,gx); Kx=Ux@Ux.T; sgx=np.sqrt(gx)
    tot=0.0; s01=0.0
    for ell in range(0,21):
        Pl=np.array([[float(sp.legendre(ell,t).subs(t,sp.nsimplify(Kx[i,j],rational=True,maxdigits=14))) for j in range(len(gx))] for i in range(len(gx))])
        A=np.outer(sgx,sgx)*Pl; cf=float(C[ell]) if ell in C else float(cl(ell)); tot+=cf*float(xx@A@xx)
        if ell<=1: s01+=cf*float(xx@A@xx)
    print('   %-12s rho=%.12f  ell<=1 part=%.12f (ceiling 1.25)  ell>=2 rem=%+.12f  recon=%.12f'%(tag,rxx,s01,tot-s01,tot))
