# -*- coding: utf-8 -*-
# _r11b.py  R19: Legendre budget split of rho (numpy evaluation, exact c_ell from sympy)
import numpy as np, math, sympy as sp
from numpy.polynomial.legendre import legvander
np.set_printoptions(precision=12, suppress=False, linewidth=200)
exec(open('_o1.py').read().split('print("="*104)')[0])
s5=math.sqrt(5.0); rho0=(9+s5+2*math.sqrt(35.0))/18; phi=(1+s5)/2
tt=sp.symbols('t')
Cex={l:sp.nsimplify((2*l+1)/2*sp.integrate(tt*sp.legendre(l,tt),(tt,0,1)),rational=True) for l in range(0,25)}
C={l:float(Cex[l]) for l in Cex}
print('exact c_ell (ell<=12):',[ (l,str(Cex[l])) for l in range(13)])
def data(U,g):
    V=U*np.sqrt(g)[:,None]; Q=V@V.T; Qp=np.maximum(Q,0.0)
    w,ev=np.linalg.eigh(Qp); k=int(np.argmax(w)); rho=w[k]; x=ev[:,k].copy()
    if x.sum()<0: x=-x
    return V,Q,Qp,rho,x
def budget(U,g,L=24):
    V,Q,Qp,rho,x=data(U,g); K=np.clip(U@U.T,-1.0,1.0); sg=np.sqrt(g); N=len(g)
    Vand=legvander(K,L)                      # Vand[i,j,l]=P_l(K_ij)
    out=[]
    for l in range(L+1):
        A=np.outer(sg,sg)*Vand[:,:,l]; A=0.5*(A+A.T)
        out.append(dict(l=l,c=C[l],q=float(x@A@x),lam=float(np.linalg.eigvalsh(A).max()),term=C[l]*float(x@A@x)))
    return rho,x,out,Q
print(); print('='*104); print('champion budget:  rho = sum_ell c_ell x^T A_ell x ,  A_ell = sqrt(g g^T) odot P_ell(Kf)  (PSD kernel => A_ell >= 0 in Loewner order)'); print('='*104)
rho,x,out,Q=budget(Uc,gc)
tot=sum(o['term'] for o in out)
for o in out[:13]:
    print('  ell=%-2d c=%+.10f  x^T A x=%+.12f  lam_max(A)=%+.9f  term=%+.12f  cum=%+.12f'%(
        o['l'],o['c'],o['q'],o['lam'],o['term'],sum(p['term'] for p in out[:o['l']+1])))
print('  recon(ell<=24) = %.15f   rho = %.15f   diff = %.2e'%(tot,rho,tot-rho))
s01=out[0]['term']+out[1]['term']; srest=tot-s01
print()
print('  ell<=1 actual   = %.12f   (termwise ceiling 3/4+1/2 = 5/4 = %.12f)'%(s01,1.25))
print('  ell>=2 remainder= %+.12f'%(srest))
print('  ell=2 %+.9f | ell=4 %+.9f | ell=6 %+.9f | ell=8 %+.9f | ell=10 %+.9f | ell=12 %+.9f'%(
    out[2]['term'],out[4]['term'],out[6]['term'],out[8]['term'],out[10]['term'],out[12]['term']))
print('  sum POSITIVE terms = %+.12f ; sum NEGATIVE = %+.12f'%(sum(o['term'] for o in out if o['term']>0),sum(o['term'] for o in out if o['term']<0)))
print('  termwise ceiling sum_{c>0} c*lam_max(A) = %.9f  >> rho0=%.9f  (USELESS)'%(sum(o['c']*o['lam'] for o in out if o['c']>0),rho0))
print()
print('='*104); print('Is the ell>=2 remainder controlled?  compare with  ||(P+)_ij - 1/4 sqrt(gigj) - 1/2 P_ij|| spectral'); print('='*104)
sg=np.sqrt(gc); K=np.clip(Uc@Uc.T,-1,1); A0=np.outer(sg,sg); A1=Q
R=Qp-0.25*A0-0.5*A1
evR=np.linalg.eigvalsh(R)[::-1]
print('  R = P+ - 1/4 A0 - 1/2 A1 : eig =',np.array2string(evR,precision=9))
print('  lam_max(R) = %.12f ;  rho0 - 5/4 = %.12f ;  lam_max(R)-(rho0-5/4) = %+.12f'%(evR[0],rho0-1.25,evR[0]-(rho0-1.25)))
print('  x^T R x = %.12f  (= ell>=2 remainder, must match %.12f)'%(float(x@R@x),srest))
print()
print('='*104); print('same split on other geometries: is  lam_max(P+) <= 1/4*lam(A0)+1/2*lam(A1)+lam_max(R)  ever > rho0 ?'); print('='*104)
def randframe(Nn,seed):
    r=np.random.default_rng(seed); A=r.normal(size=(3,Nn)); M=A@A.T; ev,EV=np.linalg.eigh(M)
    V=(EV@np.diag(ev**-0.5)@EV.T@A).T; gg=(V**2).sum(1); return V/np.sqrt(gg)[:,None],gg
def icosahedral():
    gr=(1+math.sqrt(5))/2; pts=[]
    for s1 in (1,-1):
        for s2 in (1,-1):
            pts += [np.array([0.0,s1,s2*gr]),np.array([s1,s2*gr,0.0]),np.array([s2*gr,0.0,s1])]
    P=np.array(pts); P=P/np.linalg.norm(P,axis=1)[:,None]
    M=P.T@P; ev,EV=np.linalg.eigh(M); V=(EV@np.diag(ev**-0.5)@EV.T@P.T).T
    gg=(V**2).sum(1); return V/np.sqrt(gg)[:,None],gg
Ui,gi=icosahedral(); Vi,Qi,Qpi,ri,xi=data(Ui,gi)
print('  icosahedral ETF: off-diag |<u,u>| mean = %.12f (1/sqrt5=%.12f)  lam_max(|P|)=%.12f (phi=%.12f) rho=%.12f'%(
    np.abs(Qi[np.triu_indices(6,1)]).mean(),1/math.sqrt(5),np.linalg.eigvalsh(np.abs(Qi)).max(),phi,ri))
worst=-1e9
for tag,(Ux,gx) in [('champion',(Uc,gc)),('icosahedron',(Ui,gi)),('tetra',randframe(4,3)),('rand6',randframe(6,11)),('rand9',randframe(9,12)),('rand12',randframe(12,5))]:
    rr,xx,oo,QQ=budget(Ux,gx); sgx=np.sqrt(gx); Kx=Ux@Ux.T
    Rx=np.maximum(QQ,0)-0.25*np.outer(sgx,sgx)-0.5*QQ
    lmx=np.linalg.eigvalsh(Rx).max()
    print('  %-11s rho=%.9f  ell<=1=%.9f  ell>=2 rem=%+.9f  lam_max(R)=%+.9f  ceiling 1.25+lam_max(R)=%.9f'%(
        tag,rr,oo[0]['term']+oo[1]['term'],sum(o['term'] for o in oo[2:]),lmx,1.25+lmx))
    worst=max(worst,1.25+lmx)
print('  worst ceiling over the tested geometries = %.9f  (rho0 = %.9f)'%(worst,rho0))
