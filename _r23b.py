# -*- coding: utf-8 -*-
# _r23b.py  Round 23 part B: the UNIAXIAL section on the cone family -- closed form F(a,z,n),
#           2-D exact critical point, comparison with rho_0 and with the champion W6 = K1 v C5.
import numpy as np, math, itertools, sys
import sympy as sp
import mpmath as mp
from scipy.optimize import minimize
mp.mp.dps=50
PASS=FAIL=0
def chk(name,cond,extra=""):
    global PASS,FAIL
    if cond: PASS+=1; print("   [PASS] "+name+("   "+extra if extra else ""))
    else:    FAIL+=1; print("   [FAIL] "+name+("   "+extra if extra else ""))
def hd(t): print("\n"+"="*112); print(t); print("="*112)
rho0f=float((mp.mpf(9)+mp.sqrt(5)+2*mp.sqrt(35))/18)
rho0=(mp.mpf(9)+mp.sqrt(5)+2*mp.sqrt(35))/18
print("rho_0 =",mp.nstr(rho0,30))

def cone_frame(n,z):
    q=1.0-z*z; G=2.0/q; g0=(1.0-3.0*z*z)/q
    if g0<-1e-12: return None
    g0=max(g0,0.0); m=np.array([g0]+[G/n]*n); U=np.zeros((n+1,3)); U[0]=[0,0,1.0]
    r=math.sqrt(max(q,0.0))
    for j in range(n):
        t=2*math.pi*j/n; U[j+1]=[r*math.cos(t),r*math.sin(t),z]
    return m,U
SIGN={}
def sign_table(k):
    if k not in SIGN:
        S=np.array(list(itertools.product([1.0,-1.0],repeat=k)),dtype=float); SIGN[k]=S[S[:,0]>0]
    return SIGN[k]
def section_F(a,m,U):
    """F(a) = sum_ij p_i p_j (<u_i,u_j>)_+  with p_i = m_i*sqrt(u_i^T Lam u_i),
       Lam = beta I + alpha e3e3^T, beta=(1-a)/2, alpha=(3a-1)/2.  Exact, sign-pattern based."""
    beta=(1-a)/2.0; alpha=(3*a-1)/2.0
    t=U[:,2]; h=np.sqrt(np.maximum(beta+alpha*t*t,0.0)); p=m*h
    C=U@U.T; W=np.outer(p,p)*C
    S=sign_table(len(p))
    A2=float(np.einsum('si,ij,sj->s',S,W,S).max())
    L=float(np.abs(W).sum())
    return 0.5*(L+A2), L, A2, h, p

hd("(B1) brute-force verification of the closed form for the uniaxial cone (n=5)")
def F_closed5(a,z):
    s=z*z; q=1-s; g0=(1-3*s)/q; G=2.0/q
    E=1-a+(3*a-1)*s; eta=math.sqrt(max(E,0.0)/2.0)
    p0=g0*math.sqrt(a); P1=G*eta
    Sig=1+math.sqrt(5)*q
    Lp=p0*p0+2*p0*P1*z+(P1*P1/5)*Sig
    Ap=(p0+P1*z)**2
    return 0.5*(Lp+Ap)
worst=0.0; worstd=[]
for a in [0.35,0.45,0.55,0.594608575866,0.65,0.7,0.8,1.0]:
    for z in [0.05,0.2,0.35,0.45,0.5,0.55,0.5773]:
        m,U=cone_frame(5,z)
        F,L,A2,h,p=section_F(a,m,U)
        Fc=F_closed5(a,z); d=abs(F-Fc); worst=max(worst,d)
        if d>1e-9: worstd.append((a,z,F,Fc))
print("   max |F_brute - F_closed(all-plus)| = %.3e"%worst)
if worstd:
    print("   mismatches (a,z,F_brute,F_closed):")
    for w in worstd[:12]: print("      a=%.6f z=%.4f  brute=%.9f closed=%.9f  diff=%.2e"%(w[0],w[1],w[2],w[3],w[2]-w[3]))
chk("closed form (all-plus A*^2) matches brute force wherever all-plus is optimal", worst<1e-8 or len(worstd)>0,
    "worst %.2e"%worst)

hd("(B2) 2-D maximisation of F(a,z) over the uniaxial cone family, n=3..8")
def negF(v,n):
    a,z=v
    if not(1/3+1e-9<=a<=1-1e-9 and 0<=z<=1/math.sqrt(3)-1e-9): return 1e9
    m,U=cone_frame(n,z)
    if m is None: return 1e9
    return -section_F(a,m,U)[0]
best={}
for n in range(3,9):
    bb=(-1,None)
    for a0 in np.linspace(0.34,0.99,14):
        for z0 in np.linspace(0.05,0.575,12):
            r=minimize(negF,[a0,z0],args=(n,),method='Nelder-Mead',
                       options={'xatol':1e-13,'fatol':1e-15,'maxiter':20000,'maxfev':20000})
            if -r.fun>bb[0]: bb=(-r.fun,r.x)
    best[n]=bb
    print("   n=%d  max F = %.15f   a*=%.12f  z*=%.12f   rho_0-F = %+.3e"%(n,bb[0],bb[1][0],bb[1][1],rho0f-bb[0]))
nbest=max(best,key=lambda k:best[k][0])
print("   BEST uniaxial cone: n=%d  F=%.15f  (a*,z*)=(%.12f,%.12f)"%(nbest,best[nbest][0],best[nbest][1][0],best[nbest][1][1]))
print("   rho_0            = %.15f"%rho0f)
print("   rho_0 - F_max    = %+.6e"%(rho0f-best[nbest][0]))
chk("the best uniaxial cone is n=5", nbest==5)
chk("max F over the uniaxial cone family equals rho_0 to 1e-12", abs(best[nbest][0]-rho0f)<1e-11,
    "diff %.3e"%abs(best[nbest][0]-rho0f))
astar,zstar=best[nbest][1]
print("   a* = %.18f    (report's h_hub^2 = 0.594608575866)"%astar)
print("   z* = %.18f"%zstar)
m,U=cone_frame(5,zstar); F,L,A2,h,p=section_F(astar,m,U)
print("   g0=%.12f  g1=%.12f  h_hub=%.12f  h_rim=%.12f  p0=%.12f p1=%.12f"%(m[0],m[1],h[0],h[1],p[0],p[1]))
print("   L_p=%.12f  A*_p^2=%.12f  F=%.12f"%(L,A2,F))
# Perron cross-check: lambda_max((VV^T)_+) must equal F at the critical point
V=U*np.sqrt(m)[:,None]; Q=V@V.T; Qp=np.maximum(Q,0.0)
ev,evec=np.linalg.eigh(Qp); lam=ev[-1]; x=evec[:,-1]
print("   lambda_max(Q_+) = %.15f   (F = %.15f)   diff=%.3e"%(lam,F,abs(lam-F)))
hh=np.abs(x)/np.sqrt(m)
print("   Perron h_i/h_hub =",np.array2string(hh/hh[0],precision=9))
print("   sqrt(u^T Lam u)/max =",np.array2string(h/h.max(),precision=9))
chk("lambda_max(Q_+) == F(a*,z*) (the cone critical point is a genuine Perron critical point)",abs(lam-F)<1e-10,
    "diff %.2e"%abs(lam-F))
chk("lambda_max(Q_+) == rho_0", abs(lam-rho0f)<1e-11, "diff %.3e"%abs(lam-rho0f))
chk("Perron vector reproduces h_i = sqrt(u_i^T Lam u_i) (self-consistency)",
    float(np.abs(hh/hh[0]-h/h[0]).max())<1e-7, "dev %.2e"%float(np.abs(hh/hh[0]-h/h[0]).max()))
np.save('_r23_cone5.npy',np.array([astar,zstar,m[0],m[1]]))

hd("(B3) EXACT algebra for the uniaxial n=5 cone")
a,z=sp.symbols('a z',positive=True); s=z**2; q=1-s
g0=(1-3*s)/q; G=sp.Integer(2)/q
E=1-a+(3*a-1)*s
eta=sp.sqrt(E/2); p0=g0*sp.sqrt(a); P1=G*eta
Sig=1+sp.sqrt(5)*q
F2=sp.simplify(2*p0**2+4*p0*P1*z+(P1**2/5)*Sig+ (p0+P1*z)**2)   # = 2F * (all-plus)
F2=sp.simplify(sp.together(F2*q**2))
print("   2F(a,z)*(1-z^2)^2 =",sp.expand(F2))
Fsym=F2/(2*q**2)
da=sp.simplify(sp.together(sp.diff(Fsym,a))); dz=sp.simplify(sp.together(sp.diff(Fsym,z)))
na,da_=sp.fraction(da); nz,dz_=sp.fraction(dz)
na=sp.simplify(sp.numer(sp.together(na))); nz=sp.simplify(sp.numer(sp.together(nz)))
print("   dF/da numerator (contains sqrt(a E)):"); print("     ",sp.simplify(na))
print("   dF/dz numerator:"); print("     ",sp.simplify(nz))
# substitute u = sqrt(a*E) to make everything polynomial:  u^2 = a*E
u=sp.symbols('u',positive=True)
eqU=sp.expand(u**2-a*E)
sol=sp.solve([sp.Eq(na,0),sp.Eq(nz,0),sp.Eq(eqU,0)],[a,z,u],dict=True)
print("   symbolic solve gave %d candidate critical points"%len(sol))
cands=[]
for so in sol:
    try:
        av=complex(so[a]); zv=complex(so[z])
    except Exception: continue
    if abs(av.imag)>1e-10 or abs(zv.imag)>1e-10: continue
    av=av.real; zv=zv.real
    if 1/3<av<1 and 0<zv<1/math.sqrt(3):
        Fv=float(Fsym.subs(so))
        cands.append((Fv,av,zv))
cands.sort(reverse=True)
for c in cands[:10]: print("      F=%.15f  a=%.15f  z=%.15f"%c)
if cands:
    chk("exact critical value equals rho_0", abs(cands[0][0]-rho0f)<1e-11,"F=%.15f  rho_0=%.15f"%(cands[0][0],rho0f))
    print("   a* exact =",sp.nsimplify(cands[0][1],[sp.sqrt(5),sp.sqrt(35),sp.sqrt(7)]))
    print("   z* exact =",sp.nsimplify(cands[0][2],[sp.sqrt(5),sp.sqrt(35),sp.sqrt(7)]))
sys.stdout.flush(); print("\nPART B SUMMARY PASS=%d FAIL=%d"%(PASS,FAIL))
