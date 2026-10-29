# -*- coding: utf-8 -*-
# _r23a.py  Round 23 part A: CONE FAMILIES, exact closed forms (centre + uniaxial sections).
import numpy as np, math, itertools, sys
import sympy as sp
import mpmath as mp
mp.mp.dps = 50
PASS=FAIL=0
def chk(name, cond, extra=""):
    global PASS, FAIL
    if cond: PASS+=1; print("   [PASS] "+name+("   "+extra if extra else ""))
    else:    FAIL+=1; print("   [FAIL] "+name+("   "+extra if extra else ""))
def hd(t): print("\n"+"="*112); print(t); print("="*112)
rho0=(mp.mpf(9)+mp.sqrt(5)+2*mp.sqrt(35))/18; phi=(1+mp.sqrt(5))/2
print("rho_0 =",mp.nstr(rho0,25)," 6rho_0 =",mp.nstr(6*rho0,25)," 3rho_0 =",mp.nstr(3*rho0,25))
print("3phi  =",mp.nstr(3*phi,25),"  19/81 =",mp.nstr(mp.mpf(19)/81,25))

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
def Wmat(m,U):
    """W_ij = m_i m_j <u_i,u_j>  -- the correct weight for L and A*^2."""
    return np.outer(m,m)*(U@U.T)
def Lof(m,U): return float(np.abs(Wmat(m,U)).sum())
def Astar2_brute(m,U):
    W=Wmat(m,U); S=sign_table(W.shape[0]); return float(np.einsum('si,ij,sj->s',S,W,S).max())
def Astar2_w(m,U,nz=801,npp=121):
    th=np.arctan2(U[:,1],U[:,0]); rr=np.hypot(U[:,0],U[:,1]); zz=U[:,2]
    psi=np.linspace(0,2*math.pi,npp,endpoint=False); best=-1e18
    C=np.cos(th[None,:]-psi[:,None])                     # npp x n
    for zeta in np.linspace(-1,1,nz):
        b=math.sqrt(max(1-zeta*zeta,0.0))
        A=np.abs(rr*b*C+zz*zeta)@m
        best=max(best,float((A*A).max()))
    return best
def isodef(m,U): return float(abs((m[:,None,None]*(U[:,:,None]*U[:,None,:])).sum(0)-np.eye(3)).max())

hd("(A) CENTRE SECTION, CONE FAMILY: exact L_n(z), exact A*^2, sup over n")
def L_closed(n,z):
    s=z*z; q=1-s; G=2.0/q; g0=(1-3*s)/q
    tot=sum(abs(q*math.cos(2*math.pi*d/n)+s) for d in range(n))
    return g0*g0+2*g0*G*z+(G*G/n)*tot
worst=0.0; worstA=0.0; worstI=0.0
for n in range(3,9):
    for z in [0.05,0.2,0.35,0.45,0.52,0.55,0.5773]:
        m,U=cone_frame(n,z); worstI=max(worstI,isodef(m,U))
        worst=max(worst,abs(Lof(m,U)-L_closed(n,z)))
        worstA=max(worstA,abs(Astar2_brute(m,U)-Astar2_w(m,U,nz=1201,npp=181)))
chk("cone frames are isotropic: sum m_i u_i u_i^T = I_3", worstI<1e-12, "max dev %.2e"%worstI)
chk("L closed form == sum_ij m_i m_j |<u_i,u_j>| (brute force)", worst<1e-11, "max dev %.2e"%worst)
chk("A*^2 : sign enumeration == max over w=(zeta,psi) of (sum m|<u,w>|)^2", worstA<1e-7, "max dev %.2e"%worstA)

zg=np.linspace(0.0,1/math.sqrt(3),601)
print("      n   max_z [L+A*^2]      z*     kappa3^c   all-plus optimal?   6rho_0-J    L      A*^2")
bestn=None
for n in range(3,21):
    vals=np.empty(len(zg)); plus=np.empty(len(zg)); Lv=np.empty(len(zg)); Av=np.empty(len(zg))
    for i,zz in enumerate(zg):
        m,U=cone_frame(n,zz); L=L_closed(n,zz)
        A=Astar2_brute(m,U) if n<=10 else Astar2_w(m,U,nz=601,npp=91)
        Lv[i]=L; Av[i]=A; vals[i]=L+A; plus[i]=L+(m[0]+m[1:].sum()*zz)**2
    k=int(np.argmax(vals)); flag="yes" if abs(vals[k]-plus[k])<1e-8 else "NO"
    print("     %2d   %.12f  %.6f  %.9f     %-4s   %+.6f   %.6f %.6f"%(n,vals[k],zg[k],vals[k]/6,flag,float(6*rho0)-vals[k],Lv[k],Av[k]))
    if bestn is None or vals[k]>bestn[0]: bestn=(vals[k],n,zg[k],Lv[k],Av[k])
print("   BEST cone (centre section): n=%d  J=%.12f  z=%.6f  L=%.9f  A*^2=%.9f"%(bestn[1],bestn[0],bestn[2],bestn[3],bestn[4]))
chk("n=5 cone is the best cone in the centre section", bestn[1]==5)
chk("best cone J < 6 rho_0", bestn[0]<float(6*rho0), "margin %.9f"%(float(6*rho0)-bestn[0]))

hd("(A2) EXACT algebra: n=5 cone, centre section")
z=sp.symbols('z',positive=True); s=z**2; q=1-s
Sig5=1+sp.sqrt(5)*q; g0=(1-3*s)/q; G=sp.Integer(2)/q
Lsym=sp.simplify(g0**2+2*g0*G*z+4*Sig5/(5*q**2)); Asym=sp.simplify((g0+G*z)**2)
A5=sp.expand(sp.together((Lsym+Asym)*5*q**2))
print("   A(z) := 5(1-z^2)^2 J_5(z) =",sp.factor(A5))
J5sym=A5/(5*q**2); crit=sp.simplify(sp.factor(sp.numer(sp.together(sp.diff(J5sym,z)))))
P=15*z**4-(35-sp.sqrt(5))*z**3+30*z**2+(3-sp.sqrt(5))*z-5
chk("dJ_5/dz = 0  <=>  P(z)=15z^4-(35-sqrt5)z^3+30z^2+(3-sqrt5)z-5=0",sp.rem(crit,P,z)==0,
    "crit = (%s)*P"%sp.simplify(crit/sp.expand(P)))
rts=sp.nroots(P,n=40); print("   roots of P:",[str(t) for t in rts])
zin=[t for t in rts if abs(sp.im(t))<1e-30 and 0<sp.re(t)<1/mp.sqrt(3)]
chk("P has exactly one root in (0,1/sqrt(3))",len(zin)==1)
zstar=mp.mpf(str(sp.re(zin[0]))); print("   z* =",mp.nstr(zstar,40))
def J5mp(zz): return (2*(1-3*zz**2)**2+8*zz*(1-3*zz**2)+4*zz**2+(mp.mpf(4)/5)*(1+mp.sqrt(5)*(1-zz**2)))/(1-zz**2)**2
J5s=J5mp(zstar)
print("   J_5(z*) =",mp.nstr(J5s,40)); print("   kappa3^c(cone)=J_5(z*)/6 =",mp.nstr(J5s/6,40))
print("   margin 6rho_0-J_5(z*) =",mp.nstr(6*rho0-J5s,30))
chk("J_5(z*) == R22 numerical centre sup 7.595950307218",abs(float(J5s)-7.595950307218)<1e-9,"diff %.2e"%abs(float(J5s)-7.595950307218))
mz,Uz=cone_frame(5,float(zstar))
print("   at z*: g0=%.12f  g1=%.12f  L=%.12f  A*^2=%.12f  (all-plus ||M1||^2=%.12f)"%(
      mz[0],mz[1],Lof(mz,Uz),Astar2_brute(mz,Uz),(mz[0]+mz[1:].sum()*float(zstar))**2))
chk("at z* the all-plus sign pattern is optimal, A*^2=||M1||^2=(g0+Gz)^2",
    abs(Astar2_brute(mz,Uz)-(mz[0]+mz[1:].sum()*float(zstar))**2)<1e-10)
JJ=sp.Symbol('JJ')
res=sp.together(sp.resultant(P,A5-JJ*5*q**2,z)); resn=sp.factor(sp.numer(res))
fl=sp.factor_list(resn)
print("   resultant_z(P, A(z)-5JJ(1-z^2)^2) :  const",sp.nsimplify(fl[0]))
good=None
for f,e in fl[1]:
    lam=sp.lambdify(JJ,f,'mpmath')
    try: v=abs(lam(J5s))
    except Exception: v=mp.mpf('nan')
    print("      ^%d deg=%s  |f(J5*)|=%s"%(e,sp.degree(sp.Poly(f,JJ)) if f.is_polynomial(JJ) else "?",mp.nstr(v,4)))
    if good is None and v==v and v<mp.mpf('1e-15'): good=sp.expand(f)
chk("J_5(z*) is algebraic of degree <=4 over Q(sqrt(5)); minimal polynomial found",good is not None,
    ("\n        "+str(good)) if good is not None else "")
sys.stdout.flush(); print("\nPART A SUMMARY PASS=%d FAIL=%d"%(PASS,FAIL))
