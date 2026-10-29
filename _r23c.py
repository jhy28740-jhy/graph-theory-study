# -*- coding: utf-8 -*-
# _r23c.py  Round 23 parts B3(exact)/C/D/E
import numpy as np, math, itertools, sys
import sympy as sp
import mpmath as mp
from scipy.optimize import minimize
mp.mp.dps=70
PASS=FAIL=0
def chk(n,c,e=""):
    global PASS,FAIL
    if c: PASS+=1; print("   [PASS] "+n+("   "+e if e else ""))
    else:    FAIL+=1; print("   [FAIL] "+n+("   "+e if e else ""))
def hd(t): print("\n"+"="*112); print(t); print("="*112)
S5=mp.sqrt(5); rho0=(mp.mpf(9)+S5+2*mp.sqrt(35))/18
def cone_frame(n,z):
    q=1.0-z*z; G=2.0/q; g0=(1.0-3.0*z*z)/q
    m=np.array([g0]+[G/n]*n); U=np.zeros((n+1,3)); U[0]=[0,0,1.0]; r=math.sqrt(max(q,0.0))
    for j in range(n):
        t=2*math.pi*j/n; U[j+1]=[r*math.cos(t),r*math.sin(t),z]
    return m,U
SIGN={}
def sign_table(k):
    if k not in SIGN:
        S=np.array(list(itertools.product([1.0,-1.0],repeat=k)),dtype=float); SIGN[k]=S[S[:,0]>0]
    return SIGN[k]
def section_F(a,m,U):
    beta=(1-a)/2.0; alpha=(3*a-1)/2.0; t=U[:,2]
    h=np.sqrt(np.maximum(beta+alpha*t*t,0.0)); p=m*h; C=U@U.T; W=np.outer(p,p)*C
    S=sign_table(len(p)); A2=float(np.einsum('si,ij,sj->s',S,W,S).max())
    return 0.5*(float(np.abs(W).sum())+A2)
def Fc(a,z):
    s=z*z; R=1-3*s; q=1-s; E=1-a+(3*a-1)*s
    Phi=a*R*R+2*mp.sqrt(2)*z*R*mp.sqrt(a*E)+(E/5)*((5-S5)*s+1+S5)
    return Phi/(q*q)

hd("(B3) EXACT closed form and exact critical point of the uniaxial 5-cone")
m_,U_=cone_frame(5,0.500440049768212)
dv=abs(section_F(0.594608575865948,m_,U_)-float(Fc(mp.mpf('0.594608575865948'),mp.mpf('0.500440049768212'))))
chk("closed form F(a,z)=Phi/(1-z^2)^2 agrees with brute force at the maximiser",dv<1e-12,"dev %.2e"%dv)
da=lambda a,z: mp.diff(lambda x: Fc(x,z), a); dz=lambda a,z: mp.diff(lambda y: Fc(a,y), z)
sol=mp.findroot(lambda A,Z:(da(A,Z),dz(A,Z)),(mp.mpf('0.594608575842171083'),mp.mpf('0.500440048218085076')),
                solver='mnewton',tol=1e-58,maxsteps=300)
as_,zs_=mp.mpf(sol[0]),mp.mpf(sol[1]); Fs=Fc(as_,zs_)
print("   a* =",mp.nstr(as_,55)); print("   z* =",mp.nstr(zs_,55)); print("   F(a*,z*) =",mp.nstr(Fs,55))
chk("F(a*,z*) = rho_0  (|diff| < 1e-50)",abs(Fs-rho0)<mp.mpf('1e-50'),"|diff|=%.2e"%float(abs(Fs-rho0)))
print("   g_hub = %.30f   g_cone = %.30f"%(float((1-3*zs_**2)/(1-zs_**2)),float(2/(5*(1-zs_**2)))))
found={}
for nm,val in [("a*",as_),("z*",zs_),("g_hub",(1-3*zs_**2)/(1-zs_**2)),("z*^2",zs_**2),("F*",Fs)]:
    got=None
    for d in (2,4,6,8,10,12,16,20,24):
        p=mp.findpoly(val,d)
        if p is not None and abs(mp.polyval(p,val))<mp.mpf('1e-40'): got=(d,p); break
    found[nm]=got
    print("   %-6s : %s"%(nm,("deg %d   %s"%(got[0],got[1])) if got else "no integer polynomial, deg<=24"))
aexact=(3*sp.sqrt(35)+5*sp.sqrt(7)-7-7*sp.sqrt(5))/14
print("   conjectured radical form a* = (3sqrt35+5sqrt7-7-7sqrt5)/14 =",mp.nstr(mp.mpf(sp.N(aexact,60)),40))
chk("a* = (3sqrt(35)+5sqrt(7)-7-7sqrt(5))/14  EXACTLY",abs(mp.mpf(sp.N(aexact,60))-as_)<mp.mpf('1e-50'),
    "|diff|=%.2e"%float(abs(mp.mpf(sp.N(aexact,60))-as_)))
Pr=sp.Poly(81*sp.Symbol('x')**4-162*sp.Symbol('x')**3+49*sp.Symbol('x')**2+32*sp.Symbol('x')+1,sp.Symbol('x'))
chk("rho_0 is a root of 81x^4-162x^3+49x^2+32x+1",abs(mp.polyval([81,-162,49,32,1],rho0))<mp.mpf('1e-60'),
    "|p(rho_0)|=%.2e"%float(abs(mp.polyval([81,-162,49,32,1],rho0))))
Pa=sp.Poly(found["a*"][1],sp.Symbol('x'))
print("   -> a* is a root of",sp.expand(Pa.as_expr()))
print("      factorisation over Q(sqrt5):",sp.factor(Pa.as_expr(),extension=[sp.sqrt(5)]))
print("      factorisation over Q(sqrt5,sqrt7):",sp.factor(Pa.as_expr(),extension=[sp.sqrt(5),sp.sqrt(7)]))
print("      F* polynomial:",found["F*"][1] if found["F*"] else None,"   rho_0 check |p(rho0)|=",
      mp.nstr(abs(mp.polyval(found["F*"][1],rho0)),4) if found["F*"] else "-")
# symbolic elimination with the CORRECTED relation  u(R-K) = sqrt2 z (aR - E)
a,z,u=sp.symbols('a z u',positive=True); s=z**2; R=1-3*s; q=1-s; E=1-a+(3*a-1)*s
M=(5-sp.sqrt(5))*s+1+sp.sqrt(5); K=M/5
Phi=a*R**2+2*sp.sqrt(2)*z*R*u+E*M/5
usub=sp.sqrt(2)*z*(a*R-E)/(R-K)
Poly1=sp.expand(a*E*(5*R-M)**2-50*z**2*(a*R-E)**2)
tot=2*u*(sp.diff(Phi,z)*q+4*z*Phi)+sp.diff(Phi,u)*a*sp.diff(E,z)*q
Poly2=sp.expand(sp.together(sp.simplify(tot.subs(u,usub)*(R-K))).as_numer_denom()[0])
f1=sp.lambdify((a,z),Poly1,'mpmath'); f2=sp.lambdify((a,z),Poly2,'mpmath')
v1=f1(as_,zs_); v2=f2(as_,zs_)
print("   Poly1 = aE(5R-M)^2 - 50z^2(aR-E)^2 ;  |Poly1(a*,z*)| =",mp.nstr(abs(v1),4),
      " deg_a=%d deg_z=%d"%(sp.degree(sp.Poly(Poly1,a)),sp.degree(sp.Poly(Poly1,z))))
print("   Poly2 (dF/dz=0, u eliminated)      ;  |Poly2(a*,z*)| =",mp.nstr(abs(v2),4),
      " deg_a=%d deg_z=%d"%(sp.Poly(sp.expand(Poly2),a).degree(),sp.Poly(sp.expand(Poly2),z).degree()))
chk("Poly1(a*,z*) = 0",abs(v1)<mp.mpf('1e-30'),"%.2e"%float(abs(v1)))
chk("Poly2(a*,z*) = 0",abs(v2)<mp.mpf('1e-25'),"%.2e"%float(abs(v2)))
print("   Poly2 =",str(sp.expand(Poly2))[:300],"...")
# minimal polynomial of z* : resultant over Q(sqrt5) then norm to Q
P1=sp.Poly(Poly1,a,domain='EX'); P2=sp.Poly(sp.expand(Poly2/sp.sqrt(2)),a,domain='EX')
rr=sp.resultant(P1,P2,a).as_expr()
rr=sp.expand(rr); rrc=sp.expand(rr.subs(sp.sqrt(5),-sp.sqrt(5)))
Rz=sp.expand(sp.together(rr*rrc).as_numer_denom()[0])
Rz=sp.expand(sp.Poly(Rz,z).primitive()[1].as_expr())
print("   norm-to-Q resultant in z : degree",sp.degree(sp.Poly(Rz,z)))
lam=sp.lambdify(z,Rz,'mpmath'); vr=lam(zs_)
chk("the Q-resultant in z vanishes at z*",abs(vr)<mp.mpf('1e-8'),"val %.3e"%float(abs(vr)))
fl=sp.factor_list(Rz)
for f,e in fl[1]:
    dg=sp.degree(sp.Poly(f,z))
    lf=sp.lambdify(z,f,'mpmath')
    try: vv=abs(lf(zs_))
    except Exception: vv=mp.mpf('nan')
    print("      factor ^%d deg_z=%d |f(z*)|=%s"%(e,dg,mp.nstr(vv,4)))
    if vv==vv and vv<mp.mpf('1e-25') and dg<=12:
        print("      ==> minimal polynomial of z* over Q:")
        print("         ",sp.expand(f))
sys.stdout.flush()

hd("(C) HARMONIC-4 DECOMPOSITION  m4 = 1/5 + ||H4||_F^2/9,  H4 = A4 - 3I4,  A4 = sum m u^(x4)")
I4=(np.einsum('ij,kl->ijkl',np.eye(3),np.eye(3))+np.einsum('ik,jl->ijkl',np.eye(3),np.eye(3))
    +np.einsum('il,jk->ijkl',np.eye(3),np.eye(3)))/15.0
rng=np.random.default_rng(7); w1=0.0; w2=0.0
for k in range(300):
    N=int(rng.integers(4,13)); Q,_=np.linalg.qr(rng.normal(size=(N,3))); g=(Q*Q).sum(1)
    if g.min()<1e-8: continue
    U=Q/np.sqrt(g)[:,None]; C=U@U.T
    m4=float(np.einsum('i,j,ij->',g,g,C**4))/9.0
    A4=np.einsum('n,ni,nj,nk,nl->ijkl',g,U,U,U,U); H4=A4-3*I4
    w1=max(w1,abs(m4-(0.2+(H4*H4).sum()/9.0)))
    w2=max(w2,float(np.abs(np.einsum('iijk->k',A4)-np.ones(3)).max()))
chk("m4 = 1/5 + ||H4||_F^2/9 (300 random Parseval frames)",w1<1e-12,"max dev %.2e"%w1)
chk("tr_12(A4)=I_3 => H4 harmonic (totally traceless)",w2<1e-12,"max dev %.2e"%w2)
chk("||3I4||_F^2 = 9/5",abs(9*(I4*I4).sum()-1.8)<1e-14,"%.2e"%abs(9*(I4*I4).sum()-1.8))

hd("(D) COUPLED INEQUALITY  m4 >= c1(A*^2-c0)")
c1=(20+12*S5)/225; c0=6*rho0-(51*S5-21)/16; c1f,c0f=float(c1),float(c0)
print("   c1 =",mp.nstr(c1,18),"  c0 =",mp.nstr(c0,18),"  threshold A*^2 =",mp.nstr(c0+9*c1/5,18))
print("   at A*^2=3 need m4 >=",mp.nstr(c1*(3-c0),18)," ; equi-latitude floor 19/81 =",mp.nstr(mp.mpf(19)/81,18))
def stats(m,U):
    C=U@U.T; W=np.outer(m,m)*C; S=sign_table(len(m))
    return float(np.abs(W).sum()),float(np.einsum('si,ij,sj->s',S,W,S).max()),float(np.einsum('i,j,ij->',m,m,C**4))/9.0
libs=[]
for n in range(3,9):
    for z in np.linspace(0.02,1/math.sqrt(3)-1e-6,60): libs.append(cone_frame(n,float(z)))
for n in range(3,13):
    t=1/math.sqrt(3); r=math.sqrt(1-t*t)
    U=np.array([[r*math.cos(2*math.pi*j/n),r*math.sin(2*math.pi*j/n),t] for j in range(n)]); libs.append((np.full(n,3.0/n),U))
ph=(1+math.sqrt(5))/2
ic=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1],[0,1/ph,ph],[0,-1/ph,ph]],dtype=float)
ic=ic/np.linalg.norm(ic,axis=1,keepdims=True); libs.append((np.full(6,0.5),ic))
rng=np.random.default_rng(11)
for k in range(800):
    N=int(rng.integers(4,9)); Q,_=np.linalg.qr(rng.normal(size=(N,3))); g=(Q*Q).sum(1)
    if g.min()<1e-6: continue
    libs.append((g,Q/np.sqrt(g)[:,None]))
rows=[]
for m,U in libs:
    iso=abs((m[:,None,None]*(U[:,:,None]*U[:,None,:])).sum(0)-np.eye(3)).max()
    if iso>1e-8: continue
    L,A2,m4=stats(m,U); rows.append((float(6*rho0)-(L+A2), m4-c1f*(A2-c0f), A2,L,m4))
mins=min(rows,key=lambda r:r[1]); mint=min(rows,key=lambda r:r[0])
print("   library size =",len(rows))
print("   WORST coupled slack m4-c1(A*^2-c0) = %+.6e  (A*^2=%.9f L=%.9f m4=%.9f)"%(mins[1],mins[2],mins[3],mins[4]))
print("   WORST true slack 6rho0-(L+A*^2)    = %+.6e  (A*^2=%.9f L=%.9f m4=%.9f coupled=%+.3e)"%(mint[0],mint[2],mint[3],mint[4],mint[1]))
chk("no library configuration violates m4 >= c1(A*^2-c0)",mins[1]>-1e-9,"min slack %.3e"%mins[1])
print("   --- constrained scan: minimise m4 at fixed V=6-2A*^2 (6 atoms) ---")
def unpack(x,n):
    m=np.exp(x[:n]); m=3*m/m.sum(); th=x[n:2*n]; tt=np.clip(np.tanh(x[2*n:3*n]),-0.999,0.999)
    r=np.sqrt(1-tt**2); return m,np.stack([r*np.cos(th),r*np.sin(th),tt],1)
def obj(x,n,Vt,mu):
    m,U=unpack(x,n)
    iso=(m[:,None,None]*(U[:,:,None]*U[:,None,:])).sum(0)-np.eye(3); pen=float((iso**2).sum())
    C=U@U.T; W=np.outer(m,m)*C; S=sign_table(n)
    A2=float(np.einsum('si,ij,sj->s',S,W,S).max()); m4=float(np.einsum('i,j,ij->',m,m,C**4))/9.0
    return m4+1e4*pen+mu*(A2-(3-Vt/2))**2
wg=-1e9; wV=None
for V in [0.0,0.02,0.05,0.10,0.15,0.20,0.25,0.30,0.32889]:
    bb=(1e9,None)
    for sd in range(3):
        r2=np.random.default_rng(200+sd)
        x0=np.concatenate([np.zeros(6)+math.log(0.5),r2.uniform(0,2*math.pi,6),r2.uniform(-0.4,0.9,6)])
        rr=minimize(obj,x0,args=(6,V,500.0),method='Nelder-Mead',options={'maxiter':20000,'maxfev':20000,'xatol':1e-11,'fatol':1e-13})
        if rr.fun<bb[0]: bb=(rr.fun,rr.x)
    m,U=unpack(bb[1],6); L,A2,m4=stats(m,U)
    iso=abs((m[:,None,None]*(U[:,:,None]*U[:,None,:])).sum(0)-np.eye(3)).max()
    need=c1f*((3-V/2)-c0f); gap=m4-need
    if iso<1e-5 and gap>wg: wg,wV=gap,V
    print("      V=%.5f : m4=%.9f need>=%.9f gap=%+.6e  A*^2=%.6f L=%.6f iso=%.1e"%(V,m4,need,gap,A2,L,iso))
chk("constrained V-scan satisfies the coupled inequality",wg>-1e-6,"worst gap %.3e at V=%s"%(wg,wV))
sys.stdout.flush()

hd("(E) EQUI-LATITUDE CIRCLE (A*^2 = 3): maximise L subject to nu^1=nu^2=0")
def lat_cfg(x,n):
    m=np.exp(x[:n]); m=3*m/m.sum(); th=x[n:]
    t=1/math.sqrt(3); r=math.sqrt(1-t*t)
    U=np.stack([r*np.cos(th),r*np.sin(th),np.full(n,t)],1); return m,U
def lat_pen(x,n):
    m,U=lat_cfg(x,n); th=x[n:]
    c1v=abs((m*np.cos(th)).sum())+abs((m*np.sin(th)).sum())
    c2v=abs((m*np.cos(2*th)).sum())+abs((m*np.sin(2*th)).sum())
    return c1v+c2v
def negL(x,n):
    m,U=lat_cfg(x,n); C=U@U.T
    return -float(np.einsum('i,j,ij->',m,m,np.abs(C)))+1e4*lat_pen(x,n)
best=(-1,None,None)
for n in range(3,11):
    bn=(-1,None)
    for sd in range(30):
        r2=np.random.default_rng(500+sd+n)
        x0=np.concatenate([np.zeros(n)+math.log(3.0/n),r2.uniform(0,2*math.pi,n)])
        rr=minimize(negL,x0,args=(n,),method='L-BFGS-B',options={'maxiter':20000,'maxfun':40000,'ftol':1e-16,'gtol':1e-14})
        rr=minimize(negL,rr.x,args=(n,),method='Nelder-Mead',options={'maxiter':40000,'maxfev':40000,'xatol':1e-13,'fatol':1e-15})
        m,U=lat_cfg(rr.x,n); pen=lat_pen(rr.x,n)
        L=float(np.einsum('i,j,ij->',m,m,np.abs(U@U.T)))
        if pen<1e-8 and L>bn[0]: bn=(L,rr.x)
    if bn[1] is not None:
        m,U=lat_cfg(bn[1],n); print("   n=%2d  L_max = %.12f   pen=%.1e   T+=(L+3)/2=%.9f  kappa3^c=%.9f"%(n,bn[0],lat_pen(bn[1],n),(bn[0]+3)/2,(bn[0]+3)/6))
        if bn[0]>best[0]: best=(bn[0],n,bn[1])
print("   BEST equi-latitude L = %.12f  (n=%d)"%(best[0],best[1]))
m,U=lat_cfg(best[2],best[1])
print("      masses =",np.array2string(m,precision=8)); print("      angles/2pi =",np.array2string(best[2][best[1]:]/(2*math.pi),precision=8))
need=6*float(rho0)-3
print("   required for the centre bound: L <= 6rho_0-3 = %.12f ; margin = %+.9f"%(need,need-best[0]))
chk("equi-latitude L_max < 6rho_0-3 (the A*^2=3 end of the coupled inequality closes)",best[0]<need,
    "margin %.6f"%(need-best[0]))
print("\nSUMMARY PASS=%d FAIL=%d"%(PASS,FAIL))
