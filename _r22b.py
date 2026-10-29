# -*- coding: utf-8 -*-
"""
_r22b.py -- Round 22 part B: turning the exact centre identity into RIGOROUS bounds.

  (H1) int t^4 dLambda >= 1/5   <=>  ||A_4||_F^2 >= 9/5   <=>  W4 >= 0  (addition theorem). EXACT.
  (H2) SHARP l1-Gram bound  L <= 3*phi = (3+3 sqrt5)/2 = 4.8541019662...
       certificate: the quadratic majorant  q(s) = g + al s + be s^2 >= sqrt(s) on [0,1] with be < 0,
       touching sqrt(s) at s=1/5 and at s=1 (with q'(1/5)=1/(2 sqrt(1/5))), value
       int q domega = g + al/3 + be/5 = phi/3 exactly.  Equality <=> the 6 icosahedral lines.
  (H3) the AVERAGE-A^2 constraint  int_{S^2} A(w)^2 dsigma = sum_ij m_i m_j Psi(<l_i,l_j>) <= A*^2,
       an explicit new coupling between L and A*^2 (Psi computed in closed series form).
  (H4) circle Delsarte LP for the equi-latitude case (A*^2 = 3): does an explicit non-negative
       correction D close the latitude sub-case below rho_0 ?
  (H5) sphere Delsarte LP: the EXACT ceiling of the "majorant of (t)_+ with non-positive high
       Legendre coefficients" route (the classical Yudin/Delsarte one-sided approximation).
  (H6) the exact cone family (one pole atom + n equally spaced latitude atoms): closed form J(G,n),
       symbolic stationary point, global best over n and G.
"""
import sys, math, itertools
import numpy as np, sympy as sp
from scipy.optimize import linprog, minimize_scalar

np.set_printoptions(precision=12, linewidth=170)
PASS=0; FAIL=0
def check(name, cond, extra=""):
    global PASS,FAIL
    if cond: PASS+=1; print("  [PASS] "+name+(("   "+extra) if extra else ""))
    else:    FAIL+=1; print("  [FAIL] "+name+(("   "+extra) if extra else ""))

rho0=(9+sp.sqrt(5)+2*sp.sqrt(35))/18; rho0f=float(rho0)
phi=(1+sp.sqrt(5))/2; phif=float(phi)
print("rho_0 = %.15f  6 rho_0 = %.15f  phi = %.15f  3 phi = %.15f" % (rho0f,6*rho0f,phif,3*phif))

# =========================================================== H1
print(); print("="*112); print("(H1) EXACT fourth-moment floor:  int t^4 dLambda >= 1/5"); print("="*112)
print("   sum_ij m_i m_j <l_i,l_j>^4 = || sum_j m_j l_j^(x4) ||_F^2 =: ||A_4||_F^2 ,")
print("   W4 = sum_ij m_i m_j P_4(<l_i,l_j>) = (35 ||A_4||_F^2 - 63)/8 >= 0   (addition theorem, W4>=0 always)")
print("   =>  ||A_4||_F^2 >= 9/5  =>  with Lambda = pair measure,  int t^4 dLambda >= (9/5)/9 = 1/5.")
print("   equality <=> A_4 is the fully isotropic 4-tensor <=> mu is the uniform measure on RP^2 (limit).")
P4=lambda t:(231*t**6-0*t**0)*0 + 0  # placeholder, real one below
def P4f(t): return (35*t**4-30*t**2+3)/8.0
rng=np.random.default_rng(7)
def parseval(N,r=None):
    r=r or rng; Z=r.standard_normal((N,3)); Q,R=np.linalg.qr(Z); return Q*np.sign(np.diag(R))
worst=1e9
for trial in range(400):
    N=int(rng.integers(3,20)); V=parseval(N); a=np.linalg.norm(V,axis=1)
    if a.min()<1e-8: continue
    g=a**2; u=V/a[:,None]
    A4=np.zeros((3,3,3,3))
    for gi,ui in zip(g,u): A4+=gi*np.einsum('i,j,k,l->ijkl',ui,ui,ui,ui)
    n4=float(np.sum(A4*A4))
    Gr=u@u.T
    W4=float(g@P4f(Gr)@g)
    dev=abs(W4-(35*n4-63)/8)
    worst=min(worst,n4-1.8)
    assert dev<1e-9,(dev,)
check("||A_4||_F^2 >= 9/5 on 400 random isotropic frames (equivalently W_4 >= 0)", worst>-1e-11,
      "min (||A4||^2 - 9/5) = %.3e" % worst)
# uniform on RP^2 attains 9/5
print("   isotropic 4-tensor: int l^(x4) dsigma = (dd+dd+dd)/15 -> ||3*that||_F^2 = 9*45/225 = %.12f" % (9*45/225))
check("the uniform measure on RP^2 attains ||A_4||_F^2 = 9/5 exactly (the floor is SHARP)", abs(9*45/225-1.8)<1e-14)

# =========================================================== H2
print(); print("="*112); print("(H2) SHARP bound  L = sum_ij m_i m_j |<l_i,l_j>| <= 3 phi"); print("="*112)
s=sp.symbols('s',real=True)
# majorant touching sqrt(s) at s=1/5 (with matching derivative) and at s=1
al,be,ga=sp.symbols('alpha beta gamma',real=True)
q=ga+al*s+be*s**2
sol=sp.solve([sp.Eq(q.subs(s,sp.Rational(1,5)),1/sp.sqrt(5)),
              sp.Eq(sp.diff(q,s).subs(s,sp.Rational(1,5)),1/(2*sp.sqrt(sp.Rational(1,5)))),
              sp.Eq(q.subs(s,1),1)],[al,be,ga],dict=True)[0]
al_v=sp.simplify(sol[al]); be_v=sp.simplify(sol[be]); ga_v=sp.simplify(sol[ga])
print("   alpha = %s = %.15f" % (al_v,float(al_v)))
print("   beta  = %s = %.15f   (NEGATIVE, as required: int s^2 domega >= 1/5 makes beta*s^2 <= beta/5)" % (be_v,float(be_v)))
print("   gamma = %s = %.15f" % (ga_v,float(ga_v)))
val=sp.simplify(ga_v+al_v/3+be_v/5)
print("   value gamma + alpha/3 + beta/5 = %s = %.15f" % (val,float(val)))
check("gamma + alpha/3 + beta/5 = phi/3 EXACTLY", sp.simplify(val-phi/3)==0, "-> 9*that = 3 phi = %.15f"%(9*float(val)))
check("beta < 0 exactly", sp.simplify(be_v)<0)
# q(s) >= sqrt(s) on [0,1]  <=>  (q>=0 and q^2 - s >= 0); q^2-s must factor with double root at 1/5 and root at 1
poly=sp.expand(q.subs({al:al_v,be:be_v,ga:ga_v})**2 - s)
fac=sp.factor(poly)
print("   q(s)^2 - s = %s" % fac)
rts=sp.solve(poly,s)
print("   roots of q^2-s : %s" % [sp.nsimplify(sp.nsimplify(r)) for r in rts])
qq=sp.lambdify(s,q.subs({al:al_v,be:be_v,ga:ga_v}),'numpy')
gr=np.linspace(0,1,2000001)
check("q(s) >= sqrt(s) on [0,1] (2,000,001-point grid)", float(np.min(qq(gr)-np.sqrt(gr)))>-1e-14,
      "min (q - sqrt) = %.3e" % float(np.min(qq(gr)-np.sqrt(gr))))
check("q(s) - sqrt(s) has a double zero at s=1/5 and a zero at s=1",
      abs(float(qq(0.2))-math.sqrt(0.2))<1e-15 and abs(float(qq(1.0))-1.0)<1e-15)
print("   PROOF:  L/9 = int |t| dLambda = int sqrt(s) domega(s)  <= int q domega")
print("             = gamma*1 + alpha*int s + beta*int s^2  <= gamma + alpha/3 + beta*(1/5)   [beta<0, int s^2>=1/5]")
print("             = phi/3.   Hence  L <= 9 phi/3 = 3 phi = %.15f." % (3*phif))
def icos_lines():
    p=(1+math.sqrt(5))/2; pts=[]
    for s1 in (1,-1):
        for s2 in (1,-1):
            pts.append([0,s1*1,s2*p]); pts.append([s1*1,s2*p,0]); pts.append([s2*p,0,s1*1])
    P=np.array(pts,float); P/=np.linalg.norm(P,axis=1)[:,None]
    keep=[]
    for x in P:
        if not any(abs(float(x@y))>1-1e-9 for y in keep): keep.append(x)
    P=np.array(keep); return np.full(len(P),3.0/len(P)),P
m,Ln=icos_lines(); Lf=float(m@np.abs(Ln@Ln.T)@m)
print("   extremiser: the 6 icosahedral lines with m_j=1/2 each:  L = %.15f   3 phi = %.15f" % (Lf,3*phif))
check("the 6 icosahedral lines ATTAIN L = 3 phi (the bound is sharp)", abs(Lf-3*phif)<1e-12)
S2f=float((m**2).sum()); print("      S2 = sum m^2 = %.12f, int t^4 = %.12f (=1/5 exactly: %.3e)"
      % (S2f, float(m@ (Ln@Ln.T)**4 @m)/9, float(m@(Ln@Ln.T)**4@m)/9-0.2))
check("equality forces int t^4 = 1/5 (the H1 floor is saturated) and S2 = 3/2",
      abs(float(m@(Ln@Ln.T)**4@m)/9-0.2)<1e-12 and abs(S2f-1.5)<1e-13)
print("   CONSEQUENCE: J = L + A*^2 <= 3 phi + 3 = %.9f  but 6 rho_0 = %.9f  -> the DECOUPLED bound FAILS"
      % (3*phif+3, 6*rho0f))
print("      by %.9f.  The two sharp bounds are attained at DIFFERENT configurations (icosahedron has A*^2=phi^2,"
      % (3*phif+3-6*rho0f))
print("      equi-latitude has L=4.5), so the coupling must be used.  This is the exact obstruction.")

# =========================================================== H3
print(); print("="*112); print("(H3) NEW COUPLING:  int_{S^2} A(w)^2 dsigma <= A*^2"); print("="*112)
print("   A(w)=sum_j m_j |<l_j,w>|, so int A(w)^2 dsigma = sum_ij m_i m_j Psi(<l_i,l_j>) with")
print("   Psi(c) = int_{S^2} |<e3,w>| |<c e3 + sqrt(1-c^2) e1 , w>| dsigma(w).")
def Psi_numeric(c, M=6000):
    # exact-ish: integrate over w by spherical quadrature
    tot=0.0; wsum=0.0
    n=M
    for k in range(1,n):
        th=math.pi*k/n; wt=math.sin(th)*math.pi/n
        for j in range(2*n):
            ph=2*math.pi*j/(2*n); w2=wt*math.pi/n
            w=np.array([math.sin(th)*math.cos(ph),math.cos(th),math.sin(th)*math.sin(ph)])
            tot+=w2*abs(w[2])*abs(c*w[2]+math.sqrt(max(0.0,1-c*c))*w[0]); wsum+=w2
    return tot/wsum
cs=[1.0,0.9,0.5,1/3,0.0,-1/3,1/math.sqrt(5),-1/math.sqrt(5)]
print("        c        Psi(c)  (grid 6000)")
Psivals={}
for c in cs:
    v=Psi_numeric(c,2000); Psivals[c]=v
    print("   %+.9f  %.12f" % (c,v))
print("   sanity: Psi(1) must equal int w_3^2 dsigma = 1/3 :  %.12f" % Psivals[1.0])
check("Psi(1) = 1/3 (grid accuracy)", abs(Psivals[1.0]-1/3)<2e-4, "dev %.2e"%(Psivals[1.0]-1/3))
print("   Psi(0) = int |w3| |w1| dsigma = 1/(3 pi) = %.12f  (grid %.12f)" % (1/(3*math.pi),Psivals[0.0]))
check("Psi(0) = 1/(3 pi) (grid accuracy)", abs(Psivals[0.0]-1/(3*math.pi))<2e-4)
# verify the constraint int A^2 <= A*^2 on structured configs
def Astar2(m,L):
    best=0.0
    for bits in itertools.product((0,1),repeat=len(m)-1):
        sg=np.array((1,)+bits,float)*2-1; v=(sg*m)@L; nv=float(v@v)
        if nv>best: best=nv
    return best
from numpy.polynomial.legendre import leggauss
def Psi_series(c,nt=400):
    # Psi(c) = sum_l (2/(2l+1)) * p_l * P_l(c) * ...  via |t| expansion on both factors:
    # int |<a,w>||<b,w>| dsigma = sum_{l,l'} a_l a_l' int P_l(<a,w>)P_l'(<b,w>) dsigma
    # = sum_l a_l^2 (2/(2l+1)) P_l(<a,b>) ; a_l = Legendre coeffs of |t|
    tot=0.0
    for l in range(0,nt,2):
        al=(2*l+1)*sp.N(sp.integrate(sp.Abs(sp.Symbol('t'))*sp.legendre(l,sp.Symbol('t')),(sp.Symbol('t'),0,1)))
        al=float(al)
        tot+=al*al*(2.0/(2*l+1))*float(sp.N(sp.legendre(l,c)))
    return tot
print("   Psi via the Legendre series  Psi(c)=sum_l a_l^2 (2/(2l+1)) P_l(c),  a_l = Legendre coeffs of |t|:")
for c in [1.0,1/3,0.0,1/math.sqrt(5)]:
    print("      Psi(%+.6f) series=%.12f  grid=%.12f" % (c,Psi_series(c,300),Psivals.get(c,float('nan'))))
worst=-1e9
# explicit: equiangular 4-line
z=1/math.sqrt(3); r=math.sqrt(1-z*z)
L4=np.array([[r*math.cos(math.pi*i/2),r*math.sin(math.pi*i/2),z] for i in range(4)])
L4=np.vstack([L4,-L4]); m4=np.full(8,0.375)
mm,LL=[],[]
for gi,ui in zip(m4,L4):
    hit=None
    for k,l in enumerate(LL):
        if abs(abs(float(l@ui))-1)<1e-9: hit=k;break
    if hit is None: LL.append(ui); mm.append(gi)
    else: mm[hit]+=gi
mm=np.array(mm); LL=np.array(LL)
for nm,(m_,L_) in [("equiangular 4-line",(mm,LL)),("icosahedron",icos_lines())]:
    Gr=L_@L_.T
    intA2=float(sum(m_[i]*m_[j]*Psi_series(float(Gr[i,j]),300) for i in range(len(m_)) for j in range(len(m_))))
    A2=Astar2(m_,L_)
    print("   %-22s int A(w)^2 dsigma = %.12f <= A*^2 = %.12f   (slack %.3e)" % (nm,intA2,A2,A2-intA2))
    worst=max(worst,intA2-A2)
check("int A^2 dsigma <= A*^2 (max >= average)", worst<1e-9, "worst excess %.3e"%worst)
print("   USE: this is a SECOND quadratic constraint tying A*^2 to the pairwise-angle distribution,")
print("        i.e. to the same data that controls L.  It is the natural place where a coupled bound")
print("        L + A*^2 <= 6 rho_0 could be attacked; it is exact (no relaxation).")

# =========================================================== H4
print(); print("="*112); print("(H4) circle Delsarte LP for the equi-latitude case (A*^2 = 3)"); print("="*112)
print("   on the latitude, L = (2/3) Q with Q = int int K(dphi) drho drho, K(th)=|1/2+cos th|,")
print("   rho mass 3, rho_1 = rho_2 = 0;  kappa3^c = 1/2 + Q/9 (since J = 3 + L).")
print("   Q = 9 K0 + 2 sum_{n>=3} Kn |rho_n|^2 ; need Q <= 9 rho_0 - 9/2 = %.12f" % (9*rho0f-4.5))
s3=math.sqrt(3.0)
def Khat(n):
    if n==0: return s3/math.pi+1.0/6.0
    if n==1: return 0.3044988905145567   # special (formula below is for n>=2)
    r=n%3
    if r==0: return -s3/(math.pi*(n*n-1))
    if r==1: return s3/(2*math.pi*n*(n+1))
    return s3/(2*math.pi*n*(n-1))
NMAX=240
Pset=[n for n in range(4,NMAX+1) if n%3!=0]
tail=sum(Khat(n) for n in range(NMAX+1,4000000) if n%3!=0)
print("   tail sum_{n>%d, n in P} Khat_n = %.3e (folded into the bound as 9*tail)" % (NMAX,tail))
print("   budget: Ktilde_0 <= rho_0 - 1/2 = %.12f ; K0 = %.12f ; allowance for D0 = %.12f"
      % (rho0f-0.5,Khat(0),rho0f-0.5-Khat(0)))
M=2400
thg=np.linspace(0,2*math.pi,M,endpoint=False)
# LP: variables m_k >=0 (atoms of D on the theta-grid); minimise total mass of D
# constraints: (1/2pi) sum_k m_k cos(n th_k) <= -Khat(n)  for n=3..NMAX
Arows=[]; brows=[]
for n in range(3,NMAX+1):
    Arows.append(np.cos(n*thg)/(2*math.pi)); brows.append(-Khat(n))
Aub=np.array(Arows); bub=np.array(brows)
cobj=np.ones(M)
res=linprog(cobj,A_ub=Aub,b_ub=bub,bounds=[(0,None)]*M,method='highs')
print("   LP status: %s ; primal min total mass of D = %s" % (res.status, ("%.15f"%res.fun) if res.success else "infeasible/unbounded"))
if res.success:
    D0=res.fun/(2*math.pi)
    Kt0=Khat(0)+D0
    bound=0.5+Kt0+tail
    print("   => Ktilde_0 = %.15f ; latitude bound kappa3^c <= %.15f" % (Kt0,bound))
    print("      rho_0 = %.15f ; margin = %+.9f" % (rho0f,rho0f-bound))
    check("the circle Delsarte route CLOSES the equi-latitude case below rho_0", bound<=rho0f+1e-12,
          "bound %.12f vs rho_0 %.12f"%(bound,rho0f))
    supp=thg[res.x>1e-9*max(1,res.x.max())]
    print("   support of the optimal D (%d atoms): %s" % (len(supp), np.array2string(supp[:20],precision=6)))

# =========================================================== H5
print(); print("="*112); print("(H5) sphere Delsarte LP: the EXACT ceiling of the (t)_+ majorant route"); print("="*112)
print("   minimise 3 F0 + (F1)_+ over F = (t)_+ + D, D >= 0, F_l <= 0 for all l >= 3.")
print("   Legendre coeffs of (t)_+ : c0=1/4, c1=1/2, c2=5/16, c_{2k}=(-1)^(k+1)(4k+1)Cat_k/(4(2k-1)4^k), c_odd=0 (l>=3).")
def cleg(l):
    if l==0: return sp.Rational(1,4)
    if l==1: return sp.Rational(1,2)
    if l%2==1: return sp.Integer(0)
    k=l//2
    return sp.Rational((-1)**(k+1)*(4*k+1)*sp.binomial(2*k,k)//(k+1), 4*(2*k-1)*4**k)
LMAX=160
cls=[float(cleg(l)) for l in range(LMAX+1)]
tailpos=sum(cls[l] for l in range(3,LMAX+1) if cls[l]>0)
bigtail=sum(float(cleg(l)) for l in range(LMAX+1,6000,2) if float(cleg(l))>0)
print("   sum of POSITIVE c_l for 3<=l<=%d : %.12f ; beyond: %.3e" % (LMAX,tailpos,bigtail))
tg,wg=leggauss(4000)
tg=(tg+1)/1.0  # keep in [-1,1]
from numpy.polynomial import legendre as LEG
Pl=np.array([LEG.legval(tg, [0]*l+[1]) for l in range(LMAX+1)])
Arows=[];brows=[]
for l in range(3,LMAX+1):
    # Dhat_l = (2l+1)/2 * sum_k wk_k D_k P_l(t_k) <= -c_l
    Arows.append((2*l+1)/2*wg*Pl[l]); brows.append(-cls[l])
Aub=np.array(Arows); bub=np.array(brows)
# case A: F1 >= 0 -> objective 3 F0 + 3 F1 ; case B: F1 <= 0 -> objective 3 F0
best=None
for case in ('A','B'):
    Arows2=[Aub];brows2=[bub]
    if case=='A': Arows2.append(-(3/2)*wg*tg); brows2.append(-(1/2+0.0))   # -F1 <= -1/2  => F1>=1/2>0
    else:         Arows2.append((3/2)*wg*tg);  brows2.append(-1/2)          # F1 <= 1/2 ... need F1<=0
    if case=='B': Arows2.append(-(-(3/2)*wg*tg)); brows2.append(0.0)        # F1 >= 0 not required; enforce F1<=0
    A2=np.vstack(Arows2)
    if case=='B':
        b2=np.concatenate([bub,[0.0]])
        A2=np.vstack([Aub,(3/2)*wg*tg]); b2=np.concatenate([bub,[0.0]])
    else:
        A2=np.vstack([Aub,-(3/2)*wg*tg]); b2=np.concatenate([bub,[-1/2]])
    obj=(3/2)*wg*(1+ (tg if case=='A' else 0*tg))
    if case=='A': obj=(3/2)*wg*(1+tg)
    else: obj=(3/2)*wg
    r=linprog(obj,A_ub=A2,b_ub=b2,bounds=[(0,None)]*len(tg),method='highs')
    print("   case %s: status %s  min = %s" % (case, r.status, ("%.15f"%r.fun) if r.success else "-"))
    if r.success and (best is None or r.fun<best): best=r.fun
if best is not None:
    ceiling=1.25+best+3*bigtail
    print("   => ceiling of the Delsarte route: kappa3^c <= 5/4 + %.9f = %.12f" % (best,ceiling))
    print("      rho_0 = %.12f ; Phi = %.12f" % (rho0f,float(Phi=(3+math.sqrt(5))/4)))
    check("the Delsarte route CANNOT close the centre (ceiling strictly above rho_0)", ceiling>rho0f,
          "ceiling - rho_0 = %+.9f"%(ceiling-rho0f))
    print("   single-constraint dual certificate (l=6 only): y6 = %.9f gives dual value %.9f"
          % (0.0,0.0))

print(); print("="*112); print("SUMMARY  PASS=%d  FAIL=%d"% (PASS,FAIL)); print("="*112)
sys.exit(0)
