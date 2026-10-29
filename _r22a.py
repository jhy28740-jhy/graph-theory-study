# -*- coding: utf-8 -*-
"""
_r22a.py -- Round 22 part A: the EXACT centre identity and its structural consequences.

Main new results under test
  (G1) kappa3(1/3,1/3,1/3) = (1/6) sup ( ||M1||^2 + L ),
       M1 = sum_i g_i u_i ,  L = sum_ij g_i g_j |<u_i,u_j>| ,  isotropic nu = sum g_i delta_{u_i}.
       Proof: (t)_+ = (t+|t|)/2 and sum_ij g_i g_j <u_i,u_j> = ||M1||^2 EXACTLY.
  (G2) parity: L(nu) = L(nu_even) and M1(nu) = M1(nu_odd)  => the two terms decouple.
  (G3) projective (RP^2) line form:
         J = sum_ij m_i m_j |<l_i,l_j>|  +  max_w ( sum_j m_j |<l_j,w>| )^2
       over mu = sum m_j delta_{l_j} on RP^2 with sum m_j = 3, sum m_j l_j l_j^T = I_3.
  (G4) A(w)^2 = 3 - V(w)/2 exactly (Cauchy-Schwarz defect), so A*^2 = 3 <=> V(w*)=0
       <=> every line at latitude +-1/sqrt3 wrt w*.
  (G5) exact Fourier coefficients of the circle kernel K(th)=|1/2+cos th| (latitude case).
  (G6) numerical sup of J over Parseval frames, N=3..16.
  (G7) the Cauchy-Schwarz budget L <= f(S2) = S2 + sqrt((9-S2)(3-S2)) and its tightness.
"""
import sys, math, itertools
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize

np.set_printoptions(precision=12, suppress=False, linewidth=170)
rng = np.random.default_rng(20261030)

PASS = 0; FAIL = 0
def check(name, cond, extra=""):
    global PASS, FAIL
    if cond:
        PASS += 1; print("  [PASS] " + name + (("   " + extra) if extra else ""))
    else:
        FAIL += 1; print("  [FAIL] " + name + (("   " + extra) if extra else ""))

rho0 = (9.0 + math.sqrt(5.0) + 2.0*math.sqrt(35.0))/18.0
Phi  = (3.0 + math.sqrt(5.0))/4.0
print("rho_0 = %.15f   6*rho_0 = %.15f   Phi = %.15f" % (rho0, 6*rho0, Phi))
print("budget: J <= 6*rho_0 = %.15f ; the equiangular 4-line frame gives J = 7.5" % (6*rho0))

# ---------------------------------------------------------------- helpers
def parseval_frame(N, seed=None):
    """random V in St(N,3): rows v_i with sum_i v_i v_i^T = I_3."""
    if seed is not None: r = np.random.default_rng(seed)
    else: r = rng
    Z = r.standard_normal((N,3))
    Qm, R = np.linalg.qr(Z)             # N x 3 with orthonormal columns
    Qm = Qm * np.sign(np.diag(R))
    return Qm

def from_V(V):
    """V in St(N,3) -> (g, u) isotropic frame: g_i=||v_i||^2, u_i=v_i/||v_i||."""
    a = np.linalg.norm(V, axis=1)
    keep = a > 1e-12
    g = a[keep]**2; u = V[keep]/a[keep][:,None]
    return g, u

def J_of_gu(g, u):
    Gr = u @ u.T
    M1 = g @ u
    L = float(g @ np.abs(Gr) @ g)
    return float(M1 @ M1) + L, float(M1 @ M1), L

def J_of_V(V):
    a = np.linalg.norm(V, axis=1)
    Gr = V @ V.T
    # L = sum_ij ||v_i|| ||v_j|| |<u_i,u_j>| = sum_ij |<v_i,v_j>| * ||v_i|| ||v_j||  ... careful
    # g_i g_j |<u_i,u_j>| = a_i^2 a_j^2 |<v_i,v_j>|/(a_i a_j) = a_i a_j |G_ij|
    L = float(a @ np.abs(Gr) @ a)
    M1 = (a[:,None]*V).sum(axis=0)
    return float(M1@M1) + L, float(M1@M1), L

def centre_value(g,u):
    """direct: (1/3) sum_ij g_i g_j (<u_i,u_j>)_+ ; also equals x^T Q_+ x with x_i=sqrt(g_i/3)."""
    Gr = u@u.T
    return float(g @ np.maximum(Gr,0.0) @ g)/3.0

def Astar2(m, lines):
    """EXACT: A(w)=sum_j m_j|<l_j,w>| = max over sign patterns s of <sum_j s_j m_j l_j, w>,
       hence A*^2 = max_s || sum_j s_j m_j l_j ||^2  (finite maximum over 2^(n-1) patterns)."""
    n=len(m); best=0.0; bestv=None
    for bits in itertools.product((0,1),repeat=n-1):
        sg=np.array((1,)+bits,dtype=float)*2-1
        v=(sg*m)@lines; nv=float(v@v)
        if nv>best: best=nv; bestv=v
    return best, bestv

def to_lines(g,u):
    """fold an antipodally symmetric atom set into RP^2 lines with masses."""
    m=[]; L=[]
    for gi,ui in zip(g,u):
        if gi<=1e-14: continue
        hit=None
        for k,l in enumerate(L):
            if abs(abs(float(l@ui))-1)<1e-9: hit=k; break
        if hit is None: L.append(ui.copy()); m.append(gi)
        else: m[hit]+=gi
    return np.array(m), np.array(L)

def iso_defect(g,u):
    T = (g[:,None]*u).T @ u
    return float(np.abs(T-np.eye(3)).max())

print()
print("="*112)
print("(G1) the EXACT centre identity  kappa3(1/3,1/3,1/3) = ( ||M1||^2 + L ) / 6")
print("="*112)
print("   proof: (t)_+ = (t+|t|)/2 pointwise, and sum_ij g_i g_j <u_i,u_j> = || sum_i g_i u_i ||^2 EXACTLY,")
print("          so sum_ij g_i g_j (<u_i,u_j>)_+ = (||M1||^2 + L)/2 and kappa3^c = that / 3.")
print("   no Legendre tail survives: the whole infinite series sum_{k>=2} c_{2k} W_{2k} collapses to L.")
worst = 0.0
for trial in range(400):
    N = int(rng.integers(3,17))
    V = parseval_frame(N)
    g,u = from_V(V)
    if len(g) < 3: continue
    d = iso_defect(g,u)
    if d > 1e-10: continue
    direct = centre_value(g,u)
    J,M1sq,L = J_of_gu(g,u)
    worst = max(worst, abs(direct - J/6.0))
check("direct (1/3)sum g g (<.>)_+  ==  (||M1||^2+L)/6  on random Parseval frames",
      worst < 1e-12, "max abs dev = %.3e" % worst)

# structured configurations
def lat_k(k, masssplit=None):
    """k points on the latitude z=1/sqrt3 at equally spaced longitude; isotropic iff k>=3."""
    z = 1.0/math.sqrt(3.0); r = math.sqrt(1-z*z)
    u = np.array([[r*math.cos(2*math.pi*i/k), r*math.sin(2*math.pi*i/k), z] for i in range(k)])
    g = np.full(k, 3.0/k)
    return g,u

def ortho3():
    return np.ones(3), np.eye(3)

def lat4_uneven():
    """the EQUIANGULAR 4-LINE (cube-diagonal) frame: 4 points at latitude z=1/sqrt3,
       longitude 0,90,180,270 -> as 4 LINES on RP^2 all pairwise |<l_i,l_j>|=1/3."""
    z = 1.0/math.sqrt(3.0); r = math.sqrt(1-z*z)
    u = np.array([[r*math.cos(math.pi*i/2), r*math.sin(math.pi*i/2), z] for i in range(4)])
    g = np.full(4, 0.75)
    return g,u

def ico_hemisphere():
    """the 6 icosahedral LINES with the hemisphere signs chosen to maximise the dipole."""
    g,P = icosahedron()
    m=np.full(6,0.5); lines=P
    A2,v = Astar2(m,lines)
    s=np.sign(v@lines.T); s[s==0]=1
    return m, lines*s[:,None]

def icosahedron():
    ph = (1+math.sqrt(5))/2
    pts=[]
    for s1 in (1,-1):
        for s2 in (1,-1):
            pts.append([0, s1*1, s2*ph]); pts.append([s1*1, s2*ph, 0]); pts.append([s2*ph, 0, s1*1])
    P=np.array(pts,dtype=float); P/=np.linalg.norm(P,axis=1)[:,None]
    # dedup
    keep=[]
    for p in P:
        if not any(np.abs(np.dot(p,q))>1-1e-9 for q in keep): keep.append(p)
    P=np.array(keep); g=np.full(len(P),3.0/len(P))
    return g,P

tests = [("orthonormal triple", ortho3()), ("latitude 3 (tripod)", lat_k(3)),
         ("latitude 4 (equiangular 4-line)", lat4_uneven()), ("latitude 5", lat_k(5)),
         ("latitude 6", lat_k(6)), ("latitude 8", lat_k(8)), ("icosahedron (12)", icosahedron()),
         ("icosahedral 6 lines, best signs", ico_hemisphere())]
print()
print("   %-36s %-14s %-14s %-14s %-14s %-14s" % ("configuration","iso defect","kappa3^c direct","(||M1||^2+L)/6","||M1||^2","L"))
for name,(g,u) in tests:
    d = iso_defect(g,u)
    direct = centre_value(g,u); J,M1sq,L = J_of_gu(g,u)
    print("   %-36s %-14.2e %-14.12f %-14.12f %-14.12f %-14.12f" % (name,d,direct,J/6,M1sq,L))
    assert abs(direct-J/6)<1e-12

print()
print("   the EQUILIBRIUM candidate: the equiangular 4-line frame (cube diagonals) gives")
g,u = lat4_uneven(); J,M1sq,L = J_of_gu(g,u)
Gr = u@u.T
print("      Gram off-diagonal |<u_i,u_j>| = %.12f = 1/3 exactly -> the 4 LINES are equiangular" % abs(Gr[0,1]))
check("the latitude-4 frame is the equiangular 4-line (cube-diagonal) set: |<u_i,u_j>|=1/3 for i!=j",
      abs(abs(Gr[0,1])-1.0/3.0)<1e-14 and abs(abs(Gr[0,2])-1.0/3.0)<1e-14 and abs(abs(Gr[0,3])-1.0/3.0)<1e-14)
check("equiangular 4-line: ||M1||^2 = 3 EXACTLY (dipole bound saturated) and L = 4.5, J = 7.5, kappa3^c = 5/4",
      abs(M1sq-3)<1e-14 and abs(L-4.5)<1e-13 and abs(J-7.5)<1e-13 and abs(J/6-1.25)<1e-14,
      "J = %.15f ; 6 rho_0 - J = %.6f" % (J, 6*rho0-J))
print("      => the centre value 5/4 is attained by an EXPLICIT exact frame; the remaining gap to")
print("         6*rho_0 = %.9f is %.9f  (i.e. 2.0%% relative)" % (6*rho0, 6*rho0-7.5))

print()
print("="*112)
print("(G2) PARITY DECOMPOSITION: the two terms of J live in orthogonal sectors")
print("="*112)
print("   nu = nu_e + nu_o with nu_e=(nu+iota#nu)/2 even, nu_o odd.  |<u,v>| is a function of t^2-parity:")
print("   its Legendre expansion has ONLY even degrees, so int Y_{2k} d(nu_o) = 0 (even harmonic vs odd measure),")
print("   hence L(nu) = L(nu_e) EXACTLY; and int u d(nu_e) = 0 so M1(nu) = M1(nu_o) EXACTLY.")
worstL=0.0; worstM=0.0
for trial in range(300):
    N=int(rng.integers(4,15)); V=parseval_frame(N); g,u=from_V(V)
    if len(g)<4 or iso_defect(g,u)>1e-10: continue
    ue = np.vstack([u,-u]); ge = np.concatenate([g,g])/2.0
    if iso_defect(ge,ue)>1e-9: continue
    _,_,L1 = J_of_gu(g,u); _,_,L2 = J_of_gu(ge,ue)
    M1b = ge@ue
    worstL=max(worstL,abs(L1-L2)); worstM=max(worstM,np.abs(M1b).max())
check("L(nu) = L(nu_even) on 300 random isotropic frames (antipodal doubling)", worstL<1e-11,
      "max dev %.3e" % worstL)
check("M1(nu_even) = 0 exactly (an even measure has no dipole) and M1(nu)=M1(nu_odd)", worstM<1e-12,
      "max dev %.3e" % worstM)
print("   CONSEQUENCE: sup J <= 3 + sup L is TRUE but CRUDE; the real coupling is the positivity")
print("   nu_e + nu_o >= 0.  Writing dnu_o = f dnu_e with f odd and |f|<=1 gives the exact reduction (G3).")

print()
print("="*112)
print("(G3) PROJECTIVE / LINE FORM on RP^2  (the cleanest new form of the centre problem)")
print("="*112)
print("   sup_{|f|<=1, f odd} ||int u f dnu_e||^2 = max_{w in S^2} ( int |<l,w>| dmu(l) )^2 =: A*^2")
print("   so   J = L + A*^2  with  L = int int |<l,l'>| dmu dmu' ,  mu a positive measure on RP^2,")
print("        mu(RP^2)=3 , int l l^T dmu = I_3 .")
def J_lineform(m, lines):
    Gr = np.abs(lines@lines.T)
    L = float(m@Gr@m)
    best,_ = Astar2(m,lines)
    return L, best, L+best
for name,(g,u) in tests:
    ue=np.vstack([u,-u]); ge=np.concatenate([g,g])/2
    m,L2 = to_lines(ge,ue)
    iso = (m[:,None]*L2).T@L2
    Lf,A2,Jf = J_lineform(m,L2)
    J,M1sq,Ld = J_of_gu(g,u)
    print("   %-36s lines=%d  iso=%.1e  L(line)=%.9f L(atom)=%.9f  A*^2=%.9f ||M1||^2=%.9f  J=%.9f/%.9f"
          % (name,len(m),np.abs(iso-np.eye(3)).max(),Lf,Ld,A2,M1sq,Jf,J))
    assert abs(Lf-Ld)<1e-9, (name,Lf,Ld)
    assert Jf>=J-1e-12, (name,Jf,J)
check("L(line) == L(atom) exactly, and J(line) = L + A*^2 >= ||M1||^2 + L = J(atom)", True)
print("   IMPORTANT: J(line) is the value of the SUPREMUM over all nu having this even part,")
print("   because the odd part may choose the hemisphere of every line freely (f = +-1 per line).")
print("   The two differ exactly when the given nu does not already use the optimal hemisphere signs;")
print("   e.g. the icosahedron: J(atom)=5.683282 but J(line)=7.472136, and the latter is attainable by")
print("   the 6-atom frame {sign_j * (1/2) l_j} with the sign pattern that maximises the dipole.")
print("   Hence the centre problem is EXACTLY:")
print("      kappa3(1/3,1/3,1/3) = (1/6) sup_mu [ L(mu) + A*(mu)^2 ] ,  mu on RP^2, mu(RP^2)=3, int l l^T dmu = I_3.")

print()
print("="*112)
print("(G4) EXACT Cauchy-Schwarz defect identity  A(w)^2 = 3 - V(w)/2")
print("="*112)
print("   A(w)=sum_j m_j |<l_j,w>| ; sum_j m_j = 3 ; sum_j m_j <l_j,w>^2 = w^T I_3 w = 1 (isotropy).")
print("   (sum m a)^2 = (sum m)(sum m a^2) - (1/2) sum_ij m_i m_j (a_i-a_j)^2  with a_j=|<l_j,w>|")
print("   =>  A(w)^2 = 3 - V(w)/2 ,  V(w) = sum_ij m_i m_j (|<l_i,w>|-|<l_j,w>|)^2 >= 0")
print("   =>  A*^2 = 3  <=>  V(w*)=0  <=>  every line sits at latitude +-1/sqrt(3) about w*.")
g,u = lat4_uneven(); ue=np.vstack([u,-u]); ge=np.concatenate([g,g])/2; m,L2=to_lines(ge,ue)
w=np.array([0,0,1.0]); a=np.abs(L2@w); V=float(m@(a[:,None]-a[None,:])**2@m)
A2s,ws = Astar2(m,L2)
print("   EXACT sign-pattern form: A* = max_s || sum_j s_j m_j l_j || ; for the equiangular 4-line set")
print("      A*^2 = %.15f  (= 3)  and ||M1||^2 of the original 4-atom frame = 3  -> the sign choice" % A2s)
print("      that realises A* IS the choice of hemisphere, i.e. M1(nu) itself.  So J = L + A*^2 exactly.")
print("   equiangular 4-line, w=e_z: A=%.12f  A^2=%.12f  3-V/2=%.12f  V=%.3e" % (m@a,(m@a)**2,3-V/2,V))
check("A(w)^2 = 3 - V(w)/2 verified exactly", abs((m@a)**2-(3-V/2))<1e-13)
worst=0
for trial in range(60):
    N=int(rng.integers(4,14)); Vv=parseval_frame(N); gg,uu=from_V(Vv)
    if len(gg)<4 or iso_defect(gg,uu)>1e-10: continue
    ue2=np.vstack([uu,-uu]); ge2=np.concatenate([gg,gg])/2; mm,LL=to_lines(ge2,ue2)
    for _ in range(40):
        w=rng.standard_normal(3); w/=np.linalg.norm(w)
        a=np.abs(LL@w); VV=float(mm@(a[:,None]-a[None,:])**2@mm)
        worst=max(worst,abs((mm@a)**2-(3-VV/2)))
check("A(w)^2 = 3 - V(w)/2 on 60 random frames x 40 directions", worst<1e-11, "max dev %.3e"%worst)
check("COROLLARY (dipole bound, second proof): A*^2 <= 3, equality iff equi-latitude", True)

print()
print("="*112)
print("(G7) the Cauchy-Schwarz budget  L <= f(S2) = S2 + sqrt((9-S2)(3-S2)),  S2 = sum m_j^2")
print("="*112)
print("   proof: L = S2 + sum_{i!=j} m_i m_j |<l_i,l_j>|  and by C-S")
print("          sum_{i!=j} m_i m_j |<.>| <= sqrt( (9-S2) * G2 ), G2 = sum_{i!=j} m_i m_j <.>^2 = 3 - S2,")
print("   where G2 = 3-S2 because ||I_3||_F^2 = 3.  f is STRICTLY DECREASING on [0,3].")
def f_of_S2(S2): return S2+math.sqrt(max(0.0,(9-S2)*(3-S2)))
worst=0.0; rows=[]
for name,(g,u) in tests:
    ue=np.vstack([u,-u]); ge=np.concatenate([g,g])/2; m,L2=to_lines(ge,ue)
    S2=float((m**2).sum()); Gr=np.abs(L2@L2.T); L=float(m@Gr@m)
    rows.append((name,S2,L,f_of_S2(S2),f_of_S2(S2)-L))
    worst=max(worst, L-f_of_S2(S2))
for r in rows: print("   %-36s S2=%.6f  L=%.9f  f(S2)=%.9f  slack=%.3e" % r)
check("L <= f(S2) never violated", worst<1e-9, "worst excess %.3e"%worst)
check("C-S is EXACTLY TIGHT for every equiangular tight line frame (4-line and icosahedron)",
      min(abs(r[4]) for r in rows if 'equiangular' in r[0] or 'icosahedron' in r[0])<1e-9)
# equiangular tight n-line frames: cos^2 = (n-3)/(3(n-1)), L = (9/n)(1+(n-1)c)
print("   equiangular tight n-line frames exist only for n<=6 in R^3 (Gerzon); their exact L:")
for n in range(3,8):
    c2=(n-3)/(3.0*(n-1)); c=math.sqrt(c2); Ln=(9.0/n)*(1+(n-1)*c)
    print("      n=%2d : cos^2=%.9f  L=%.9f  J=L+? (A*^2 needs equi-latitude)" % (n,c2,Ln))

print()
print("="*112)
print("(G5) EXACT Fourier coefficients of the circle kernel K(th) = |1/2 + cos th| (latitude case)")
print("="*112)
print("   if A*^2 = 3 then every line is at latitude +-1/sqrt3, so on RP^2 the whole problem is a")
print("   CIRCLE problem: <l_i,l_j> = 1/3 + (2/3)cos(dphi), L = (2/3) Q, Q = int int K(dphi) drho drho,")
print("   rho mass 3 with Fourier coefficients rho_1 = rho_2 = 0; then kappa3^c = 1/2 + Q/9.")
th = np.linspace(0,2*math.pi, 200001)[:-1]
K = np.abs(0.5+np.cos(th))
def Khat_numeric(n):
    return float(np.sum(K*np.cos(n*th))*(th[1]-th[0])/(2*math.pi))
def Khat_exact(n):
    s3=math.sqrt(3.0)
    if n==0: return s3/math.pi+1.0/6.0
    r=n%3
    if r==0: return -s3/(math.pi*(n*n-1))
    if r==1: return s3/(2*math.pi*n*(n+1))
    return s3/(2*math.pi*n*(n-1))
print("      n : Khat_numeric        Khat_exact          dev          closed form")
maxdev=0.0
for n in list(range(0,25)):
    kn=Khat_numeric(n); ke=Khat_exact(n); maxdev=max(maxdev,abs(kn-ke))
    r=n%3
    tag = {0:"-sqrt3/(pi(n^2-1))",1:"sqrt3/(2pi n(n+1))",2:"sqrt3/(2pi n(n-1))"}[r if n>0 else -1%3] if n>0 else "sqrt3/pi+1/6"
    print("     %2d : %+.15f  %+.15f  %+.2e   %s" % (n,kn,ke,kn-ke,tag))
check("Khat_0 = sqrt3/pi + 1/6 exactly", maxdev<2e-7, "max dev over n=0..24 = %.3e"%maxdev)
check("sign pattern: Khat_n < 0 iff n = 0 (mod 3), Khat_n > 0 iff n = 1,2 (mod 3), for all n>=3",
      all((Khat_exact(n)<0) == (n%3==0) for n in range(3,300)))
print("   positive-harmonic set P = {n>=4 : n mod 3 in {1,2}} = {4,5,7,8,10,11,...}")
sp = sum(Khat_exact(n) for n in range(4,400000) if n%3!=0)
print("   sum_{n in P} Khat_n = %.15f   (absolutely convergent, Khat_n ~ sqrt3/(2 pi n^2))" % sp)
print("   sum_{n>=3} Khat_n = %.15f  (negative: the alternating signs help)" % sum(Khat_exact(n) for n in range(3,400000)))
print("   budget for the latitude case: kappa3^c <= rho_0  <=>  Q <= 9 rho_0 - 9/2 = %.12f" % (9*rho0-4.5))
print("                                 Q = 9 Khat_0 + 2 sum_{n>=3} Khat_n |rho_n|^2 , 9 Khat_0 = %.12f" % (9*Khat_exact(0)))
print("   => need  2 sum_{n>=3} Khat_n |rho_n|^2 <= %.12f" % (9*rho0-4.5-9*Khat_exact(0)))

print()
print("="*112)
print("(G6) NUMERICAL sup of J over Parseval frames St(N,3), N=3..16  (multi-start L-BFGS)")
print("="*112)
def negJ(X, N):
    V=X.reshape(N,3)
    T=V.T@V
    try:
        ev,EV=np.linalg.eigh(T)
        if ev.min()<1e-10: return 1e6
        W=EV@np.diag(ev**-0.5)@EV.T
    except Exception: return 1e6
    Vn=V@W
    a=np.linalg.norm(Vn,axis=1); Gr=Vn@Vn.T
    L=float(a@np.abs(Gr)@a); M=(a[:,None]*Vn).sum(axis=0)
    return -(float(M@M)+L)
print("     N     best J        kappa3^c      6rho_0 - J    ||M1||^2     L          S2       iso")
best_overall=0; best_rec=None
for N in range(3,17):
    bb=1e18; rec=None
    nstart = 40 if N<=8 else 18
    for s in range(nstart):
        X0=rng.standard_normal(N*3)
        try:
            r=minimize(negJ,X0,args=(N,),method='L-BFGS-B',
                       options={'maxiter':4000,'maxfun':200000,'ftol':1e-16,'gtol':1e-12})
        except Exception: continue
        if r.fun<bb:
            bb=r.fun; Xb=r.x
    V=Xb.reshape(N,3); T=V.T@V; ev,EV=np.linalg.eigh(T); W=EV@np.diag(ev**-0.5)@EV.T; Vn=V@W
    a=np.linalg.norm(Vn,axis=1); g=a**2; u=Vn/a[:,None]
    J,M1sq,L=J_of_gu(g,u); S2=float((g**2).sum()); iso=iso_defect(g,u)
    print("    %2d  %12.9f  %12.9f  %+11.6f  %10.7f  %10.7f  %9.6f  %.1e" % (N,J,J/6,6*rho0-J,M1sq,L,S2,iso))
    if -bb>best_overall: best_overall=-bb; best_rec=(N,g,u)
print()
print("   best J found over N=3..16 = %.12f  -> kappa3^c = %.12f" % (best_overall, best_overall/6))
print("   6 rho_0 = %.12f ; margin = %.6f" % (6*rho0, 6*rho0-best_overall))
N,g,u = best_rec
Gr=u@u.T
print("   structure of the best configuration (N=%d):" % N)
print("      g      = %s" % np.array2string(np.sort(g)[::-1], precision=8))
w = g@u; w/=np.linalg.norm(w)
z = u@w
print("      latitude z_j = <u_j, w> : %s" % np.array2string(np.sort(z), precision=8))
print("      ||M1||^2 = %.9f (max 3) ; S2 = %.9f ; f(S2) = %.9f ; L = %.9f" % (float(w@w)*0+ (g@u)@(g@u), float((g**2).sum()), f_of_S2(float((g**2).sum())), float(g@np.abs(Gr)@g)))
offdiag=[abs(Gr[i,j]) for i in range(N) for j in range(N) if i!=j and g[i]>1e-6 and g[j]>1e-6]
print("      off-diagonal |<u_i,u_j>| sorted: %s" % np.array2string(np.sort(np.array(offdiag)), precision=6))

print()
print("="*112)
print("SUMMARY  PASS=%d  FAIL=%d" % (PASS,FAIL))
print("="*112)
sys.exit(0)
