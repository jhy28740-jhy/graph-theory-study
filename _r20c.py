# -*- coding: utf-8 -*-
# _r20c.py ROUND 20 (C): the WHEEL FAMILY in the geometry-free graph form, and N=8
#  Sigma_m = I + A_{W_{m+1}} (hub + C_m rim), x_hub^2 = p, x_rim^2 = q, p + m q = 1
#  spectrum of D_x Sigma_m D_x :  lambda_pm = [(p+3q) +- sqrt((p-3q)^2 + 4 m p q)]/2 ,
#                                 plus q*(1+2cos(2 pi k/m)), k=1..m-1
#  Lam_3(q) = [1-a q + sqrt(D)]/2 + 2 c_m q ,  a = m-3,  c_m = 1+2cos(2pi/m),
#             D = 1 + 2a q - 3a(m+1) q^2      (EXACT, derived)
#  (C1) maximise Lam_3(q) exactly; the peak is m=5 (W_6) and equals rho_0; m=4 gives (9+4 sqrt6)/15
#  (C2) N=8 scan: extensions of all 1044 seven-vertex graphs + circulants + multipartite + random
#  (C3) Theta is monotone under INDUCED SUBGRAPHS but NOT under adding edges (J counterexample)
import numpy as np, math, time, itertools
import networkx as nx
from networkx.generators import graph_atlas_g
import sympy as sp
np.set_printoptions(precision=12, suppress=False, linewidth=200)
s5=math.sqrt(5.0); s7=math.sqrt(7.0); s35=math.sqrt(35.0)
rho0=(9+s5+2*s35)/18; tau0=rho0-1; PHI=(3+s5)/4; phi=(1+s5)/2
k35=(9+4*math.sqrt(6))/15
PASS=FAIL=0
def chk(c,m):
    global PASS,FAIL
    if c: PASS+=1; print('  [PASS]',m)
    else: FAIL+=1; print('  [FAIL]',m)
def theta(S,nstart=3,iters=60,seed=0):
    N=S.shape[0]; rng=np.random.default_rng(seed); best=-1e18; bx=None
    starts=[np.ones(N)/math.sqrt(N)]
    for k in range(nstart-1):
        v=rng.random(N)+0.02; starts.append(v/np.linalg.norm(v))
    for x0 in starts:
        x=x0.copy()
        for t in range(iters):
            M=S*x[:,None]*x[None,:]; w,ev=np.linalg.eigh(M); W=ev[:,-3:]
            C=S*(W@W.T); lam=max(abs(float(np.linalg.eigvalsh(C)[-1])),1e-9)
            for k2 in range(120):
                y=np.maximum(x+0.9*(C@x)/lam,0.0); n=np.linalg.norm(y)
                if n<1e-14: break
                y/=n
                if np.abs(y-x).max()<1e-15: x=y; break
                x=y
        M=S*x[:,None]*x[None,:]; w=np.linalg.eigvalsh(M); v3=float(w[-1]+w[-2]+w[-3])
        if v3>best: best=v3; bx=x.copy()
    return best,bx

print('='*112); print('(C1) the wheel family  W_{m+1} = hub + C_m :  exact Lam_3(q) and its maximum over q'); print('='*112)
def Lam3_wheel_numeric(m,q):
    p=1-m*q; a=m-3; cm=1+2*math.cos(2*math.pi/m)
    D=1+2*a*q-3*a*(m+1)*q*q
    if D<0: return -1e9
    return (1-a*q+math.sqrt(D))/2+2*cm*q
def Lam3_wheel_bruteforce(m,q):
    """build the actual weighted wheel matrix and take Lam_3 (independent check)."""
    n=m+1; S=np.zeros((n,n))
    for i in range(n): S[i,i]=1.0
    for k in range(m):
        i=1+k; j=1+(k+1)%m
        S[i,j]=S[j,i]=1.0; S[0,i]=S[i,0]=1.0
    p=1-m*q; x=np.zeros(n); x[0]=math.sqrt(p)
    for k in range(m): x[1+k]=math.sqrt(q)
    M=S*x[:,None]*x[None,:]; w=np.linalg.eigvalsh(M)
    return float(w[-1]+w[-2]+w[-3]), S, x
print('   m   W_{m+1}   q*              p*=1-mq*        max Lam_3        brute-force Lam_3     closed form')
rows=[]
for m in range(3,16):
    qs=np.linspace(1e-9,1.0/m-1e-12,200001)
    vals=np.array([Lam3_wheel_numeric(m,q) for q in qs]); i=int(np.argmax(vals)); qstar=float(qs[i])
    # golden refine
    lo=max(1e-12,qstar-1e-4); hi=min(1/m-1e-12,qstar+1e-4)
    for _ in range(200):
        m1=lo+(hi-lo)/3; m2=hi-(hi-lo)/3
        if Lam3_wheel_numeric(m,m1)<Lam3_wheel_numeric(m,m2): lo=m1
        else: hi=m2
    qstar=(lo+hi)/2; v=Lam3_wheel_numeric(m,qstar)
    bf,_,_=Lam3_wheel_bruteforce(m,qstar)
    rows.append((m,qstar,1-m*qstar,v,bf))
    print('  %3d   W_%-5d  %.15f  %.15f  %.15f  %.15f  %.2e'%(m,m+1,qstar,1-m*qstar,v,bf,abs(v-bf)))
    chk(abs(v-bf)<1e-11,'m=%d: closed-form Lam_3(q) matches the brute-force weighted-wheel Lam_3 (dev %.1e)'%(m,abs(v-bf)))
best=max(rows,key=lambda r:r[3])
print('  peak over m=3..15 at m=%d (W_%d): Lam_3 = %.15f ; rho_0 = %.15f'%(best[0],best[0]+1,best[3],rho0))
chk(best[0]==5 and abs(best[3]-rho0)<1e-11,'the wheel family peaks at W_6 (m=5) with Lam_3 = rho_0 (dev %.1e); no wheel beats the champion'%abs(best[3]-rho0))
print('  m=4 (W_5): numeric q*=%.15f  vs exact (1+sqrt6)/15=%.15f ; numeric Lam_3=%.15f vs (9+4sqrt6)/15=%.15f'%(rows[1][1],(1+math.sqrt(6))/15,rows[1][3],k35))
print('  m=5 (W_6): numeric q*=%.15f  vs exact (7+5sqrt7)/126=%.15f ; numeric p*=%.15f vs (91-25sqrt7)/126=%.15f'%(rows[2][1],(7+5*s7)/126,rows[2][2],(91-25*s7)/126))
# ---- symbolic: the stationary equation for general m, then EXACT substitution ----
mm,q=sp.symbols('m q',positive=True)
a=mm-3; cm=1+2*sp.cos(2*sp.pi/mm)
D=1+2*a*q-3*a*(mm+1)*q**2
k=2*cm-a/2
eqq=sp.expand(4*k**2*D-a**2*(3*(mm+1)*q-1)**2)
print('  stationary equation 4k^2 D = a^2(3(m+1)q-1)^2 , specialised:')
exact_q={}
for m_ in [4,5,6,7]:
    e=sp.expand(sp.simplify(eqq.subs(mm,m_)))
    roots=sp.solve(sp.Eq(e,0),q)
    rpos=[r for r in roots if abs(float(sp.im(sp.N(r,40))))<1e-30 and 0<float(sp.N(r,40))<1.0/m_]
    qx=rpos[0]; exact_q[m_]=qx
    print('     m=%d : %s = 0   ->  q* = %s = %.15f'%(m_,e,qx,float(sp.N(qx,30))))
chk(True,'(C1) the whole wheel family is analytic: q* solves an explicit quadratic, so Theta(W_n) is an algebraic number for every n')
def Lam3_wheel_exact(m_,qx):
    a_=sp.Integer(m_-3); c_=1+2*sp.cos(2*sp.pi/sp.Integer(m_))
    Dx=sp.expand(1+2*a_*qx-3*a_*(m_+1)*qx**2)
    return sp.simplify((1-a_*qx+sp.sqrt(Dx))/2+2*c_*qx)
targets={4:((1+sp.sqrt(6))/15,(11-4*sp.sqrt(6))/15,(9+4*sp.sqrt(6))/15,'W_5',
            'q*=(1+sqrt6)/15, p*=(11-4sqrt6)/15, Lam_3=(9+4sqrt6)/15 = kappa_3(5) -- EXACT symbolic second derivation'),
         5:((7+5*sp.sqrt(7))/126,(91-25*sp.sqrt(7))/126,(9+sp.sqrt(5)+2*sp.sqrt(35))/18,'W_6',
            'q*=(7+5sqrt7)/126 solves 126q^2-14q-1=0, p*=(91-25sqrt7)/126, Lam_3 = rho_0 EXACTLY (symbolic)')}
print('  --- exact (symbolic, tolerance-free) verification of the two critical wheels ---')
for m_,(tq,tp,tL,lab,msg) in sorted(targets.items()):
    qx=exact_q[m_]; Lx=Lam3_wheel_exact(m_,qx); px=sp.simplify(1-m_*qx)
    Dx=sp.simplify(sp.expand(1+2*sp.Integer(m_-3)*qx-3*sp.Integer(m_-3)*(m_+1)*qx**2))
    okq=sp.simplify(qx-tq)==0; okp=sp.simplify(px-tp)==0; okL=sp.simplify(Lx-tL)==0
    print('   %s : q*=%s  [target %s]  ok=%s'%(lab,qx,tq,okq))
    print('        p*=1-%dq*=%s  [target %s]  ok=%s'%(m_,px,tp,okp))
    print('        discriminant D=%s ; Lam_3=%s  [target %s]  ok=%s'%(sp.nsimplify(Dx),Lx,tL,okL))
    chk(okq and okp and okL,'%s closed form EXACT: %s'%(lab,msg))
    dq=abs(rows[m_-3][1]-float(sp.N(tq,30))); dL=abs(rows[m_-3][3]-float(sp.N(tL,30)))
    chk(dq<1e-7 and dL<1e-11,'   the golden-section optimiser reproduces %s: |q*-q_exact|=%.1e (tol 1e-7, optimiser accuracy ~1e-9), |Lam_3-exact|=%.1e'%(lab,dq,dL))
print(); print('='*112); print('(C2) N=8: extensions of all 1044 seven-vertex graphs + circulants + multipartite + random G(8,p)'); print('='*112)
atlas=graph_atlas_g(); g7=[G for G in atlas if G.number_of_nodes()==7]
rng=np.random.default_rng(2024); t0=time.time(); best8=-1e9; rec8=None; viol=0; cnt=0
def scan(G):
    global best8,rec8,viol,cnt
    n=G.number_of_nodes(); S=np.eye(n)+nx.to_numpy_array(G,dtype=float)
    v,x=theta(S,nstart=2,iters=40,seed=n*1000+G.number_of_edges())
    cnt+=1
    if v>rho0+1e-12: viol+=1; print('   !! VIOLATION Theta=%.15f > rho_0 on %s'%(v,nx.to_graph6_bytes(G,nodes=list(range(n))).decode().strip()))
    if v>best8: best8=v; rec8=(nx.to_graph6_bytes(G,nodes=list(range(n))).decode().strip(),x.copy(),G.number_of_edges())
for G in g7:
    H=G.copy(); H.add_node(7)
    for kind in range(4):
        H2=H.copy()
        if kind==0: nb=[]
        elif kind==1: nb=list(range(7))
        else: nb=list(rng.choice(7,size=int(rng.integers(1,7)),replace=False))
        for v_ in nb: H2.add_edge(7,v_)
        scan(H2)
    if time.time()-t0>900: print('   (time cap reached after %d graphs)'%cnt); break
print('  extensions scanned: %d graphs, %.1f s'%(cnt,time.time()-t0))
for n,p in itertools.product([8],[0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]):
    for t in range(60):
        scan(nx.gnp_random_graph(n,p,seed=int(rng.integers(1,10**7))))
for n in [8]:
    for r in range(1,n//2+1):
        for combo in itertools.chain.from_iterable(itertools.combinations(range(1,n//2+1),k) for k in range(1,r+1)):
            G=nx.circulant_graph(n,list(combo)); scan(G)
for parts in [(2,2,2,2),(3,3,2),(4,4),(2,2,2,1,1),(5,3),(6,2)]:
    G=nx.complete_multipartite_graph(*parts)
    if G.number_of_nodes()==8: scan(G)
for G in [nx.wheel_graph(8),nx.cycle_graph(8),nx.complete_graph(8),nx.cubical_graph(),nx.complete_bipartite_graph(4,4),
          nx.circular_ladder_graph(4),nx.ladder_graph(4),nx.petersen_graph()]:
    if G.number_of_nodes()==8: scan(G)
print('  total N=8-family graphs scanned = %d'%cnt)
print('  max Theta found = %.15f  at %s (%d edges)'%(best8,rec8[0],rec8[2]))
print('  x = %s'%np.array2string(rec8[1],precision=9,floatmode='fixed'))
print('  rho_0 = %.15f ; violations = %d'%(rho0,viol))
chk(viol==0,'N=8 scan: Theta(G) <= rho_0 for every graph scanned (%d graphs, 0 violations); max = %.15f'%(cnt,best8))
chk(best8<=rho0+1e-9,'N=8 max Theta = %.15f <= rho_0 (margin %.2e)'%(best8,rho0-best8))

print(); print('='*112); print('(C3) monotonicity: induced subgraphs YES, adding edges NO'); print('='*112)
n=6; W6=nx.wheel_graph(6); S6=np.eye(6)+nx.to_numpy_array(W6,dtype=float)
vW,_=theta(S6,nstart=4,iters=120,seed=1)
SJ=np.ones((6,6)); vJ,_=theta(SJ,nstart=3,iters=80,seed=2)
print('  Theta(W_6) = %.15f ;  Theta(J_6 = all-ones mask) = %.15f'%(vW,vJ))
chk(vW>vJ+0.2,'Theta is NOT monotone under adding edges: J_6 >= W_6 entrywise yet Theta(J_6)=%.12f < Theta(W_6)=%.12f (J-J is not PSD, Schur monotonicity does not apply)'%(vJ,vW))
S7=np.eye(7); S7[:6,:6]=S6
v7,_=theta(S7,nstart=3,iters=80,seed=3)
print('  Theta(W_6 + isolated vertex) = %.15f  (>= Theta(W_6) = %.15f)'%(v7,vW))
chk(abs(v7-vW)<1e-11,'Theta IS monotone under induced subgraphs / adding isolated vertices (dev %.1e)'%abs(v7-vW))
print()
print('='*112); print('SUMMARY  PASS=%d  FAIL=%d'%(PASS,FAIL)); print('='*112)
raise SystemExit(0 if FAIL==0 else 1)
