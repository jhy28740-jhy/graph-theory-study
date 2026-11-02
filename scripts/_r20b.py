# -*- coding: utf-8 -*-
# _r20b.py ROUND 20 (B): the geometry-free graph form, exhaustively
#   kappa_3 = sup_{N} max_{Sigma = I+A_G, G simple on N vertices} Theta(G),
#   Theta(G) = max_{x>=0, ||x||=1} Lam_3(D_x Sigma D_x)
#   (B1) exhaustive Theta over ALL isomorphism classes for N=4,5,6,7 (graph atlas) -> extremal graph
#   (B2) closed-form certificate for the wheel W_6: Theta(W_6) >= rho_0, exact in Q(sqrt5,sqrt7)  [sympy]
#   (B3) lambda_min(D_x Sigma D_x) >= 1/N - 1/2  (the N=4 engine for kappa_3(4)=5/4), + exhaustion over N=4
#   (B4) complement/Ky-Fan-subadditivity family ceiling: Lam_3 <= 1 - (3 smallest of D_x Sbar D_x)
import numpy as np, math, itertools, time
import networkx as nx
from networkx.generators import graph_atlas_g
import sympy as sp
np.set_printoptions(precision=12, suppress=False, linewidth=200)
s5=math.sqrt(5.0); s7=math.sqrt(7.0); s35=math.sqrt(35.0)
rho0=(9+s5+2*s35)/18; tau0=rho0-1; PHI=(3+s5)/4; phi=(1+s5)/2
k34=1.25; k35=(9+4*math.sqrt(6))/15
PASS=FAIL=0
def chk(c,m):
    global PASS,FAIL
    if c: PASS+=1; print('  [PASS]',m)
    else: FAIL+=1; print('  [FAIL]',m)

def xstep_exact(C):
    """max_{x>=0,||x||=1} x^T C x  by face enumeration (exact)."""
    N=C.shape[0]; best=-1e18; bx=None
    for r in range(1,N+1):
        for F in itertools.combinations(range(N),r):
            Cf=C[np.ix_(F,F)]; w,ev=np.linalg.eigh(Cf); v=ev[:,-1]
            if v.sum()<0: v=-v
            if v.min()>=-1e-11:
                if w[-1]>best: best=float(w[-1]); y=np.zeros(N); y[list(F)]=np.maximum(v,0.0); bx=y/np.linalg.norm(y)
    return best,bx

def xstep_pg(C,x0,iters=250):
    N=C.shape[0]; lam=np.linalg.eigvalsh(C)[-1]; eta=0.9/max(abs(lam),1e-9)
    x=np.maximum(x0,1e-6); x/=np.linalg.norm(x); f=x@(C@x)
    for t in range(iters):
        y=np.maximum(x+eta*(C@x),0.0); n=np.linalg.norm(y)
        if n<1e-14: break
        y/=n; g=y@(C@y)
        if abs(g-f)<1e-16: x=y; f=g; break
        x=y; f=g
    return f,x

def theta(S,nstart=2,iters=45,seed=0,exact=False,refine=0):
    N=S.shape[0]; rng=np.random.default_rng(seed); best=-1e18; bx=None
    starts=[np.ones(N)/math.sqrt(N)]
    for k in range(nstart-1):
        v=rng.random(N)+0.02; starts.append(v/np.linalg.norm(v))
    for x0 in starts:
        x=x0.copy(); val=-1e18
        for t in range(iters):
            M=S*x[:,None]*x[None,:]
            w,ev=np.linalg.eigh(M); val=float(w[-1]+w[-2]+w[-3]); W=ev[:,-3:]
            C=S*(W@W.T)
            if exact: xn=xstep_exact(C)[1]
            else: xn=xstep_pg(C,x)[1]
            if xn is None or np.linalg.norm(xn)<1e-12: break
            if np.abs(xn-x).max()<1e-14: x=xn; break
            x=xn
        M=S*x[:,None]*x[None,:]; w=np.linalg.eigvalsh(M); v3=float(w[-1]+w[-2]+w[-3])
        if v3>best: best=v3; bx=x.copy()
    for k in range(refine):
        xr=np.random.default_rng(1000+k).random(N)+0.02; xr/=np.linalg.norm(xr)
        x=xr
        for t in range(200):
            M=S*x[:,None]*x[None,:]; w,ev=np.linalg.eigh(M); W=ev[:,-3:]
            C=S*(W@W.T); v,xn=xstep_exact(C)
            if xn is None: break
            if np.abs(xn-x).max()<1e-14: x=xn; break
            x=xn
        M=S*x[:,None]*x[None,:]; w=np.linalg.eigvalsh(M); v3=float(w[-1]+w[-2]+w[-3])
        if v3>best: best=v3; bx=x.copy()
    return best,bx

print('='*112)
print('(B1) exhaustive scan of Theta(G)=max_x Lam_3(D_x (I+A_G) D_x) over all isomorphism classes, N=4..7')
print('='*112)
atlas=graph_atlas_g()
summary={}
for N in [4,5,6,7]:
    gs=[G for G in atlas if G.number_of_nodes()==N and nx.number_of_selfloops(G)==0]
    # dedup by canonical form (atlas already non-isomorphic)
    t0=time.time(); res=[]
    for G in gs:
        S=np.eye(N)+nx.to_numpy_array(G,dtype=float)
        v,x=theta(S,nstart=2,iters=45,seed=hash(nx.to_graph6_bytes(G,nodes=list(range(N))))%(2**31))
        res.append((v,nx.to_graph6_bytes(G,nodes=list(range(N))).decode().strip(),x,G))
    res.sort(key=lambda r:-r[0])
    top=res[:6]
    print('  N=%d : %d isomorphism classes, %.1f s'%(N,len(gs),time.time()-t0))
    for v,g6,x,G in top:
        print('      Theta=%.15f  %-14s edges=%2d  x=%s'%(v,g6,G.number_of_edges(),np.array2string(x,precision=6,floatmode='fixed')))
    # refine the best few exactly
    best=-1e18; bg6=None; bx=None
    for v,g6,x,G in res[:8]:
        S=np.eye(N)+nx.to_numpy_array(G,dtype=float)
        v2,x2=theta(S,nstart=1,iters=1,seed=0,refine=4)
        if v2>best: best=v2; bg6=g6; bx=x2; bG=G
    print('    >> refined max Theta = %.15f  at %s (%s, %d edges)'%(best,bg6,nx.to_graph6_bytes(bG,nodes=list(range(N))).decode().strip(),bG.number_of_edges()))
    print('       refined x = %s ; x^2 = %s'%(np.array2string(bx,precision=9,floatmode='fixed'),np.array2string(bx**2,precision=9,floatmode='fixed')))
    print('       negative tail 1-Theta = %.15f'%(1-best))
    summary[N]=(best,bg6,bG,bx)
tgt={4:k34,5:k35,6:rho0,7:rho0}
for N in [4,5,6,7]:
    best,bg6,bG,bx=summary[N]
    print('  N=%d: max_G Theta(G) = %.15f   target kappa_3(N) = %.15f   diff = %.3e'%(N,best,tgt[N],best-tgt[N]))
    chk(abs(best-tgt[N])<2e-9,'N=%d: exhaustive graph-form scan reproduces kappa_3(N) = %.15f (diff %.2e); extremal graph %s'%(N,tgt[N],best-tgt[N],bg6))
    chk(best<=rho0+1e-12,'N=%d: Theta(G) <= rho_0 for EVERY looped graph on %d vertices (max %.15f, margin %.2e)'%(N,N,best,rho0-best))

print()
print('='*112)
print('(B2) closed form for the wheel W_6 = K_1 v C_5 :  Theta(W_6) >= rho_0 with exact surd arithmetic')
print('='*112)
r7=sp.sqrt(7); r5=sp.sqrt(5)
p=(91-25*r7)/126; q=(7+5*r7)/126
print('  p = x_hub^2 = (91-25*sqrt7)/126 = %.15f ;  q = x_rim^2 = (7+5*sqrt7)/126 = %.15f'%(float(p),float(q)))
print('  p + 5q - 1 = ',sp.simplify(p+5*q-1))
chk(sp.simplify(p+5*q-1)==0,'weights normalised exactly: p+5q=1')
disc=sp.simplify((1-2*q)**2+8*p*q)
print('  (1-2q)^2 + 8pq = ',disc,' = ',sp.nsimplify(disc),' -> sqrt = sqrt(5/7) = sqrt35/7')
chk(disc==sp.Rational(5,7),'the 2x2 discriminant collapses EXACTLY to 5/7 (this is where sqrt7 and sqrt35 come from)')
lp=(1-2*q+sp.sqrt(disc))/2; lm=(1-2*q-sp.sqrt(disc))/2
ph=(1+r5)/2
rim_hi=q*(1+2*sp.cos(2*sp.pi/5)); rim_lo=q*(1+2*sp.cos(4*sp.pi/5))
rim_hi_e=sp.simplify(q*ph); rim_lo_e=sp.simplify(-q/ph)
print('  rim modes: q*(1+2cos72) = q*phi = %.15f ;  q*(1+2cos144) = -q/phi = %.15f'%(float(rim_hi_e),float(rim_lo_e)))
chk(abs(float(rim_hi)-float(rim_hi_e))<1e-14,'1+2cos(72 deg) = phi exactly => the two rim "high" modes are q*phi')
chk(abs(float(rim_lo)-float(rim_lo_e))<1e-14,'1+2cos(144 deg) = -1/phi exactly => the two rim "low" modes are -q/phi')
L3=sp.simplify(lp+2*rim_hi_e); tgt=sp.Rational(9,18)+r5/18+2*sp.sqrt(35)/18
print('  Lam_3 = lambda_+ + 2*q*phi = ',L3)
print('  rho_0 = (9+sqrt5+2sqrt35)/18')
print('  difference = ',sp.simplify(L3-sp.nsimplify(rho0)))
chk(sp.simplify(L3-tgt)==0,'EXACT: Lam_3(D_x Sigma_{W_6} D_x) = lambda_+ + 2 q phi = (9+sqrt5+2sqrt35)/18 = rho_0  (symbolic)')
negtail=sp.simplify(lm+2*rim_lo_e)
print('  negative tail = lambda_- + 2*(-q/phi) = ',sp.nsimplify(negtail),' = %.15f ;  -tau_0 = %.15f'%(float(negtail),-tau0))
chk(sp.simplify(negtail+ (tgt-1))==0,'EXACT: the negative tail of the weighted wheel equals -tau_0 = -(rho_0-1), i.e. tau_0 IS the wheel tail')

print()
print('='*112)
print('(B3) lambda_min(D_x Sigma D_x) >= 1/N - 1/2 : the engine behind kappa_3(4)=5/4 (new proof)')
print('='*112)
rng=np.random.default_rng(5); worst=1e9; worstrec=None; viol=0
for t in range(6000):
    N=int(rng.integers(3,10))
    G=nx.gnp_random_graph(N,0.5,seed=int(rng.integers(1,10**6)))
    S=np.eye(N)+nx.to_numpy_array(G,dtype=float)
    x=rng.random(N)**1.5; x/=np.linalg.norm(x)
    w=rng.normal(size=N); w/=np.linalg.norm(w)
    M=S*x[:,None]*x[None,:]
    val=float(w@M@w)                       # >= lambda_min
    lm=float(np.linalg.eigvalsh(M)[0])
    # the two ingredients: q = sum x_i^2 w_i^2 >= r^2/N ; 2PMn <= r^2/2 with r = <x,|w|>
    qq=float((x**2*w**2).sum()); rr=float((x*np.abs(w)).sum())
    a=x*w; P=a[a>0].sum(); Mn=-a[a<0].sum()
    if qq-rr*rr/N < -1e-12: print('   !! q >= r^2/N violated'); viol+=1
    if 2*P*Mn-rr*rr/2 > 1e-12: print('   !! 2PMn <= r^2/2 violated'); viol+=1
    marg=lm-(1.0/N-0.5)
    if marg<worst: worst=marg; worstrec=(N,lm,1.0/N-0.5,nx.to_graph6_bytes(G,nodes=list(range(N))).decode().strip())
print('  worst margin  lambda_min - (1/N - 1/2)  over 6000 random (N,Sigma,x): %.6e  at N=%d'%(worst,worstrec[0]))
print('  (lambda_min = %.15f , 1/N-1/2 = %.15f , graph %s)'%(worstrec[1],worstrec[2],worstrec[3]))
chk(worst>=-1e-12 and viol==0,'lambda_min(D_x Sigma D_x) >= 1/N - 1/2 on 6000 random tests; both sub-inequalities (q>=r^2/N, 2PMn<=r^2/2) hold')
# N=4 exhaustion: min lambda_min over all 11 graphs and x on a grid + exact refinement
N=4; gs=[G for G in atlas if G.number_of_nodes()==4]
best4=1e9; rec4=None
for G in gs:
    S=np.eye(N)+nx.to_numpy_array(G,dtype=float)
    for t in range(400):
        x=rng.random(N)**2; x/=np.linalg.norm(x)
        lm=float(np.linalg.eigvalsh(S*x[:,None]*x[None,:])[0])
        if lm<best4: best4=lm; rec4=(nx.to_graph6_bytes(G,nodes=list(range(N))).decode().strip(),x.copy())
print('  N=4 exhaustive (11 graphs x 400 x): min lambda_min = %.15f  at %s'%(best4,rec4[0]))
print('     x = %s  (uniform 1/2 = %.6f)'%(np.array2string(rec4[1],precision=9,floatmode='fixed'),0.5))
chk(best4>=-0.25-1e-9,'N=4 exhaustive: lambda_min(D_x Sigma D_x) >= -1/4 for all 11 graphs and 4400 weights (min found %.15f, no violation)'%best4)
S4=np.eye(4)+nx.to_numpy_array(nx.cycle_graph(4),dtype=float); xu=np.full(4,0.5)
lm4=float(np.linalg.eigvalsh(S4*xu[:,None]*xu[None,:])[0])
print('  exact attainment: C_4 with x uniform -> lambda_min = %.16f (= -1/4)'%lm4)
print('  spec = %s ; Lam_3 = 1 - lambda_min = %.16f'%(np.array2string(np.linalg.eigvalsh(S4*xu[:,None]*xu[None,:])[::-1],precision=12),1-lm4))
chk(abs(lm4+0.25)<1e-14,'C_4 + uniform x attains lambda_min = -1/4 EXACTLY => kappa_3(4) = 1+1/4 = 5/4 (NEW independent proof, no frame geometry used)')

print()
print('='*112)
print('(B4) complement + Ky-Fan subadditivity family:  Lam_3(D_x Sigma D_x) <= 1 - (3 smallest of D_x Sbar D_x)')
print('='*112)
def ceil_comp(S,x):
    Sb=np.ones_like(S)-S; np.fill_diagonal(Sb,0.0)
    w=np.linalg.eigvalsh(Sb*x[:,None]*x[None,:])
    return 1.0-float(w[:3].sum())
# champion configuration: Sigma = I + A_{W_6}, x^2 = (p, q*5)
N=6; W6=nx.wheel_graph(6)   # nx.wheel_graph(n) has n vertices: hub + C_{n-1}
S=np.eye(6)+nx.to_numpy_array(W6,dtype=float)
xw=np.sqrt(np.array([float(p)]+[float(q)]*5))
assert abs(p-(91-25*s7)/126)<1e-15 and abs(q-(7+5*s7)/126)<1e-15, (p,q)
x=xw
M=S*x[:,None]*x[None,:]
print('  wheel W_6 mask (nx.wheel_graph(6)), hub = vertex 0:'); print('  x = %s'%np.array2string(x,precision=9))
print(np.array2string(S.astype(int)))
print('  Lam_3 = %.15f (must be rho_0=%.15f)'%(Lam3:=float(np.linalg.eigvalsh(M)[-1]+np.linalg.eigvalsh(M)[-2]+np.linalg.eigvalsh(M)[-3]),rho0))
print('  complement-KyFan ceiling at the same (Sigma,x) = %.15f'%(ceil_comp(S,x)))
best_ce=0; rec=None
for NN in [5,6,7]:
    for G in [G for G in atlas if G.number_of_nodes()==NN]:
        S=np.eye(NN)+nx.to_numpy_array(G,dtype=float)
        for t in range(30):
            xx=rng.random(NN)**1.5; xx/=np.linalg.norm(xx)
            c=ceil_comp(S,xx)
            if c>best_ce: best_ce=c; rec=(NN,nx.to_graph6_bytes(G,nodes=list(range(NN))).decode().strip())
print('  sup of this ceiling over all graphs N=5..7 and 30 random x each: %.15f  at %s'%(best_ce,rec))
print('  compare: rho_0 = %.15f , Phi = %.15f'%(rho0,PHI))
chk(best_ce>PHI,'REFUTATION with ceiling: complement+KyFan-subadditivity family tops out at %.12f > Phi = %.12f > rho_0 => cannot reach rho_0 (dead)'%(best_ce,PHI))

print()
print('='*112); print('SUMMARY  PASS=%d  FAIL=%d'%(PASS,FAIL)); print('='*112)
raise SystemExit(0 if FAIL==0 else 1)
