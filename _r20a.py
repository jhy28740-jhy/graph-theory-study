# -*- coding: utf-8 -*-
# _r20a.py  ROUND 20 (A): the geometry-free Ky Fan form
#   kappa_3 = sup_{N, Sigma in {0,1}^{N x N} sym diag-1, x>=0, ||x||=1}  Lam_3(D_x Sigma D_x)
#   (self-consistency in thm:kappa(ii) is REDUNDANT -- new proof: Ky Fan duality + convexity in Sigma)
#   (A1) champion ALIGNMENT test: is Lam_3(D_x Sigma D_x) == rho_0 exactly, or strictly larger?
#        (strictly larger  =>  kappa_3 > rho_0  =>  the champion is NOT extremal; huge)
#   (A2) the mask/KyFan fixed-point iteration Q -> Pi_3(D_x Sigma_Q D_x) from random starts
#   (A3) the complementary-dimension dual form: E(x,W) = sum x_i^2||w_i||^2 - sum_{i!=j} x_i x_j (-<w_i,w_j>)_+
#   (A4) Gram-energy budget  sum_{i!=j} <w_i,w_j>^2 <= 3d/(d+3)  for d+3 vectors with sum w_i w_i^T = I_d
#   (A5) B_1 <= 3 (Cauchy-Schwarz) and B_2 == 0 (isotropy) re-check on random isotropic frames
import numpy as np, math
np.set_printoptions(precision=12, suppress=False, linewidth=200)
s5=math.sqrt(5.0); s7=math.sqrt(7.0); s35=math.sqrt(35.0)
rho0=(9+s5+2*s35)/18; tau0=rho0-1; PHI=(3+s5)/4
PASS=FAIL=0
def chk(c,m):
    global PASS,FAIL
    if c: PASS+=1; print('  [PASS]',m)
    else: FAIL+=1; print('  [FAIL]',m)

def champion():
    Vhh=14.0/(-7*s5-7+5*s7+3*s35); Vhc=math.sqrt(7.0/5.0); Vcc=-3.5-7*s5/10+3*s35/10+5*s7/2
    Ah=math.sqrt(Vhh); Ac=Vhc/Ah; Rc=math.sqrt(Vcc-Ac**2)
    uh=(91-25*s7)/126; uc=(7+5*s7)/126; wv=np.array([uh]+[uc]*5)
    cols=[np.array([Ah,0.0,0.0])]
    for j in range(5):
        th=2*math.pi*j/5; cols.append(np.array([Ac,Rc*math.cos(th),Rc*math.sin(th)]))
    A0=np.array([math.sqrt(wv[k])*cols[k][rr] for rr in range(3) for k in range(6)]).reshape(3,6)
    P=A0.T@np.linalg.inv(A0@A0.T)@A0
    w,Vs=np.linalg.eigh(P); idx=np.argsort(-w); Vs=Vs[:,idx]; w=w[idx]
    V=Vs[:,:3]*np.sqrt(np.maximum(w[:3],0))
    g=np.linalg.norm(V,axis=1)**2; U=V/np.sqrt(g)[:,None]
    return U,g

def perron_of_Qp(Q):
    Qp=np.maximum(Q,0.0); w,ev=np.linalg.eigh(Qp); k=int(np.argmax(w)); rho=float(w[k]); x=ev[:,k].copy()
    if x.sum()<0: x=-x
    if np.abs(x).min()<1e-14: pass
    return rho,x,Qp

def Lam3(M):
    w=np.linalg.eigvalsh(M); return float(w[-1]+w[-2]+w[-3])

def Pi3(M):
    w,ev=np.linalg.eigh(M); return ev[:,-3:]@ev[:,-3:].T

print('='*108)
print('(A1) CHAMPION ALIGNMENT:  Lam_3(D_x Sigma D_x)  vs  rho(Q_+) = rho_0')
print('='*108)
Uc,gc=champion()
Vc=Uc*np.sqrt(gc)[:,None]; Qc=Vc@Vc.T
rho_c,x_c,Qp_c=perron_of_Qp(Qc)
Sig_c=(Qc>=0).astype(float)
print('  champion: rho(Q_+)           = %.16f'%rho_c)
print('            rho_0              = %.16f   (diff %.3e)'%(rho0,rho_c-rho0))
print('  x (Perron) =',np.array2string(x_c,precision=12))
print('  Sigma (sign mask of Q, loops on) =')
print(np.array2string(Sig_c.astype(int)))
M_c=Sig_c*x_c[:,None]*x_c[None,:]
wM=np.linalg.eigvalsh(M_c)
print('  spec(D_x Sigma D_x) =',np.array2string(wM[::-1],precision=12))
L3=Lam3(M_c)
print('  Lam_3(D_x Sigma D_x) = %.16f'%(L3))
print('  Lam_3 - rho_0        = %.6e'%(L3-rho0))
print('  Lam_3 - rho(Q_+)     = %.6e'%(L3-rho_c))
print('  tr(M) = %.15f (must be 1)   negative tail = %.15f'%(np.trace(M_c),float(wM[:-3].sum())))
Qpi=Pi3(M_c)
d_sub=np.abs(Qpi-Qc).max()
print('  || Pi_3(D_x Sigma D_x) - Q_champion ||_max = %.6e'%d_sub)
chk(abs(L3-rho_c)<1e-12,'ALIGNMENT: Lam_3(D_x Sigma D_x) == rho(Q_+) at the champion (gap %.2e) => the champion frame IS the Ky Fan maximizer for its own mask'%abs(L3-rho_c))
chk(d_sub<1e-11,'Pi_3(D_x Sigma D_x) == Q_champion (dev %.2e): the report thm:kappa(ii) self-consistency claim holds at the champion'%d_sub)
chk(abs(L3-rho0)<1e-12,'therefore Lam_3(D_x Sigma D_x) = rho_0 exactly: NO kappa_3 > rho_0 witness from the champion (gap %.2e)'%(L3-rho0))

print()
print('='*108)
print('(A2) mask / KyFan fixed-point iteration  Q -> Pi_3(D_x Sigma_Q D_x)   (a NEW, non-gradient ascent)')
print('='*108)
def mask_iter(N,seed,iters=400,verbose=False):
    rng=np.random.default_rng(seed)
    A=rng.normal(size=(N,3)); V,_=np.linalg.qr(A); V=V.T  # N x 3 random Stiefel
    Q=V@V.T; best=rho_c*0
    hist=[]
    for t in range(iters):
        rho,x,_=perron_of_Qp(Q)
        Sig=(Q>=0).astype(float)
        M=Sig*x[:,None]*x[None,:]
        L3=Lam3(M)
        hist.append((rho,L3))
        Qn=Pi3(M)
        if np.abs(Qn-Q).max()<1e-15: Q=Qn; break
        Q=Qn
    rho,x,_=perron_of_Qp(Q)
    return rho,Q,hist
print('   N  seed   rho(Q_+) after iteration     max Lam_3 seen      >= rho_0 ?')
worst=0.0; worstinfo=None
for N in [5,6,7,8,9,10,12]:
    for seed in range(6):
        rho,Q,hist=mask_iter(N,seed*97+N)
        mx=max(h[1] for h in hist)
        if mx>worst: worst=mx; worstinfo=(N,seed,rho,mx)
        flag='*** ABOVE rho_0 ***' if mx>rho0+1e-12 else ''
        if seed<2 or mx>rho0-1e-9:
            print('  %3d %4d   %.15f   %.15f   %s'%(N,seed,rho,mx,flag))
print('  global max Lam_3 seen over all runs = %.16f  (rho_0 = %.16f)'%(worst,rho0))
chk(worst<=rho0+1e-12,'mask-iteration never exceeded rho_0 (max %.15f, excess %.2e) => no counterexample to kappa_3<=rho_0 found by this new ascent'%(worst,worst-rho0))

print()
print('='*108)
print('(A3) complementary-dimension dual form at the champion:  E(x,W) >= -tau_0 with equality?')
print('     W = any orthonormal (N-3)-frame; E = sum_i x_i^2 ||w_i||^2 - sum_{i!=j} x_i x_j (-<w_i,w_j>)_+')
print('='*108)
def Efun(x,W,Sig=None):
    G=W@W.T; neg=np.maximum(-G,0.0); np.fill_diagonal(neg,0.0)
    if Sig is None:
        cross=(x[:,None]*x[None,:])*neg
    else:
        cross=(x[:,None]*x[None,:])*neg*Sig
    d=(x**2)*np.diag(G)
    return float(d.sum()-cross.sum())
wM2,evM=np.linalg.eigh(M_c)
Wbot=evM[:,:3]   # the 3 smallest eigenvectors = the (N-3)-frame
print('  spec(M) ascending =',np.array2string(wM2,precision=12))
print('  sum of 3 smallest eigenvalues of D_x Sigma D_x = %.16f'%float(wM2[:3].sum()))
print('  -tau_0 = %.16f'%(-tau0))
Ebot=Efun(x_c,Wbot)
print('  E(x_champion, W_bottom3) = %.16f   (should equal sum of the 3 smallest = 1 - Lam_3)'%(Ebot))
chk(abs(Ebot-float(wM2[:3].sum()))<1e-12,'dual functional E reproduces tr(W^T M W) exactly (dev %.2e)'%abs(Ebot-float(wM2[:3].sum())))
chk(abs(Ebot+tau0)<1e-12,'champion dual attains E = -tau_0 exactly (dev %.2e) => the dual form is TIGHT at the champion'%abs(Ebot+tau0))
# E minimized over the mask (worst-case Sigma) for a random frame:
rng=np.random.default_rng(7)
mins=[]
for trial in range(400):
    Nn=rng.integers(4,11); d=Nn-3
    A=rng.normal(size=(Nn,d)); W,_=np.linalg.qr(A); W=W  # Nn x d, orthonormal cols
    xx=rng.random(Nn)**2; xx=xx/np.linalg.norm(xx)
    Efree=Efun(xx,W)      # Sigma = mask of negativity (worst case, loops on)
    mins.append(Efree)
mins=np.array(mins)
print('  random scan (N=4..10, 400 trials): min E = %.15f ; -tau_0 = %.15f ; violations = %d'%(mins.min(),-tau0,int((mins< -tau0-1e-12).sum())))
chk(mins.min()>=-tau0-1e-9,'no random violation of E >= -tau_0 (worst %.15f, margin %.3e)'%(mins.min(),mins.min()+tau0))

print()
print('='*108)
print('(A4) Gram-energy budget: for d+3 vectors with sum w_i w_i^T = I_d,  sum_{i!=j}<w_i,w_j>^2 <= 3d/(d+3)')
print('='*108)
rng=np.random.default_rng(11); worstb=0; worstrec=None; per={}
for d in range(1,15):
    Nn=d+3
    for t in range(300):
        A=rng.normal(size=(Nn,d))
        # project to make sum w_i w_i^T = I_d exactly: A (A^T A)^{-1/2}
        ev,EV=np.linalg.eigh(A.T@A); A=A@(EV@np.diag(ev**-0.5)@EV.T)
        W=A; G=W@W.T; off=G**2; np.fill_diagonal(off,0.0)
        b=float(off.sum()); bound=3*d/(d+3)
        if b-worstb>0: worstb=b; worstrec=(d,b,bound)
        per[d]=max(per.get(d,-1.0),b)
        if b>bound+1e-12: print('   VIOLATION d=%d b=%.15f bound=%.15f'%(d,b,bound))
print('  max off-diagonal Gram energy over d=1..14, 300 random frames each: %.15f'%worstb)
print('  per-d: max off-diag Gram energy vs bound 3d/(d+3):')
for d in sorted(per): print('     d=%2d   %.15f   <=  %.15f   %s'%(d,per[d],3*d/(d+3),'OK' if per[d]<=3*d/(d+3)+1e-12 else 'VIOLATION'))
print('  tight frame (equal norms) attains 3d/(d+3): d=3 -> 9/2=4.5? no: 3*3/6=1.5 ; champion Q: sum_{i!=j}<v_i,v_j>^2 = 3 - sum g_i^2')
print('  equality case = equal-norm tight frames (||w_i||^2 = d/(d+3)); built by row/column Sinkhorn:')
for d in [1,2,3,4,6,10]:
    Nn=d+3; rng2=np.random.default_rng(100+d); A=rng2.normal(size=(Nn,d))
    for it in range(400):
        A=A*(math.sqrt(d/Nn)/np.linalg.norm(A,axis=1))[:,None]
        ev,EV=np.linalg.eigh(A.T@A); A=A@(EV@np.diag(np.maximum(ev,1e-14)**-0.5)@EV.T)
    G=A@A.T; off=G**2; np.fill_diagonal(off,0.0)
    iso=np.abs(A.T@A-np.eye(d)).max(); nr=np.abs((A**2).sum(1)-d/Nn).max()
    print('   d=%2d  sum_{i!=j}<.>^2 = %.15f   bound 3d/(d+3) = %.15f   gap=%.2e  (iso %.1e, norm-dev %.1e)'%(d,off.sum(),3*d/(d+3),3*d/(d+3)-off.sum(),iso,nr))
    chk(abs(off.sum()-3*d/(d+3))<1e-8,'d=%d: equal-norm tight frame ATTAINS the budget 3d/(d+3) = %.12f (gap %.1e)'%(d,3*d/(d+3),3*d/(d+3)-off.sum()))
chk(worstb<=3+1e-12,'dimension-free form: sum_{i!=j}<w_i,w_j>^2 <= 3d/(d+3) <= 3 on 4200 random frames (max %.12f)'%worstb)

print()
print('='*108)
print('(A5) B_1 <= 3 and B_2 == 0 on random isotropic frames (Cauchy-Schwarz proof: lem:b1three)')
print('='*108)
rng=np.random.default_rng(3); worstB1=0; worstB2=0; rec=None; worstexc=-1e9
for t in range(3000):
    Nn=int(rng.integers(3,13))
    A=rng.normal(size=(3,Nn)); M=A@A.T
    if np.linalg.cond(M)>1e10: continue
    ev,EV=np.linalg.eigh(M); V=(EV@np.diag(ev**-0.5)@EV.T@A).T
    g=(V**2).sum(1); U=V/np.sqrt(g)[:,None]
    K=U@U.T
    B1=float(np.linalg.norm((g[:,None]*U).sum(0))**2)
    P2=1.5*K**2-0.5
    B2=float((np.outer(g,g)*P2).sum())
    Gm=(U*g[:,None]).T@U*g[0]*0+np.einsum('i,ij,ik->jk',g,U,U)
    lam=float(np.linalg.eigvalsh(Gm)[-1]); excess=B1-3*lam
    if excess>worstexc: worstexc=excess; rec=(Nn,lam)
    if B1>worstB1: worstB1=B1
    worstB2=max(worstB2,abs(B2))
print('  max B_1 = %.15f ; max |B_2| = %.3e'%(worstB1,worstB2))
print('  C-S-with-respect-to-mu form:  B_1 <= 3*lambda_max(sum g_i u_i u_i^T) ; worst excess = %.3e (at N=%d)'%(worstexc,rec[0]))
# equality case: 3 points on the latitude circle u.e = 1/sqrt(3), 120 deg apart, equal weights 1
c=1.0/math.sqrt(3.0); s=math.sqrt(1-c*c)
Ueq=np.array([[s,0.0,c],[s*math.cos(2*math.pi/3),s*math.sin(2*math.pi/3),c],[s*math.cos(4*math.pi/3),s*math.sin(4*math.pi/3),c]])
geq=np.ones(3)
iso=np.abs((geq[:,None]*Ueq).T@ (geq[:,None]*Ueq)*0+ (Ueq.T*geq)@Ueq - np.eye(3)).max()
B1eq=float(np.linalg.norm((geq[:,None]*Ueq).sum(0))**2)
print('  latitude tripod (3 pts at u.e=1/sqrt3, 120 deg, g=1): isotropy dev %.2e ; B_1 = %.15f'%(iso,B1eq))
Uor=np.eye(3); gor=np.ones(3)
print('  orthonormal triple: B_1 = %.15f'%float(np.linalg.norm((gor[:,None]*Uor).sum(0))**2))
chk(worstexc<=1e-12,'B_1 <= 3*lambda_max(G) on 3000 random frames (worst excess %.2e) => B_1<=3 exactly when G=I_3'%worstexc)
chk(abs(B1eq-3)<1e-12 and abs(iso)<1e-12,'B_1 = 3 attained by the latitude tripod (isotropy dev %.1e) => sup B_1 = 3 EXACT'%iso)
chk(worstB2<1e-12,'B_2 == 0 on every isotropic frame (max |B_2| = %.2e)'%worstB2)

print()
print('='*108); print('SUMMARY  PASS=%d  FAIL=%d'%(PASS,FAIL)); print('='*108)
raise SystemExit(0 if FAIL==0 else 1)
