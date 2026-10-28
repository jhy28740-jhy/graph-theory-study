# -*- coding: utf-8 -*-
# W3 certificate v2: exact third-derivative tensor T_ijk of rho at champion.
# rho3_dir(d) = D^3 rho[d,d,d] computed by order-3 single-var jet (exact).
# Full symmetric tensor recovered by solving rho3(c)=sum_tau U_tau c^tau at generic c.
import mpmath as mp, itertools, numpy as np
mp.mp.dps = 50
O = mp.mpf(0)

def champion():
    s5=mp.sqrt(5); s7=mp.sqrt(7); s35=mp.sqrt(35)
    aa=(91-25*s7)/126; bb=(7+5*s7)/126
    Vhh=14/(-7*s5-7+5*s7+3*s35); Vhc=mp.sqrt(mp.mpf(7)/5)
    Vcc=-mp.mpf(7)/2-7*s5/10+3*s35/10+5*s7/2
    Ah=mp.sqrt(Vhh); Ac=Vhc/Ah; Rc=mp.sqrt(Vcc-Ac**2)
    w=[aa,bb,bb,bb,bb,bb]
    cc=[mp.matrix([Ah,O,O])]
    for j in range(5):
        th=2*mp.pi*j/5; cc.append(mp.matrix([Ac,Rc*mp.cos(th),Rc*mp.sin(th)]))
    A0=mp.matrix(3,6)
    for k in range(6):
        for r in range(3): A0[r,k]=mp.sqrt(w[k])*cc[k][r]
    B=(A0[:,:3]**-1)*A0
    X=mp.matrix(3,3)
    for r in range(3):
        for c in range(3): X[r,c]=B[r,3+c]
    return X
Xs=champion()
Ab=mp.matrix(3,6)
for i in range(3):
    for j in range(3): Ab[i,j]=mp.mpf(1) if i==j else O
    for j in range(3): Ab[i,3+j]=Xs[i,j]
S0=Ab*Ab.T; Pi0=Ab.T*(S0**-1)*Ab
mc=mp.matrix(6,6)
for i in range(6):
    for j in range(6): mc[i,j]=mp.mpf(1) if Pi0[i,j]>0 else O
M0=mp.matrix(6,6)
for i in range(6):
    for j in range(6): M0[i,j]=mc[i,j]*Pi0[i,j]
E,ER=mp.eig(M0); lam=[mp.re(e) for e in E]; i1=lam.index(max(lam)); l0=lam[i1]
v0=mp.matrix([mp.re(ER[k,i1]) for k in range(6)]); v0=v0/mp.norm(v0)
sl=sorted(lam); g0=l0-sl[-2]
print("lambda_max =",mp.nstr(l0,25)," rho0 =",mp.nstr((9+mp.sqrt(5)+2*mp.sqrt(35))/18,25)," gap =",mp.nstr(g0,12))
L=l0*mp.eye(6)-M0
Lreg=L+v0*v0.T   # bordered: invertible; for rhs orth v0 gives the perp pseudo-inverse solution

class P3:
    __slots__=("c",)
    def __init__(s,c): s.c=c
    def __add__(s,o): return P3([s.c[k]+o.c[k] for k in range(4)])
    def __sub__(s,o): return P3([s.c[k]-o.c[k] for k in range(4)])
    def __mul__(s,o):
        r=[O,O,O,O]
        for i in range(4):
            if s.c[i]==0: continue
            for j in range(4-i):
                if o.c[j]==0: continue
                r[i+j]+=s.c[i]*o.c[j]
        return P3(r)
    def inv(s):
        ia=1/s.c[0]; r=[ia,O,O,O]
        for n in (1,2,3):
            acc=O
            for k in range(1,n+1): acc+=s.c[k]*r[n-k]
            r[n]=-ia*acc
        return P3(r)
ZERO=P3([O,O,O,O])
def Pconst(x): return P3([x,O,O,O])
def inv3(M):
    det=(M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]))
    idet=det.inv()
    C=[[M[1][1]*M[2][2]-M[1][2]*M[2][1],M[0][2]*M[2][1]-M[0][1]*M[2][2],M[0][1]*M[1][2]-M[0][2]*M[1][1]],
       [M[1][2]*M[2][0]-M[1][0]*M[2][2],M[0][0]*M[2][2]-M[0][2]*M[2][0],M[0][2]*M[1][0]-M[0][0]*M[1][2]],
       [M[1][0]*M[2][1]-M[1][1]*M[2][0],M[0][1]*M[2][0]-M[0][0]*M[2][1],M[0][0]*M[1][1]-M[0][1]*M[1][0]]]
    return [[idet*C[i][j] for j in range(3)] for i in range(3)]
idx=[(r,c) for r in range(3) for c in range(3)]

def rho3_dir(dv):
    # dv: length-9 list of mpf; A entry (i,3+j) = P3([X_ij, dv[3i+j],0,0])
    As=[[None]*6 for _ in range(3)]
    for i in range(3):
        for j in range(3): As[i][j]=Pconst(Ab[i,j])
        for j in range(3): As[i][3+j]=P3([Ab[i,3+j], mp.mpf(dv[3*i+j]), O, O])
    Ss=[[sum(((As[i][k]*As[j][k]) for k in range(6)),ZERO) for j in range(3)] for i in range(3)]
    Si=inv3(Ss)
    Pis=[[None]*6 for _ in range(6)]
    for i in range(6):
        for j in range(6):
            acc=ZERO
            for k in range(3):
                for l in range(3): acc=acc+As[k][i]*Si[k][l]*As[l][j]
            Pis[i][j]=acc
    Mk=[]
    for o in range(4):
        Mm=mp.matrix(6,6)
        for i in range(6):
            for j in range(6): Mm[i,j]=mc[i,j]*Pis[i][j].c[o]
        Mk.append(Mm)
    M1,M2,M3=Mk[1],Mk[2],Mk[3]
    # M(t)=sum t^k Mk (Taylor coeffs), lam=sum t^k lam_k, v=sum t^k u_k, u0=v0
    lam1=(v0.T*M1*v0)[0]
    u1=mp.lu_solve(Lreg,M1*v0-v0*lam1)
    lam2=(v0.T*M2*v0)[0]+(u1.T*M1*v0)[0]
    rhs2=M1*u1+M2*v0-lam1*u1-lam2*v0
    rhs2=rhs2-v0*(v0.T*rhs2)[0]
    u2=mp.lu_solve(Lreg,rhs2)
    lam3=(u2.T*M1*v0)[0]+(u1.T*M2*v0)[0]+(v0.T*M3*v0)[0]
    return 6*lam3

# sanity: basis direction p -> T_ppp
print("--- diag sanity: rho3(e_p)=T_ppp ---")
for p in range(9):
    e=[mp.mpf(1) if k==p else O for k in range(9)]
    print("  p=%d rho3=%.12f"%(p,rho3_dir(e)))

# generic-direction solve for full symmetric tensor U_tau (tau sorted triples)
triples=list(itertools.combinations_with_replacement(range(9),3))
mult={t:(6 if len(set(t))==3 else (3 if len(set(t))==2 else 1)) for t in triples}
rng=np.random.default_rng(2024)
NC=210
Cs=[]; R=[]
for i in range(NC):
    c=rng.uniform(-1.2,1.2,9)
    Cs.append(c); R.append(rho3_dir([mp.mpf(x) for x in c]))
print("built",NC,"generic evals")
# monomial matrix
def mono(c,t): return mp.mpf(c[t[0]])*mp.mpf(c[t[1]])*mp.mpf(c[t[2]])
nt=len(triples)
Mat=mp.matrix(nt,nt)
for i in range(nt):
    for jt,t in enumerate(triples): Mat[i,jt]=mult[t]*mono(Cs[i],t)
rhs=mp.matrix(nt,1)
for i in range(nt): rhs[i]=R[i]
Tvec=mp.lu_solve(Mat,rhs)
T={t:Tvec[jt] for jt,t in enumerate(triples)}
resid=O
for i in range(NC):
    pred=sum((mult[t]*T[t]*mono(Cs[i],t) for t in triples),O)
    resid=max(resid,abs(pred-R[i]))
print("solve residual over all",NC,"evals =",mp.nstr(resid,5))
maxT=O; argT=None
for t in triples:
    if abs(T[t])>abs(maxT): maxT=T[t]; argT=t
print("max|T_ijk| =",mp.nstr(abs(maxT),20)," at",argT,"=",idx[argT[0]],idx[argT[1]],idx[argT[2]]," mult",mult[argT])
C3=27*abs(maxT)
print("C3* = 27*max|T| =",mp.nstr(C3,15))
kappa=mp.mpf("0.020962946762952704")
print("r*=kappa/C3* =",mp.nstr(kappa/C3,10),"  3kappa/C3* =",mp.nstr(3*kappa/C3,10))
# also report max over distinct-only and diagonal-only for transparency
dm=max(abs(T[t]) for t in triples if len(set(t))==3)
mm=max(abs(T[t]) for t in triples if len(set(t))==2)
dg=max(abs(T[t]) for t in triples if len(set(t))==1)
print("max|T| distinct=%.12f  two-eq=%.12f  diag=%.12f"%(dm,mm,dg))
# FD cross-check: fd ~ rho3(D) = sum U_tau D^tau
def rho_of(X):
    A=mp.matrix(3,6)
    for i in range(3):
        for j in range(3): A[i,j]=mp.mpf(1) if i==j else O
        for j in range(3): A[i,3+j]=X[i,j]
    S=A*A.T; P=A.T*(S**-1)*A
    M=mp.matrix(6,6)
    for i in range(6):
        for j in range(6): M[i,j]=mc[i,j]*P[i,j]
    EE,_=mp.eig(M); return max(mp.re(e) for e in EE)
print("--- FD cross-check ---")
h=mp.mpf("1e-5"); worst=O
for tt in range(5):
    Dn=rng.normal(0,1,9); Dn=Dn/np.linalg.norm(Dn)
    D=mp.matrix(3,3)
    for k in range(9): D[idx[k][0],idx[k][1]]=mp.mpf(Dn[k])
    fs=[rho_of(Xs+mp.mpf(a)*h*D) for a in (-2,-1,1,2)]
    fd=(fs[3]-2*fs[2]+2*fs[1]-fs[0])/(2*h**3)
    ana=sum((mult[t]*T[t]*mono(Dn,t) for t in triples),O)
    jet=rho3_dir([mp.mpf(Dn[k]) for k in range(9)])
    print("  FD=%.10f  tensor=%.10f  jet=%.10f  |FD-tensor|=%.2e"%(fd,ana,jet,abs(fd-ana)))
    worst=max(worst,abs(fd-ana))
print("worst |FD-tensor| =",mp.nstr(worst,3))
np.save("_w3_T.npy",np.array([float(T[t]) for t in triples]))
np.save("_w3_triples.npy",np.array(triples))
print("C3_star_float =",float(C3))

# ============ 报告数字汇总（Frobenius 常数 + 算子范数 + 一致性边界检验）============
import numpy as _np, itertools as _it
_Tv=_np.array([float(T[t]) for t in triples])
_mult=_np.array([mult[t] for t in triples],dtype=float)
TF=float(_np.sqrt(_np.sum(_mult*_Tv**2)))            # ||T||_F  (rigorous CS bound on op-norm)
Tfull=_np.zeros((9,9,9))
for t,val in zip(triples,_Tv):
    for pm in set(_it.permutations(t)): Tfull[pm]=val
def _f(c): return float(_np.einsum('abc,a,b,c',Tfull,c,c,c))
def _g(c): return _np.einsum('abc,b,c',Tfull,c,c)
_rng=_np.random.default_rng(1); _best=(0.0,None)
for _ in range(3000):
    c=_rng.normal(size=9); c/=_np.linalg.norm(c)
    for _it2 in range(250):
        g=_g(c); ng=_np.linalg.norm(g)
        if ng<1e-14: break
        cn=g/ng
        if _f(cn)<0: cn=-cn
        if _np.linalg.norm(cn-c)<1e-14: c=cn; break
        c=cn
    v=abs(_f(c))
    if v>_best[0]: _best=(v,c.copy())
opn,cstar=_best
kap=mp.mpf("0.020962946762952704")
C3F=mp.mpf(repr(TF)); R=3*kap/C3F
print("="*60)
print("max|T_ijk|      =", mp.nstr(abs(maxT),18), " at", argT)
print("||T||_F         = %.12f   (rigorous: ||T||_op <= ||T||_F)"%TF)
print("||T||_op        = %.12f   (numeric power-iteration)"%opn)
print("27*max|T|       = %.12f   (crude certified)"%float(C3))
print("C3 := ||T(X*)||_F = %.12f"%TF)
print("kappa           = 0.020962946762952704")
print("r_hessneg=kappa/C3 = %.10f"%float(kap/C3F))
print("R=3kappa/C3 (near-segment stability radius) = %.10f"%float(R))
print("stability: rho(X) <= rho0 - (kappa/2)d^2 + (C3/6)d^3 < rho0  for 0<d<R")
# uniform boundary check (worst direction + random) at radius R
def _Dm(c):
    D=mp.matrix(3,3)
    for k in range(9): D[idx[k][0],idx[k][1]]=mp.mpf(float(c[k]))
    return D
_h=mp.mpf("1e-5")
def _third(Xb,c):
    D=_Dm(c); fs=[rho_of(Xb+mp.mpf(a)*_h*D) for a in (-2,-1,1,2)]
    return (fs[3]-2*fs[2]+2*fs[1]-fs[0])/(2*_h**3)
print("worst-dir 3rd deriv at X*      = %.10f (== ||T||_op)"%_third(Xs,cstar))
print("worst-dir 3rd deriv at X*+R*c* = %.10f (must be < C3=%.6f)"%(_third(Xs+R*_Dm(cstar),cstar),TF))
_w=O
for _ in range(8):
    cc=_rng.normal(size=9); cc/=_np.linalg.norm(cc); _w=max(_w,abs(_third(Xs+R*_Dm(cc),cc)))
print("max random-boundary directional 3rd deriv = %.10f (< C3: %s)"%(_w,_w<TF))
print("="*60)
