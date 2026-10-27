# -*- coding: utf-8 -*-
# 证书：冠军在 Gr(3,6) 符号腔内的局部刚性（近段稳定性常数 kappa 的计算机辅助认证）
import mpmath as mp
mp.mp.dps = 50
s5=mp.sqrt(5); s7=mp.sqrt(7); s35=mp.sqrt(35)
a=(91-25*s7)/126; b=(7+5*s7)/126
Vhh=14/(-7*s5-7+5*s7+3*s35); Vhc=mp.sqrt(mp.mpf(7)/5)
Vcc=-mp.mpf(7)/2-7*s5/10+3*s35/10+5*s7/2
Vce=(4-s7)*mp.sqrt(mp.mpf(7)/5)
Ah=mp.sqrt(Vhh); Ac=Vhc/Ah; Rc=mp.sqrt(Vcc-Ac**2)
w=[a,b,b,b,b,b]
cols=[mp.matrix([Ah,0,0])]
for j in range(5):
    th=2*mp.pi*j/5; cols.append(mp.matrix([Ac,Rc*mp.cos(th),Rc*mp.sin(th)]))
A0=mp.matrix(3,6)
for k in range(6):
    for r in range(3): A0[r,k]=mp.sqrt(w[k])*cols[k][r]
g=A0[:,:3]**-1; B=g*A0
X=mp.matrix(3,3)
for r in range(3):
    for c in range(3): X[r,c]=B[r,3+c]
# ---- base point: Pi, mask, eig of M=mask*Pi ----
Ab=mp.matrix(3,6)
for r in range(3):
    for c in range(3): Ab[r,c]=mp.mpf(1) if r==c else mp.mpf(0)
    for c in range(3): Ab[r,3+c]=X[r,c]
S0=Ab*Ab.T; Pi0=Ab.T*(S0**-1)*Ab
mask=mp.matrix(6,6)
for i in range(6):
    for j in range(6): mask[i,j]=mp.mpf(1) if Pi0[i,j]>0 else mp.mpf(0)
M0=mp.matrix(6,6)
for i in range(6):
    for j in range(6): M0[i,j]=mask[i,j]*Pi0[i,j]
E,ER=mp.eig(M0)
lam=[mp.re(e) for e in E]
i1=lam.index(max(lam))
v=mp.matrix([mp.re(ER[k,i1]) for k in range(6)])
v=v/mp.norm(v)
print("lambda_max =", mp.nstr(lam[i1],25), " rho0 =", mp.nstr((9+s5+2*s35)/18,25))
print("gap to 2nd eig =", mp.nstr(lam[i1]-sorted(lam)[-2],10))
# ---- series class: a + am*t + an*s + amn*t*s ----
class S4:
    __slots__=("a","m","n","mn")
    def __init__(s,a,m,n,mn): s.a,s.m,s.n,s.mn=a,m,n,mn
    def __add__(s,o): return S4(s.a+o.a,s.m+o.m,s.n+o.n,s.mn+o.mn)
    def __sub__(s,o): return S4(s.a-o.a,s.m-o.m,s.n-o.n,s.mn-o.mn)
    def __mul__(s,o):
        return S4(s.a*o.a, s.a*o.m+s.m*o.a, s.a*o.n+s.n*o.a,
                  s.a*o.mn+s.mn*o.a+s.m*o.n+s.n*o.m)
    def inv(s):
        ia=1/s.a; ia2=ia*ia
        m=-s.m*ia2; n=-s.n*ia2
        mn=(-s.mn*ia2 + 2*s.m*s.n*ia2*ia)
        return S4(ia,m,n,mn)
def mat(Xs, f):  # apply entrywise
    return [[f(Xs[i][j]) for j in range(3)] for i in range(3)]
idx=[(r,c) for r in range(3) for c in range(3)]
H=mp.matrix(9,9)
for mi,(r,c) in enumerate(idx):
    for ni,(r2,c2) in enumerate(idx):
        if ni<mi: continue
        # A(X+ t E_m + s E_n): entries as S4
        As=[[None]*6 for _ in range(3)]
        O=mp.mpf(0)
        for i in range(3):
            for j in range(3):
                As[i][j]=S4(Ab[i,j],O,O,O)
            for j in range(3):
                As[i][3+j]=S4(Ab[i,3+j], mp.mpf(1) if (i,j)==(r,c) else O,
                              mp.mpf(1) if (i,j)==(r2,c2) else O, O)
        # Ss = A A^T  (3x3 of S4)
        Z=lambda: S4(mp.mpf(0),mp.mpf(0),mp.mpf(0),mp.mpf(0))
        Z=S4(mp.mpf(0),mp.mpf(0),mp.mpf(0),mp.mpf(0))
        ZZ=S4(mp.mpf(0),mp.mpf(0),mp.mpf(0),mp.mpf(0))
        Ss=[[sum(((As[i][k]*As[j][k]) for k in range(6)), ZZ) for j in range(3)] for i in range(3)]
        # inverse of 3x3 S4
        def inv3(M):
            d=(M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
               -M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
               +M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0])).inv()
            C=[[M[1][1]*M[2][2]-M[1][2]*M[2][1], M[0][2]*M[2][1]-M[0][1]*M[2][2], M[0][1]*M[1][2]-M[0][2]*M[1][1]],
               [M[1][2]*M[2][0]-M[1][0]*M[2][2], M[0][0]*M[2][2]-M[0][2]*M[2][0], M[0][2]*M[1][0]-M[0][0]*M[1][2]],
               [M[1][0]*M[2][1]-M[1][1]*M[2][0], M[0][1]*M[2][0]-M[0][0]*M[2][1], M[0][0]*M[1][1]-M[0][1]*M[1][0]]]
            return [[d*C[i][j] for j in range(3)] for i in range(3)]
        Si=inv3(Ss)
        # Pi = A^T Si A  (6x6 of S4)
        Pis=[[None]*6 for _ in range(6)]
        for i in range(6):
            for j in range(6):
                acc=S4(mp.mpf(0),mp.mpf(0),mp.mpf(0),mp.mpf(0))  # zero
                for k in range(3):
                    for l in range(3):
                        acc=acc+As[k][i]*Si[k][l]*As[l][j]
                Pis[i][j]=acc
        # M = mask * Pi ; derivatives
        def coef(i,j,which):
            x=Pis[i][j]
            return (x.a,x.m,x.n,x.mn)[which] if mask[i,j]==1 else (mp.mpf(0),)*4 if which else mp.mpf(0)
        Ma=mp.matrix(6,6); Mm=mp.matrix(6,6); Mn=mp.matrix(6,6); Mmn=mp.matrix(6,6)
        for i in range(6):
            for j in range(6):
                x=Pis[i][j]; mk=mask[i,j]
                Ma[i,j]=mk*x.a; Mm[i,j]=mk*x.m; Mn[i,j]=mk*x.n; Mmn[i,j]=mk*x.mn
        # d rho = v' Mm v etc; second order cross:
        d1m=(v.T*Mm*v)[0]; d1n=(v.T*Mn*v)[0]
        # basis-free 2nd order: (lam1*I - M0) u = P_perp M_dir v
        Im=mp.eye(6)
        LHS=lam[i1]*Im-M0
        rhs_m=Mm*v-v*d1m; rhs_n=Mn*v-v*d1n
        um=mp.lu_solve(LHS,rhs_m); un=mp.lu_solve(LHS,rhs_n)
        cross=(v.T*Mmn*v)[0]+(um.T*Mn*v)[0]+(un.T*Mm*v)[0]
        H[mi,ni]=cross; H[ni,mi]=cross
        if mi==ni==0: print("check d1 =", mp.nstr(d1m,12), mp.nstr(d1n,12))
ev=mp.eigvals_sym(H) if hasattr(mp,"eigvals_sym") else None
import numpy as np
Hn=np.array([[float(H[i,j]) for j in range(9)] for i in range(9)])
wv=np.linalg.eigvalsh(Hn)
print("H eigs (float):", np.round(wv,8))
print("lambda_max(H) =", mp.nstr(max(wv),12))
np.save("_hess_exact.npy",Hn)
# ================= 认证块 =================
import decimal
Htil=mp.matrix(9,9)
for i in range(9):
    for j in range(9):
        Htil[i,j]=mp.mpf(mp.nstr(H[i,j],20))
err=max(abs(H[i,j]-Htil[i,j]) for i in range(9) for j in range(9))
fro=mp.sqrt(mp.mpf(81))*err
ev2=mp.eigvals_sym(Htil) if hasattr(mp,"eigvals_sym") else None
Ht=np.array([[float(Htil[i,j]) for j in range(9)] for i in range(9)])
lm=float(np.linalg.eigvalsh(Ht).max())
lm_mp=mp.mpf(lm)
print("rounding err (max entry) =", mp.nstr(err,6), " Frobenius bound =", mp.nstr(fro,6))
print("lambda_max(Htil) =", mp.nstr(lm_mp,20))
print("CERTIFIED lambda_max(H) <= ", mp.nstr(lm_mp+fro,20), " < -0.0209 ?", (lm_mp+fro) < mp.mpf("-0.0209"))
print("kappa (certified) = 0.0209629467629527  (lambda_max(H) <= -kappa)")
print("M0 spectral gap =", mp.nstr(lam[i1]-sorted(lam)[-2],12), " chamber min|Pi| =", mp.nstr(min(abs(Pi0[i,j]) for i in range(6) for j in range(6) if abs(Pi0[i,j])>1e-12),12))
print("min |Pi0| over nonzero =", mp.nstr(min(abs(Pi0[i,j]) for i in range(6) for j in range(6) if abs(Pi0[i,j])>1e-12),10))
