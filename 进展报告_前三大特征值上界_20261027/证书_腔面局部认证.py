# -*- coding: utf-8 -*-
# 1b-cert: 三类单条目腔面上界函数 F_m 在其 argmax 处的严格局部最大认证（50 位 + Weyl）
import mpmath as mp, numpy as np
from scipy.optimize import minimize
mp.mp.dps=50
s5=mp.sqrt(5); s7=mp.sqrt(7); s35=mp.sqrt(35)
a=(91-25*s7)/126; b=(7+5*s7)/126
Vhh=14/(-7*s5-7+5*s7+3*s35); Vhc=mp.sqrt(mp.mpf(7)/5)
Vcc=-mp.mpf(7)/2-7*s5/10+3*s35/10+5*s7/2
Ah=mp.sqrt(Vhh); Ac=Vhc/Ah; Rc=mp.sqrt(Vcc-Ac**2)
w=[a,b,b,b,b,b]
cols=[mp.matrix([Ah,0,0])]
for j in range(5):
    th=2*mp.pi*j/5; cols.append(mp.matrix([Ac,Rc*mp.cos(th),Rc*mp.sin(th)]))
A0=mp.matrix(3,6)
for k in range(6):
    for r in range(3): A0[r,k]=mp.sqrt(w[k])*cols[k][r]
Xs0=(A0[:,:3]**-1)*A0
Xf=np.array([[float(Xs0[r,3+c]) for c in range(3)] for r in range(3)])
def Pi_f(X):
    Am=np.hstack([np.eye(3),X.reshape(3,3)])
    return Am.T@np.linalg.solve(Am@Am.T,Am)
mc=(Pi_f(Xf)>0).astype(float)
def mk(i,j,v):
    m=mc.copy(); m[i,j]=m[j,i]=v; return m
masks={"class1 hc":mk(0,1,0.0),"class2 d1":mk(1,2,0.0),"class3 d2+":mk(1,3,1.0)}
# ---- S-jet order 2 in 2 vars (a, am, an, amn) ----
class S4:
    __slots__=("a","m","n","mn")
    def __init__(s,a,m,n,mn): s.a,s.m,s.n,s.mn=a,m,n,mn
    def __add__(s,o): return S4(s.a+o.a,s.m+o.m,s.n+o.n,s.mn+o.mn)
    def __sub__(s,o): return S4(s.a-o.a,s.m-o.m,s.n-o.n,s.mn-o.mn)
    def __mul__(s,o): return S4(s.a*o.a,s.a*o.m+s.m*o.a,s.a*o.n+s.n*o.a,s.a*o.mn+s.mn*o.a+s.m*o.n+s.n*o.m)
    def inv(s):
        ia=1/s.a; ia2=ia*ia
        return S4(ia,-s.m*ia2,-s.n*ia2,-s.mn*ia2+2*s.m*s.n*ia2*ia)
def certify(tag, m, Xstart):
    F=lambda x: -np.linalg.eigvalsh(m*Pi_f(x))[-1]
    best=(-9,None)
    rng=np.random.default_rng(5)
    for sd in [Xstart]+[Xstart+rng.normal(0,0.25,9) for _ in range(4)]+[rng.normal(0,1.0,9) for _ in range(4)]:
        r=minimize(F,sd,method="Nelder-Mead",options=dict(maxiter=25000,maxfev=50000,xatol=1e-13,fatol=1e-15))
        r2=minimize(F,r.x,method="Powell",options=dict(maxiter=25000,xtol=1e-13,ftol=1e-15))
        if -r2.fun>best[0]: best=(-r2.fun,r2.x)
    val,Xh=best
    X=mp.matrix(3,3)
    for r in range(3):
        for c in range(3): X[r,c]=mp.mpf(Xh[3*r+c])
    Ab=mp.matrix(3,6)
    O=mp.mpf(0)
    for i in range(3):
        for j in range(3): Ab[i,j]=mp.mpf(1) if i==j else O
        for j in range(3): Ab[i,3+j]=X[i,j]
    S0=Ab*Ab.T; Pi0=Ab.T*(S0**-1)*Ab
    M0=mp.matrix(6,6)
    for i in range(6):
        for j in range(6): M0[i,j]=mp.mpf(m[i,j])*Pi0[i,j]
    E,ER=mp.eig(M0); lam=[mp.re(e) for e in E]; i1=lam.index(max(lam))
    v=mp.matrix([mp.re(ER[k,i1]) for k in range(6)]); v=v/mp.norm(v)
    gap=lam[i1]-sorted(lam)[-2]
    idx=[(r,c) for r in range(3) for c in range(3)]
    H=mp.matrix(9,9)
    for mi,(r,c) in enumerate(idx):
        for ni,(r2c,c2) in enumerate(idx):
            if ni<mi: continue
            As=[[None]*6 for _ in range(3)]
            for i in range(3):
                for j in range(3): As[i][j]=S4(Ab[i,j],O,O,O)
                for j in range(3):
                    As[i][3+j]=S4(Ab[i,3+j], mp.mpf(1) if (i,j)==(r,c) else O, mp.mpf(1) if (i,j)==(r2c,c2) else O, O)
            Ss=[[sum(((As[i][k]*As[j][k]) for k in range(6)), S4(O,O,O,O)) for j in range(3)] for i in range(3)]
            def inv3(M):
                d=(M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])-M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])+M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0])).inv()
                C=[[M[1][1]*M[2][2]-M[1][2]*M[2][1],M[0][2]*M[2][1]-M[0][1]*M[2][2],M[0][1]*M[1][2]-M[0][2]*M[1][1]],
                   [M[1][2]*M[2][0]-M[1][0]*M[2][2],M[0][0]*M[2][2]-M[0][2]*M[2][0],M[0][2]*M[1][0]-M[0][0]*M[1][2]],
                   [M[1][0]*M[2][1]-M[1][1]*M[2][0],M[0][1]*M[2][0]-M[0][0]*M[2][1],M[0][0]*M[1][1]-M[0][1]*M[1][0]]]
                return [[d*C[i][j] for j in range(3)] for i in range(3)]
            Si=inv3(Ss)
            Pis=[[None]*6 for _ in range(6)]
            for i in range(6):
                for j in range(6):
                    acc=S4(O,O,O,O)
                    for k in range(3):
                        for l in range(3): acc=acc+As[k][i]*Si[k][l]*As[l][j]
                    Pis[i][j]=acc
            Mm=mp.matrix(6,6); Mn=mp.matrix(6,6); Mmn=mp.matrix(6,6)
            for i in range(6):
                for j in range(6):
                    mm=mp.mpf(m[i,j]); x=Pis[i][j]
                    Mm[i,j]=mm*x.m; Mn[i,j]=mm*x.n; Mmn[i,j]=mm*x.mn
            d1m=(v.T*Mm*v)[0]; d1n=(v.T*Mn*v)[0]
            LHS=lam[i1]*mp.eye(6)-M0
            um=mp.lu_solve(LHS,Mm*v-v*d1m); un=mp.lu_solve(LHS,Mn*v-v*d1n)
            cr=(v.T*Mmn*v)[0]+(um.T*Mn*v)[0]+(un.T*Mm*v)[0]
            H[mi,ni]=cr; H[ni,mi]=cr
    Ht=np.array([[float(mp.nstr(H[i,j],20)) for j in range(9)] for i in range(9)])
    err=max(abs(H[i,j]-mp.mpf(Ht[i,j])) for i in range(9) for j in range(9))
    lm=mp.mpf(float(np.linalg.eigvalsh(Ht).max()))
    cert=lm+mp.sqrt(mp.mpf(81))*err
    print("%s: Fmax=%.12f  gap=%.6f  |grad|~%.1e  lambda_max(H)<=%s  NEG=%s"%(tag,val,float(gap),float(max(abs(d1m),abs(d1n))),mp.nstr(cert,15),cert<0))
    return cert
Xst=Xf.reshape(9)
for tag,m in masks.items():
    certify(tag,m,Xst)
