# -*- coding: utf-8 -*-
import numpy as np, math
from scipy.optimize import minimize
s5=math.sqrt(5); s7=math.sqrt(7); s35=math.sqrt(35)
a=(91-25*s7)/126; b=(7+5*s7)/126
Vhh=14/(-7*s5-7+5*s7+3*s35); Vhc=math.sqrt(7/5)
Vcc=-7/2-7*s5/10+3*s35/10+5*s7/2
Vce=(4-s7)*math.sqrt(7/5)
Ah=math.sqrt(Vhh); Ac=Vhc/Ah; Rc=math.sqrt(Vcc-Ac**2)
w=np.array([a,b,b,b,b,b])
cols=[np.array([Ah,0,0])]
for j in range(5):
    th=2*math.pi*j/5; cols.append(np.array([Ac,Rc*math.cos(th),Rc*math.sin(th)]))
A0=np.stack([math.sqrt(w[k])*cols[k] for k in range(6)],axis=1)
Xs=(np.linalg.inv(A0[:,:3])@A0)[:,3:6]
def Pi_of(X):
    Am=np.hstack([np.eye(3),X.reshape(3,3)])
    return Am.T@np.linalg.solve(Am@Am.T,Am)
P0=Pi_of(Xs); rho0=(9+s5+2*s35)/18
mc=(P0>0).astype(float)
def sym(M,i,j,v):
    M[i,j]=M[j,i]=v
def masks():
    out={}
    m=mc.copy(); sym(m,0,1,0); sym(m,0,2,0); out["P1 hc+hc"]=m
    m=mc.copy(); sym(m,0,1,0); sym(m,1,2,0); out["P2 hc+d1"]=m
    m=mc.copy(); sym(m,1,2,0); sym(m,2,3,0); out["P3 d1+d1 adj"]=m
    m=mc.copy(); sym(m,1,2,0); sym(m,3,4,0); out["P4 d1+d1 dis"]=m
    m=mc.copy(); sym(m,1,3,1); sym(m,2,4,1); out["P5 d2+ d2+"]=m
    m=mc.copy(); sym(m,0,1,0); sym(m,1,3,1); out["P6 hc- d2+"]=m
    return out
rng=np.random.default_rng(31)
seeds=[Xs.reshape(9)]+[ (Xs.reshape(9)+rng.normal(0,0.3,9)) for _ in range(4)]+[rng.normal(0,1.0,9) for _ in range(6)]
for tag,m in masks().items():
    F=lambda x: -np.linalg.eigvalsh(m*Pi_of(x))[-1]
    best=(-9,None)
    for sd in seeds:
        r=minimize(F,sd,method="Nelder-Mead",options=dict(maxiter=25000,maxfev=50000,xatol=1e-13,fatol=1e-15))
        r2=minimize(F,r.x,method="Powell",options=dict(maxiter=25000,xtol=1e-13,ftol=1e-15))
        if -r2.fun>best[0]: best=(-r2.fun,r2.x)
    print("%-14s max=%.12f margin=%+.6f"%(tag,best[0],rho0-best[0]))
    np.save("_w1p_%s.npy"%tag.split()[0], best[1])
