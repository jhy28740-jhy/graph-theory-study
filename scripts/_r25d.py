# -*- coding: utf-8 -*-
"""R25-d: Phi 泛函的极小值搜索 + 锥族剖面 + 与 thm:r23uni 的 F(a,z) 对照."""
import numpy as np
from numpy.linalg import norm, eigh
from scipy.optimize import minimize
exec(open('_r25c.py',encoding='utf-8-sig').read().split('# ---- champion ----')[0])

print("delta0/3 = %.15f   2*rho0 = %.15f"%(del0/3,2*rho0))
astar=(3*S35+5*np.sqrt(7)-7-7*S5)/14; zstar=0.5004400497682123079777569143187991903444
rho_st=(3*astar-1)/2

def cone(n,z):
    G=2/(n*(1-z*z)); g0=1-n*G*z*z
    U=[np.array([0,0,1.])]
    for k in range(n):
        p=2*np.pi*k/n
        U.append(np.array([np.sqrt(1-z*z)*np.cos(p),np.sqrt(1-z*z)*np.sin(p),z]))
    return np.array([g0]+[G]*n), np.array(U)

# ---- A) 验证 thm:r23uni 的闭式 F(a,z) 就是 J_w/2 ----
print("\n=== A) 闭式 F(a,z) 与 J_w/2 的对照 (n=5 锥) ===")
def Fclosed(a,z):
    E=1-a+(3*a-1)*z*z
    Psi=a*(1-3*z*z)**2+2*np.sqrt(2)*z*(1-3*z*z)*np.sqrt(a*E)+E/5*((5-S5)*z*z+1+S5)
    return Psi/(1-z*z)**2
worst=0
for a in [0.45,0.5,0.55,astar,0.65,0.7]:
    for z in [0.3,0.4,zstar,0.6]:
        g,U=cone(5,z); r=full(g,U,uniax((3*a-1)/2))
        if abs(r['half']-Fclosed(a,z))>1e-9: continue
        worst=max(worst,abs(r['half']-Fclosed(a,z)))
print("  在闭式适用区域内 max|J_w/2 - F(a,z)| = %.3e"%worst)
g,U=cone(5,zstar); rc=full(g,U,uniax(rho_st))
print("  冠军: J_w/2=%.15f  F(a*,z*)=%.15f  差=%.2e"%(rc['half'],Fclosed(astar,zstar),rc['half']-Fclosed(astar,zstar)))
print("  冠军: Phi=%.15f  delta0/3=%.15f  差=%.2e"%(rc['Phi'],del0/3,rc['Phi']-del0/3))
print("  冠军 A_w=%.15f ; 1-E/2 形式核对: E=%.12f, A_w-(1-E/2)=%.2e"%(
    rc['A'], 2-2*rc['A'], rc['A']-(1-(2-2*rc['A'])/2)))

# ---- B) 锥族 Phi 剖面: min over (n,a,z) ----
print("\n=== B) 锥族 Phi 的最小值 (应 = delta0/3, 在 n=5,a=a*,z=z*) ===")
best=(1e9,None)
for n in range(3,11):
    bn=(1e9,None)
    for a in np.linspace(1/3+1e-4,0.9999,120):
        rhoa=(3*a-1)/2
        for z in np.linspace(0.05,0.95,120):
            G=2/(n*(1-z*z)); g0=1-n*G*z*z
            if g0<=1e-8: continue
            g,U=cone(n,z); r=full(g,U,uniax(rhoa),seed=1)
            if r['Phi']<bn[0]: bn=(r['Phi'],(a,z,r['half']))
    print("  n=%2d  min Phi=%.12f at a=%.5f z=%.5f (J/2=%.9f)"%(n,bn[0],bn[1][0],bn[1][1],bn[1][2]))
    if bn[0]<best[0]: best=bn
print("  锥族总最小 Phi=%.15f ; 与 delta0/3 之差=%.3e"%(best[0],best[0]-del0/3))

# ---- C) 随机帧 + Stiefel 局部下降, 搜索 Phi < delta0/3 ----
print("\n=== C) 随机帧上 Phi 的局部极小化 (Stiefel 投影梯度/Nelder-Mead) ===")
def Phi_of_V(V,rhoa,seed=3):
    g=(V*V).sum(1); U=V/np.sqrt(g)[:,None]
    return full(g,U,uniax(rhoa),seed=seed)['Phi']
def obj(x,N,rhoa):
    M=x.reshape(N,3); Q,_=np.linalg.qr(M)
    R=np.linalg.cholesky(np.linalg.inv(Q.T@(M.T@M)@np.linalg.inv(Q.T@M).T)) if False else None
    V=Q@np.linalg.cholesky(Q.T@Q)  # Q already orthonormal
    return Phi_of_V(V,rhoa)
rng=np.random.default_rng(2024)
gbest=(1e9,None)
for N in [5,6,7,8]:
    for rhoa in [0.25,0.35,rho_st,0.45,0.6,0.8]:
        loc=(1e9,None)
        for trial in range(6):
            M=rng.standard_normal((N,3)); Q,_=np.linalg.qr(M); x0=Q.reshape(-1)
            f=lambda x: obj(x,N,rhoa)
            r=minimize(f,x0,method='Powell',options=dict(maxiter=4000,xtol=1e-8,ftol=1e-11))
            if r.fun<loc[0]: loc=(r.fun,r.x)
        if loc[0]<gbest[0]: gbest=(loc[0],(N,rhoa))
        print("  N=%d rho=%.4f : min Phi=%.12f  (delta0/3=%.12f, 差=%+.3e)"%(N,rhoa,loc[0],del0/3,loc[0]-del0/3))
print("  全局搜索最小 Phi=%.12f at %s"%(gbest[0],gbest[1]))
