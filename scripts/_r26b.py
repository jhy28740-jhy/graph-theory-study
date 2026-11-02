# -*- coding: utf-8 -*-
"""秩一截面 sup L_w 的一般帧搜索 (非锥), 与锥族闭式 (6+phi)/5 对比"""
import numpy as np
from numpy.linalg import norm, qr, svd
from scipy.optimize import minimize
rng=np.random.default_rng(20261030)
S5=np.sqrt(5.0); phi=(1+S5)/2
rho0=(9+S5+2*np.sqrt(35))/18
LIM=2*rho0-1; CONE=(6+phi)/5
print("目标上界 2rho0-1 = %.15f"%LIM)
print("锥族 sup  (6+phi)/5 = %.15f   余量 %.6e"%(CONE,LIM-CONE))

def frame(V):
    g=(V*V).sum(1); return g, V/np.sqrt(g)[:,None]

def Lw(g,U,z):
    a=np.abs(U@z); w=g*a; G=U@U.T
    return float((np.outer(w,w)*np.abs(G)).sum())

def best_z(g,U):
    """max over z of L_w : 用候选 z (原子方向/随机) + 局部细化"""
    cand=list(U)+[rng.standard_normal(3) for _ in range(6)]
    bs=(-1,None)
    for c in cand:
        z=c/norm(c); v=Lw(g,U,z)
        if v>bs[0]: bs=(v,z)
    z=bs[1]
    for _ in range(60):   # 幂迭代式细化
        gstep=1e-3
        improved=False
        for _ in range(30):
            p=rng.standard_normal(3); p-= (p@z)*z; p/=norm(p)
            for st in (gstep,-gstep,gstep*3,-gstep*3):
                z2=z+st*p; z2/=norm(z2); v=Lw(g,U,z2)
                if v>bs[0]+1e-15: bs=(v,z2); z=z2; improved=True; break
            if improved: break
        if not improved: gstep*=0.5
        if gstep<1e-9: break
    return bs

def negLw_from_V(V,zpar):
    g,U=frame(V@V.T+1e-12*np.eye(V.shape[0])) if False else frame(V)
    z=zpar/norm(zpar)
    return -Lw(g,U,z)

print("\n=== B) 一般帧随机 + 局部极大化 (N=4..12) ===")
glob=(-1,None)
for N in range(4,13):
    best=(-1,None)
    for trial in range(400):
        M=rng.standard_normal((N,3))
        Q,_=qr(M); V=Q[:,:3] if N>=3 else None
        V=Q[:, :3]
        # 让 V 是 St(N,3): V^T V = I_3
        g,U=frame(V)
        v,z=best_z(g,U)
        if v>best[0]: best=(v,(N,V.copy(),z.copy()))
    # 局部极大化 (Nelder-Mead on V entries + z), 保持 Stiefel 用 QR 重投影
    V0=best[1][1].copy(); z0=best[1][2].copy()
    def obj(x):
        Vm=x[:N*3].reshape(N,3); zz=x[N*3:]
        Qm,_=qr(Vm); Vs=Qm[:,:3]
        g,U=frame(Vs)
        return -Lw(g,U,zz/norm(zz))
    x0=np.concatenate([V0.ravel(),z0])
    r=minimize(obj,x0,method="Powell",options=dict(maxiter=60000,maxfev=120000,xtol=1e-13,ftol=1e-14))
    vloc=-r.fun
    if vloc>best[0]: best=(vloc,("loc",r.x,z0))
    print("  N=%2d  best L_w = %.12f   (锥族 sup = %.12f, 差 %+.3e)"%(N,best[0],CONE,best[0]-CONE))
    if best[0]>glob[0]: glob=best
print("\n  一般帧全局最好 L_w = %.12f"%glob[0])
print("  与锥族 sup 之差 = %+.6e ; 与 2rho0-1 之差 = %+.6e"%(glob[0]-CONE, LIM-glob[0]))
