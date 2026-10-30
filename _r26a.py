# -*- coding: utf-8 -*-
"""秩一截面: L_w = 1 + 2 Neg 的精确研究 (锥族闭式 + 一般帧随机/局部搜索)"""
import numpy as np
from numpy.linalg import norm
from scipy.optimize import minimize_scalar, minimize
from mpmath import mp, mpf, sqrt as msqrt
mp.dps=50
S5=np.sqrt(5.0)
rho0=(9+S5+2*np.sqrt(35))/18; phi=(1+S5)/2
del0=3*phi+3-6*rho0; TGT=del0/3
LIM=2*rho0-1
print("rho0   = %.18f"%rho0)
print("2rho0-1= phi - delta0/3 = %.18f   (= %.18f)"%(LIM, phi-TGT))
print("rho0-1 = %.18f   (Neg 的等价上界)"%(rho0-1))

def ring(n,s):
    """纯环: n 个等权原子, 纬度 t^2=s. 各向同性 => g=3/n, 且 s=1/3 才可能(纯环)."""
    return None

def cone(n,s):
    """极点 + 正 n 边形, 纬度 t^2=s (s<=1/3). 返回 (g, U)."""
    s=mpf(s) if not isinstance(s,float) else s
    g1=2.0/(n*(1.0-s)); g0=(1.0-3.0*s)/(1.0-s)
    t=np.sqrt(s); b=np.sqrt(1.0-s)
    U=[np.array([0,0,1.0])]
    for k in range(n):
        p=2*np.pi*k/n; U.append(np.array([b*np.cos(p),b*np.sin(p),t]))
    return np.array([g0]+[g1]*n), np.array(U), g0, g1

def Lw_rank1(g,U,z=np.array([0,0,1.0])):
    a=np.abs(U@z); w=g*a
    G=U@U.T
    return float((np.outer(w,w)*np.abs(G)).sum())

def Lw_closed(n,s):
    """闭式: L_w = [(1-3s)(1+s) + 4s(1+R_n(s))/n]/(1-s)^2"""
    R=sum(abs((1-s)*np.cos(2*np.pi*d/n)+s) for d in range(1,n))
    return ((1-3*s)*(1+s)+4*s*(1+R)/n)/(1-s)**2

print("\n=== A) 锥族 (极点 + 正 n 边形) 在秩一截面 z=轴 上的 L_w ===")
print("    s=t^2 在 [0,1/3] 上扫描; 需要 L_w <= %.10f"%LIM)
best=(0,None,None)
for n in range(3,25):
    ss=np.linspace(0,1/3,20001)
    vals=np.array([Lw_closed(n,s) for s in ss])
    k=int(np.argmax(vals)); smax=ss[k]; vmax=vals[k]
    r=minimize_scalar(lambda s:-Lw_closed(n,s),bounds=(0,1/3),method="bounded",
                      options=dict(xatol=1e-14))
    s2=-r.fun; sm2=float(r.x)
    if s2>vmax: vmax,smax=s2,sm2
    # 交叉核对: 用显式帧算
    gg,UU,g0,g1=cone(n,smax)
    direct=Lw_rank1(gg,UU)
    err=abs(direct-vmax)
    ok = "OK " if err<2e-9 else "CHK"
    print("  n=%2d  max L_w=%.12f at s=%.9f  (g0=%.6f)  闭式-显式差=%.2e %s  余量=%+.6f"
          %(n,vmax,smax,g0,err,ok,LIM-vmax))
    if vmax>best[0]: best=(vmax,n,smax)
print("  锥族总 sup L_w = %.12f  at n=%d, s=%.9f ;  余量 = %.6e"
      %(best[0],best[1],best[2],LIM-best[0]))
print("  Neg = (L_w-1)/2 = %.12f  vs  rho0-1 = %.12f"%( (best[0]-1)/2, rho0-1))
