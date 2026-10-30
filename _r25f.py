# -*- coding: utf-8 -*-
"""R25-f: Phi 的全局极小值搜索 (一般帧 + 一般 Lambda), 寻找反例."""
import numpy as np
from numpy.linalg import norm, eigh
from scipy.optimize import minimize
src=open('_r25e.py',encoding='utf-8-sig').read().split('astar=(3*S35')[0]
exec(src)
astar=(3*S35+5*np.sqrt(7)-7-7*S5)/14; zstar=0.5004400497682123079777569143187991903444
rho_st=(3*astar-1)/2
print("target delta0/3 = %.15f"%TGT)
rng=np.random.default_rng(987654321)

def randLam(kind,rng):
    if kind==0: r=rng.random(); return uniax(r)
    if kind==1: ev=np.sort(rng.random(3))[::-1]; ev=ev/ev.sum()
    else:
        ev=np.array([1.0,0,0]); 
    M=rng.standard_normal((3,3)); O=eigh(M@M.T)[1]
    return O@np.diag(ev)@O.T

print("\n=== 1) 大规模随机扫描 (一般帧 x 一般 Lambda) ===")
gmin=(1e9,None); cnt=0
for N in range(4,13):
    bn=(1e9,None)
    for _ in range(1200):
        M=rng.standard_normal((N,3)); Q,_=np.linalg.qr(M); V=Q
        g=(V*V).sum(1); U=V/np.sqrt(g)[:,None]
        Lam=randLam(int(rng.integers(0,3)),rng)
        d=fast(g,U,Lam); cnt+=1
        if d['Phi']<bn[0]: bn=(d['Phi'],(N,))
        if d['Phi']<gmin[0]: gmin=(d['Phi'],(N,Lam.copy(),d['half']))
    print("  N=%2d  min Phi over 1200 = %.12f  (diff to target %+.3e)"%(N,bn[0],bn[0]-TGT))
print("  总计 %d 个随机构型, 全局最小 Phi = %.12f (diff %+.3e)"%(cnt,gmin[0],gmin[0]-TGT))

print("\n=== 2) 局部极小化 (Powell on Stiefel x Lambda谱), 从冠军与随机点出发 ===")
def pack(V,ev): return np.concatenate([V.reshape(-1),ev])
def unpack(x,N):
    M=x[:3*N].reshape(N,3); Q,_=np.linalg.qr(M); V=Q
    ev=np.abs(x[3*N:3*N+3]); s=ev.sum(); ev=ev/s if s>1e-9 else np.array([1/3]*3)
    return V,ev
def makeLam(ev,M):
    O=eigh(M@M.T)[1]; return O@np.diag(ev)@O.T
def Phifun(x,N,M):
    V,ev=unpack(x,N); g=(V*V).sum(1); U=V/np.sqrt(g)[:,None]
    return fast(g,U,makeLam(ev,M))['Phi']
gb=(1e9,None)
for N in [5,6,7,8]:
    for tr in range(8):
        if tr==0 and N==6:
            g0,U0=cone(5,zstar); V0=(np.sqrt(g0)[:,None]*U0); ev0=np.sort(np.linalg.eigvalsh(uniax(rho_st)))[::-1]; x0=pack(V0,ev0)
        else:
            M0=rng.standard_normal((N,3)); Q,_=np.linalg.qr(M0); x0=pack(Q,rng.random(3))
        M=rng.standard_normal((3,3))
        r=minimize(Phifun,x0,args=(N,M),method='Powell',options=dict(maxiter=6000,maxfev=60000,xtol=1e-10,ftol=1e-13))
        if r.fun<gb[0]: gb=(r.fun,(N,tr,r.x))
    print("  N=%d 最好局部极小 Phi = %.12f  (diff %+.3e)"%(N,gb[0],gb[0]-TGT))
print("  全局最好 = %.15f  at N=%d trial=%d"%(gb[0],gb[1][0],gb[1][1]))
V,ev=unpack(gb[1][2],gb[1][0]); print("   谱 =",np.round(np.sort(ev)[::-1],9))
g=(V*V).sum(1); U=V/np.sqrt(g)[:,None]; d=fast(g,U,makeLam(ev,rng.standard_normal((3,3))))
print("   该点 J_w/2 = %.12f (rho0=%.12f)"%(d['half'],rho0))
