# -*- coding: utf-8 -*-
"""R25-b: 单轴/全局截面极值剖面探测 —— ½[L_w+A_w^2] 的上确界."""
import numpy as np
from numpy.linalg import eigh, norm, qr, svd

S5=np.sqrt(5.0); S35=np.sqrt(35.0)
rho0=(9+S5+2*S35)/18

def Aw(U,w,nstart=60,seed=0):
    """max_{||z||=1} sum w_i |<u_i,z>|  via sign fixed-point, multi-start."""
    n=len(w); rng=np.random.default_rng(seed); best=(0.0,None)
    starts=[U[i]/norm(U[i]) for i in range(n)]
    starts+= [rng.standard_normal(3) for _ in range(nstart)]
    starts.append(np.array([0,0,1.0])); starts.append(np.array([1.0,0,0]))
    for z0 in starts:
        z=z0/norm(z0)
        for _ in range(300):
            s=U@z; e=np.where(s>=0,1.0,-1.0)
            y=(e*w)@U; ny=norm(y)
            if ny<1e-14: break
            zn=y/ny
            if norm(zn-z)<1e-14: z=zn; break
            z=zn
        val=float(np.abs(U@z)@w)
        if val>best[0]: best=(val,z)
    return best[0],best[1]

def halfJ(g,U,Lam):
    """Lam: 3x3 psd trace1.  returns ½[L_w+A_w^2] and parts."""
    h=np.sqrt(np.einsum('ij,jk,ik->i',U,Lam,U))
    w=g*h; G=U@U.T
    Lw=(np.outer(w,w)*np.abs(G)).sum()
    A,z=Aw(U,w)
    return 0.5*(Lw+A*A), Lw, A*A, w, h

def uniax(rhoa):
    b=(1-rhoa)/3; return b*np.eye(3)+(rhoa-b)*np.diag([0,0,1.0])

def frame_from_V(V):
    gg=(V*V).sum(1); U=V/np.sqrt(gg)[:,None]; return gg,U

def max_over_rho(g,U,ngrid=41):
    f=lambda r: halfJ(g,U,uniax(r))[0]
    rs=np.linspace(0,1,ngrid); vals=[f(r) for r in rs]
    i=int(np.argmax(vals))
    lo=rs[max(0,i-1)]; hi=rs[min(len(rs)-1,i+1)]
    for _ in range(60):
        m1=lo+(hi-lo)*0.382; m2=lo+(hi-lo)*0.618
        if f(m1)<f(m2): lo=m1
        else: hi=m2
    r=0.5*(lo+hi); return r,f(r)

# ---------- 1) cone family scan ----------
print("=== 1) 锥族 (极点 + 正 n 边形) 扫描 ===")
best=(0,None)
for n in range(3,13):
    bn=(0,None)
    for a in np.linspace(1/3+1e-6,1-1e-6,140):
        rhoa=(3*a-1)/2
        for zz in np.linspace(-0.999,0.999,200):
            Gc=2/(n*(1-zz*zz)); gh=1-n*Gc*zz*zz
            if gh<=1e-9: continue
            U=[np.array([0,0,1.0])]
            for k in range(n):
                ph=2*np.pi*k/n
                U.append(np.array([np.sqrt(1-zz*zz)*np.cos(ph),np.sqrt(1-zz*zz)*np.sin(ph),zz]))
            U=np.array(U); g=np.array([gh]+[Gc]*n)
            v=halfJ(g,U,uniax(rhoa))[0]
            if v>bn[0]: bn=(v,(a,zz))
    print("  n=%2d  max=%.12f  at a=%.6f z=%.6f   (rho0=%.12f, diff=%.3e)"%(n,bn[0],bn[1][0],bn[1][1],rho0,bn[0]-rho0))
    if bn[0]>best[0]: best=bn
print("  锥族总最大 = %.15f ; 与 rho0 之差 = %.3e"%(best[0],best[0]-rho0))

# ---------- 2) random frames ----------
print("\n=== 2) 随机帧: max over rho ===")
rng=np.random.default_rng(12345)
worst=(0,None)
for N in [4,5,6,7,8,9,10,12,14]:
    bm=(0,None)
    for trial in range(40):
        M=rng.standard_normal((N,3)); Q,_=np.linalg.qr(M); V=Q[:,:3]
        g,U=frame_from_V(V)
        r,v=max_over_rho(g,U,ngrid=25)
        if v>bm[0]: bm=(v,r)
    print("  N=%2d  best over 40 random frames: %.12f at rho=%.4f  (diff to rho0 = %+.3e)"%(N,bm[0],bm[1],bm[0]-rho0))
    if bm[0]>worst[0]: worst=bm
print("  随机帧全局最好 = %.15f"%worst[0])
