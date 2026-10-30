# -*- coding: utf-8 -*-
"""R25-e: 快速版 —— Phi 泛函; 符号枚举向量化; 锥族剖面 + 全局搜索."""
import numpy as np
from numpy.linalg import norm, eigh
S5=np.sqrt(5.0); S35=np.sqrt(35.0)
rho0=(9+S5+2*S35)/18; phi=(1+S5)/2
gam=(1+S5)/16; alp=(7*S5-5)/8; bet=(25-15*S5)/16
cW=(2*S5-1)/14; cS=(10+S5)/14; cZ=(5-3*S5)/14
del0=3*phi+3-6*rho0; TGT=del0/3
def P4(t): return (35*t**4-30*t*t+3)/8
def qf(s): return gam+alp*s+bet*s*s
def uniax(r): 
    L=((1-r)/3)*np.eye(3); L[2,2]+=r; return L

_SIG={}
def signpat(n):
    if n not in _SIG:
        E=np.array([[(i>>k)&1 for k in range(n-1)] for i in range(2**(n-1))],dtype=float)
        E=1-2*E; E=np.hstack([E,np.ones((E.shape[0],1))])
        _SIG[n]=E
    return _SIG[n]

def fast(g,U,Lam):
    h=np.sqrt(np.einsum('ij,jk,ik->i',U,Lam,U)); w=g*h; n=len(g)
    G=U@U.T; W=w.sum()
    Sw=(w[:,None,None]*(U[:,:,None]*U[:,None,:])).sum(0); SF=(Sw*Sw).sum()
    Lw=(np.outer(w,w)*np.abs(G)).sum()
    M=(np.outer(w,w)*G); E=signpat(n)
    nrm2=((E@M)*E).sum(1); k=int(np.argmax(nrm2)); A=np.sqrt(max(nrm2[k],0))
    sg=E[k]; m=sg*w@U; z=m/norm(m) if norm(m)>0 else np.array([0,0,1.])
    a=np.abs(U@z); Vw=(np.outer(w,w)*(a[:,None]-a[None,:])**2).sum()
    M4=(np.outer(w,w)*G**4).sum(); Z4=(np.outer(w,w)*P4(G)).sum()
    Dw=(np.outer(w,w)*(qf(G**2)-np.abs(G))).sum()
    dh=h[:,None]-h[None,:]; gg=np.outer(g,g)
    H=(gg*dh**2).sum(); Et=(gg*dh**2*G**2).sum()
    Ea=float((g*(h-a)**2).sum())
    Phi=(1-A*A)+0.5*(cW*H+cS*Et)+(-cZ)*Z4+Dw
    r_ID=(2*rho0-(Lw+A*A))-((1-A*A)-TGT+0.5*(cW*H+cS*Et)+(-cZ)*Z4+Dw)
    r_A=A-(1-Ea/2); r_H=H-(6-2*W*W); r_E=SF-(1-Et/2)
    r_M=M4-(-3*W*W/35+6*SF/7+8*Z4/35)
    return dict(h=h,w=w,W=W,SF=SF,A=A,Aw2=A*A,z=z,Lw=Lw,M4=M4,Z4=Z4,Dw=Dw,H=H,Et=Et,
                E=Ea,Jw=Lw+A*A,half=(Lw+A*A)/2,Phi=Phi,Vw=Vw,
                t1=1-A*A,t2=0.5*cW*H,t3=0.5*cS*Et,t4=-cZ*Z4,t5=Dw,
                r_ID=r_ID,r_A=r_A,r_H=r_H,r_E=r_E,r_M=r_M)

def cone(n,z):
    Gc=2/(n*(1-z*z)); g0=1-n*Gc*z*z
    U=[np.array([0,0,1.])]
    for k in range(n):
        p=2*np.pi*k/n
        U.append(np.array([np.sqrt(1-z*z)*np.cos(p),np.sqrt(1-z*z)*np.sin(p),z]))
    return np.array([g0]+[Gc]*n), np.array(U)

astar=(3*S35+5*np.sqrt(7)-7-7*S5)/14; zstar=0.5004400497682123079777569143187991903444
rho_st=(3*astar-1)/2
print("=== 冠军核对 ===")
g6,U6=cone(5,zstar); d=fast(g6,U6,uniax(rho_st))
print(" J_w/2=%.15f rho0=%.15f diff=%.2e"%(d['half'],rho0,d['half']-rho0))
print(" Phi=%.15f  delta0/3=%.15f  diff=%.2e"%(d['Phi'],TGT,d['Phi']-TGT))
print(" 五项: (1-A^2)=%.9f  cW*H/2=%.9f  cS*Et/2=%.9f  |cZ|Z4=%.9f  Dw=%.9f"%(d['t1'],d['t2'],d['t3'],d['t4'],d['t5']))
print(" A_w=%.12f  1-E/2=%.12f  残差=%.2e ; E=%.12f"%(d['A'],1-d['E']/2,d['r_A'],d['E']))
print(" 残差: 主恒等=%.1e  H=6-2W^2:%.1e  SF=1-Et/2:%.1e  M4分解:%.1e"%(d['r_ID'],d['r_H'],d['r_E'],d['r_M']))

print("\n=== 锥族剖面: max J_w/2 与 min Phi (n=3..12) ===")
for n in range(3,13):
    bj=(-1,None); bp=(1e9,None)
    for a in np.linspace(1/3+1e-4,0.9999,240):
        r=(3*a-1)/2; L=uniax(r)
        for z in np.linspace(0.02,0.97,240):
            Gc=2/(n*(1-z*z)); g0=1-n*Gc*z*z
            if g0<=1e-8: continue
            g,U=cone(n,z); dd=fast(g,U,L)
            if dd['half']>bj[0]: bj=(dd['half'],(a,z))
            if dd['Phi']<bp[0]: bp=(dd['Phi'],(a,z))
    print("  n=%2d  max J/2=%.12f (a=%.4f z=%.4f, diff=%+.3e) | min Phi=%.12f (a=%.4f z=%.4f, diff=%+.3e)"%(
        n,bj[0],bj[1][0],bj[1][1],bj[0]-rho0,bp[0],bp[1][0],bp[1][1],bp[0]-TGT))
np.save('_r25_astar.npy',np.array([astar,zstar,rho_st]))

