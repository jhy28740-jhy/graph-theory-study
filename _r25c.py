# -*- coding: utf-8 -*-
"""R25-c: 通用五项恒等式 + A_w^2<=1 + 锥族/随机帧极值剖面 (修正 Lambda 构造)."""
import numpy as np
from numpy.linalg import eigh, norm

S5=np.sqrt(5.0); S35=np.sqrt(35.0)
rho0=(9+S5+2*S35)/18; phi=(1+S5)/2
gam=(1+S5)/16; alp=(7*S5-5)/8; bet=(25-15*S5)/16
cW=(2*S5-1)/14; cS=(10+S5)/14; cZ=(5-3*S5)/14
del0=3*phi+3-6*rho0
assert abs(2*rho0-phi-(1-del0/3))<1e-14
def P4(t): return (35*t**4-30*t*t+3)/8
def qfun(s): return gam+alp*s+bet*s*s

def uniax(rhoa):
    b=(1-rhoa)/3; L=b*np.eye(3); L[2,2]=b+rhoa; return L   # eigen (b,b,a), a=b+rho

def Aw(U,w,seed=0,nrand=40):
    n=len(w); rng=np.random.default_rng(seed); best=(0.0,None)
    starts=[U[i] for i in range(n)]+[rng.standard_normal(3) for _ in range(nrand)]
    starts += [np.array([0,0,1.]),np.array([1.,0,0]),np.array([0,1.,0])]
    for z0 in starts:
        z=z0/norm(z0)
        for _ in range(400):
            e=np.where(U@z>=0,1.0,-1.0); y=(e*w)@U; ny=norm(y)
            if ny<1e-14: break
            zn=y/ny
            if norm(zn-z)<1e-15: z=zn; break
            z=zn
        v=float(np.abs(U@z)@w)
        if v>best[0]: best=(v,z)
    return best

def full(g,U,Lam,seed=0):
    h=np.sqrt(np.einsum('ij,jk,ik->i',U,Lam,U)); w=g*h
    n=len(g); G=U@U.T; Wo=w.sum()
    Sw=(w[:,None,None]*(U[:,:,None]*U[:,None,:])).sum(0)
    Lw=(np.outer(w,w)*np.abs(G)).sum()
    A,z=Aw(U,w,seed=seed)
    a=np.abs(U@z)
    Vw=(np.outer(w,w)*(a[:,None]-a[None,:])**2).sum()
    zet=z@Sw@z
    M4=(np.outer(w,w)*G**4).sum(); Z4=(np.outer(w,w)*P4(G)).sum()
    Dw=(np.outer(w,w)*(qfun(G**2)-np.abs(G))).sum()
    H=(np.outer(g,g)*(h[:,None]-h[None,:])**2).sum()
    Et=(np.outer(g,g)*(h[:,None]-h[None,:])**2*G**2).sum()
    Jw=Lw+A*A
    term1=1-A*A; term2=0.5*(cW*H+cS*Et); term3=-cZ*Z4; term4=Dw
    Phi=term1+term2+term3+term4
    return dict(h=h,w=w,W=Wo,Sw=Sw,nrmSF=(Sw*Sw).sum(),lam=eigh(Sw)[0],zet=zet,z=z,
                Lw=Lw,A=A,Aw2=A*A,Vw=Vw,M4=M4,Z4=Z4,Dw=Dw,H=H,Et=Et,Jw=Jw,half=Jw/2,
                t1=term1,t2=term2,t3=term3,t4=term4,Phi=Phi,
                r_M4=M4-(-3*Wo*Wo/35+6*(Sw*Sw).sum()/7+8*Z4/35),
                r_W2=Wo*Wo-(3-H/2), r_SF=(Sw*Sw).sum()-(1-Et/2),
                r_Aw2=A*A-(Wo*zet-0.5*Vw),
                r_Lw=Lw-(gam*Wo*Wo+alp*(Sw*Sw).sum()+bet*M4-Dw),
                r_ID=(2*rho0-Jw)-((1-A*A)-del0/3+0.5*(cW*H+cS*Et)+(-cZ)*Z4+Dw))

# ---- champion ----
astar=(3*S35+5*np.sqrt(7)-7-7*S5)/14; zstar=0.5004400497682123079777569143187991903444
rho_st=(3*astar-1)/2; Gc=2/(5*(1-zstar**2)); gh=1-5*Gc*zstar**2
def cone(n,z):
    G=2/(n*(1-z*z)); g0=1-n*G*z*z
    U=[np.array([0,0,1.])]
    for k in range(n):
        p=2*np.pi*k/n
        U.append(np.array([np.sqrt(1-z*z)*np.cos(p),np.sqrt(1-z*z)*np.sin(p),z]))
    return np.array([g0]+[G]*n), np.array(U)
g6,U6=cone(5,zstar)
d=full(g6,U6,uniax(rho_st))
print("=== 冠军 (5-锥, a=a*) ===")
print(" J_w/2 = %.15f   rho0 = %.15f   diff = %.3e"%(d['half'],rho0,d['half']-rho0))
print(" 五项: 1-A_w^2=%.12f  (cW*H+cS*Et)/2=%.12f  |cZ|Z4=%.12f  D_w=%.12f"%(d['t1'],d['t2'],d['t3'],d['t4']))
print(" Phi = %.15f    delta0/3 = %.15f    Phi-delta0/3 = %.3e"%(d['Phi'],del0/3,d['Phi']-del0/3))
print(" A_w^2=%.12f  W=%.12f  H=%.12f  Et=%.12f  V_w=%.12f"%(d['Aw2'],d['W'],d['H'],d['Et'],d['Vw']))
print(" 恒等式残差: M4=%.1e W2=%.1e SF=%.1e Aw2=%.1e Lw=%.1e 主恒等=%.1e"%(
    d['r_M4'],d['r_W2'],d['r_SF'],d['r_Aw2'],d['r_Lw'],d['r_ID']))
print(" eig(S_w)=",np.round(d['lam'],9))

# ---- identity stress test on random frames & random Lambda ----
print("\n=== 通用恒等式压力测试 (随机帧 x 随机 Lambda, 含非单轴) ===")
rng=np.random.default_rng(7)
mx={k:0 for k in ['r_M4','r_W2','r_SF','r_Aw2','r_Lw','r_ID']}
maxAw2=0; minPhi=1e9; argmin=None
for it in range(400):
    N=int(rng.integers(4,11)); V,_=np.linalg.qr(rng.standard_normal((N,3)))
    g=(V*V).sum(1); U=V/np.sqrt(g)[:,None]
    M=rng.standard_normal((3,3)); M=M@M.T; ev=np.sort(rng.random(3))[::-1]; ev=ev/ev.sum()
    Om=eigh(M)[1]; Lam=Om@np.diag(ev)@Om.T
    dd=full(g,U,Lam,seed=it)
    for k in mx: mx[k]=max(mx[k],abs(dd[k]))
    maxAw2=max(maxAw2,dd['Aw2'])
    if dd['Phi']<minPhi: minPhi=dd['Phi']; argmin=(N,ev,dd['half'])
print(" 最大残差:",{k:"%.2e"%v for k,v in mx.items()})
print(" max A_w^2 over 400 random = %.12f  (应 <=1)"%maxAw2)
print(" min Phi = %.12f  vs delta0/3 = %.12f ; at %s"%(minPhi,del0/3,argmin))
