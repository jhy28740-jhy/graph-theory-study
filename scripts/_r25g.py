# -*- coding: utf-8 -*-
"""R25-g: 冠军处 Phi 的严格局部极小性 + 高精度核对 + 窗口定理常数."""
import numpy as np
from numpy.linalg import norm
from mpmath import mp, mpf, sqrt as msqrt, polyroots
mp.dps=60
src=open('_r25e.py',encoding='utf-8-sig').read().split('astar=(3*S35')[0]
exec(src)

# ---- exact constants at 60 digits ----
r5=msqrt(5); r7=msqrt(7); r35=msqrt(35)
rho0_mp=(9+r5+2*r35)/18; phi_mp=(1+r5)/2
del0_mp=3*phi_mp+3-6*rho0_mp; TGT_mp=del0_mp/3
cW_mp=(2*r5-1)/14
astar_mp=(3*r35+5*r7-7-7*r5)/14
print("rho0   =",mp.nstr(rho0_mp,30))
print("delta0 =",mp.nstr(del0_mp,30)," delta0/3 =",mp.nstr(TGT_mp,30))
print("a*     =",mp.nstr(astar_mp,30))
# z* root of 11225 z^8 -6700 z^6 +430 z^4 +260 z^2 -31
co=[mpf(11225),0,-6700,0,430,0,260,0,-31]
rts=polyroots(co)
cand=[r for r in rts if abs(r.imag)<mpf('1e-30') and mpf('0.49')<r.real<mpf('0.51')]
zstar_mp=cand[0].real if cand else None
print("z*     =",mp.nstr(zstar_mp,30))
print("F(a*,z*)-rho0 check via 8次多项式: %.3e"%float(abs(11225*zstar_mp**8-6700*zstar_mp**6+430*zstar_mp**4+260*zstar_mp**2-31)))
# W_c 窗口常数
Wc2=3-TGT_mp/(3*cW_mp); Wc=msqrt(Wc2)
print("\n窗口定理常数: delta0/(3cW) =",mp.nstr(TGT_mp/(3*cW_mp),25))
print("  W_c =",mp.nstr(Wc,25),"  H_c = 2*(3-W_c^2) =",mp.nstr(2*(3-Wc2),25))
print("  cos(theta1)_c = W_c/sqrt3 =",mp.nstr(Wc/msqrt(3),25),"  theta1_c(deg) =",mp.nstr(float(np.degrees(float(np.arccos(float(Wc/msqrt(3)))))),12))
print("  sin^2(theta1) 窗口上界 =",mp.nstr((3-Wc2)/3,25))
astar_f=float(astar_mp); zstar_f=float(zstar_mp); rho_st=(3*astar_f-1)/2
print("\n冠军 a*=%.18f  z*=%.18f  rho*=%.18f"%(astar_f,zstar_f,rho_st))
g6,U6=cone(5,zstar_f); d=fast(g6,U6,uniax(rho_st))
print("Phi(champion)=%.18e   delta0/3=%.18e   diff=%.3e"%(d['Phi'],float(TGT_mp),d['Phi']-float(TGT_mp)))
print("J_w/2=%.18f  rho0=%.18f  diff=%.3e"%(d['half'],float(rho0_mp),d['half']-float(rho0_mp)))

# ---- strict local minimality of Phi at the champion ----
print("\n=== 冠军处 Phi 的严格局部极小性 (随机方向扰动) ===")
rng=np.random.default_rng(4242)
def Phi_of(V,Lam):
    g=(V*V).sum(1); U=V/np.sqrt(g)[:,None]; return fast(g,U,Lam)['Phi']
V6=np.sqrt(g6)[:,None]*U6
base=d['Phi']; bad=0; mins=[]
for eps in [1e-1,3e-2,1e-2,3e-3,1e-3,3e-4,1e-4]:
    mn=1e9; nbad=0; ntr=400
    for _ in range(ntr):
        P=rng.standard_normal((6,3)); Q,_=np.linalg.qr(np.eye(6)+eps*P)
        Vp=Q@V6
        M=rng.standard_normal((3,3)); O=np.linalg.qr(np.eye(3)+eps*M)[0]
        Lp=O@uniax(rho_st)@O.T
        v=Phi_of(Vp,Lp); mn=min(mn,v)
        if v<base-1e-12: nbad+=1
    bad+=nbad; mins.append((eps,mn,mn-base,nbad))
    print("  eps=%.0e : min Phi=%.15f  (Phi-base = %+.3e)  低于base的个数=%d/%d"%(eps,mn,mn-base,nbad,ntr))
print("  合计低于 base 的扰动: %d / %d"%(bad,len(mins)*400))
print("\n二阶增长估计: (Phi(eps)-Phi(0))/eps^2 =", ["%.4f"%((m[1]-base)/m[0]**2) for m in mins])
