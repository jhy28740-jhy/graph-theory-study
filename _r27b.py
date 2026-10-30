# -*- coding: utf-8 -*-
# _r27b.py  第 27 轮 B 段：钝角余弦平方质量 P 的 Fourier 上界（这是本轮的主定理）
#
#   P(mu) = \iint (cos(phi-phi'))_-^2 dmu dmu'
#   约束:  mu>=0,  mu(S^1)=2,  \int e^{2i phi} dmu = 0   (二维各向同性 sum rho_i p_i p_i^T = I_2)
#
#   K(u)=(cos u)_-^2 的 Fourier 系数（解析推导）:
#       Khat_0 = 1/4 ;  Khat_2 = 1/8 ;  Khat_k = 0 (k 偶, k>=4) ;
#       Khat_k = 2 sin(k pi/2) / (pi k (k^2-4))   (k 奇)
#   =>  P = 4 Khat_0 + 2 sum_{k>=1} Khat_k |muhat(k)|^2
#         = 1 + 2 sum_{k>=1} Khat_k |muhat(k)|^2          (Khat_2 项被 muhat(2)=0 杀掉)
#   用 |muhat(k)| <= muhat(0) = 2 且丢掉 Khat_k<0 的项:
#       P <= 1 + 8 * sum_{k>=1 odd, Khat_k>0} Khat_k
#         = 1 + (16/pi) * S,     S = sum_{m>=1} 1/((4m+1)((4m+1)^2-4))
import numpy as np
from mpmath import mp, mpf, pi, sqrt, nsum, inf, pslq
mp.dps = 60
rng = np.random.default_rng(7)

print("="*100); print("B1  Fourier 系数 Khat_k 的解析式 vs 数值积分"); print("="*100)
def K(u): return np.maximum(np.cos(u),0.0)**0 * np.where(np.cos(u)<0, np.cos(u)**2, 0.0)
def Khat_num(k, M=4000000):
    u=np.linspace(0,2*np.pi,M,endpoint=False); du=2*np.pi/M
    return float((K(u)*np.exp(-1j*k*u)).sum()*du/(2*np.pi))
def Khat_for(k):
    if k==0: return mpf(1)/4
    if k==2: return mpf(1)/8
    if k%2==0: return mpf(0)
    return 2*mp.sin(k*pi/2)/(pi*k*(k**2-4))
for k in [0,1,2,3,4,5,6,7,8,9]:
    a=Khat_num(k); b=float(Khat_for(k))
    print("  k=%d  数值=%+.12f  解析=%+.12f  差=%.2e"%(k,a,b,abs(a-b)))

print()
print("="*100); print("B2  Fourier 恒等式 P = sum_k Khat_k |muhat(k)|^2 的数值核验"); print("="*100)
def P_direct(rho,ang):
    c=np.cos(ang[:,None]-ang[None,:])
    M=(rho[:,None]*rho[None,:])*np.where(c<0,c**2,0.0)
    return float(M.sum())
def P_fourier(rho,ang,K=400):
    mh=np.array([np.sum(rho*np.exp(-1j*k*ang)) for k in range(-K,K+1)])
    ks=np.arange(-K,K+1)
    kh=np.array([float(Khat_for(abs(int(k)))) for k in ks])
    return float(np.sum(kh*np.abs(mh)**2))
for trial in range(12):
    N=int(rng.integers(3,13))
    rho=rng.random(N)*0.8; ang=rng.random(N)*2*np.pi
    print("  N=%2d  P_direct=%.12f  P_fourier=%.12f  差=%.2e"%(N,P_direct(rho,ang),P_fourier(rho,ang),abs(P_direct(rho,ang)-P_fourier(rho,ang))))

print()
print("="*100); print("B3  常数 S 的闭式与最终上界"); print("="*100)
S_num = nsum(lambda m: 1/((4*m+1)*((4*m+1)**2-4)), [1, inf])
print("  S = sum_{m>=1} 1/((4m+1)((4m+1)^2-4)) = %s"%mp.nstr(S_num,40))
S_part = sum(mpf(1)/((4*m+1)*((4*m+1)**2-4)) for m in range(1,200000))
print("  部分和(m<2e5) = %s   与 nsum 差 = %s"%(mp.nstr(S_part,30),mp.nstr(abs(S_part-S_num),3)))
cand = pslq([S_num, mpf(1), pi, sqrt(2), sqrt(3), mp.ln(2), mp.ln(3), mpf(1)/mpf(2)])
print("  pslq 猜测:",cand)
Pbound = 1 + 16*S_num/pi
print("  P <= 1 + 16S/pi = %s"%mp.nstr(Pbound,30))
rho0 = (9+sqrt(5)+2*sqrt(35))/18
Neg_bound = Pbound/4
print("  Neg <= P/4      = %s"%mp.nstr(Neg_bound,30))
print("  rho_0 - 1       = %s"%mp.nstr(rho0-1,30))
print("  余量 (rho_0-1) - Neg_bound = %s"%mp.nstr((rho0-1)-Neg_bound,6))
print("  kappa_3(秩一) <= 1 + Neg_bound = %s   vs rho_0 = %s"%(mp.nstr(1+Neg_bound,20),mp.nstr(rho0,20)))
print("  L_w <= 1 + 2 Neg_bound = %s   vs (13+sqrt5)/10 = %s"%(mp.nstr(1+2*Neg_bound,20),mp.nstr((13+sqrt(5))/10,20)))
print("  是否关闭秩一截面: ", 1+Neg_bound < rho0)
print("  对比 STW: Phi_STW = %s ; 本界 = %s ; 改进 = %s"
      %(mp.nstr((3+sqrt(5))/4,20), mp.nstr(1+Neg_bound,20), mp.nstr((3+sqrt(5))/4-(1+Neg_bound),6)))
print("  对比锥族(条件性): (11+phi)/10 = %s ; 本界 - 锥族 = %s"
      %(mp.nstr((11+(1+sqrt(5))/2)/10,20), mp.nstr((1+Neg_bound)-(11+(1+sqrt(5))/2)/10,6)))

print()
print("B1--B3 DONE (B4 moved to _r27c.py)")
