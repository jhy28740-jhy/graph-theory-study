# -*- coding: utf-8 -*-
# _r27a.py  第 27 轮 A 段：秩一截面的新恒等式链与「钝角余弦平方质量」P 的上确界探索
#   链条：  Neg = sum_{N}(-sigma_i sigma_j - r_i r_j c_ij)
#         <= (1/4) sum_{N} rho_i rho_j c_ij^2 =: Theta/4        [Young, 等号 iff kappa=2]
#         <= (1/4) sum_{i!=j, c_ij<0} rho_i rho_j c_ij^2 =: P/4  [N ⊂ {c<0}, 等号 iff N={c<0}]
#   刚性恒等式： sum_ij rho_i rho_j c_ij^2 = ||I_2||_F^2 = 2
import numpy as np, itertools
np.set_printoptions(precision=12, suppress=False)
sqrt5=np.sqrt(5.0); phi=(1+sqrt5)/2
rho0=(9+sqrt5+2*np.sqrt(35))/18
rng=np.random.default_rng(20261030)

def champ():
    """冠军 W_6 = K_1 v C_5 的秩一数据：极点 + 正五边形纬度圈，s=s0=phi/(4+phi)。"""
    s0=phi/(4+phi); g0=(1-3*s0)/(1-s0); g1=2/(5*(1-s0))
    sig=[g0*1.0]+[g1*s0]*5
    rh =[g0*0.0]+[g1*(1-s0)]*5
    ang=[0.0]+[2*np.pi*k/5 for k in range(5)]
    return np.array(sig),np.array(rh),np.array(ang)

def quantities(sig,rh,ang):
    sig=np.asarray(sig,float); rh=np.asarray(rh,float); ang=np.asarray(ang,float)
    r=np.sqrt(sig*rh); lam=np.where(rh>0, np.sqrt(np.divide(sig,rh,out=np.zeros_like(sig),where=rh>0)), np.inf)
    c=np.cos(ang[:,None]-ang[None,:])
    t=sig[:,None]*sig[None,:]+r[:,None]*r[None,:]*c
    Nmask=(t<0)
    Neg=float((-t[Nmask]).sum())
    Lw=float(np.abs(t).sum())
    Theta=float((rh[:,None]*rh[None,:]*c*c)[Nmask].sum())
    P=float((rh[:,None]*rh[None,:]*c*c)[(c<0)&(~np.eye(len(sig),dtype=bool))].sum())
    kappa=np.where(Nmask&(sig[:,None]*sig[None,:]>0), np.abs(c)/np.maximum(lam[:,None]*lam[None,:],1e-300), np.nan)
    defect=float(np.nansum((sig[:,None]*sig[None,:])*(kappa-2)**2))
    return dict(Neg=Neg,Lw=Lw,Theta=Theta,P=P,defect=defect,t=t,c=c,r=r,lam=lam,sig=sig,rh=rh,Nmask=Nmask)

print("="*100)
print("A1  冠军处的恒等式核验（N=6：极点 + 正五边形，s0=phi/(4+phi)）")
print("="*100)
sig,rh,ang=champ()
q=quantities(sig,rh,ang)
print("  sigma =",np.array2string(sig,precision=12))
print("  rho   =",np.array2string(rh ,precision=12))
print("  sum sigma = %.15f (应为 1)     sum rho = %.15f (应为 2)"%(sig.sum(),rh.sum()))
print("  sum rho e^{2i phi} = %.3e (应为 0)"%abs((rh*np.exp(2j*ang)).sum()))
print("  sum r e^{i phi}    = %.3e (应为 0)"%abs((q['r']*np.exp(1j*ang)).sum()))
print("  Neg              = %.15f   phi^2/10 = %.15f   差 %.3e"%(q['Neg'],phi**2/10,q['Neg']-phi**2/10))
print("  L_w = sum|t|     = %.15f   (13+sqrt5)/10 = %.15f   差 %.3e"%(q['Lw'],(13+sqrt5)/10,q['Lw']-(13+sqrt5)/10))
print("  L_w - 1 - 2Neg   = %.3e"%abs(q['Lw']-1-2*q['Neg']))
print("  sum t            = %.15f (应为 1)"%q['t'].sum())
print("  刚性恒等式 sum_ij rho_i rho_j c_ij^2 = %.15f (应为 2)"%(rh[:,None]*rh[None,:]*q['c']**2).sum())
print("  Theta            = %.15f   Theta/4 = %.15f"%(q['Theta'],q['Theta']/4))
print("  P (全部钝角对)   = %.15f   P/4     = %.15f"%(q['P'],q['P']/4))
print("  Theta == P ?     差 %.3e   （冠军处 N = {c<0}，故应相等）"%abs(q['Theta']-q['P']))
print("  kappa-2 亏格项   = %.3e   （冠军处 kappa==2，故应为 0）"%q['defect'])
print("  恒等式 Neg = Theta/4 - defect/4 : 残差 %.3e"%abs(q['Neg']-(q['Theta']-q['defect'])/4))
print("  Young 链:  Neg <= Theta/4 <= P/4  ->  %.15f <= %.15f <= %.15f  (两处均取等)"%(q['Neg'],q['Theta']/4,q['P']/4))
print("  目标 rho_0-1 = %.15f ; 冠军余量 = %.6e"%(rho0-1,rho0-1-q['Neg']))
print("  关闭秩一截面所需:  P <= 4(rho_0-1) = %.15f"%(4*(rho0-1)))
print("  冠军 P 与所需上界的余量 = %.6e"%(4*(rho0-1)-q['P']))

print()
print("="*100)
print("A2  随机核验：Neg = Theta/4 - (1/4)sum_N sigma_i sigma_j (kappa_ij-2)^2  与  sum rho rho c^2 = 2")
print("="*100)
mx1=mx2=mx3=0.0; worst=None
for trial in range(4000):
    N=rng.integers(4,11)
    # 随机 Stiefel 帧 V in St(N,3)（V^T V = I_3），半球翻转后 v_{i,1}>=0
    A=rng.standard_normal((N,3)); Qm,_=np.linalg.qr(A); V=Qm[:,:3].copy()
    V[V[:,0]<0]*=-1
    g=(V**2).sum(1); sig=V[:,0]**2; rh=g-sig
    ok=(rh>1e-12)&(sig>1e-12)
    if ok.sum()<3: continue
    p=V[:,1:3]/np.sqrt(np.maximum((V[:,1:3]**2).sum(1,keepdims=True),1e-300))
    ang=np.arctan2(p[:,1],p[:,0])
    q=quantities(sig,rh,ang)
    mx1=max(mx1,abs(q['Neg']-(q['Theta']-q['defect'])/4))
    mx2=max(mx2,abs((rh[:,None]*rh[None,:]*q['c']**2).sum()-2.0))
    mx3=max(mx3,q['Neg']-q['P']/4+1e-15)
print("  kappa-2 恒等式最大残差      = %.3e"%mx1)
print("  刚性恒等式 sum rho rho c^2=2 最大残差 = %.3e"%mx2)
print("  Neg <= P/4 的最大违例        = %.3e  (应 <=0)"%mx3)

print()
print("="*100)
print("A3  新引理 ||Z_S||^2 <= sigma_S 对所有子集 S（=> A_w = 1 的一行证明）")
print("="*100)
def check_ZS(N=8,trials=300):
    bad=0; mxmargin=-1e9
    for _ in range(trials):
        A=rng.standard_normal((N,3)); Qm,_=np.linalg.qr(A); V=Qm[:,:3]
        if V[:,0].sum()<0: V=-V
        sig=V[:,0]**2; z=V[:,0][:,None]*V
        assert np.allclose(z.sum(0),[1,0,0],atol=1e-12), z.sum(0)
        Mi=np.diag(1/sig)
        assert np.allclose(z.T@Mi@z,np.eye(3),atol=1e-10)
        for S in itertools.product([0,1],repeat=N):
            idx=np.array(S,dtype=bool)
            if idx.sum()==0: continue
            ZS=z[idx].sum(0); sS=sig[idx].sum()
            m=(ZS**2).sum()-sS
            if m>1e-12: bad+=1
            mxmargin=max(mxmargin,m)
    return bad,mxmargin
bad,mx=check_ZS(7,120)
print("  N=7, 120 个随机 Stiefel 帧, 全部 2^7-1 个子集: 违例 %d 个, 最大 margin = %.3e (应 <=0)"%(bad,mx))
print("  推论: max_eps ||sum eps_i z_i||^2 = 4||Z_S||^2-4sigma_S+1 <= 1, 而全 + 号给出 =1,")
print("        故 A_w = 1（秩一刚性 prop:r25rank1 的新证明，只用 Cauchy-Schwarz 与 sum_S v_i v_i^T <= I_3）。")
print("ALL A-SECTION CHECKS DONE")
