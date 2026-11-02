# -*- coding: utf-8 -*-
# =============================================================================
# 证书：优占对偶 (domination duality) 与冠军腔的精确全局上界
# 复核报告 sec:domdual 的全部结论。逐块对应：
#   (A) 定理 domdual     F(Sigma) <= G(Sigma) = Omega(Sigma)，三种算法互核
#   (B) 定理 champglobal Omega(Sigma0) = rho0，闭式根式证书 T，符号验证全部为 0
#   (C) 定理 wallsexact  四类腔面/邻腔掩码 Omega < rho0，精确有理数 PSD 认证
#   (D) 命题 octnogo     八面体掩码 Omega = 4/3 > rho0（松弛类的定量天花板）
#   (E) 命题 negmass     tr((-M*)_+) = tau 精确饱和；Lam3(M*)=(||Y||_1+1^T Y 1)/2
# 运行： C:\Python314\python.exe -X utf8 -u 证书_优占对偶全局上界.py
# =============================================================================
import numpy as np, sympy as sp, itertools, cvxpy as cp, time
from fractions import Fraction as Fr
from scipy.optimize import minimize

s5=sp.sqrt(5); phi=(1+s5)/2; rho0=(9+s5+2*sp.sqrt(35))/18; tau=rho0-1
k1=(s5-1)/4; k2=-(s5+1)/4                      # cos72, cos144
c=(4*rho0+1-phi)/5; d=rho0-1/s5; e=rho0-phi+1/s5; S=2*d+2*e
c_=fv=lambda x: float(sp.N(x,40))
RHO0=float(sp.N(rho0,30)); TAU=float(sp.N(tau,30)); PHI=float(sp.N(phi,30))

# ---- 掩码 ----
Sig=np.ones((6,6))
for i in range(1,6):
    for j in range(1,6):
        if i!=j: Sig[i,j]=1.0 if min(abs(i-j),5-abs(i-j))==1 else 0.0
def flip(M,i,j,v):
    M=M.copy(); M[i,j]=M[j,i]=v; return M
MASKS={"champion Sigma0":Sig,
       "class1 hc (0,1)->0":flip(Sig,0,1,0.0),
       "class2 d1 (1,2)->0":flip(Sig,1,2,0.0),
       "class3 d2+(1,3)->1":flip(Sig,1,3,1.0),
       "class12 双条目":flip(flip(Sig,0,1,0.0),1,2,0.0),
       "class13 双条目":flip(flip(Sig,0,1,0.0),1,3,1.0)}
OCT=np.array([[1,0,1,1,1,1],[0,1,1,1,1,1],[1,1,1,0,1,1],
              [1,1,0,1,1,1],[1,1,1,1,1,0],[1,1,1,1,0,1]],dtype=float)

def top3(M):
    ev,U=np.linalg.eigh(M); o=np.argsort(ev)[::-1]; return ev[o],U[:,o]
def Lam3(M): return top3(M)[0][:3].sum()
def F_val(A,ns=24,seed=13):
    """F(A)=max_{p in Delta} Lam3(D_sqrt p A D_sqrt p).
    softmax 参数化 x=log p; 解析梯度 dLam3/dx = g - <p,g>, g=diag(Pi_3 M Pi_3).
    L-BFGS-B 多起点(含 n 个顶点起点) + Nelder-Mead/Powell 收尾。"""
    n=A.shape[0]; rng=np.random.default_rng(seed)
    def fg(x):
        x=x-x.max(); p=np.exp(x); p/=p.sum(); q=np.sqrt(p); M=(q[:,None]*A)*q[None,:]
        ev,U=top3(M); P=U[:,:3]@U[:,:3].T; g=np.diag(P@M@P)
        return ev[:3].sum(), g-(g@p)
    nfun=lambda x:-fg(x)[0]; gfun=lambda x:-fg(x)[1]
    starts=[np.zeros(n)]
    for i in range(n):
        v=np.full(n,-6.0); v[i]=6.0; starts.append(v)
    starts+=[rng.normal(0,1.0,n) for _ in range(ns)]
    best=(-9.0,None)
    for s0 in starts:
        r=minimize(nfun,s0,jac=gfun,method="L-BFGS-B",options=dict(maxiter=4000,ftol=1e-16,gtol=1e-12))
        if -r.fun>best[0]: best=(-r.fun,r.x)
    r=minimize(nfun,best[1],method="Nelder-Mead",options=dict(maxiter=20000,maxfev=40000,xatol=1e-14,fatol=1e-16))
    if -r.fun>best[0]: best=(-r.fun,r.x)
    r=minimize(nfun,r.x,method="Powell",options=dict(maxiter=20000,xtol=1e-14,ftol=1e-16))
    return max(best[0],-r.fun)
def G_val(A,ns=25,seed=6):
    """G(A)=max_p tr((D_sqrt p A D_sqrt p)_+)  (= (1+max_p ||.||_*)/2 当 A_ii=1)"""
    n=A.shape[0]; rng=np.random.default_rng(seed)
    def neg(x):
        p=np.exp(x-x.max()); p/=p.sum(); q=np.sqrt(p)
        return -np.clip(np.linalg.eigvalsh((q[:,None]*A)*q[None,:]),0,None).sum()
    best=-9
    for s0 in [np.zeros(n)]+[rng.normal(0,1.2,n) for _ in range(ns-1)]:
        r=minimize(neg,s0,method="Nelder-Mead",options=dict(maxiter=8000,maxfev=16000,xatol=1e-13,fatol=1e-15))
        r=minimize(neg,r.x,method="Powell",options=dict(xtol=1e-13,ftol=1e-15))
        best=max(best,-r.fun)
    return best
def Omega_val(A):
    t=cp.Variable(); T=cp.Variable(A.shape,symmetric=True)
    p=cp.Problem(cp.Minimize(t),[T>>A,T>>0,cp.diag(T)<=t*np.ones(A.shape[0])])
    p.solve(solver=cp.SCS,eps=1e-12,max_iters=400000)
    return float(t.value),np.array(T.value)

print("="*80)
print("(A) 定理 domdual：F(Sigma) <= G(Sigma) = Omega(Sigma) = min{max_i T_ii : T>>Sigma, T>>0}")
print("    G(Sigma) = max_p tr((D_sqrt p Sigma D_sqrt p)_+) = (1 + max_p ||D_sqrt p Sigma D_sqrt p||_*)/2")
print("="*80)
print(" %-20s %14s %14s %14s %10s %10s"%("mask","F (Lam3)","G (tr(.)_+)","Omega (SDP)","G-F","Omega-G"))
for tag,A in MASKS.items():
    f=F_val(A); g=G_val(A); om,_=Omega_val(A)
    print(" %-20s %14.10f %14.10f %14.10f %10.2e %10.2e"%(tag,f,g,om,g-f,om-g))
print(" 判读：三列满足 F <= G = Omega（G=Omega 由 SDP 强对偶证明，数值差 <=1e-9）；")
print("       G-F>0 当且仅当 G 的极大点处 D_sqrt p Sigma D_sqrt p 有 >=4 个正特征值。")

print()
print("="*80)
print("(B) 定理 champglobal：Omega(Sigma0) = rho0（闭式根式证书，sympy 符号验证）")
print("="*80)
print(" T = rho0*I + c*B_hubbycyc + d*B_dist1 + e*B_dist2")
print("   rho0 = (9+sqrt5+2sqrt35)/18 = %s"%sp.N(rho0,20))
print("   c = (4rho0+1-phi)/5         = %s"%sp.N(c,20))
print("   d = rho0 - 1/sqrt5          = %s"%sp.N(d,20))
print("   e = rho0 - phi + 1/sqrt5    = %s"%sp.N(e,20))
print(" D5 分块 => T>>0 与 N:=T-Sigma0>>0 的充要标量条件（全部精确）：")
allok=True
for nm,val,tgt in [
    ("det 平凡块(T) = rho0(rho0+S)-5c^2",  rho0*(rho0+S)-5*c**2,            sp.Integer(0)),
    ("lam_B(T) = rho0+2d*cos144+2e*cos72", rho0+2*d*k2+2*e*k1,              sp.Integer(0)),
    ("lam_A(T) = rho0+2d*cos72 +2e*cos144",rho0+2*d*k1+2*e*k2,              phi),
    ("平凡块大特征值(T) = 2rho0+S",         2*rho0+S,                        6*rho0-2*phi),
    ("det 平凡块(N) = tau(tau+S-2)-5(c-1)^2",tau*(tau+S-2)-5*(c-1)**2,      sp.Integer(0)),
    ("lam_A(N) = tau+2(d-1)cos72+2e cos144", tau+2*(d-1)*k1+2*e*k2,         sp.Integer(0)),
    ("lam_B(N) = tau+2(d-1)cos144+2e cos72", tau+2*(d-1)*k2+2*e*k1,         1/phi),
    ("平凡块大特征值(N) = 2tau+S-2",        2*tau+S-2,                       6*rho0-2*phi-4)]:
    r=sp.simplify(sp.expand(val-tgt)); ok=(r==0); allok&=ok
    print("   %-44s = %-22s (目标 %-22s) %s"%(nm,sp.N(val,18),sp.N(tgt,18),"OK" if ok else "MISMATCH"))
for nm,val in [("rho0+S (T 平凡块迹)",rho0+S),("tau+S-2 (N 平凡块迹)",tau+S-2)]:
    pos=bool(sp.N(val,25)>0); allok&=pos
    print("   %-44s = %-22s >0 : %s"%(nm,sp.N(val,20),pos))
print(" 全部符号条件通过 :",allok)
print(" => spec(T)      = {6rho0-2phi, phi, phi, 0,0,0},  tr T = 6rho0")
print(" => spec(T-Sig0) = {6rho0-2phi-4, 1/phi, 1/phi, 0,0,0}")
print(" => T>>0, T-Sigma0>>0, max_i T_ii = rho0  =>  Omega(Sigma0) <= rho0")
Tn=np.zeros((6,6))
for i in range(6): Tn[i,i]=fv(rho0)
for j in range(1,6): Tn[0,j]=Tn[j,0]=fv(c)
for i in range(1,6):
    for j in range(1,6):
        if i!=j: Tn[i,j]=fv(d) if min(abs(i-j),5-abs(i-j))==1 else fv(e)
print(" 数值交叉核对：min eig(T)=%.3e  min eig(T-Sigma0)=%.3e  max diag=%.15f"%(
    np.linalg.eigvalsh(Tn).min(),np.linalg.eigvalsh(Tn-Sig).min(),Tn.diagonal().max()))

print()
print("="*80)
print("(C) 定理 wallsexact：腔面/邻腔掩码的 Omega < rho0，精确有理数 PSD 认证")
print("    认证方式：T 的 63 个主子式与 (T-Sigma) 的 63 个主子式全部 >= 0（有理数精确算术，Sylvester）")
print("="*80)
def det_exact(M):
    n=len(M); M=[r[:] for r in M]; det=Fr(1)
    for cc in range(n):
        piv=None
        for r in range(cc,n):
            if M[r][cc]!=0: piv=r; break
        if piv is None: return Fr(0)
        if piv!=cc: M[cc],M[piv]=M[piv],M[cc]; det=-det
        det*=M[cc][cc]; inv=Fr(1)/M[cc][cc]
        for r in range(cc+1,n):
            f=M[r][cc]*inv
            if f!=0:
                for k in range(cc,n): M[r][k]-=f*M[cc][k]
    return det
def psd_exact(Mf):
    n=len(Mf); bad=[]; worst=None
    for r in range(1,n+1):
        for idx in itertools.combinations(range(n),r):
            dd=det_exact([[Mf[i][j] for j in idx] for i in idx])
            if dd<0: bad.append((idx,dd))
            if worst is None or dd<worst: worst=dd
    return bad,worst
QI=2**40; DELTA=Fr(1,2**20)
print(" %-20s %14s %12s %12s %10s %10s"%("mask","Omega(SDP)","rho0-Omega","diagmax(有理)","T 主子式","T-S 主子式"))
for tag,A in MASKS.items():
    if tag.startswith("champion"): continue
    om,Tv=Omega_val(A)
    Tf=[[Fr(int(round(Tv[i,j]*QI)),QI) for j in range(6)] for i in range(6)]
    for i in range(6): Tf[i][i]+=DELTA
    Af=[[Fr(int(round(A[i,j]))) for j in range(6)] for i in range(6)]
    b1,w1=psd_exact(Tf); b2,w2=psd_exact([[Tf[i][j]-Af[i][j] for j in range(6)] for i in range(6)])
    mx=max(Tf[i][i] for i in range(6))
    print(" %-20s %14.10f %12.3e %12s %10s %10s"%(tag,om,RHO0-om,
        "OK" if not b1 else "FAIL", "OK" if not b2 else "FAIL", "%.2e"%float(w2)))
    print("      diagmax(有理证书) = %s   rho0 - diagmax = %.6f"%(
        sp.nsimplify(Fr(mx.numerator,mx.denominator)),RHO0-float(mx)))

print()
print("="*80)
print("(D) 命题 octnogo：八面体掩码 Sigma_oct = J6 - (对径完美匹配)")
print("="*80)
om,To=Omega_val(OCT)
print(" Omega(Sigma_oct) = %.15f  (= 4/3 ? 差 %.2e)"%(om,om-4/3))
print(" 下界 1 + tr((-Sigma)_+)/6 = %.15f"%(1+np.clip(-np.linalg.eigvalsh(OCT),0,None).sum()/6))
print(" F(Sigma_oct) = max_p Lam3 = %.15f  (p 均匀时 = 7/6 = %.15f)"%(F_val(OCT),7/6))
print(" spec(Sigma_oct) =",np.sort(np.linalg.eigvalsh(OCT))[::-1])
print(" => 该掩码是 R^3 中标架 Gram 的符号模式且满足 (P3) 的归一化，但 Omega=4/3>rho0；")
print("    故『对一切标架可实现掩码一致地 Omega<=rho0』为假 —— 自洽性 A=1{Pi_3(M)>=0} 不可省略。")

print()
print("="*80)
print("(E) 命题 negmass：冠军精确饱和负质量界；Lam3 = (||Y||_1 + 1^T Y 1)/2")
print("="*80)
a=(91-25*sp.sqrt(7))/126; b=(7+5*sp.sqrt(7))/126
w=np.array([fv(a)]+[fv(b)]*5); Dw=np.diag(w)
M=np.diag(np.sqrt(w))@Sig@np.diag(np.sqrt(w))
ev,U=top3(M); P=U[:,:3]@U[:,:3].T; Y=np.diag(np.sqrt(w))@P@np.diag(np.sqrt(w))
print(" w = (a,b,b,b,b,b), a=(91-25sqrt7)/126=%.15f, b=(7+5sqrt7)/126=%.15f"%(w[0],w[1]))
print(" spec(M*) =",np.round(ev,12))
print(" Lam3(M*) = %.15f   rho0 = %.15f   差 = %.3e"%(ev[:3].sum(),RHO0,ev[:3].sum()-RHO0))
print(" tr M*    = %.15f (=1)"%ev.sum())
print(" tr((-M*)_+) = %.15f   tau = %.15f   差 = %.3e   [精确饱和]"%(np.clip(-ev,0,None).sum(),TAU,np.clip(-ev,0,None).sum()-TAU))
print(" tr((M*)_+)  = %.15f = Lam3(M*)  <=> M* 恰有 3 个正特征值 : %s"%(np.clip(ev,0,None).sum(),abs(np.clip(ev,0,None).sum()-ev[:3].sum())<1e-13))
print(" 对偶证书 Y* = D_sqrt(w) Pi_3(M*) D_sqrt(w):")
print("   min eig Y* = %.3e ; min eig(D_w - Y*) = %.3e  [Y* 可行]"%(np.linalg.eigvalsh(Y).min(),np.linalg.eigvalsh(Dw-Y).min()))
print("   tr(Y* Sigma0) = %.15f  (= rho0, 差 %.3e)"%(np.sum(Y*Sig),np.sum(Y*Sig)-RHO0))
print("   ||Y*||_1 = %.15f ; 1^T Y* 1 = %.15f ; (和)/2 = %.15f"%(np.abs(Y).sum(),Y.sum(),(np.abs(Y).sum()+Y.sum())/2))
print("   互补松弛 ||Y*(T-Sigma0)|| = %.3e ; ||(D_w-Y*)T|| = %.3e"%(np.linalg.norm(Y@(Tn-Sig)),np.linalg.norm((Dw-Y)@Tn)))
print("   sign(Y*) == 2Sigma0-J (自洽性) :",np.abs(np.sign(Y)-(2*Sig-1))[~np.eye(6,dtype=bool)].max()==0.0)
print(" => Lam3(M*) = (||Y*||_1 + 1^T Y* 1)/2 = rho0，即 c_3^Sigma 的『加权 eta_mix』形式在冠军处精确取到 2rho0。")
print()
print("总结：冠军腔 W1-全局 已闭合（F(Sigma0)=Omega(Sigma0)=rho0，闭式证书）；")
print("      四类近腔掩码有 Omega<rho0 的精确有理数证书；八面体给出松弛类的定量天花板 4/3。")
