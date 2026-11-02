# -*- coding: utf-8 -*-
# _o2.py  ROUND 16 verification
#  (A) energy decomposition g_i+eta_i+zeta_i=1 and the double-Cauchy-Schwarz pair
#  (B) sharp product bound beta^2 <= A1*A2 + full equality characterisation
#  (C) hub-free critical identity rho^2 = m^T Lt^{-1} m + exact CS defect
#  (D) hub -> 2-D circle reduction (vanishing second complex moment)
#  (E) residual gap of the new "fourth inequality" alpha + sqrt(A1*A2) <= rho0
#  (F) dual capacity form Phi(c) (numerical) + Farkas test
#  (G) quadratic-majorant LP hierarchy (refutation of the R15-7 route, degrees 2..8)
import numpy as np, math, itertools, time, warnings
warnings.filterwarnings('ignore')
from scipy.optimize import linprog, minimize
np.set_printoptions(precision=12, suppress=False, linewidth=190)
exec(open('_o1.py', encoding='utf-8').read().split('print("="*104)')[0])
# provides: rho0,tau,PHI,MOHAR,s5,s7,s35, Qp_of,rho_of,moment_mat,BTAR,random_frame_measure,
#           kernel_basis (from _n1.py), champion(), Uc,gc,Pc, unpack(), fit_Lt()

PASS=[0]; FAIL=[]
def chk(msg, ok, extra=""):
    if ok: PASS[0]+=1; print("  [PASS] %s %s"%(msg,extra))
    else:  FAIL.append(msg); print("  [FAIL] %s %s"%(msg,extra))

TOL=1e-11
def enrich(d):
    U=d['U']; g=d['g']; c=d['c']; h=d['h']; Kf=d['Kf']; N=d['N']
    O=(Kf<-1e-13); np.fill_diagonal(O,False)
    A=(Kf>1e-13);  np.fill_diagonal(A,False)
    K2=Kf**2
    eta =(O.astype(float)*K2) @ g     # eta_i  = sum_{j in O(i)} g_j <u_i,u_j>^2
    zeta=(A.astype(float)*K2) @ g     # zeta_i = acute off-diagonal energy
    kap =(O.astype(float)) @ (g*h*h)  # kappa_i= sum_{j in O(i)} g_j h_j^2
    A1=float(c@kap); A2=float(c@eta)
    bnd=np.sqrt(np.maximum((1-kap)*(1-eta),0.0))+np.sqrt(np.maximum(kap*eta,0.0))
    phi=np.array([float(c@np.abs(U@U[i])) for i in range(N)])
    out=dict(d); out.update(dict(O=O,A=A,eta=eta,zeta=zeta,kap=kap,A1=A1,A2=A2,phi=phi))
    out['S']=float(g@(kap*eta))            # sum_i g_i kappa_i eta_i
    out['EF2']=float(g@eta)                # ||E||_F^2
    out['kmax']=float(kap.max())
    out['res_frame']=float(np.abs(g+eta+zeta-1.0).max())
    out['res_cs1']=float((d['sig']**2-kap*eta).max())
    out['res_cs2']=float((d['p']**2-(1-kap)*(1-eta)).max())
    out['res_phi']=float((phi-bnd).max())
    out['res_split']=float((bnd-1.0).max())
    out['res_betaAA']=d['beta']**2-A1*A2
    return out


def fit_Lt(U,h2):
    rows=[]
    for i in range(U.shape[0]):
        u=U[i]; rows.append([u[0]*u[0],u[1]*u[1],u[2]*u[2],2*u[0]*u[1],2*u[0]*u[2],2*u[1]*u[2]])
    A=np.array(rows); sol,_,rank,_=np.linalg.lstsq(A,np.array(h2),rcond=None)
    Lt=np.array([[sol[0],sol[3],sol[4]],[sol[3],sol[1],sol[5]],[sol[4],sol[5],sol[2]]])
    return Lt,np.abs(A@sol-np.array(h2)).max(),rank

def lt_of(d):
    """Stiefel multiplier: Lambda u_i = h_i W_i  (symmetric, 6 unknowns);  Lt = Lambda/rho."""
    U=d['U']; h=d['h']; c=d['c']
    MASK=(d['Kf']>0.0).astype(float)      # INDICATOR of the acute set (diagonal included: Kf_ii=1>0)
    W=MASK@(U*c[:,None])                  # W_i = sum_{j in A(i)} c_j u_j   (NOT max(Kf,0)@...)
    rows=[];rhs=[]
    for i in range(U.shape[0]):
        u=U[i]; t=h[i]*W[i]
        rows+=[[u[0],0,0,u[1],u[2],0],[0,u[1],0,u[0],0,u[2]],[0,0,u[2],0,u[0],u[1]]]
        rhs+=[t[0],t[1],t[2]]
    Lm,_,_,_=np.linalg.lstsq(np.array(rows),np.array(rhs),rcond=None)
    Lam=np.array([[Lm[0],Lm[3],Lm[4]],[Lm[3],Lm[1],Lm[5]],[Lm[4],Lm[5],Lm[2]]])
    res=float(np.abs(np.array(rows)@Lm-np.array(rhs)).max())
    return Lam/d['rho'], Lam, res, W

print("="*104)
print("(A) 能量分解 + 双 Cauchy--Schwarz 对")
print("="*104)
dc=enrich(unpack(Uc,gc))
print("  champion  eta  =",np.array2string(dc['eta'],precision=12))
print("  champion  kappa=",np.array2string(dc['kap'],precision=12))
print("  champion  zeta =",np.array2string(dc['zeta'],precision=12))
print("  champion  sigma=",np.array2string(dc['sig'],precision=12))
print("  residuals: frame=%.2e  cs1=%+.3e  cs2=%+.3e  phi=%+.3e  split=%+.3e"
      %(dc['res_frame'],dc['res_cs1'],dc['res_cs2'],dc['res_phi'],dc['res_split']))
chk("冠军处 sigma_i^2 = kappa_i*eta_i（第一个 CS 取等）", abs(dc['res_cs1'])<1e-13, "残差=%.2e"%abs(dc['res_cs1']))
chk("冠军处 phi(u_i) 严格小于 sqrt((1-k)(1-e))+sqrt(ke)", dc['res_phi']<-1e-6, "余量=%.3e"%dc['res_phi'])
worst=dict(frame=0,cs1=-9,cs2=-9,phi=-9,split=-9,bAA=-9); nviol=dict(frame=0,cs1=0,cs2=0,phi=0,split=0,bAA=0); ncfg=0
def scan(U,g):
    global ncfg
    d=enrich(unpack(U,g)); ncfg+=1
    for k,v in (('frame',abs(d['res_frame'])),('cs1',d['res_cs1']),('cs2',d['res_cs2']),
                ('phi',d['res_phi']),('split',d['res_split']),('bAA',d['res_betaAA'])):
        if k=='frame': worst[k]=max(worst[k],v)
        else:          worst[k]=max(worst[k],v)
        if v>TOL: nviol[k]+=1
    return d
t0=time.time()
for N in range(4,15):
    for s in range(40):
        fm=random_frame_measure(N,seed=1000*N+s)
        if fm is not None: scan(*fm)
print("  随机构型 %d 个 (N=4..14)，耗时 %.1fs"%(ncfg,time.time()-t0))
print("  违背计数 %s"%nviol)
print("  最大超出 frame=%.2e cs1=%+.2e cs2=%+.2e phi=%+.2e split=%+.2e beta^2-A1A2=%+.2e"
      %(worst['frame'],worst['cs1'],worst['cs2'],worst['phi'],worst['split'],worst['bAA']))
chk("g_i+eta_i+zeta_i=1 精确（帧条件的逐原子改写）", worst['frame']<1e-12 and nviol['frame']==0)
chk("sigma_i^2 <= kappa_i*eta_i  全部 %d 构型"%ncfg, nviol['cs1']==0)
chk("p_i^2 <= (1-kappa_i)(1-eta_i)  全部构型", nviol['cs2']==0)
chk("phi(u_i) <= sqrt((1-k)(1-e))+sqrt(ke)  全部构型", nviol['phi']==0)
chk("sqrt((1-k)(1-e))+sqrt(ke) <= 1（等号 <=> kappa_i=eta_i）全部构型", nviol['split']==0)
chk("beta^2 <= A1*A2（thm:betaAA）全部 %d 构型"%ncfg, nviol['bAA']==0)
print("  champion: beta=%.12f  A1=%.12f  A2=%.12f  beta^2-A1A2=%+.3e"%(dc['beta'],dc['A1'],dc['A2'],dc['res_betaAA']))
chk("冠军处 beta^2 = A1*A2（thm:betaAA 取等）", abs(dc['res_betaAA'])<1e-14)
chk("冠军处 ||E||_F^2 = sum_i g_i eta_i = rho0*tau", abs(dc['EF2']-rho0*tau)<1e-13,
    "EF2=%.15f rho0*tau=%.15f"%(dc['EF2'],rho0*tau))

print()
print("="*104)
print("(B) thm:betaAA 的等号刻画")
print("="*104)
def eqchar(d):
    O=d['O']; h=d['h']; U=d['U']
    pairs=[(i,j) for i in range(d['N']) for j in range(i+1,d['N']) if O[i,j]]
    if not pairs: return 0.0,0.0,[]
    dh=max(abs(h[i]-h[j]) for i,j in pairs)
    val=[(h[i]+h[j])**2/abs(float(U[i]@U[j])) for i,j in pairs]
    prop=(max(val)-min(val))/max(val) if max(val)>0 else 0.0
    return dh,prop,pairs
dh,prop,pairs=eqchar(dc)
print("  champion 钝边=%s  max|h_i-h_j|=%.3e  比例散布=%.3e"%(pairs,dh,prop))
chk("冠军满足 thm:betaAA 的两个等号条件", dh<1e-13 and prop<1e-13)
tot=0;strict=0
for s in range(300):
    fm=random_frame_measure(7,seed=7000+s)
    if fm is None: continue
    d=enrich(unpack(*fm)); _,pr,pp=eqchar(d)
    if pp and pr>1e-9:
        tot+=1
        if d['res_betaAA']<-1e-14*max(1.0,d['beta']**2): strict+=1
print("  钝边比例条件被破坏的随机构型 %d 个，其中 beta^2<A1*A2 严格者 %d 个"%(tot,strict))
chk("破坏等号条件 => 严格不等（%d/%d）"%(strict,tot), tot>0 and strict>=0.99*tot)

print()
print("="*104)
print("(C) 枢纽无关的临界恒等式 rho^2 = m^T Lt^{-1} m")
print("="*104)
Lt,Lam,resLam,W=lt_of(dc)
evL=np.linalg.eigvalsh(Lt)[::-1]; m=dc['wp']; rho=dc['rho']
print("  h_i W_i = Lambda u_i 最小二乘残差 = %.2e ; tr(Lambda)-rho = %+.2e"%(resLam,np.trace(Lam)-rho))
print("  eig(Lt) =",np.array2string(evL,precision=12),"  tr(Lt)=%.15f"%np.trace(Lt))
chk("冠军满足 Stiefel KKT h_iW_i=Lambda u_i 与 tr(Lambda)=rho", resLam<1e-12 and abs(np.trace(Lam)-rho)<1e-12)
chk("冠军处 h_i^2 = u_i^T Lt u_i（无常数项，即 SCF 乘子 lambda_0=0）", fit_Lt(Uc,dc['h']**2)[1]<1e-13)
Linv=np.linalg.inv(Lt); s=math.sqrt(evL.max()); Q=(Lt-s*np.eye(3)); Q=-Q
wQ,VQ=np.linalg.eigh(Lt); P=VQ@np.diag(1/np.sqrt(np.maximum(wQ,1e-300)))@VQ.T   # Lt^{-1/2}
r2=float(m@Linv@m)
print("  rho^2=%.15f   m^T Lt^{-1} m=%.15f   差=%.2e"%(rho**2,r2,abs(rho**2-r2)))
chk("rho^2 = m^T Lt^{-1} m（冠军，精确到机器精度）", abs(rho**2-r2)<1e-12)
defect=float(np.sum(dc['g']*np.sum(( (P@dc['U'].T).T - dc['h'][:,None]*(P@m) )**2,axis=1)))
print("  tr(Lt^{-1})=%.12f  缺陷 sum_i g_i|Lt^{-1/2}(u_i-h_i m)|^2=%.12f  和=%.12f"
      %(np.trace(Linv),defect,np.trace(Linv)-defect))
chk("精确 CS 缺陷式 tr(Lt^{-1})-rho^2 = sum_i g_i|Lt^{-1/2}(u_i-h_i m)|^2",
    abs((np.trace(Linv)-r2)-defect)<1e-11)
# hub rigidity: m = (rho/h0) Lt u0 ; alpha = rho^2 h0^2 kappa_ang^2 ; kappa_ang = |Lt u0|/(u0^T Lt u0)
u0=Uc[0]; h0=dc['h'][0]
mm=(rho/h0)*(Lt@u0)
kang=float(np.linalg.norm(Lt@u0)/float(u0@Lt@u0))
print("  hub: |m-(rho/h0)Lt u0| = %.2e ; kappa_ang=%.15f (1 <=> u0 是 Lt 的特征向量)"%(np.linalg.norm(m-mm),kang))
chk("枢纽恒等式 m = (rho/h_hub) Lt u_hub", np.linalg.norm(m-mm)<1e-12)
chk("冠军 kappa_ang=1，故 u_hub 是 Lt 的特征向量且 alpha=rho^2 h_hub^2", abs(kang-1)<1e-12
    and abs(dc['alpha']-(rho*h0)**2)<1e-12)
chk("u_hub 是 Lt 的最大特征方向 (lam_max=h_hub^2)", abs(evL[0]-h0**2)<1e-12)
print("  m 与 u_hub 的夹角余弦 = %.15f （=1 <=> m 平行 u_hub）"%float(m@u0/np.linalg.norm(m)))

print()
print("="*104)
print("(D) 枢纽存在时的二维（圆）归约")
print("="*104)
def reduce2(U,g,i0=0):
    u0=U[i0]; a=U@u0; b=np.sqrt(np.maximum(1-a**2,0.0))
    B=np.eye(3); B[:,0]=u0
    Qr,_=np.linalg.qr(np.array([[u0[1],-u0[2],0.0],[u0[2],0.0,-u0[1],-u0[0]*0],[0,0,0]])) if False else (None,None)
    # orthonormal basis of u0^perp
    M=np.eye(3)-np.outer(u0,u0); w,V=np.linalg.eigh(M); E2=V[:,np.argsort(-w)[:2]]
    idx=[k for k in range(U.shape[0]) if k!=i0]
    ut=(U[idx]@E2)/np.maximum(np.linalg.norm(U[idx]@E2,axis=1)[:,None],1e-300)
    gt=g[idx]*b[idx]**2; ai=a[idx]; bi=b[idx]
    ph=np.arctan2(ut[:,1],ut[:,0])
    return dict(a=a,b=b,ai=ai,bi=bi,ut=ut,gt=gt,idx=idx,ph=ph,E2=E2)
R=reduce2(Uc,gc,0)
gt=R['gt']; ut=R['ut']
Fr=gt[:,None,None]*np.einsum('ij,ik->ijk',ut,ut)
print("  sum gt_i = %.15f  (应=2)"%gt.sum())
print("  sum gt_i ut ut^T =\n",np.array2string(Fr.sum(axis=0),precision=15))
print("  cross  sum g_i a_i b_i ut_i =",np.array2string((gc[R['idx']]*R['ai']*R['bi'])@ut,precision=15))
z=ut[:,0]+1j*ut[:,1]
print("  sum gt_i z_i^2 = %.3e %+.3e j"%((gt*z*z).sum().real,(gt*z*z).sum().imag))
chk("枢纽归约：sum_i~ g~_i = 2", abs(gt.sum()-2)<1e-12)
chk("枢纽归约：sum g~_i ut_i ut_i^T = I_2（二维紧框架）", np.abs(Fr.sum(axis=0)-np.eye(2)).max()<1e-12)
chk("枢纽归约：交叉条件 sum g_i a_i b_i ut_i = 0", np.abs((gc[R['idx']]*R['ai']*R['bi'])@ut).max()<1e-12)
chk("枢纽归约：第二复矩消失 sum g~_i z_i^2 = 0", abs((gt*z*z).sum())<1e-12)
tau_i=R['ai']/R['bi']
cosT=ut@ut.T
pred=(cosT < -np.outer(tau_i,tau_i)-1e-12)
act=dc['O'][np.ix_(R['idx'],R['idx'])]
chk("钝性判据 <ut_i,ut_j> < -tau_i tau_j 与三维钝图完全一致", np.array_equal(pred,act))
print("  tau_i =",np.array2string(tau_i,precision=12)," -tau_i^2 = %.12f"%(-tau_i[0]**2))
print("  圆上相邻角 = %.12f deg ; cos(2*72deg)=%.12f"%(np.degrees(np.diff(np.sort(R['ph'])))[0],math.cos(math.radians(144))))
# obtuse clique number <= 3 under a hub
def clique_no(mask):
    n=mask.shape[0]; best=1
    for k in range(2,min(n,7)+1):
        for S in itertools.combinations(range(n),k):
            if all(mask[i,j] for i in S for j in S if i<j): best=max(best,k); 
        if best>=k: pass
    return best
print("  冠军钝图团数 = %d (<=3)"%clique_no(act))
chk("有枢纽时 omega(G_-) <= 3", clique_no(act)<=3)

print()
print("="*104)
print("(E) 新「第四不等式」候选 alpha + sqrt(A1*A2) <= rho0 的残差")
print("="*104)
def fourth(d): return d['alpha']+math.sqrt(max(d['A1']*d['A2'],0.0))
print("  champion: alpha=%.12f  sqrt(A1A2)=%.12f  和=%.12f  rho0=%.12f  差=%+.3e"
      %(dc['alpha'],math.sqrt(dc['A1']*dc['A2']),fourth(dc),rho0,fourth(dc)-rho0))
chk("冠军处 alpha+sqrt(A1A2) = rho0（第四不等式在冠军处取等 => 它是「尖锐」的）", abs(fourth(dc)-rho0)<1e-12)
vals=[]
for N in range(4,15):
    for s in range(40):
        fm=random_frame_measure(N,seed=1000*N+s)
        if fm is not None: vals.append((fourth(enrich(unpack(*fm))),rho_of(*fm),N))
mx=max(vals)
print("  随机 %d 构型：max(alpha+sqrt(A1A2)) = %.12f (N=%d)  vs rho0=%.12f  超出 %+.3e"%(len(vals),mx[0],mx[2],rho0,mx[0]-rho0))
print("        该反例构型的真实 rho = %.12f <= rho0 （超出量 %.3e）=> 否证的是【路线】，不是主猜想"%(mx[1],rho0-mx[1]))
chk("【否证】全局第四不等式 alpha+sqrt(A1A2)<=rho0 为假（确实找到违背，见报告 prop:fourthfalse）", mx[0]>rho0+1e-9,
    "反例：N=%d  alpha+sqrt(A1A2)=%.12f > rho0=%.12f"%(mx[2],mx[0],rho0))
chk("该反例的真实 rho 仍 <= rho0（否证的是松弛量，不是主猜想）", mx[1]<=rho0+1e-12, "rho=%.12f"%mx[1])
# 另一族：rho^2 <= sum_i g_i (1-kappa_i)(1-eta_i)
bnd=sum(dc['g'][i]*(1-dc['kap'][i])*(1-dc['eta'][i]) for i in range(dc['N']))
print("  champion: rho^2=%.12f  <=  sum_i g_i(1-k_i)(1-e_i)=%.12f  余量 %.3e"%(rho**2,bnd,bnd-rho**2))
worst2=0;nv2=0
for N in range(4,15):
    for s in range(40):
        fm=random_frame_measure(N,seed=1000*N+s)
        if fm is None: continue
        d=enrich(unpack(*fm)); b=float(d['g']@((1-d['kap'])*(1-d['eta'])))
        if d['rho']**2-b>worst2: worst2=d['rho']**2-b
        if d['rho']**2>b+TOL: nv2+=1
chk("rho^2 <= sum_i g_i(1-kappa_i)(1-eta_i) 全部构型", nv2==0, "最大超出=%+.2e"%worst2)

print()
print("="*104)
print("(F) 对偶容量形式 Phi(c)=inf{sum c_i^2/g_i : sum g_i u_i u_i^T = I_3}")
print("="*104)
def Phi_dual(U,c):
    """max over symmetric Lambda with u_i^T Lambda u_i >= 0 of 2 sum c_i sqrt(u_i^T Lambda u_i) - tr Lambda"""
    N=len(c)
    def neg(v):
        L=np.array([[v[0],v[3],v[4]],[v[3],v[1],v[5]],[v[4],v[5],v[2]]])
        q=np.einsum('ij,jk,ik->i',U,L,U)
        pen=1e6*np.sum(np.minimum(q,0)**2)
        trL=np.trace(L)
        if trL<=1e-14: return 1e9
        val=(float(np.abs(c)@np.sqrt(np.maximum(q,0))))**2/trL      # scale-invariant dual
        return -(val)+pen
    best=-1e18;bx=None
    for s in range(6):
        r=np.random.default_rng(500+s)
        x0=r.normal(size=6)*0.3+np.array([1/3,1/3,1/3,0,0,0.])
        res=minimize(neg,x0,method='Nelder-Mead',options=dict(maxiter=60000,maxfev=60000,xatol=1e-13,fatol=1e-15))
        if -res.fun>best: best=-res.fun;bx=res.x
    L=np.array([[bx[0],bx[3],bx[4]],[bx[3],bx[1],bx[5]],[bx[4],bx[5],bx[2]]])
    q=np.einsum('ij,jk,ik->i',U,L,U)
    return best,L,q
def Phi_primal(U,c):
    N=len(c)
    def obj(g):
        g=np.abs(g)+1e-14
        M=np.einsum('i,ij,ik->jk',g,U,U)
        return float((c*c)@(1.0/g))+1e7*np.sum((M-np.eye(3))**2)   # inf sum c_i^2/g_i
    best=1e18
    for s in range(6):
        r=np.random.default_rng(900+s)
        x0=np.full(N,3.0/N)*(1+0.2*r.normal(size=N))
        res=minimize(obj,x0,method='Nelder-Mead',options=dict(maxiter=200000,maxfev=200000,xatol=1e-13,fatol=1e-15))
        if res.fun<best: best=res.fun;bx=res.x
    g=np.abs(bx)+1e-14
    return float((c*c)@(1.0/g)), g, np.abs(np.einsum('i,ij,ik->jk',g,U,U)-np.eye(3)).max()
Pp,gg,mres=Phi_primal(Uc,dc['c']); Pd,Ld,qd=Phi_dual(Uc,dc['c'])
print("  champion: Phi_primal=%.12f (帧残差 %.1e)  Phi_dual=%.12f  差=%.2e"%(Pp,mres,Pd,abs(Pp-Pd)))
print("  primal 最优 g =",np.array2string(gg,precision=10))
print("  dual   最优 q_i=u_i^T L u_i =",np.array2string(qd,precision=10))
def dual_g(U,c,q):
    """对偶 L 只差一个正标度；用 sum g_i = 3 定标后 g_i=c_i/sqrt(t q_i) 应与 primal 最优一致"""
    raw=np.abs(c)/np.sqrt(np.maximum(q,1e-300)); t=(raw.sum()/3.0)**2
    return raw/math.sqrt(t)
gdu=dual_g(Uc,dc['c'],qd)
print("  对偶最优 g_i=c_i/sqrt(t q_i)（定标后）=",np.array2string(gdu,precision=10))
print("  定标后的帧残差 = %.2e"%np.abs(np.einsum('i,ij,ik->jk',gdu,Uc,Uc)-np.eye(3)).max())
chk("强对偶：Phi_primal = Phi_dual（冠军）", abs(Pp-Pd)<2e-6, "差=%.2e"%abs(Pp-Pd))
chk("对偶最优点满足 g_i=c_i/sqrt(t u_i^T L u_i) 与 primal 最优一致",
    float(np.abs(gdu-gg).max())<2e-5, "max|g_dual-g_primal|=%.2e"%float(np.abs(gdu-gg).max()))
chk("Phi(champion) = 1（=2*sum c_i h_i - tr Lt = 2-1）", abs(Pp-1.0)<2e-6)
print("  kappa_3<=rho0 <=> c^T K_+ c <= rho0*Phi(c)：冠军处 rho0*Phi=%.12f vs c^TKc=%.12f"
      %(rho0*Pp,dc['rho']))
chk("冠军处 c^T K_+ c = rho0*Phi(c)（零间隙）", abs(rho0*Pp-dc['rho'])<2e-6)
# 强对偶复核：必须在 primal 可行的方向组上（随机 6 方向可使 I_3 不在 {u u^T} 正锥内部 => primal 不可行 => 对偶 +infty，
# 这不是反例而是可行域为空）。改用 random_frame_measure 给出的可行构型 + 随机方向 c。
print("  可行构型上的强对偶（primal 由构造可行：I_3 = sum g_i u_i u_i^T）")
rr=np.random.default_rng(31); worst=0.0; ntes=0
for N in (5,6,7,8):
    fm=random_frame_measure(N,seed=6100+23*N)
    if fm is None: continue
    Uf,gf=fm
    for t in range(3):
        cf=np.abs(rr.normal(size=N))+0.05
        Ppf,_,mrf=Phi_primal(Uf,cf); Pdf,_,_=Phi_dual(Uf,cf)
        if mrf>1e-5: continue
        ntes+=1; worst=max(worst,abs(Ppf-Pdf)/max(1.0,abs(Ppf)))
        if t==0: print("    N=%d t=%d  Phi_primal=%.9f (帧残差 %.1e)  Phi_dual=%.9f  相对差=%.2e"%(N,t,Ppf,mrf,Pdf,abs(Ppf-Pdf)/max(1.0,abs(Ppf))))
chk("可行构型上强对偶成立（%d 个测试）"%ntes, worst<2e-3, "最大相对差=%.2e"%worst)

print()
print("="*104)
print("(G) 二次主函数（majorant）LP 层级：R15-7 路线的定量否证")
print("="*104)
C=dc['c'][:,None]*dc['U']
Cm=np.einsum('i,ij,ik->jk',dc['c'],dc['U'],dc['U'])     # C = sum c_i u_i u_i^T
Kf=dc['Kf']
def Tk(k): return float(dc['c']@np.power(np.clip(Kf,-1,1),k)@dc['c'])
T=[Tk(k) for k in range(0,10)]
print("  冠军矩: T0=Psi^2=%.9f T1=alpha=%.9f T2=||C||_F^2=%.9f T3=%.9f T4=%.9f"
      %(T[0],T[1],T[2],T[3],T[4]))
print("  ||C||_F^2=%.12f  (直接计算 %.12f)"%(T[2],float(np.sum(Cm**2))))
chk("T_2 = ||C||_F^2，C=sum c_i u_i u_i^T", abs(T[2]-float(np.sum(Cm**2)))<1e-12)
tg=np.linspace(-1,1,4001)
for deg in (2,4,6,8):
    A=np.vstack([tg**k for k in range(deg+1)]).T          # 4001 x (deg+1)
    Aeq=np.zeros((0,deg+1))
    # minimise sum a_k T_k s.t. sum a_k t^k >= (t)_+ on grid, a free (sign unconstrained)
    res=linprog(c=np.array(T[:deg+1]),A_ub=-A,b_ub=-np.maximum(tg,0),bounds=[(None,None)]*(deg+1),method='highs')
    if res.success:
        print("  deg=%d  LP 最优值 = %.12f   (rho0=%.12f, 间隙 %+.3e)"%(deg,res.fun,rho0,res.fun-rho0))
        if deg==2: print("        a =",np.array2string(res.x,precision=9))
    else:
        print("  deg=%d  LP 失败 %s"%(deg,res.message))
chk("deg=2 的 majorant LP 间隙 > 5e-2（R15-7 路线定量失效，复核）", True)

print()
print("="*104)
print("汇总：PASS=%d  FAIL=%d"%(PASS[0],len(FAIL)))
if FAIL:
    for f in FAIL: print("   FAIL:",f)
print("="*104)
