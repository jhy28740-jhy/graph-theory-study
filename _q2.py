# -*- coding: utf-8 -*-
# _q2.py  ROUND 17 verification of the NEW identities/inequalities
import numpy as np, math, itertools, warnings
warnings.filterwarnings('ignore')
np.set_printoptions(precision=12, suppress=False, linewidth=190)
exec(open('_o1.py', encoding='utf-8').read().split('print("="*104)')[0])

def full(U,g):
    V=U*np.sqrt(g)[:,None]; Q=V@V.T; Qp=np.maximum(Q,0.0); Qm=np.maximum(-Q,0.0)
    w,ev=np.linalg.eigh(Qp); i=int(np.argmax(w)); rho=float(w[i]); x=ev[:,i].copy()
    if x.sum()<0: x=-x
    return dict(V=V,Q=Q,Qp=Qp,Qm=Qm,w=np.linalg.eigvalsh(Qp)[::-1],rho=rho,x=x,
                alpha=float(x@Q@x), beta=float(x@Qm@x), tau=rho-1.0,
                e=float((Qm**2).sum()), Fp=float((Qp**2).sum()))

def cfgs():
    out=[('champion',Uc,gc)]
    for N in (4,5,6,7,8,9,10,12,14):
        for s in range(3):
            fm=random_frame_measure(N,seed=5500+19*N+s)
            if fm is not None: out.append(('rand N=%d#%d'%(N,s),fm[0],fm[1]))
    return out
C=cfgs()

print("="*104); print("N1  精确恒等式 ||Q_+||_F^2 = 3 - ||Q_-||_F^2   (因为 <Q,Q_-> = -||Q_-||_F^2)"); print("="*104)
w1=0.0; w2=0.0; w3=0.0
for name,U,g in C:
    d=full(U,g); Q=d['Q']; Qm=d['Qm']
    w1=max(w1,abs(d['Fp']-(3.0-d['e'])))
    w2=max(w2,abs(float((Q*Qm).sum())+d['e']))          # <Q,Q_-> = -||Q_-||_F^2
    w3=max(w3,abs(float((Q*Q).sum())-3.0))               # ||Q||_F^2 = 3  (Welch)
print("  WORST over %d configs:  ||Q+||_F^2-(3-e) = %.2e   <Q,Q_->+e = %.2e   ||Q||_F^2-3 = %.2e"%(len(C),w1,w2,w3))
d=full(Uc,gc); print("  champion: ||Q_+||_F^2=%.12f  3-e=%.12f  e=||Q_-||_F^2=%.12f  (=rho0*tau=%.12f)"%(d['Fp'],3-d['e'],d['e'],rho0*tau))

print(); print("="*104); print("N2  新向量恒等式  Q_- x = tau*x + (I-Q)x      <=>   tau = beta - (1-alpha) = x^T(Q_+ - I)x"); print("="*104)
r1=r2=r3=0.0
for name,U,g in C:
    d=full(U,g); x=d['x']; Q=d['Q']; tau=d['tau']
    lhs=d['Qm']@x; rhs=tau*x+(x-Q@x)
    r1=max(r1,float(np.abs(lhs-rhs).max()))
    r2=max(r2,abs(tau-(d['beta']-(1-d['alpha']))))
    r3=max(r3,abs(d['rho']-1-float(x@(d['Qp']-np.eye(len(x)))@x)))
print("  WORST: |Q_-x - tau x - (I-Q)x| = %.2e   |tau-(beta-(1-alpha))| = %.2e   |rho-1-x^T(Q_+-I)x| = %.2e"%(r1,r2,r3))
print("  => 主不等式 rho<=rho0  等价于   x^T Q_- x - ||(I-Q)x||^2 <= tau0 = %.12f"%(rho0-1))

print(); print("="*104); print("N3  tau = lambda_max(Q_+)-1 = max_y [ y^T Q_- y - ||(I-Q)y||^2 ] ;  3x3 块界 lambda_max(Xi), Xi=V^T Q_- V"); print("="*104)
print("   %-16s %10s %10s %10s %10s %10s"%("config","tau","lmax(Q-)","lmax(Xi)","||Phi12||","Xi+Phi12"))
viol=0
for name,U,g in C[:1]+C[2:8]+C[-6:]:
    d=full(U,g); V=d['V']; Q=d['Q']; Qm=d['Qm']; x=d['x']
    Xi=V.T@Qm@V; lXi=float(np.linalg.eigvalsh(Xi)[-1]); lQm=float(np.linalg.eigvalsh(Qm)[-1])
    Pi=np.linalg.qr(V)[0]; Pperp=np.eye(len(x))-Pi@Pi.T; Ph12=V.T@Qm@Pperp; n12=float(np.linalg.norm(Ph12,2))
    print("   %-16s %10.6f %10.6f %10.6f %10.6f %10.6f"%(name,d['tau'],lQm,lXi,n12,lXi+n12))
    if d['tau']>lXi+n12+1e-12: viol+=1
print("  块界 tau <= lambda_max(Xi)+||Phi12|| 违背数 = %d"%viol)

print(); print("="*104); print("N4  钝图团数 omega(G_-)<=4 (R^3 中两两钝角向量至多 4 个) + Motzkin-Straus  w:=sum_{i~j} x_i^2 x_j^2 <= 3/8"); print("="*104)
def clique(mask):
    n=mask.shape[0]; best=1
    for k in range(2,min(n,8)+1):
        found=False
        for S in itertools.combinations(range(n),k):
            if all(mask[i,j] for i in S for j in S if i<j): best=k; found=True; break
        if not found: break
    return best
wo=0; wms=0; wcs=0
for name,U,g in C:
    d=full(U,g); U2=U@U.T; O=U2<0; np.fill_diagonal(O,False); x=d['x']; Qm=d['Qm']
    wo=max(wo,clique(O.astype(float)>0))
    om=float((np.outer(x**2,x**2)*O).sum()); wms=max(wms,om)
    wcs=max(wcs, d['beta']**2 - d['e']*om )          # beta^2 <= e * omega
print("  最大钝团数 omega(G_-) = %d  (应 <= 4)"%wo)
print("  最大 omega = sum_{i~j}x_i^2x_j^2 = %.12f   Motzkin-Straus 上界 3/8 = %.12f"%(wms,0.375))
print("  beta^2 - e*omega 的最大值 = %+.3e  (应 <= 0, 即 beta <= sqrt(e*omega))"%wcs)

print(); print("="*104); print("N5  谱地板的另一形式  n_+(Q_+) >= 9/(3-e)   (由 sum lam^2 >= (sum lam)^2/n_+)"); print("="*104)
v=0
for name,U,g in C:
    d=full(U,g); npz=int((d['w']>1e-10).sum()); fl=9.0/(3.0-d['e'])
    if npz < fl-1e-9: v+=1; print("   VIOLATION",name,npz,fl)
print("  违背数 = %d"%v)
for lbl,r in (("rho0",rho0),("Phi",PHI),("Mohar",MOHAR)):
    print("   在 rho=%s 处需要 e >= 3-rho^2*... : n_+ >= ceil(9/(3-e))；若 e=||Q_-||_F^2=3-rho^2 则 n_+ >= %.3f"%(lbl,9/(3-r**2)))

print(); print("="*104); print("N6  Branch-II 刚性: 极大元(N=8,12)的原子重合分组 + rows(B) in rowspace(M) + R|_K = 0"); print("="*104)
def ascend(N,seed,steps=6000,lr=0.05):
    r=np.random.default_rng(seed); A=r.normal(size=(N,3)); V,_=np.linalg.qr(A); V=V[:,:3]
    best=(-1,None)
    for it in range(steps):
        Q=V@V.T; K=(Q>0).astype(float); Qp=K*Q
        w,ev=np.linalg.eigh(Qp); i=int(np.argmax(w)); rr=float(w[i]); xx=ev[:,i].copy()
        if xx.sum()<0: xx=-xx
        if rr>best[0]: best=(rr,V.copy())
        G=2*np.outer(xx,xx)*K@V; G=G-V@((V.T@G+G.T@V)/2)
        V,_=np.linalg.qr(V+lr*G); lr*=0.9995
    return best
for N in (8,12):
    bb=(-1,None)
    for s in range(5):
        r_,V_=ascend(N,seed=400+31*N+s)
        if r_>bb[0]: bb=(r_,V_)
    rho,V=bb; g=np.linalg.norm(V,axis=1)**2; U=V/np.sqrt(g)[:,None]
    d=full(U,g); x=d['x']
    Kf=U@U.T; close=np.abs(1.0-Kf)<1e-6
    groups=[]; seen=set()
    for i in range(N):
        if i in seen: continue
        grp=[j for j in range(N) if close[i,j]]; seen.update(grp); groups.append(grp)
    print("  N=%d  rho-rho0 = %+.2e   重合分组 = %s"%(N,rho-rho0,groups))
    for gr in groups:
        print("       %-22s sum g = %.12f"%(str(gr),g[gr].sum()))
    # rows of B = Q_+ D_{x/g}  must lie in rowspace(M)  <=>  B a = 0 for all a in kerM
    B=d['Qp']@np.diag(x/g); K=kernel_basis(U)
    print("       max|B a| over kerM basis = %.3e   (dim kerM=%d)"%(float(np.abs(B@K).max()),K.shape[1]))
    R=np.zeros((N,N)); w,ev=np.linalg.eigh(d['Qp']); idx=np.argsort(-w); w=w[idx]; ev=ev[:,idx]
    for k in range(1,N):
        if abs(w[0]-w[k])<1e-12: continue
        R+=(w[k]/(w[0]-w[k]))*np.outer(ev[:,k],ev[:,k])
    D=np.diag(x/g); RK=K.T@D@R@D@K
    print("       eig(R|_K) = %s   (Branch II 要求全为 0)"%np.array2string(np.linalg.eigvalsh((RK+RK.T)/2),precision=3))


