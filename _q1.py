# -*- coding: utf-8 -*-
# _q1.py  ROUND 17 (GLOBAL): moment-kernel second variation -> atom-count reduction
#   T1  re-verify  rho''(a) = 2 rho sum_{k>=2} rho_k (q.w_k)^2/(rho-rho_k),  q_i=(a_i/g_i)x_i
#   T2  if Q_+ >= 0 then rho''(a) > 0 for every 0 != a in ker M  (=> no local max with N>=7)
#   T3  necessary condition at a local max:  R = sum_{k>=2} rho_k/(rho-rho_k) w_k w_k^T  <= 0 on K
#       => n_+(Q_+) <= 6, and  N-6 <= nu_-(Q_+) + dim ker Q_+
#   T4  Stiefel gradient ascent for N=5..12: is the maximizer a zero-row padding, or nu_->0 ?
import numpy as np, math, time, warnings
warnings.filterwarnings('ignore')
np.set_printoptions(precision=12, suppress=False, linewidth=190)
exec(open('_o1.py', encoding='utf-8').read().split('print("="*104)')[0])

def fit_Lt(U,h2):
    rows=[]
    for i in range(U.shape[0]):
        u=U[i]; rows.append([u[0]*u[0],u[1]*u[1],u[2]*u[2],2*u[0]*u[1],2*u[0]*u[2],2*u[1]*u[2]])
    A=np.array(rows); sol,_,rank,_=np.linalg.lstsq(A,np.array(h2),rcond=None)
    Lt=np.array([[sol[0],sol[3],sol[4]],[sol[3],sol[1],sol[5]],[sol[4],sol[5],sol[2]]])
    return Lt,np.abs(A@sol-np.array(h2)).max(),rank

def spec(U,g):
    V=U*np.sqrt(g)[:,None]; Q=V@V.T; Qp=np.maximum(Q,0.0)
    w,ev=np.linalg.eigh(Qp); idx=np.argsort(-w); return w[idx],ev[:,idx],Q,Qp

print("="*104); print("T1  二阶变分公式 rho''(a)=2*rho*sum_{k>=2} rho_k (q.w_k)^2/(rho-rho_k)   [沿矩核]"); print("="*104)
worst1=worst2=0.0
for N in (7,8,9,10):
    for s in range(3):
        fm=random_frame_measure(N,seed=700+13*N+s)
        if fm is None: continue
        U,g=fm; w,ev,Q,Qp=spec(U,g); rho=w[0]; x=ev[:,0].copy()
        if x.sum()<0: x=-x
        K=kernel_basis(U)
        if K.shape[1]==0: continue
        for j in range(min(4,K.shape[1])):
            a=K[:,j]/np.linalg.norm(K[:,j])
            h=2e-4; f=lambda t: rho_of(U,g+t*a)
            n2=(f(h)-2*f(0)+f(-h))/h**2
            q=(a/g)*x
            f2=sum(2*rho*w[k]*float(q@ev[:,k])**2/(rho-w[k]) for k in range(1,N) if abs(rho-w[k])>1e-12)
            worst1=max(worst1,abs(n2-f2))
            # also: first-order condition  <a,h^2> = 0  for a in ker M
            worst2=max(worst2,abs(float(a@(x**2/g))))
        print("  N=%2d s=%d dim kerM=%d  rho=%.9f  nu_-(Q+)=%d  n_+(Q+)=%d"%(N,s,K.shape[1],rho,int((w<-1e-10).sum()),int((w>1e-10).sum())))
print("  WORST |rho''_num - formula| = %.3e    WORST |<a,h^2>| = %.3e"%(worst1,worst2))

print(); print("="*104); print("T2/T3  Q_+ >= 0  =>  rho''(a) > 0 对所有 0 != a in ker M   (等价: R|_K 的最小特征值 > 0)"); print("="*104)
npsd=0; ntot=0; minlam=1e9; worstF=0.0
for N in (7,8,9,10,11,12):
    for s in range(6):
        fm=random_frame_measure(N,seed=3100+29*N+s)
        if fm is None: continue
        U,g=fm; w,ev,Q,Qp=spec(U,g); rho=w[0]; x=ev[:,0].copy()
        if x.sum()<0: x=-x
        K=kernel_basis(U)
        if K.shape[1]==0: continue
        ntot+=1
        psd = bool((w>=-1e-11).all())
        # R = sum_{k>=2} rho_k/(rho-rho_k) w_k w_k^T   ;  rho''(a) = 2 rho * q^T R q
        R=np.zeros((N,N))
        for k in range(1,N):
            if abs(rho-w[k])<1e-12: continue
            R += (w[k]/(rho-w[k]))*np.outer(ev[:,k],ev[:,k])
        Qb=(K*(x/g)[:,None]).T          # columns span cal-K in q-space:  q = Qb a-coeffs
        # pull R back to ker M coordinates:  q = D a, D=diag(x/g)
        D=np.diag(x/g); RK=K.T@D@R@D@K
        evRK=np.linalg.eigvalsh((RK+RK.T)/2)
        if psd:
            npsd+=1; minlam=min(minlam,evRK[0])
        worstF=min(worstF,evRK[0])
        if s==0:
            print("  N=%2d  Q_+ psd=%-5s nu_-=%d n_+=%d  dim kerM=%d  eig(R|_K) min=%+.6e max=%+.6e"
                  %(N,psd,int((w<-1e-10).sum()),int((w>1e-10).sum()),K.shape[1],evRK[0],evRK[-1]))
print("  随机构型 %d 个，其中 Q_+ >= 0 的 %d 个；这些上 min eig(R|_K) = %+.4e  (应 > 0)"%(ntot,npsd,minlam))
print("  所有构型上 min eig(R|_K) 的最小值 = %+.4e  (<0 <=> Q_+ 有负特征值时允许局部极大)"%worstF)

print(); print("="*104); print("T4  Stiefel 梯度上升 (已验证的解析梯度 grad_V rho = 2 (x x^T o K) V)：N=5..12 的极大元结构"); print("="*104)
def ascend(N,seed,steps=4000,lr=0.05):
    r=np.random.default_rng(seed)
    A=r.normal(size=(N,3)); V,_=np.linalg.qr(A); V=V[:,:3]
    best=(-1,None)
    for it in range(steps):
        Q=V@V.T; K=(Q>0).astype(float); Qp=K*Q
        w,ev=np.linalg.eigh(Qp); i=int(np.argmax(w)); rho=float(w[i]); x=ev[:,i].copy()
        if x.sum()<0: x=-x
        if rho>best[0]: best=(rho,V.copy())
        G=2*np.outer(x,x)*K@V
        G=G-V@((V.T@G+G.T@V)/2)              # project onto Stiefel tangent
        Vn=V+lr*G; V,_=np.linalg.qr(Vn)
        lr*=0.9995
    return best
for N in (5,6,7,8,9,10,12):
    bb=(-1,None)
    for s in range(4):
        r_,V_=ascend(N,seed=90+N*7+s)
        if r_>bb[0]: bb=(r_,V_)
    rho,V=bb; g=np.linalg.norm(V,axis=1)**2; w,ev,Q,Qp=spec(V/np.sqrt(g)[:,None],g)
    print("  N=%2d  max rho=%.12f   rho0=%.12f  gap=%+.3e | min g_i=%.3e  nu_-(Q_+)=%d  n_+(Q_+)=%d"
          %(N,rho,rho0,rho-rho0,g.min(),int((w<-1e-9).sum()),int((w>1e-9).sum())))
    print("        g = %s"%np.array2string(np.sort(g),precision=6))
