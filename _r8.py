import numpy as np, math
np.set_printoptions(precision=12, suppress=False, linewidth=200)
exec(open('_o1.py').read().split('print("="*104)')[0])
def core(U,g):
    N=len(g); Kf=U@U.T; V=U*np.sqrt(g)[:,None]; Q=V@V.T; Qp=np.maximum(Q,0.0)
    w,ev=np.linalg.eigh(Qp); i=int(np.argmax(w)); rho=w[i]; x=ev[:,i].copy()
    if x.sum()<0: x=-x
    h=x/np.sqrt(g); c=g*h; m=c@U
    return rho,Kf,float(m@m),np.abs((U*g[:,None]).T@U-np.eye(3)).max()
def whiten(U,g):
    M=(U*g[:,None]).T@U; ev,EV=np.linalg.eigh(M)
    if ev.min()<=1e-9: return None,None
    S=EV@np.diag(ev**-0.5)@EV.T; W=U@S; nr=np.linalg.norm(W,axis=1)
    return W/nr[:,None], g*nr**2
def transitive(Kf,tol=0.0):
    N=Kf.shape[0]
    best=None
    for m in range(1<<(N-1)):
        sg=np.array([1.0]+[1.0 if (m>>k)&1 else -1.0 for k in range(N-1)])
        r=(np.outer(sg,sg)*Kf).min()
        if best is None or r>best[0]: best=(r,sg)
    return best
rng=np.random.default_rng(888); mx=-1e9; ntr=0; worst=None; skipframe=0
for t in range(6000):
    N=int(rng.integers(3,9)); k=int(rng.integers(2,5)); lab=rng.integers(0,k,size=N)
    dirs=rng.normal(size=(k,3)); dirs/=np.linalg.norm(dirs,axis=1)[:,None]
    U=np.array([dirs[l]+1e-3*rng.normal(size=3) for l in lab]); U/=np.linalg.norm(U,axis=1)[:,None]
    g=rng.dirichlet(np.ones(N))*3.0; U,g=whiten(U,g)
    if U is None: continue
    rho,Kf,alpha,fr=core(U,g)
    if fr>1e-10: skipframe+=1; continue
    r,sg=transitive(Kf)
    if r<0.0: continue           # 严格可传递：sigma_i sigma_j <u_i,u_j> >= 0
    ntr+=1
    if rho>mx: mx=rho; worst=(N,alpha,fr,r)
print('严格可传递构型 %d 个（帧残差>1e-10 而剔除 %d 个）：max rho = %.15f'%(ntr,skipframe,mx))
print('该最坏构型: N=%d alpha=%.15f 帧残差=%.2e min(sigma_i sigma_j Kf_ij)=%.2e'%(worst[0],worst[1],worst[2],worst[3]))
print('结论：rho <= 1 %s （超出 %.3e）'%('成立' if mx<=1+1e-11 else '不成立',mx-1))
