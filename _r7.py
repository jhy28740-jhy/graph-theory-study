# -*- coding: utf-8 -*-
# _r7.py  聚簇/可传递符号模式定理的正确形式与检验
import numpy as np, math
np.set_printoptions(precision=12, suppress=False, linewidth=200)
exec(open('_o1.py').read().split('print("="*104)')[0])
s5=math.sqrt(5.0); rho0=(9+s5+2*math.sqrt(35.0))/18
PASS=FAIL=0
def chk(c,m):
    global PASS,FAIL
    if c: PASS+=1; print('  [PASS]',m)
    else: FAIL+=1; print('  [FAIL]',m)
def core(U,g):
    N=len(g); Kf=U@U.T; V=U*np.sqrt(g)[:,None]; Q=V@V.T; Qp=np.maximum(Q,0.0)
    w,ev=np.linalg.eigh(Qp); idx=np.argsort(-w); w=w[idx]; ev=ev[:,idx]
    rho=w[0]; x=ev[:,0].copy()
    if x.sum()<0: x=-x
    h=x/np.sqrt(g); c=g*h; m=c@U
    return dict(N=N,rho=rho,h=h,alpha=float(m@m),Kf=Kf,g=g,U=U,x=x)
def whiten(U,g):
    M=(U*g[:,None]).T@U; ev,EV=np.linalg.eigh(M); ev=np.maximum(ev,1e-12)
    S=EV@np.diag(ev**-0.5)@EV.T; W=U@S; nr=np.linalg.norm(W,axis=1)
    return W/nr[:,None], g*nr**2
def transitive(Kf,tol=1e-9):
    """exists sigma in {+-1}^N with sigma_i sigma_j <u_i,u_j> >= 0 for all i,j ?  (brute force over 2^(N-1))"""
    N=Kf.shape[0]
    if N>18: return None
    for m in range(1<<(N-1)):
        sg=np.array([1]+[1 if (m>>k)&1 else -1 for k in range(N-1)],dtype=float)
        if np.all(np.outer(sg,sg)*Kf>=-tol): return sg
    return None
print('='*100); print('(H-fix) 可传递符号模式定理：存在 sigma 使 sigma_i sigma_j <u_i,u_j> >= 0  =>  rho <= 1'); print('='*100)
print('  证明要点：此时 sum_ij x_i x_j |<v_i,v_j>| = ||sum_i sigma_i x_i v_i||^2 <= lam_max(sum v_i v_i^T)=1，且 alpha=x^TQx<=1，')
print('            而 rho = (alpha + x^T|Q|x)/2 <= 1。（k>=3 个团的「聚簇」并不给出可传递符号模式——这是上一版检验的 bug）')
rng=np.random.default_rng(888); mx=-1e9; ntr=0; nskip=0
for t in range(4000):
    N=int(rng.integers(3,9)); k=int(rng.integers(2,5)); lab=rng.integers(0,k,size=N)
    dirs=rng.normal(size=(k,3)); dirs/=np.linalg.norm(dirs,axis=1)[:,None]
    U=np.array([dirs[l]+1e-3*rng.normal(size=3) for l in lab]); U/=np.linalg.norm(U,axis=1)[:,None]
    g=rng.dirichlet(np.ones(N))*3.0; U,g=whiten(U,g)
    Kf=U@U.T
    if transitive(Kf) is None: nskip+=1; continue
    C=core(U,g); ntr+=1; mx=max(mx,C['rho'])
print('  可传递构型 %d 个（其余 %d 个被排除）：max rho = %.12f'%(ntr,nskip,mx))
chk(ntr>=50 and mx<=1+1e-9,'可传递符号模式 => rho<=1（%d 个构型，max=%.12f）'%(ntr,mx))
print('  --- 反向：所有 rho>1 的构型都必须【不可传递】（钝图 G_- 不是完全二部） ---')
bad=0; n2=0
for t in range(3000):
    N=int(rng.integers(4,12)); Z=rng.normal(size=(N,3)); U,g=whiten(Z/np.linalg.norm(Z,axis=1)[:,None],rng.dirichlet(np.ones(N))*3.0)
    C=core(U,g)
    if C['rho']>1+1e-9:
        n2+=1
        if transitive(C['Kf']) is not None: bad+=1
print('  rho>1 的随机构型 %d 个，其中可传递者 %d 个'%(n2,bad))
chk(bad==0,'rho>1 => 符号模式不可传递（%d/%d 检验通过）'%(n2-bad,n2))
print('  --- 冠军与 Stiefel 极大元 ---')
Cc=core(Uc,gc)
print('  冠军: rho=%.12f  G_- = C5 + 孤立枢纽（非完全二部）  可传递? %s'%(Cc['rho'],transitive(Cc['Kf']) is not None))
chk(transitive(Cc['Kf']) is None,'冠军的符号模式【不可传递】——这正是 rho0>1 的结构原因；G_- 含奇圈 C5，解释了 rho0 中的 sqrt(5)')
print(); print('第七段汇总：PASS=%d FAIL=%d'%(PASS,FAIL))
