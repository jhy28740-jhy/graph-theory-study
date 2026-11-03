# -*- coding: utf-8 -*-
import re
p='证书_Lam参数化与稀释恒等式.py'
s=open(p,encoding='utf-8').read()
i=s.index("print(); print('='*100); print('(H)")
head=s[:i]
tail="""print(); print('='*100); print('(H) 可传递符号模式定理：若存在 sigma 使 sigma_i sigma_j <u_i,u_j> >= 0，则 rho <= 1'); print('='*100)
print('   证明：此时 sum_ij x_i x_j|<v_i,v_j>| = ||sum_i sigma_i x_i v_i||^2 <= lam_max(sum v_i v_i^T)=1，且 alpha=x^TQx<=1，')
print('         而 rho = (alpha + x^T|Q|x)/2 <= 1。（注意：k>=3 个团的聚簇并不给出可传递符号模式）')
def transitive(Kf):
    N=Kf.shape[0]; best=-1e9
    if N>20: return None
    for msk in range(1<<(N-1)):
        sg=np.array([1.0]+[1.0 if (msk>>k)&1 else -1.0 for k in range(N-1)])
        best=max(best,(np.outer(sg,sg)*Kf).min())
    return best
r4=np.random.default_rng(888); mx=-1e9; ntr=0; nskip=0; wfr=0.0
for t in range(6000):
    N=int(r4.integers(3,9)); k=int(r4.integers(2,5)); lab=r4.integers(0,k,size=N)
    dirs=r4.normal(size=(k,3)); dirs/=np.linalg.norm(dirs,axis=1)[:,None]
    Uu=np.array([dirs[l]+1e-3*r4.normal(size=3) for l in lab]); Uu/=np.linalg.norm(Uu,axis=1)[:,None]
    gg=r4.dirichlet(np.ones(N))*3.0; Uu,gg=whiten(Uu,gg)
    C=core(Uu,gg); fr=np.abs((Uu*gg[:,None]).T@Uu-np.eye(3)).max()
    if fr>1e-11: nskip+=1; continue
    if transitive(C['Kf'])<0.0: continue
    ntr+=1; wfr=max(wfr,fr)
    if C['rho']>mx: mx=C['rho']
print('   严格可传递构型 %d 个（帧残差>1e-11 剔除 %d 个）：max rho = %.15f （超出 1 的量 %+.2e，帧残差上界 %.1e）'%(ntr,nskip,mx,mx-1,wfr))
chk(ntr>=50 and mx<=1+1e-8,'可传递符号模式 => rho<=1（%d 个构型，max=%.12f，超出量与帧残差同阶）'%(ntr,mx))
r5=np.random.default_rng(31337); n2=0; bad=0
for t in range(1500):
    N=int(r5.integers(4,12)); Z=r5.normal(size=(N,3)); Uu,gg=whiten(Z/np.linalg.norm(Z,axis=1)[:,None],r5.dirichlet(np.ones(N))*3.0)
    if np.abs((Uu*gg[:,None]).T@Uu-np.eye(3)).max()>1e-11: continue
    C=core(Uu,gg)
    if C['rho']>1+1e-9:
        n2+=1
        if transitive(C['Kf'])>=0.0: bad+=1
print('   逆否命题：rho>1 的随机构型 %d 个，其中可传递者 %d 个'%(n2,bad))
chk(bad==0,'rho>1 => 符号模式不可传递（钝图 G_- 不是完全二部）： %d/%d 通过'%(n2-bad,n2))
chk(transitive(Cc['Kf'])<0.0,'冠军的符号模式【不可传递】：G_- = C5 + 孤立枢纽，含奇圈 C5——这正是 rho0>1 的结构原因，也解释了 rho0 中的 sqrt(5)')
print(); print('='*100)
print('证书汇总：PASS=%d  FAIL=%d'%(PASS,FAIL))
if FAILED:
    print('未通过项：'); [print('  -',m) for m in FAILED]
print('='*100)
"""
open(p,'w',encoding='utf-8').write(head+tail)
print('patched section H')
