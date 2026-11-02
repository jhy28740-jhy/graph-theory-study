import numpy as np, itertools
from numpy.linalg import eigh, norm
sq5=np.sqrt(5)
d0=-25624284539519/53655078400000+sq5/2-41*np.sqrt(3)/120
d=d0+103/9000000+1/3000000
c=(3+sq5)/4-d/2
b=(7+5*np.sqrt(7))/126; a=1-5*b
w=np.array([a]+[b]*5)
A6=np.zeros((6,6))
for i in range(5): A6[0,i+1]=A6[i+1,0]=1
for i in range(5):
    j=(i+1)%5; A6[i+1,j+1]=A6[j+1,i+1]=1
Bq=A6+np.eye(6)
M=np.diag(np.sqrt(w))@Bq@np.diag(np.sqrt(w))
ev,V=eigh(M); o=np.argsort(-ev); ev=ev[o]; V=V[:,o]
Z=V[:,:3]; Lam=ev[:3]; L=np.diag(Lam)
# sign convention: f1 >= 0
F=Z/np.sqrt(w)[:,None]
if F[:,0].sum()<0: F[:,0]*=-1; Z[:,0]*=-1
G=F@F.T
print("=== champion canonical law (structure theorem checks) ===")
print(" f1 = F[:,0] =",np.round(F[:,0],8),"  all > 0 ?",np.all(F[:,0]>0))
print(" isotropy  max dev          =",f"{np.max(np.abs(np.einsum('i,ij,ik->jk',w,F,F)-np.eye(3))):.3e}")
print(" v_i^T L v_i == M           =",f"{np.max(np.abs(np.einsum('ij,jk,ik->i',F,L,F)-Lam.sum())):.3e}")
print(" v_i . v_j != 0 (min |dot|) =",f"{np.min(np.abs(G[~np.eye(6,dtype=bool)])):.8f}")
print(" 1{v_i.v_j>0} == W_ij       =",np.array_equal((G>0).astype(int),Bq.astype(int)))
res=np.einsum('ij,j,jl->il',Bq.astype(float),w,F)-F@L
print(" sum_j p_j B_ij v_j == Lv_i max dev =",f"{np.max(np.abs(res)):.3e}")
print(" {v_i v_i^T} linearly independent? rank =",
      np.linalg.matrix_rank(np.stack([np.outer(F[i],F[i]).ravel() for i in range(6)]),tol=1e-8))
print("   (Thm six says every POSITIVE CLIQUE has <=6 independent rank-1s;")
Bm=Bq.astype(int); maxcl=0
for r in range(1,7):
    for J in itertools.combinations(range(6),r):
        if all(Bm[i,j]==1 for i in J for j in J):
            maxcl=max(maxcl,r)
print(f"   max positive clique size = {maxcl}  (|J|<=6 required: {maxcl<=6}) )")
for J in [(0,1,2),(0,1,2,3),(0,1,2,3,4),(0,1,2,3,4,5)]:
    Mx=np.stack([np.outer(F[i],F[i]).ravel() for i in J])
    print(f"   clique {J}: rank of rank-1 family = {np.linalg.matrix_rank(Mx,tol=1e-8)} / |J|={len(J)}")
print("\n=== graphon Sigma_3 vs law value ===")
AA=(w[:,None]*w[None,:]*np.abs(G)).sum(); m=(w[:,None]*F).sum(0); BB=m@m
print(f" A=E|X.X'| = {AA:.12f}   B=|EX|^2 = {BB:.12f}   (A+B)/2 = {(AA+BB)/2:.15f}")
print(f" M=Sigma_3(W) = {Lam.sum():.15f}   c_w = {(9+sq5+2*np.sqrt(35))/18:.15f}")
print(f" A+B == 2M ? dev = {abs(AA+BB-2*Lam.sum()):.3e}")
Ds=(3+sq5)/2-AA-BB
print(f" D_* = {Ds:.12f} ; 1/20 = 0.05 -> branch = {'LARGE (trivial)' if Ds>=0.05 else 'small'} ; d = {d:.12f} ; D_*>=d: {Ds>=d}")
print(f" radius ratio = {norm(F,axis=1).max()/norm(F,axis=1).min():.8f}   (need < sqrt8 = {np.sqrt(8):.6f})")
print(f" lambda_3 = {Lam[2]:.8f} > 44/375 = {44/375:.8f} ? {Lam[2]>44/375}")
print(f" lambda_1 = {Lam[0]:.8f} < 15/16 = {15/16:.8f} ? {Lam[0]<15/16}")
print(f" sqrt(lam1/lam3) = {np.sqrt(Lam[0]/Lam[2]):.8f} < sqrt(5625/704) = {np.sqrt(5625/704):.8f} ? {np.sqrt(Lam[0]/Lam[2])<np.sqrt(5625/704)}")
