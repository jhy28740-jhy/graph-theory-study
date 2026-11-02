# -*- coding: utf-8 -*-
"""_aud15.py : 修正 _aud6.py 的两处脚本 bug 后重验。
 bug1: 行方程 einsum('ij,i,jl->il') 用 p_i 加权，应为 p_j -> einsum('ij,j,jl->il')
 bug2: 穷举 worst 初始化为 0，把真正的最小余量截掉了（负值全被丢弃）
本脚本同时补齐冠军律的全部结构方程检验（含 f_1>0 的符号约定）。"""
import numpy as np, itertools
from numpy.linalg import eigh, norm
sq5=np.sqrt(5.0)
d0=-25624284539519/53655078400000+sq5/2-41*np.sqrt(3)/120
d=d0+103/9000000+1/3000000
c=(3+sq5)/4-d/2; cw=(9+sq5+2*np.sqrt(35))/18
stw=(3+sq5)/4

def champion():
    b=(7+5*np.sqrt(7))/126; a=1-5*b
    w=np.array([a]+[b]*5)
    A6=np.zeros((6,6))
    for i in range(5): A6[0,i+1]=A6[i+1,0]=1
    for i in range(5):
        j=(i+1)%5; A6[i+1,j+1]=A6[j+1,i+1]=1
    Bq=A6+np.eye(6)
    Mq=np.diag(np.sqrt(w))@Bq@np.diag(np.sqrt(w))
    ev,V=eigh(Mq); o=np.argsort(-ev); ev=ev[o]; V=V[:,o]
    F=V[:,:3]/np.sqrt(w)[:,None]
    if F[:,0].sum()<0: F[:,0]*=-1
    return w,Bq,ev,F

p,Bq,ev,F=champion()
L=np.diag(ev[:3]); M3=ev[:3].sum()
print("=== 冠军律的结构方程（修正加权后）===")
print(f"  商矩阵谱            = {np.round(ev,10)}")
print(f"  Sigma_3 = sum ev[:3] = {M3:.15f}   rho_0 = {cw:.15f}   diff = {M3-cw:+.3e}")
G=F@F.T
print(f"  (a) sum_i p_i v_i v_i^T = I_3      dev = {norm(np.einsum('i,ij,ik->jk',p,F,F)-np.eye(3)):.3e}")
print(f"  (b) v_i^T L v_i = M (all i)        dev = {np.max(np.abs(np.einsum('ij,jk,ik->i',F,L,F)-M3)):.3e}")
print(f"  (c) v_i1 > 0 (all i)               ->  {np.all(F[:,0]>0)}   col1 = {np.round(F[:,0],8)}")
print(f"  (d) 1{{v_i.v_j>0}} = W_ij            ->  {np.array_equal((G>0).astype(int),Bq.astype(int))}")
print(f"      min |v_i.v_j| offdiag          =   {np.min(np.abs(G[~np.eye(6,dtype=bool)])):.10f}  (远离 0)")
lhs=np.einsum('ij,j,jl->il',Bq.astype(float),p,F)   # 正确: 按 p_j 加权
print(f"  (e) sum_j p_j B_ij v_j = L v_i     dev = {np.max(np.abs(lhs-F@L)):.3e}   <- _aud6 此处为 1.675e-01 (bug)")
lhs_bug=np.einsum('ij,i,jl->il',Bq.astype(float),p,F)
print(f"      (复现 _aud6 的错误加权)        dev = {np.max(np.abs(lhs_bug-F@L)):.3e}")
Ri=norm(F,axis=1)
AA=(p[:,None]*p[None,:]*np.abs(G)).sum(); mm=(p[:,None]*F).sum(0); BB=mm@mm
print(f"  (f) (A+B)/2 = {((AA+BB)/2):.15f}  = Sigma_3 ? {abs((AA+BB)/2-M3)<1e-12}   => A+B = 2M 成立")
print(f"  (g) 正团线性无关性: 各正团 J 的 rank{{v_i v_i^T}}")
Pm=np.einsum('ij,ik->ijk',F,F).reshape(6,9)
for r in range(2,7):
    for J in itertools.combinations(range(6),r):
        sub=(G[np.ix_(J,J)]>0)
        if sub.all():
            rk=np.linalg.matrix_rank(Pm[list(J)],tol=1e-9)
            flag="OK" if rk==r else "!! 违反 |J|<=6 / 线性无关"
            if r<=6 and rk!=r: print(f"      J={J} |J|={r} rank={rk}  {flag}")
print("      全部正团的 rank == |J| : 已验证 (无输出即全部通过)")

print("\n=== 穷举 n=3..6 全部图 (修正 worst 初始化) ===")
def Lam3(Am):
    e=eigh(Am)[0]; return e[-1]+e[-2]+e[-3]
worst=-1e18; argn=None; argA=None; cnt=0
for n in range(3,7):
    verts=list(itertools.combinations(range(n),2))
    for mask in range(1<<len(verts)):
        Am=np.zeros((n,n))
        for idx,(i,j) in enumerate(verts):
            if mask>>idx & 1: Am[i,j]=Am[j,i]=1
        v=Lam3(Am)-(c*n-3); cnt+=1
        if v>worst: worst=v; argn=n; argA=Am.copy()
print(f"  检验图数 = {cnt}   (2^0+2^1+2^3+2^6+2^10+2^15 = {1+2+8+64+1024+32768})")
print(f"  max (Lambda_3 - (c n - 3)) = {worst:+.10f}   {'<=0 => 无违例' if worst<=0 else '!! 违例'}")
print(f"  最小余量 = {-worst:.10f}  在 n={argn}")
print(f"  最紧图的邻接矩阵:\n{argA.astype(int)}")
print(f"  该图 Lambda_3 = {Lam3(argA):.10f} ,  c*n-3 = {c*argn-3:.10f} ,  STW: (3+sqrt5)/4*n-3 = {stw*argn-3:.10f}")
print(f"  相对 STW 的额外余量 = {(stw-c)*argn:.10f}")
