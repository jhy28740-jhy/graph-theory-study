# -*- coding: utf-8 -*-
# _r28b.py  正 n 边形环的 P 值（精确 + 数值），验证 n=5 在环族中最大
import numpy as np, sympy as sp
sqrt5=np.sqrt(5.0)
Pchamp=(3+sqrt5)/5
def P_ring_exact(n):
    # rho_i = 2/n, ang_k = 2 pi k /n ; P = sum_{i!=j, cos<0} (2/n)^2 cos^2
    k=sp.symbols('k',integer=True)
    tot=sp.Integer(0)
    for d in range(1,n):   # 差 d=i-j (mod n)，每个 d 出现 n 次（有序对）
        c=sp.cos(2*sp.pi*d/n)
        if sp.simplify(c)<0:
            tot+= n*(sp.Rational(2,n)**2)*c**2
    return sp.simplify(tot), sp.nsimplify(tot)
def P_ring_num(n):
    ang=2*np.pi*np.arange(n)/n; rho=np.full(n,2.0/n)
    c=np.cos(ang[:,None]-ang[None,:])
    M=(rho[:,None]*rho[None,:])*np.where(c<0,c**2,0.0)
    return float(M.sum())
print("="*80)
print("正 n 边形环的 P 值（P = sum_{i!=j,c<0} rho_i rho_j c_ij^2, rho=2/n）")
print("  冠军 5-环 P = (3+sqrt5)/5 = %.12f"%Pchamp)
print("="*80)
best=-1;bestn=None
for n in range(3,25):
    pn=P_ring_num(n)
    ex,_=P_ring_exact(n) if n<=12 else (None,None)
    exs=str(ex) if ex is not None else ""
    mark=""
    if pn>best+1e-12: best=pn;bestn=n;mark="  <-- 环族最大"
    print("  n=%2d  P=%.12f  精确=%s%s"%(n,pn,exs,mark))
print("-"*80)
print("  环族最大 P = %.12f 在 n=%d"%(best,bestn))
print("  与 (3+sqrt5)/5 的差 = %.3e"%(best-Pchamp))
# 精确核对 n=3,4,5,6
print()
print("精确值核对:")
for n in [3,4,5,6,7,8]:
    ex,ns=P_ring_exact(n)
    print("  n=%d: P=%s = %s = %.15f"%(n,ex,ns,float(ex.evalf(30))))
print("  (3+sqrt5)/5 = %.15f"%Pchamp)
print("  n=5 是否等于 (3+sqrt5)/5:", sp.simplify(P_ring_exact(5)[0]-(3+sp.sqrt(5))/5)==0)
print("ALL R28B DONE")
