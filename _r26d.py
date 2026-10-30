# -*- coding: utf-8 -*-
import numpy as np
from scipy.optimize import minimize_scalar
def r(s):
    b=1-s
    if s>=1/3-1e-15: th=np.arccos(-1/2)
    else: th=np.arccos(-s/b)
    return (2/np.pi)*(np.sqrt(max(1-2*s,0))+s*th)-s
f=lambda s: ((1-3*s)*(1+s)+4*s*r(s))/(1-s)**2
ss=np.linspace(0,1/3,20001); v=[f(s) for s in ss]; k=int(np.argmax(v))
print("lim_{n->inf} max_s L_w = %.9f at s=%.6f ; 端点 s=1/3 值 = %.9f"%(v[k],ss[k],f(1/3-1e-15)))
print("(13+sqrt5)/10 = %.9f"%((13+np.sqrt(5))/10))
# 环族 L_w = T_n/n 复核
for n in (3,4,5,6,8):
    d=np.arange(n); T=np.abs(2*np.cos(2*np.pi*d/n)+1).sum()
    g=3/n; a=1/np.sqrt(3)
    U=np.array([[np.sqrt(2/3)*np.cos(2*np.pi*k/n),np.sqrt(2/3)*np.sin(2*np.pi*k/n),1/np.sqrt(3)] for k in range(n)])
    w=g*a*np.ones(n); Lw=float((np.outer(w,w)*np.abs(U@U.T)).sum())
    print("  n=%d: T_n=%.9f  T_n/n=%.9f  显式 L_w=%.9f"%(n,T,T/n,Lw))
