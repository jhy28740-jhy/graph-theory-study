import numpy as np, sympy as sp
from numpy.linalg import eigvalsh, eigh
rho0=(9+np.sqrt(5)+2*np.sqrt(35))/18
z=0.500440049768212308
g0=(1-3*z*z)/(1-z*z); gc=2/(5*(1-z*z))
u=[np.array([0,0,1.0])]
for j in range(5):
    c,s=np.cos(2*np.pi*j/5),np.sin(2*np.pi*j/5)
    u.append(np.array([np.sqrt(1-z*z)*c,np.sqrt(1-z*z)*s,z]))
u=np.array(u); g=np.array([g0]+[gc]*5)
V=np.array([np.sqrt(g[i])*u[i] for i in range(6)]); Q=V@V.T
Qp=np.maximum(Q,0); Qm=np.maximum(-Q,0)
w,x=eigh(Qp); xp=x[:,-1]
if xp.sum()<0: xp=-xp
rho=w.max(); Neg=xp@Qm@xp; A2=xp@Q@xp
beta= float(xp@(Qm@Qm)@xp)
print("=== beta-identity check ===")
print("rho=",rho," Neg=",Neg," A2=",A2," A2+Neg=",A2+Neg)
print("||Q- x||^2 (beta) numeric =",beta)
print("rho(1-rho)+Neg(2rho-1)     =",rho*(1-rho)+Neg*(2*rho-1))
print("match?",abs(beta-(rho*(1-rho)+Neg*(2*rho-1)))<1e-12)
print("necessary cond: Neg >= rho(rho-1)/(2rho-1) =",rho*(rho-1)/(2*rho-1),"  actual Neg=",Neg)
print("=== Q- structure: is it w*A(C5) on rim, apex isolated? ===")
print("Q- =\n",np.round(Qm,6))
print("apex row of Q- (should be 0):",np.round(Qm[0],9))
print("=== rho(|Q|) route ===")
absQ=np.abs(Q); print("rho(|Q|)=",eigvalsh(absQ).max()," vs rho0=",rho0," |Q|>=Q+ so >=rho0? ",eigvalsh(absQ).max()>=rho0-1e-12)
print("=== diag-split route: max g_i + rho(offdiag Q+) ===")
offQp=Qp-np.diag(np.diag(Qp)); print("max g_i=",g.max()," rho(offdiag Q+)=",eigvalsh(offQp).max()," sum=",g.max()+eigvalsh(offQp).max())
print("=== moment-dual gap diagnosis (P1') ===")
# nu-relaxation: max int K dnu, nu>=0, nuhat(0)=4, nuhat(2)=0.  K=(cos)_-^2
# candidate nu=2 delta_pi + 2 delta_{pi/2}: nuhat(2)=2 e^{2i pi}+2 e^{i pi}=2-2=0 OK
K=lambda th: max(np.cos(th),0.0)**0  # placeholder
def Kfun(th):
    c=np.cos(th); return (min(c,0.0))**2
pts=[np.pi,np.pi/2]; mass=[2,2]
print("nu=2d_pi+2d_(pi/2): int K dnu=",sum(m*Kfun(t) for m,t in zip(mass,pts))," nuhat2=",abs(sum(m*np.exp(2j*t) for m,t in zip(mass,pts))))
print("true sup P=(3+sqrt5)/5=",(3+np.sqrt(5))/5," => nu-LP relaxation loose by factor",2/((3+np.sqrt(5))/5))