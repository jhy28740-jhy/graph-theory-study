import numpy as np
from numpy.linalg import eigvalsh, eigh
np.set_printoptions(precision=12, suppress=False)
rho0=(9+np.sqrt(5)+2*np.sqrt(35))/18
print("rho0 =",repr(rho0), " rho0-1 =",repr(rho0-1))

# --- Sanity: champion wheel W6 = K1 v C5 ---
A=np.zeros((6,6))
for j in range(5): A[0,1+j]=A[1+j,0]=1
for j in range(5):
    k=(j+1)%5; A[1+j,1+k]=A[1+k,1+j]=1
ev=np.sort(eigvalsh(A))[::-1]; Lam3=ev[:3].sum()
print("W6 adj eigs:",np.round(ev,9))
print("Lam3(W6)=",repr(Lam3)," (Lam3+3)/6 =",repr((Lam3+3)/6)," == rho0?",abs((Lam3+3)/6-rho0)<1e-12)

def cone_frame(n,z):
    u=[np.array([0,0,1.0])]
    for j in range(n):
        c,s=np.cos(2*np.pi*j/n),np.sin(2*np.pi*j/n)
        u.append(np.array([np.sqrt(1-z*z)*c,np.sqrt(1-z*z)*s,z]))
    u=np.array(u); g0=(1-3*z*z)/(1-z*z); G=2/(1-z*z)
    g=np.array([g0]+[G/n]*n); return u,g

zstar=0.500440049768212308
u,g=cone_frame(5,zstar)
Iso=sum(g[i]*np.outer(u[i],u[i]) for i in range(6))
print("isotropy err:",np.abs(Iso-np.eye(3)).max()," sum g =",g.sum())
V=np.array([np.sqrt(g[i])*u[i] for i in range(6)])
Q=V@V.T
print("Q proj err:",np.abs(Q@Q-Q).max()," trQ=",np.trace(Q))
print("Q eigs:",np.round(np.sort(eigvalsh(Q))[::-1],9))
Qp=np.maximum(Q,0); Qm=np.maximum(-Q,0)
rhoQp=eigvalsh(Qp).max(); lamQm=eigvalsh(Qm).max()
print("rho(Q+)=",repr(rhoQp)," ==rho0?",abs(rhoQp-rho0)<1e-9)
print("lam_max(Q-)=",repr(lamQm)," rho0-1=",repr(rho0-1))
print("Weyl 1+lam_max(Q-)=",repr(1+lamQm)," >= rho0?",1+lamQm>=rho0-1e-12)
w,x=eigh(Qp); xp=x[:,-1]
if xp.sum()<0: xp=-xp
Neg=xp@Qm@xp; A2=xp@Q@xp
print("Perron x=",np.round(xp,9))
print("A2=x'Qx=",repr(A2)," Neg=x'Q-x=",repr(Neg)," A2+Neg=",repr(A2+Neg))
print("check A2+Neg==rho(Q+)?",abs(A2+Neg-rhoQp)<1e-9)