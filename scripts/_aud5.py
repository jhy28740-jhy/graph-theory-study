import numpy as np
from numpy.linalg import eigh, norm
from scipy.optimize import minimize
sq5=np.sqrt(5)
d0=-25624284539519/53655078400000+sq5/2-41*np.sqrt(3)/120
d=d0+103/9000000+1/3000000
c=(3+sq5)/4-d/2
cw=(9+sq5+2*np.sqrt(35))/18
stw=(3+sq5)/4
print(f"c={c:.15f}  c_w={cw:.15f}  STW={stw:.15f}")
print(f"d={d:.12f}  2c={2*c:.12f}")

def obj_from_atoms(p,V):
    G=V@V.T
    A=(p[:,None]*p[None,:]*np.abs(G)).sum()
    m=(p[:,None]*V).sum(0); B=m@m
    return (A+B)/2

def whitened(theta,k):
    Y=theta[:3*k].reshape(k,3); z=theta[3*k:]
    w=np.exp(z); w/=w.sum()
    S=np.einsum('i,ij,ik->jk',w,Y,Y)
    ev,U=eigh(S); ev=np.maximum(ev,1e-12)
    V=Y@(U@np.diag(1/np.sqrt(ev))@U.T)
    return w,V

def neg(theta,k):
    w,V=whitened(theta,k)
    return -obj_from_atoms(w,V)

best=(-1,None)
rng=np.random.default_rng(0)
for k in [4,5,6,7,8,9,10,12,14]:
    loc=(-1,None)
    for trial in range(60):
        th0=np.concatenate([rng.normal(size=3*k), rng.normal(size=k)*0.3])
        r=minimize(neg,th0,args=(k,),method='L-BFGS-B',options=dict(maxiter=4000,ftol=1e-16,gtol=1e-12))
        r2=minimize(neg,r.x,args=(k,),method='Nelder-Mead',options=dict(maxiter=20000,maxfev=20000,xatol=1e-12,fatol=1e-14))
        cand=r2 if r2.fun<r.fun else r
        if -cand.fun>loc[0]: loc=(-cand.fun,cand.x)
    print(f"  k={k:3d}  max (A+B)/2 = {loc[0]:.12f}   (gap to c = {c-loc[0]:+.3e}, to c_w = {loc[0]-cw:+.3e})")
    if loc[0]>best[0]: best=loc
print(f"\nBEST over k: {best[0]:.12f}   vs c={c:.12f}  vs c_w={cw:.12f}")
print("  best <= c ? ", best[0]<=c, "   best <= c_w+1e-9 ? ", best[0]<=cw+1e-9)
np.save('_aud_best.npy',best[1])
k=(len(best[1])-0)//4  # not exact; recompute
