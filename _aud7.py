import numpy as np
from numpy.linalg import eigh
from scipy.optimize import minimize
sq5=np.sqrt(5); cw=(9+sq5+2*np.sqrt(35))/18; c=(3+sq5)/4-(-25624284539519/53655078400000+sq5/2-41*np.sqrt(3)/120+103/9000000+1/3000000)/2
def unpack(theta,k):
    Y=theta[:3*k].reshape(k,3); z=theta[3*k:]-theta[3*k:].max()
    w=np.exp(z); w/=w.sum()
    S=np.einsum('i,ij,ik->jk',w,Y,Y)
    try: ev,U=eigh(S)
    except Exception: return None,None
    ev=np.maximum(ev,1e-10)
    if not np.all(np.isfinite(ev)): return None,None
    return w,Y@(U@np.diag(1/np.sqrt(ev))@U.T)
def neg(theta,k):
    w,V=unpack(theta,k)
    if V is None or not np.all(np.isfinite(V)): return 1e3
    G=V@V.T
    A=(w[:,None]*w[None,:]*np.abs(G)).sum(); m=(w[:,None]*V).sum(0)
    v=(A+m@m)/2
    return -v if np.isfinite(v) else 1e3
rng=np.random.default_rng(11)
for k in [10,12,14,16,20]:
    loc=(-1,None)
    for trial in range(70):
        th0=np.concatenate([rng.normal(size=3*k)*1.2, rng.normal(size=k)*0.2])
        try:
            r=minimize(neg,th0,args=(k,),method='L-BFGS-B',options=dict(maxiter=3000,ftol=1e-15,gtol=1e-11))
            r2=minimize(neg,r.x,args=(k,),method='Powell',options=dict(maxiter=40000,maxfev=40000,xtol=1e-12,ftol=1e-14))
            cand=r2 if r2.fun<r.fun else r
        except Exception: continue
        if np.isfinite(cand.fun) and -cand.fun>loc[0]: loc=(-cand.fun,cand.x)
    print(f"  k={k:3d}  max (A+B)/2 = {loc[0]:.12f}   gap to c_w = {loc[0]-cw:+.3e}   gap to c = {c-loc[0]:+.3e}")
print(f"\nc_w = {cw:.15f}    c = {c:.15f}")
