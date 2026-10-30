import numpy as np
from numpy.linalg import eigvalsh, eigh, qr, svd
rng=np.random.default_rng(12345)
rho0=(9+np.sqrt(5)+2*np.sqrt(35))/18

def rho_Qp(V):
    Q=V@V.T; Qp=np.maximum(Q,0.0)
    w,x=eigh(Qp); i=np.argmax(w); xp=x[:,i]
    if xp.sum()<0: xp=-xp
    return w[i],xp,Q

def proj_stiefel(V):
    U,s,Wt=svd(V,full_matrices=False)
    return U@Wt   # polar -> nearest orthonormal cols

def ascend(N,iters=4000,seed=0):
    r=np.random.default_rng(seed)
    V=proj_stiefel(r.standard_normal((N,3)))
    best=-1; bestV=None
    eta=0.05
    for t in range(iters):
        rho,x,Q=rho_Qp(V)
        Mask=(Q>0).astype(float)
        grad=2.0*(Mask*np.outer(x,x))@V
        # project gradient to tangent of Stiefel: grad - V sym(V^T grad)
        VtG=V.T@grad; grad=grad-V@(0.5*(VtG+VtG.T))
        V=proj_stiefel(V+eta*grad)
        eta*=0.9995
        if rho>best: best=rho; bestV=V.copy()
    return best,bestV

print("rho0 =",rho0)
print("="*80)
for N in [4,5,6,7,8,9]:
    bests=[]
    for s in range(40):
        b,V=ascend(N,seed=1000*s+7)
        bests.append((b,V))
    bests.sort(key=lambda t:-t[0])
    top,bV=bests[0]
    # how many starts reach within 1e-6 of top?
    nclose=sum(1 for b,_ in bests if top-b<1e-6)
    print(f"N={N}: max rho(Q+) over 40 starts = {top:.12f}  (rho0={rho0:.12f}, diff={top-rho0:+.3e})  #near-top={nclose}")
    # structure of optimizer: distinct row-norms (masses g_i) and distinct |z| (latitudes)
    g=np.array([bV[i]@bV[i] for i in range(N)])
    zn=np.array([abs(bV[i,2])/np.sqrt(g[i]) if g[i]>1e-9 else 0 for i in range(N)])
    print("   masses g_i:",np.round(np.sort(g)[::-1],6))
    print("   |cos lat| :",np.round(np.sort(zn)[::-1],6))
    print("   exceeds rho0 by >1e-7?", top>rho0+1e-7)
print("="*80)
print("DONE")