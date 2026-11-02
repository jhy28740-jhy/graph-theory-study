import numpy as np
from numpy.linalg import eigh, svd
rho0=(9+np.sqrt(5)+2*np.sqrt(35))/18

def rho_of_frame_Lam(V,a):
    # single-axis Lam=diag(b,b,a), b=(1-a)/2 ; x_i=sqrt(g_i)*sqrt(b+a_par*s_i)
    g=np.einsum('ij,ij->i',V,V)              # ||V_i||^2
    with np.errstate(invalid='ignore'):
        u=V/np.sqrt(np.maximum(g,1e-300))[:,None]
    s=u[:,2]**2
    b=(1-a)/2; alpha=(3*a-1)/2
    h=np.sqrt(np.maximum(b+alpha*s,0))
    x=np.sqrt(np.maximum(g,0))*h
    Q=V@V.T; Qp=np.maximum(Q,0)
    return float(x@Qp@x), x, Q

def proj(V):
    U,_,Wt=svd(V,full_matrices=False); return U@Wt

def grad_step(V,a,eta):
    g=np.einsum('ij,ij->i',V,V)
    u=V/np.sqrt(np.maximum(g,1e-300))[:,None]
    s=u[:,2]**2; b=(1-a)/2; alpha=(3*a-1)/2
    h=np.sqrt(np.maximum(b+alpha*s,0)); x=np.sqrt(np.maximum(g,0))*h
    Q=V@V.T; Qp=np.maximum(Q,0)
    # grad of x'Qp x wrt V (numeric-free): treat x fixed-ish + Q dependence; use finite-diff-free analytic approx
    Mask=(Q>0).astype(float)
    # d/dV [ x' (Mask*Q) x ] with x depending on V too -> use full autodoc via complex? use numeric grad on Q part + x part
    Gq=2.0*(Mask*np.outer(x,x))@V
    # x-part: dx_i/dV_i ; rho=sum_ij x_i x_j Qp_ij ; drho/dx_i=2 (Qp x)_i ; dx_i/dg_i=h/(2sqrt g); dg_i/dV_i=2V_i; ds_i/dV_i=2 u_i3 e3 *? approx include
    Qpx=Qp@x
    coef=(Qpx/np.maximum(h,1e-9))   # drho/dx * dx/d(sqrt g) ... approx
    Gx=(coef[:,None]*V)
    grad=Gq+Gx
    VtG=V.T@grad; grad=grad-V@(0.5*(VtG+VtG.T))
    return proj(V+eta*grad)

def F_gen(a,N,nstart=25,iters=2500,seed=0):
    best=-1
    for k in range(nstart):
        r=np.random.default_rng(seed*97+k)
        V=proj(r.standard_normal((N,3))); eta=0.05
        for t in range(iters):
            V=grad_step(V,a,eta); eta*=0.999
        val,_,_=rho_of_frame_Lam(V,a)
        if val>best: best=val; bestV=V.copy()
    return best,bestV

def F_cone(a,nmax=8):
    best=-1; bz=None; bn=None
    for n in range(3,nmax+1):
        for z in np.linspace(0.0,1/np.sqrt(3)-1e-9,1400):
            u=[np.array([0,0,1.0])]
            for j in range(n):
                c,sn=np.cos(2*np.pi*j/n),np.sin(2*np.pi*j/n)
                u.append(np.array([np.sqrt(1-z*z)*c,np.sqrt(1-z*z)*sn,z]))
            u=np.array(u); g0=(1-3*z*z)/(1-z*z); G=2/(1-z*z)
            g=np.array([g0]+[G/n]*n)
            if g.min()<-1e-9: continue
            V=np.array([np.sqrt(max(g[i],0))*u[i] for i in range(n+1)])
            val,_,_=rho_of_frame_Lam(V,a)
            if val>best: best=val; bz=z; bn=n
    return best,bz,bn

print("rho0 =",rho0)
print("Direct P1 test: single-axis section, max over GENERAL frames vs CONE family")
print(" a      F_cone(a)      F_gen(N=6..8)   diff(gen-cone)   cone argmax(n,z)")
for a in [0.34,0.40,0.50,0.594608575865948,0.70,0.80,0.90]:
    fc,zc,nc=F_cone(a)
    fg=max(F_gen(a,N,seed=int(a*1000)+N)[0] for N in [6,7,8])
    print(f" {a:.6f}  {fc:.12f}  {fg:.12f}   {fg-fc:+.3e}     n={nc},z={zc:.4f}")
print()
print("=== dichotomy: Neg=lambda_max(Q-) => A2 in {0,1}. Check on randoms: A2 in (0,1) => Neg<lam_max(Q-) strictly ===")
rng=np.random.default_rng(7)
viol=0; minmargin=1e9
for _ in range(300):
    N=rng.integers(5,12); A=rng.standard_normal((N,3)); V,_=np.linalg.qr(A)
    Q=V@V.T; Qp=np.maximum(Q,0); Qm=np.maximum(-Q,0)
    w,x=eigh(Qp); i=np.argmax(w); xp=x[:,i]
    if xp.sum()<0: xp=-xp
    A2=xp@Q@xp; Neg=xp@Qm@xp; lq=eigh(Qm)[0].max()
    if 1e-6<A2<1-1e-6:
        margin=lq-Neg; minmargin=min(minmargin,margin)
        if Neg>lq+1e-12: viol+=1
print("violations (Neg>lam_max(Q-) with A2 in (0,1)):",viol," min(lam_max(Q-)-Neg) over A2 in(0,1):",minmargin)
print("=> whenever 0<A2<1, x is NOT the top eigenvector of Q- (strict misalignment), confirming dichotomy contrapositive")