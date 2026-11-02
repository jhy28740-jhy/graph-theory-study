import numpy as np
from numpy.linalg import eigh, svd, matrix_rank
rng=np.random.default_rng(2026)
def rand_proj(N):
    A=rng.standard_normal((N,3)); Q,_=np.linalg.qr(A)  # N x 3 orthonormal
    V=Q; P=V@V.T; return P
print("Universal Perron-obtuse identities on random rank-3 projections Q (N=5..12):")
print("  identity1: ||Q-x||^2 =? rho(1-rho)+Neg(2rho-1)")
print("  identity2: ||Q-x||^2 =? Neg^2 + A2(1-A2)")
print("  A2+Neg =? rho ; and rho<=1+Neg (Weyl)")
maxerr=0; rows=[]
for trial in range(400):
    N=rng.integers(5,13); Q=rand_proj(N)
    Qp=np.maximum(Q,0); Qm=np.maximum(-Q,0)
    wv,x=eigh(Qp); i=np.argmax(wv); xp=x[:,i]
    if xp.sum()<0: xp=-xp
    rho=wv[i]; Neg=xp@Qm@xp; A2=xp@Q@xp
    beta=float(xp@(Qm@Qm)@xp)
    id1=rho*(1-rho)+Neg*(2*rho-1)
    id2=Neg*Neg+A2*(1-A2)
    e=max(abs(beta-id1),abs(beta-id2),abs(A2+Neg-rho))
    maxerr=max(maxerr,e)
    rows.append((rho,Neg,A2,beta))
print("  max abs error over 400 random projections:",maxerr)
rows=np.array(rows)
print("  rho range:",rows[:,0].min(),"..",rows[:,0].max())
# check necessary condition Neg >= rho(rho-1)/(2rho-1) for rho>1
sub=rows[rows[:,0]>1.0]
if len(sub):
    lhs=sub[:,1]; rhs=sub[:,0]*(sub[:,0]-1)/(2*sub[:,0]-1)
    print("  necessary cond Neg>=rho(rho-1)/(2rho-1) holds for all",len(sub),"cases with rho>1?",bool((lhs>=rhs-1e-12).all()))
    print("  min slack:",float((lhs-rhs).min()))
print()
print("Rank-1 specialization: if A2==1 then Q- x == Neg*x (shared eigenvector).")
print("  Test on champion-like single-axis: A2=%.9f (not 1) => spread=A2(1-A2)=%.6e"%(0.9765952590554369,0.9765952590554369*(1-0.9765952590554369)))
print("  For A2=1 exactly: spread=0 => Q-x=Neg x. (algebraic, see identity2)")