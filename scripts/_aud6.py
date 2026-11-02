import numpy as np, itertools
from numpy.linalg import eigh, norm, matrix_rank
sq5=np.sqrt(5)
d0=-25624284539519/53655078400000+sq5/2-41*np.sqrt(3)/120
d=d0+103/9000000+1/3000000
c=(3+sq5)/4-d/2
cw=(9+sq5+2*np.sqrt(35))/18

# ---------- 1. champion 6-part K1 v C5 ----------
b=(7+5*np.sqrt(7))/126; a=1-5*b
w=np.array([a]+[b]*5)
A6=np.zeros((6,6))
for i in range(5): A6[0,i+1]=A6[i+1,0]=1
for i in range(5):
    j=(i+1)%5; A6[i+1,j+1]=A6[j+1,i+1]=1
Bq=A6+np.eye(6)
M=np.diag(np.sqrt(w))@Bq@np.diag(np.sqrt(w))
ev,V=eigh(M); o=np.argsort(-ev); ev=ev[o]; V=V[:,o]
Z=V[:,:3]; Lam=ev[:3]
print("champion: eigs =",np.round(ev,10))
print("  Lambda_3(M) =",f"{ev[:3].sum():.15f}"," c_w =",f"{cw:.15f}"," diff=",f"{ev[:3].sum()-cw:.3e}")
F=Z/np.sqrt(w)[:,None]           # canonical frame rows F(i)
print("  isotropy ||sum p FF^T - I|| =",f"{norm(np.einsum('i,ij,ik->jk',w,F,F)-np.eye(3)):.3e}")
G=F@F.T
print("  sign pattern 1{F_i.F_j>0} equals W_ij ?", np.array_equal((G>1e-12).astype(int),Bq.astype(int)))
print("  min |F_i.F_j| offdiag =",f"{np.min(np.abs(G[~np.eye(6,dtype=bool)])):.10f}")
L=np.diag(Lam)
print("  v_i^T L v_i == M ?  max dev =",f"{np.max(np.abs(np.einsum('ij,jk,ik->i',F,L,F)-ev[:3].sum())):.3e}")
Bm=(G>0).astype(float)
print("  sum_j p_j B_ij v_j == L v_i ?  max dev =",f"{np.max(np.abs(np.einsum('ij,i,jl->il',Bm,w,F)-F@L)):.3e}")
print("  F_{i1} > 0 all?", np.all(F[:,0]>0), " F col1 =",np.round(F[:,0],8))
# A+B = 2M check for canonical law
Ri=norm(F,axis=1); p=w
AA=(p[:,None]*p[None,:]*np.abs(G)).sum(); mm=(p[:,None]*F).sum(0); BB=mm@mm
print("  (A+B)/2 =",f"{(AA+BB)/2:.15f}","  M =",f"{ev[:3].sum():.15f}","  D_* =",f"{(3+sq5)/2-AA-BB:.12f}",
      "  vs 1/20 =",0.05, " -> branch:", "LARGE" if (3+sq5)/2-AA-BB>=0.05 else "small")
print("  D_* >= d ?", (3+sq5)/2-AA-BB>=d)

# ---------- 2. finite graph bound Lambda_3(G) <= c n - 3 ----------
def Lam3(Am):
    e=eigh(Am)[0]; return e[-1]+e[-2]+e[-3]
print("\n=== exhaustive n<=6 and random/handpicked graphs: Lambda_3(G) <= c*n-3 ===")
worst=0; arg=None
for n in range(3,7):
    verts=list(itertools.combinations(range(n),2))
    for mask in range(1<<len(verts)):
        Am=np.zeros((n,n))
        for idx,(i,j) in enumerate(verts):
            if mask>>idx & 1: Am[i,j]=Am[j,i]=1
        v=Lam3(Am)-(c*n-3)
        if v>worst: worst=v; arg=(n,Am.copy())
print(f"  exhaustive n=3..6 : max (Lambda_3 - (c n - 3)) = {worst:+.6e}  (<=0 means bound holds)")
if arg is not None: print("   at n=",arg[0])
# larger structured graphs
def test(Am,name):
    n=Am.shape[0]; v=Lam3(Am); b=c*n-3; s=(3+sq5)/4*n-3
    print(f"  {name:28s} n={n:4d} Lambda_3={v:10.5f}  c n-3={b:10.5f}  margin={b-v:+9.5f}   STW margin={s-v:+9.5f}")
    return b-v
g=[]
n=6; Am=np.zeros((n,n))
for i in range(5): Am[0,i+1]=Am[i+1,0]=1
for i in range(5):
    j=(i+1)%5; Am[i+1,j+1]=Am[j+1,i+1]=1
g.append(test(Am,"K1 v C5 (champion graph)"))
for n in [10,30,60,120,300,600]:
    Am=np.triu(np.ones((n,n)),1); Am=Am+Am.T
    g.append(test(Am,f"K_{n}"))
for n in [12,30,60,121,300]:
    r=int(np.sqrt(n)); 
    if r*r!=n: continue
    Am=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if (i//r==j//r) or (i%r==j%r):
                if i!=j: Am[i,j]=1
    g.append(test(Am,f"rook {r}x{r}"))
rng=np.random.default_rng(3)
for n,pr in [(60,0.5),(200,0.5),(200,0.1),(500,0.05),(400,0.3)]:
    U=rng.random((n,n)); Am=(U<pr); Am=np.triu(Am,1); Am=(Am|Am.T).astype(float)
    g.append(test(Am,f"G(n,{pr})"))
for n in [30,60,120,240]:
    Am=np.zeros((n,n))
    h=n//2
    for i in range(h):
        for j in range(h,n): Am[i,j]=Am[j,i]=1
    g.append(test(Am,f"K_{{{h},{n-h}}}"))
# cycle/petersen/clebsch-like
C=np.zeros((10,10))
for i in range(10): C[i,(i+1)%10]=C[(i+1)%10,i]=1
pet=np.zeros((10,10))
for i in range(5):
    pet[i,(i+1)%5]=pet[(i+1)%5,i]=1; pet[i+5,(i+2)%5+5]=pet[(i+2)%5+5,i+5]=1; pet[i,i+5]=pet[i+5,i]=1
g.append(test(C+pet,"Petersen"))
print(f"\n  overall worst margin over tested graphs = {min(g):+.6e}  (>=0 means bound never violated)")
print("  is c n - 3 ever tight? min margin =",f"{min(g):.6f}")

