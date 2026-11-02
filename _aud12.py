import numpy as np
from numpy.linalg import eigh, eigvalsh
sq5=np.sqrt(5); cw=(9+sq5+2*np.sqrt(35))/18
d0=-25624284539519/53655078400000+sq5/2-41*np.sqrt(3)/120
d=d0+103/9000000+1/3000000; c=(3+sq5)/4-d/2; stw=(3+sq5)/4
a=(91-25*np.sqrt(7))/126; b=(7+5*np.sqrt(7))/126
print(f"template masses: a={a:.9f} b={b:.9f}  a+5b={a+5*b:.12f}")
def build(sizes):
    n=sum(sizes); off=[0]
    for s in sizes: off.append(off[-1]+s)
    A=np.zeros((n,n),dtype=np.int8)
    tmpl=np.zeros((6,6),dtype=int)
    for i in range(5): tmpl[0,i+1]=tmpl[i+1,0]=1
    for i in range(5):
        j=(i+1)%5; tmpl[i+1,j+1]=tmpl[j+1,i+1]=1
    for i in range(6):
        for j in range(6):
            if tmpl[i,j] or i==j:      # diagonal block = clique  (W=1 on I_i x I_i)
                sl_i=slice(off[i],off[i+1]); sl_j=slice(off[j],off[j+1])
                A[sl_i,sl_j]=1
    for i in range(n): A[i,i]=0
    return A
print("\n n    sizes                       Lambda_3      Lambda_3/n     c*n-3      margin    STW margin")
for n in [126,252,504,1008,2016]:
    raw=np.array([a]+[b]*5)*n
    sizes=np.floor(raw).astype(int); sizes[0]+= n-sizes.sum()
    A=build(sizes); ev=eigvalsh(A); L3=ev[-1]+ev[-2]+ev[-3]
    print(f"{n:5d}  {str(list(sizes)):26s} {L3:11.5f}  {L3/n:12.9f}  {c*n-3:10.4f}  {c*n-3-L3:+9.4f}  {stw*n-3-L3:+9.4f}")
print(f"\n c_w = {cw:.12f}   c = {c:.12f}   c - c_w = {c-cw:.9f}")
print(" => Lambda_3/n -> c_w, so the paper's coefficient c exceeds the true asymptotic value by", f"{c-cw:.6f}")
print(" => margin (c n - 3) - Lambda_3 grows like (c-c_w)*n =", f"{(c-cw)*2016:.3f} at n=2016")
