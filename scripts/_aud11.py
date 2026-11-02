import numpy as np, sympy as sp
from sympy import Rational as R
t=sp.Symbol('t',real=True)
a3=R(12928,41503); a4=R(10240,41503)
c0=R(8352,41503); c1=R(43806,41503); c2=R(54016,41503)
fp=c0+c1*t+c2*t**2-a3*t**3-a4*t**4
zst=t**4-R(85,64)*t**2+R(15,128)*t+R(27,128)
kap=R(3,125000)
print("=== exact symbolic factorisations (Lemma scalar / contact-square) ===")
print(" 41503*sigma(w), w>=0 :", sp.factor(sp.expand(41503*(fp-2*t))))
print(" 41503*sigma(w), w<=0 :", sp.factor(sp.expand(41503*fp)))
print(" (3/2)*41503*sigma - 41503*z_*^2 , w>=0 :", sp.factor(sp.expand(R(3,2)*41503*(fp-2*t)-41503*zst**2)))
print(" (3/2)*41503*sigma - 41503*z_*^2 , w<=0 :", sp.factor(sp.expand(R(3,2)*41503*fp-41503*zst**2)))
print(" h(w)=sigma-z_*/125+3/125000 as sum of squares (w>=0 branch) identity check:",
      sp.simplify(sp.expand((fp-2*t)-zst/125+kap - ((fp-2*t-R(2,3)*zst**2)+R(2,3)*(zst-R(3,500))**2)))==0)
# numeric scan with pure numpy
c=[8352,43806,54016,-12928,-10240]; den=41503.0
def fnum(g): return np.polyval([x/den for x in c[::-1]],g)
def signum(g): return g**4-(85/64)*g**2+(15/128)*g+27/128
grid=np.linspace(-1,1,2000001)
sg=np.where(grid>=0, fnum(grid)-2*grid, fnum(grid))
gs=(1-9*grid**2)/25
zz=signum(grid)**2
print("\n=== dense numeric scan (2,000,001 points) ===")
print(f"  min sigma(w)                 = {sg.min():+.6e}   (claim >= 0)")
print(f"  min [sigma(w)-(1-9w^2)/25]   = {(sg-gs).min():+.6e}   (claim >= 0)")
print(f"  max [z_*^2 - (3/2)sigma]     = {np.max(zz-1.5*sg):+.6e}   (claim <= 0)")
m=sg>1e-12
print(f"  max z_*^2/sigma on [0,1]     = {np.max((zz/np.where(m,sg,1))[(grid>=0)&m]):.9f}  <= 41503/27840 = {float(R(41503,27840)):.9f} ? {np.max((zz/np.where(m,sg,1))[(grid>=0)&m])<=float(R(41503,27840))+1e-12}")
print(f"  max z_*^2/sigma on [-1,0]    = {np.max((zz/np.where(m,sg,1))[(grid<0)&m]):.9f}  <= 41503/54400 = {float(R(41503,54400)):.9f} ? {np.max((zz/np.where(m,sg,1))[(grid<0)&m])<=float(R(41503,54400))+1e-12}")
print(f"  both < 3/2 ? {np.max(zz[m]/sg[m])<1.5}")
hh=sg-signum(grid)/125+float(kap)
print(f"  min h(w) = sigma - z_*/125 + 3/125000 = {hh.min():+.6e}   (claim >= 0)")
print(f"  h(1) = {np.polyval([x/den for x in c[::-1]],1.0)-2-float(signum(np.array([1.0]))[0])/125+float(kap):+.9e}   kappa = {float(kap):.9e}")
# zeros
print(f"  sigma(1)={float(fnum(np.array([1.0]))[0]-2):.3e}  sigma(1/2)={float(fnum(np.array([0.5]))[0]-1.0):.3e}  sigma(-3/8)={float(fnum(np.array([-0.375]))[0]):.3e}")
# q_* bound in contact-square
qs=(1-t)*(t+R(3,8))**2*(t+R(9,8))**2
gg=np.linspace(0,1,200001); qn=(1-gg)*(gg+3/8)**2*(gg+9/8)**2
print(f"\n  max q_*(t) on [0,1] = {qn.max():.12f}  at t={gg[qn.argmax()]:.6f} ; bound q_*(2/3)+13/13824 = {float(qs.subs(t,R(2,3))+R(13,13824)):.12f} ; 6/5 = 1.2")
print(f"  41503/(4*8352)*max q_* = {41503/(4*8352)*qn.max():.12f}  vs 41503/27840 = {float(R(41503,27840)):.12f}  vs 3/2")
ps=(1+gg)*(gg+0.5)*(9/8-gg)
print(f"  max p_*(x) on [0,1] = {ps.max():.9f} < 1 ? {ps.max()<1}")
print(f"  min 928+82x-160x^2 on [0,1] = {np.min(928+82*gg-160*gg**2):.3f} (>=850)")
# chi identity
p=sp.symbols('p0:6',positive=True); Rr=sp.symbols('R0:6',positive=True)
chiS=sum(p[i]**2*Rr[i]**2 for i in range(6)); CS=sum(p[i]*Rr[i] for i in range(6))
print("\n=== chi identity ===")
print("  chi - C^2/m == sum(p_iR_i - C/m)^2 (m=6) ?", sp.simplify(chiS-CS**2/6-sum((p[i]*Rr[i]-CS/6)**2 for i in range(6)))==0)
print("  chi - C^2/24 == sum(p_iR_i - C/m)^2 + C^2(1/m-1/24) (m=6) ?",
      sp.simplify(chiS-CS**2/24-(sum((p[i]*Rr[i]-CS/6)**2 for i in range(6))+CS**2*(R(1,6)-R(1,24))))==0)

