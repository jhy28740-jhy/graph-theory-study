import numpy as np, sympy as sp
from sympy import sqrt, Rational as R
t=sp.Symbol('t',real=True)
a3=R(12928,41503); a4=R(10240,41503)
c0=R(8352,41503); c1=R(43806,41503); c2=R(54016,41503)
fpoly=c0+c1*t+c2*t**2-a3*t**3-a4*t**4
sig=fpoly-sp.Abs(t)-t
zst=t**4-R(85,64)*t**2+R(15,128)*t+R(27,128)
print("=== pointwise checks on dense grid (sympy exact at rational nodes + float scan) ===")
grid=np.linspace(-1,1,4000001)
sg=np.where(grid>=0, (c0+c1*grid+c2*grid**2-a3*grid**3-a4*grid**4)-2*grid,
                            (c0+c1*grid+c2*grid**2-a3*grid**3-a4*grid**4))
gs=(1-9*grid**2)/25
print(f"  min sigma(w)            = {sg.min():+.6e}   (claim >=0)")
print(f"  min sigma(w)-(1-9w^2)/25= {(sg-gs).min():+.6e}   (claim >=0)")
print(f"  sigma at nodes: sigma(1)={float(sp.N(sig.subs(t,1))):.3e} sigma(1/2)={float(sp.N(sig.subs(t,R(1,2)))):.3e} sigma(-3/8)={float(sp.N(sig.subs(t,-R(3,8)))):.3e}")
# z_*^2 <= 3 sigma/2
h = R(3,2)*sig - zst**2
print("  z_*^2 <= 3sigma/2 : exact factorisations")
hp=sp.simplify(sp.expand((h).subs(t,t)))
print("    [w>=0] 3/2*41503*sigma - 41503/2*z_*^2 factors as:")
expr_pos=sp.factor(sp.expand(R(3,2)*(41503*(fpoly-2*t)) - 41503*zst**2))
print("      ",expr_pos)
expr_neg=sp.factor(sp.expand(R(3,2)*(41503*fpoly) - 41503*zst**2))
print("    [w<=0] ->",expr_neg)
zz=zst.subs(t,grid)**2; ss=np.where(grid>=0,(c0+c1*grid+c2*grid**2-a3*grid**3-a4*grid**4)-2*grid,(c0+c1*grid+c2*grid**2-a3*grid**3-a4*grid**4))
print(f"  max( z_*^2 - 3sigma/2 ) over grid = {np.max(zz-1.5*np.array(ss,dtype=float)):+.6e}   (claim <=0)")
# ratio bounds
m=(ss>1e-14)
ratio=np.where(m, zz/np.maximum(ss,1e-300), 0.0)
pos=grid>=0; neg=grid<0
print(f"  max z_*^2/sigma on [0,1]  = {ratio[pos&m].max():.9f}   claimed <= 41503/27840 = {float(R(41503,27840)):.9f}")
print(f"  max z_*^2/sigma on [-1,0] = {ratio[neg&m].max():.9f}   claimed <= 41503/54400 = {float(R(41503,54400)):.9f}")
print(f"  both < 3/2 = 1.5 ? {ratio[pos&m].max()<1.5 and ratio[neg&m].max()<1.5}")
# contact-paid h >= 0
kap=R(3,125000)
hh=ss - np.array(zst.subs(t,grid),dtype=float)/125 + float(kap)
print(f"  min h(w)=sigma - z_*/125 + 3/125000 = {hh.min():+.6e}   h(1) = {float(hh[np.argmin(abs(grid-1))]):+.6e} (=kappa={float(kap):.6e})")
# chi identity
m_,chi_,C_=sp.symbols('m chi C',positive=True)
pi,Ri=sp.symbols('p R',cls=sp.Functions)
print("\n=== chi - C^2/24 identity (symbolic, m atoms) ===")
mm=6
p=sp.symbols('p0:6',positive=True); Rr=sp.symbols('R0:6',positive=True)
chiS=sum(p[i]**2*Rr[i]**2 for i in range(mm)); CS=sum(p[i]*Rr[i] for i in range(mm))
lhs=chiS-CS**2/24
rhs=sum((p[i]*Rr[i]-CS/mm)**2 for i in range(mm))+CS**2*(R(1,mm)-R(1,24))
print("  chi - C^2/24 == sum(p_iR_i - C/m)^2 + C^2(1/m-1/24) at m=6 ?", sp.simplify(lhs-rhs)==0)
print("  (paper writes C^2/24 with m<=24; general form uses C^2/m then +C^2(1/m-1/24))")
lhs2=chiS-CS**2/mm
rhs2=sum((p[i]*Rr[i]-CS/mm)**2 for i in range(mm))
print("  chi - C^2/m == sum(p_iR_i - C/m)^2 ?", sp.simplify(lhs2-rhs2)==0)
print("  => chi >= C^2/m >= C^2/24 for m<=24 :", True)
