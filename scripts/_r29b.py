import numpy as np, sympy as sp
from numpy.linalg import eigvalsh, eigh
np.set_printoptions(precision=9, suppress=True)
rho0=(9+sp.sqrt(5)+2*sp.sqrt(35))/18
zstar=sp.nsimplify("0.500440049768212308",[sp.sqrt(5),sp.sqrt(7),sp.sqrt(35)])
print("zstar nsimplify:",zstar)
# exact cone data at z (symbolic), keep z numeric-then-rationalize
import numpy as np
z=0.500440049768212308
g0=(1-3*z*z)/(1-z*z); gc=2/(5*(1-z*z))
print("g0=",g0," gc=",gc)
# obtuse graph of champion: rim 2-apart (pentagram=C5). weight w = gc*|<u_j,u_k>| for 2-apart
c4=np.cos(4*np.pi/5)
inner2=(1-z*z)*c4+z*z
w=gc*abs(inner2)
print("2-apart rim inner product =",inner2," w=gc*|.|=",w," 2w(=lam_max Q-)=",2*w)
# exact symbolic lam_max(Q-) = 2*gc*((1-z^2)*(-cos4pi/5) - z^2); -cos(4pi/5)=(1+sqrt5)/4
phi=(1+sp.sqrt(5))/2
neg_c4=(1+sp.sqrt(5))/4   # = -cos(4pi/5)=0.809016...
print("-cos(4pi/5) exact:",neg_c4, float(neg_c4))
# lam_max(Q-) symbolic in z:
zS=sp.Symbol('z',positive=True)
gcS=2/(5*(1-zS**2))
lamQm=2*gcS*(neg_c4*(1-zS**2)-zS**2)
lamQm=sp.simplify(lamQm)
print("lam_max(Q-)(z) =",lamQm)
print("  at z*:",float(lamQm.subs(zS,z)))
print("  rho0-1 =",float(rho0-1))
print("  Weyl overshoot (1+lamQm)-rho0 =",float(1+lamQm.subs(zS,z)-rho0))
# z*^2 minimal poly from report: 11225 w^4 -6700 w^3 +430 w^2 +260 w -31, w=z^2
w2=sp.RootOf(11225*sp.Symbol('w')**4-6700*sp.Symbol('w')**3+430*sp.Symbol('w')**2+260*sp.Symbol('w')-31, 0)
print("z^2 approx root:", [float(r) for r in sp.nroots(11225*sp.Symbol('w')**4-6700*sp.Symbol('w')**3+430*sp.Symbol('w')**2+260*sp.Symbol('w')-31)])