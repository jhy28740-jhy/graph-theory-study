# -*- coding: utf-8 -*-
import sympy as sp
s=sp.Symbol('s',positive=True); r5=sp.sqrt(5); phi=(1+r5)/2
R5=r5*(1-s)
Lw=sp.together(sp.expand(((1-3*s)*(1+s)+4*s*(1+R5)/5)/(1-s)**2))
s0=(1+2*r5)/19
print("L_w(s) =",Lw)
print("s0 = (1+2sqrt5)/19 =",float(s0),"  = (1+sqrt5)/(9+sqrt5)? 差",sp.simplify(s0-(1+r5)/(9+r5)))
val=sp.simplify(Lw.subs(s,s0)); print("\nmax L_w =",val,"  nsimplify:",sp.nsimplify(val),"  数值 %.15f"%float(val))
print("  (6+phi)/5 =",sp.simplify((6+phi)/5),"  差 =",sp.simplify(val-(6+phi)/5))
print("  (5+phi^2)/5 =",sp.simplify((5+phi**2)/5),"  差 =",sp.simplify(val-(5+phi**2)/5))
rho0=(9+r5+2*sp.sqrt(35))/18; LIM=sp.simplify(2*rho0-1)
print("\n2rho0-1 =",LIM,"= %.15f"%float(LIM))
mg=sp.simplify(LIM-val); print("余量 =",mg,"= %.15e"%float(mg))
print("Neg_max =",sp.simplify((val-1)/2)," phi^2/10 =",sp.simplify(phi**2/10)," 差",sp.simplify((val-1)/2-phi**2/10))
print("rho0-1 = %.15f ; Neg 余量 = %.15e"%(float(rho0-1),float(rho0-1-(val-1)/2)))
for n,Rn,lbl in ((3,1-3*s,"1-3s"),(4,sp.Integer(1),"1"),(6,3-4*s,"3-4s")):
    L=sp.together(sp.expand(((1-3*s)*(1+s)+4*s*(1+Rn)/n)/(1-s)**2))
    d=sp.diff(L,s); cr=[c for c in sp.solve(sp.numer(sp.together(d)),s) if c.is_real]
    cand=[sp.Rational(0),sp.Rational(1,3)]+[c for c in cr if 0<c<sp.Rational(1,3)]
    vs=[(sp.simplify(L.subs(s,c)),c) for c in cand]
    mx=max(vs,key=lambda t:float(t[0]))
    print("n=%d (R=%s): max L_w = %s = %.12f at s=%s"%(n,lbl,mx[0],float(mx[0]),mx[1]))
