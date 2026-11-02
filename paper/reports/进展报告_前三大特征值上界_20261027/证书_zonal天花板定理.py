# -*- coding: utf-8 -*-
import sympy as sp
s5=sp.sqrt(5); Phi=(3+s5)/4; rho0=(9+s5+2*sp.sqrt(35))/18
# dual witness measure on t in {-1/sqrt5, +1/sqrt5, +1}
t1=-1/s5; t2=1/s5; t3=sp.Integer(1)
lam1=(5-s5)/4; lam2=(5+s5)/4; lam3=sp.Rational(1,2)
print("lam =",(lam1,lam2,lam3)," all >=0:",all(sp.simplify(x)>0 for x in (lam1,lam2,lam3)))
S0=lam1+lam2+lam3; S1=lam1*t1+lam2*t2+lam3*t3; S2=lam1*t1**2+lam2*t2**2+lam3*t3**2
OBJ=lam1*sp.Max(t1,0)+lam2*sp.Max(t2,0)+lam3*sp.Max(t3,0)
print("sum lam      =",sp.simplify(S0),"(need <= 3)")
print("sum lam*t    =",sp.simplify(S1),"(need <= 1)")
print("sum lam*t^2  =",sp.simplify(S2),"(need <= 1)  -> kills the b-term exactly")
print("sum lam*(t)+ =",sp.simplify(lam2*t2+lam3*t3),"  == Phi ?",sp.simplify(lam2*t2+lam3*t3-Phi)==0)
print("Phi =",sp.nsimplify(Phi)," rho0 =",sp.N(rho0,16)," Phi-rho0 =",sp.N(Phi-rho0,12))
# Legendre P_l ; check sum lam P_l >= 0 for all l>=3
x=sp.Symbol("x")
def P(l): return sp.legendre(l,x)
bad=[]
for l in range(3,301):
    v=sp.nsimplify(lam1*sp.legendre(l,t1)+lam2*sp.legendre(l,t2)+lam3*sp.legendre(l,t3))
    v=sp.nsimplify(v,[s5])
    fv=float(v)
    if fv < -1e-12: bad.append((l,fv))
print("sum lam*P_l < 0 for l in 3..300 :",bad if bad else "NONE  (all >= 0)")
for l in [3,4,5,6,7,8,9,10,20,50,100,200]:
    v=sp.nsimplify(lam1*sp.legendre(l,t1)+lam2*sp.legendre(l,t2)+lam3*sp.legendre(l,t3),[s5])
    print("   l=%-4d sum lam P_l = %-28s = %+.10f"%(l,str(v),float(v)))
print()
print("=== THE NO-GO, 3 lines ===")
print("Multiply the majorant  (t)_+ <= a + beta*t + b*(t^2-1/3) - sum_{l>=3} rho_l P_l(t)  by lam and sum:")
print("  Phi = sum lam (t)_+  <=  a*3 + beta*1 + b*0 - sum rho_l*(>=0)  =  3a + beta")
print("  and  Lambda_3(G) <= [max(3a,b)+beta]*n - 3  with  max(3a,b)+beta >= 3a+beta >= Phi.")
print("  => the whole zonal-majorant family (STW's method, generalized to (t)_+ and arbitrary")
print("     Legendre corrections) has ceiling EXACTLY Phi. It cannot prove c_3 <= rho_0.")
