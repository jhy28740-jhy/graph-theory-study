import sympy as sp
from sympy import sqrt, Rational as R, factorint, nsimplify
P = R(25624284539519,53655078400000)
d0 = -P + sqrt(5)/2 - 41*sqrt(3)/120
d1 = d0 + R(103,9000000)
d  = d1 + R(1,3000000)
print("denominators:", factorint(53655078400000), "|", factorint(804826176000000), "|", factorint(482895705600000))
for name,val,claim,s5u,s3l in [
  ("1/20-d0", R(1,20)-d0, R(620701786231,482895705600000), R(161,72), R(433,250)),
  ("1/20-d1", R(1,20)-d1, R(1025292188593,804826176000000), R(161,72), R(433,250)),
  ("1/20-d",  R(1,20)-d,  R(1025023913201,804826176000000), R(161,72), R(433,250)),
  ("d0-7/144",d0-R(7,144),R(940735364069,29456638041600000),None,None),
]:
    print("\n==",name,"==")
    print("  true value      =", sp.N(val,22))
    print("  paper's rational=", sp.N(claim,22), " (claim: true > rational)")
    print("  claim correct?  ", bool(val > claim), "   margin =", sp.N(val-claim,10))
    if s5u is not None:
        LB = R(1,20) + P - s5u/2 + 41*s3l/120
        if name!="1/20-d0":
            LB = LB - (d1-d0) - ((d-d1) if name=="1/20-d" else 0)
        print("  rational LB from sqrt5<161/72, sqrt3>433/250 =", sp.N(LB,22))
        print("  LB > claim ? ", bool(LB>claim), "  LB-claim =", sp.N(LB-claim,10))
        # what sqrt5 upper bound would be needed?
        need = (LB-claim)*2
        print("  slack available in sqrt5 (needs 161/72 - sqrt5 >=", sp.N(need,10), "; actual =", sp.N(s5u-sqrt(5),10),")")
        # tighter bounds
        for a,b in [(R(161,72),R(433,250)),(R(337,151),R(1351,780)),(R(377,169),R(181,105)),(R(682,305),R(1351,780))]:
            LB2 = R(1,20)+P-a/2+41*b/120
            if name=="1/20-d1": LB2 = LB2-R(103,9000000)
            if name=="1/20-d":  LB2 = LB2-R(103,9000000)-R(1,3000000)
            print("    with sqrt5<%s, sqrt3>%s -> LB=%s  > claim? %s"%(a,b,sp.N(LB2,18),bool(LB2>claim)))
# ---- fix the p_* test ----
x=sp.Symbol('x')
mid = R(15,16)+x/16-R(11,8)*x**2
print("\np_* completion: 15/16+y/16-11y^2/8 == 1321/1408-(11/8)(y-1/44)^2 ?",
      sp.expand(mid-(R(1321,1408)-R(11,8)*(x-R(1,44))**2))==0)
print("-15y^2/8 - y^3 <= -11y^2/8 on [-1/2,1/2] ?  <=> -y^3 <= y^2/2, check endpoints/critical:",
      all((-y**3 - y**2/2)<=0 for y in [R(-1,2),R(-1,4),0,R(1,4),R(1,2)]))
c = R(106106902139519,107310156800000)+41*sqrt(3)/240-R(53,9000000)
print("\nc == (3+sqrt5)/4 - d/2 ?", sp.simplify((3+sqrt(5))/4 - d/2 - c)==0)
print("c-w shift: c == (3+sqrt5)/4 - d1/2 - 1/6000000 ?", sp.simplify((3+sqrt(5))/4 - d1/2 - R(1,6000000) - c)==0)
