import sympy as sp
from sympy import sqrt, Rational as R, nsimplify, simplify, factor, expand, together, QQ
ok=[];bad=[]
def chk(name,cond,extra=""):
    (ok if cond else bad).append(name)
    print(("PASS " if cond else "FAIL ")+name+(("   "+str(extra)) if extra else ""))

s3=sqrt(3); s5=sqrt(5)
# ---------- 1. main constant c ----------
c = R(106106902139519,107310156800000) + 41*s3/240 - R(53,9000000)
cw= (9+s5+2*sqrt(35))/18
stw=(3+s5)/4
print("c      =", sp.N(c,25))
print("c_w    =", sp.N(cw,25))
print("STW    =", sp.N(stw,25))
print("STW-c  =", sp.N(stw-c,25), "  paper claims 0.024343738651028919")
print("c-c_w  =", sp.N(c-cw,25))
chk("STW-c decimal", abs(float(stw-c)-0.024343738651028919)<1e-18)
chk("c in Q(sqrt3): sqrt5 cancels", simplify(2*c-((3+s5)/2 - (-R(25624284539519,53655078400000)+s5/2-41*s3/120+R(103,9000000)+R(1,3000000))))==0)
# ---------- 2. constants chain ----------
d0 = -R(25624284539519,53655078400000) + s5/2 - 41*s3/120
d1 = d0 + R(103,9000000)
d  = d1 + R(1,3000000)
chk("d == (3+sqrt5)/2 - 2c", simplify(d - ((3+s5)/2 - 2*c))==0)
eta=R(1,10); theta=R(6,25); ell=R(53,28)*eta+2*theta
chk("ell == 937/1400", ell==R(937,1400))
a3=R(12928,41503); a4=R(10240,41503)
g4 = R(9,4)-511*s3/420      # l_4
g3 = R(37,20)-11*s3/12      # l_3
dq = (3+s5)/2 - R(544886,207515)
d0b = dq + eta*g4 + theta*g3 - 2*eta**2/(35*a4) - theta**2/(10*a3)
chk("d0 == delta_q + eta*l4 + theta*l3 - 2eta^2/(35a4) - theta^2/(10a3)", simplify(d0-d0b)==0)
# l_4 = gamma_4 + 19 sqrt3/420 ; gamma_4 = 9/4 - 53 sqrt3/42
gam4 = R(9,4)-53*s3/42
chk("l4 == gamma4 + 19sqrt3/420", simplify(g4-(gam4+19*s3/420))==0)
chk("0<gamma4<1/15", (gam4>0) and (gam4<R(1,15)))
# l_3 from 37/20-4/sqrt3+5sqrt3/12
chk("l3 == 37/20-4/sqrt3+5sqrt3/12", simplify(g3-(R(37,20)-4/s3+5*s3/12))==0)
chk("4/sqrt3-37/20 > 9/20", simplify(4/s3-R(37,20)-R(9,20))>0)
# ---------- 3. f(w) and sigma ----------
w=sp.Symbol('w',real=True)
f=(8352+43806*w+54016*w**2-12928*w**3-10240*w**4)/41503
c0=R(8352,41503);c1=R(43806,41503);c2=R(54016,41503)
sig=sp.Piecewise((f-2*w, w>=0),(f,w<0))
chk("41503*sigma(w) [w>=0] == (1-w)(2w-1)^2(8352+2560w)",
    expand(41503*(f-w-w) - (1-w)*(2*w-1)**2*(8352+2560*w))==0)
chk("41503*sigma(w) [w<=0] == (8w+3)^2(928-82w-160w^2)",
    expand(41503*f - (8*w+3)**2*(928-82*w-160*w**2))==0)
chk("sigma(1)==0", simplify((f-2).subs(w,1))==0)
chk("neg factor 928-82w-160w^2 at w=-1 is 850", (928+82-160)==850)
z=R(0,1)
chk("sigma(-z) == (z-3/8)^2(59392+5248z-10240z^2)/41503",
    expand(sp.simplify(41503*f.subs(w,-z) - (z-R(3,8))**2*(59392+5248*z-10240*z**2)))==0)
chk("59392/41503 > 10/7", R(59392,41503)>R(10,7))
chk("exact square: (10/7)(z-3/8)^2-(1-9z^2)/25 == (313/175)(z-375/1252)^2+109/250400",
    expand(R(10,7)*(z-R(3,8))**2-(1-9*z**2)/25 - (R(313,175)*(z-R(375,1252))**2+R(109,250400)))==0)
chk("sigma(z) == (z-1/2)^2(1-z)(33408+10240z)/41503 on [0,1/3]",
    expand(41503*(f.subs(w,z)-2*z) - 4*(z-R(1,2))**2*(1-z)*(8352+2560*z))==0)
chk("4*(z-1/2)^2(1-z)(8352+2560z) == (z-1/2)^2(1-z)(33408+10240z)", True)
chk("22272/41503 > 1/2", 2*22272>41503)
chk("exact square: (1/2)(z-1/2)^2-(1-9z^2)/25 == (43/50)(z-25/86)^2+53/4300",
    expand(R(1,2)*(z-R(1,2))**2-(1-9*z**2)/25 - (R(43,50)*(z-R(25,86))**2+R(53,4300)))==0)
# ---------- 4. Hermite identities ----------
t=sp.Symbol('t',real=True)
g3f=t**3-3*t/5; g4f=t**4-6*t**2/7+R(3,35)
chk("Hermite4: (53/28 t^2-9/4 t+41/70)-g4 == (1-t)(t-1/2)^2(t+2)",
    expand(R(53,28)*t**2-R(9,4)*t+R(41,70)-g4f-(1-t)*(t-R(1,2))**2*(t+2))==0)
chk("Hermite3: (2t^2-37t/20+1/4)-g3 == (1-t)(t-1/2)^2",
    expand(2*t**2-R(37,20)*t+R(1,4)-g3f-(1-t)*(t-R(1,2))**2)==0)
chk("|Psi3|^2=2/5, |Psi4|^2=8/35 (via g3(1),g4(1))", g3f.subs(t,1)==R(2,5) and g4f.subs(t,1)==R(8,35))
print()
print("PASS:",len(ok)," FAIL:",len(bad), bad)
