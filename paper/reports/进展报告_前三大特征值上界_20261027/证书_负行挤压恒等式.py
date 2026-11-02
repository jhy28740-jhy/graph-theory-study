# -*- coding: utf-8 -*-
# 证书：定理 squeeze（负行挤压与钉扎）在冠军处的精确饱和核验
# 运行：C:\Python314\python.exe -X utf8 -u 证书_负行挤压恒等式.py
import sympy as sp
s5,s7,s35 = sp.sqrt(5),sp.sqrt(7),sp.sqrt(35)
a=(91-25*s7)/126; b=(7+5*s7)/126                      # 冠军吹起权 w*
Vhh=14/(-7*s5-7+5*s7+3*s35); Vhc=sp.sqrt(sp.Rational(7,5))
Vcc=-sp.Rational(7,2)-7*s5/10+3*s35/10+5*s7/2
Vce=(4-s7)*sp.sqrt(sp.Rational(7,5)); Vc2=-sp.sqrt(sp.Rational(7,5))
PH=a*Vhh; PC=b*Vcc; H=sp.sqrt(a*b)*Vhc; A=b*Vce; B=b*Vc2
RHO=(9+s5+2*s35)/18
G=sp.zeros(6,6)
for i in range(6):
    for j in range(6):
        if i==j: G[i,j]=PH if i==0 else PC
        elif i==0 or j==0: G[i,j]=H
        else:
            d=min((i-j)%5,(j-i)%5); G[i,j]=A if d==1 else B
E=G*G-G
print("G^2=G residuals:", [float(E[p,q]) for p,q in [(0,0),(0,1),(1,1),(1,2),(1,3)]])
print("trace G =", float(sp.trace(G)))
yh=sp.Integer(1); yc=sp.simplify((RHO-PH)/(5*H))      # Perron 向量 (D5 对称, 归一 y_hub=1)
print("Perron hub residual:", float(sp.simplify(PH*yh+5*H*yc-RHO*yh)))
print("Perron cyc residual:", float(sp.simplify(H*yh+(PC+2*A)*yc-RHO*yc)))
Y2=yh**2+5*yc**2
def row(k):                                            # k=0 hub, k=1 cycle
    if k==0:
        P=5*H*yc; N=sp.Integer(0); Em=sp.Integer(0); Sk=sp.Integer(0); pk=PH; yk=yh
    else:
        P=H*yh+2*A*yc; N=-2*B*yc; Em=2*B**2; Sk=2*yc**2; pk=PC; yk=yc
    i1=sp.simplify(P-(RHO-pk)*yk)
    i2=sp.simplify((pk-pk**2-Em)*(Y2-yk**2-Sk)-((RHO-pk)*yk)**2)
    i3=sp.simplify(N**2-Em*Sk)
    return i1,i2,i3
for nm,k in [("hub",0),("cycle",1)]:
    i1,i2,i3=row(k)
    print(nm,"(i)",i1," (ii) slack",sp.nsimplify(i2,[s5,s7,s35]),float(i2)," (iii)",i3)
yGy=yh*(PH*yh+5*H*yc)+5*yc*(H*yh+PC*yc+2*A*yc+2*B*yc)
totN=5*yc*(-2*B*yc)
print("(iv) pin residual:", sp.simplify(RHO*Y2-totN-yGy))
print("Y2 =",sp.nsimplify(sp.simplify(Y2),[s5,s7,s35]),float(Y2))
print("yGy=",float(yGy)," totN=",sp.nsimplify(sp.simplify(totN),[s5,s7,s35]),float(totN))
