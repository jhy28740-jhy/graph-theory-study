# -*- coding: utf-8 -*-
# _r27c.py  第 27 轮 C 段：常数 S 的闭式（符号）+ sup P 的数值探索（检验 10/(3pi) 的紧性）
import numpy as np, sympy as sp
from scipy.optimize import minimize
rng=np.random.default_rng(2026)
pi=sp.pi

print("="*100); print("C1  S = sum_{m>=1} 1/((4m+1)((4m+1)^2-4)) 的闭式（sympy 精确）"); print("="*100)
m=sp.symbols('m',positive=True,integer=True)
term=1/((4*m+1)*((4*m+1)**2-4))
Ssym=sp.summation(term,(m,1,sp.oo))
print("  sympy summation  =",sp.simplify(Ssym))
guess=sp.Rational(5,24)-pi/16
print("  猜测 5/24 - pi/16 =",guess)
print("  差（精确化简）    =",sp.simplify(sp.nsimplify(Ssym)-guess) if Ssym.is_number else "n/a")
print("  数值: S=%.20f  5/24-pi/16=%.20f  差=%.3e"
      %(float(Ssym.evalf(40)),float(guess.evalf(40)),abs(float(Ssym.evalf(40))-float(guess.evalf(40)))))

# 独立推导：psi(1/4)-psi(3/4) = -pi
print("  独立推导: S = (1/16)(psi(1/4)-psi(3/4)) + 5/24")
Spsi=sp.Rational(1,16)*(sp.polygamma(0,sp.Rational(1,4))-sp.polygamma(0,sp.Rational(3,4)))+sp.Rational(5,24)
print("    sympy 化简 =",sp.simplify(Spsi),"  数值差 =",abs(float(Spsi.evalf(40))-float(Ssym.evalf(40))))

print()
print("="*100); print("C2  P 的上界塌缩为 10/(3pi)（精确代数）"); print("="*100)
Pb=sp.simplify(1+16*Ssym/pi)
print("  1 + 16S/pi =",Pb, " = ",sp.nsimplify(Pb))
print("  10/(3pi)   =",sp.simplify(10/(3*pi)))
print("  差         =",sp.simplify(Pb-10/(3*pi)))
print("  数值 P上界 = %.20f   10/(3pi) = %.20f"%(float(Pb.evalf(40)),float((10/(3*pi)).evalf(40))))

sqrt5=sp.sqrt(5); phi=(1+sqrt5)/2
rho0=(9+sqrt5+2*sp.sqrt(35))/18
Negb=sp.Rational(5,1)/(6*pi)
print()
print("  Neg <= 5/(6pi) = %.20f"%float(Negb.evalf(40)))
print("  kappa_3(秩一) <= 1 + 5/(6pi) = %.20f"%float((1+Negb).evalf(40)))
print("  rho_0 = %.20f"%float(rho0.evalf(40)))
gap=sp.simplify(rho0-1-Negb)
print("  余量 rho_0-1-5/(6pi) = %s = %.20f"%(gap,float(gap.evalf(40))))
print("  余量>0 的严格证明: 5/(6pi) < 5/(6*3.14) = %.10f < rho_0-1 = %.10f -> %s"
      %(float(sp.Rational(5,1)/(6*sp.Rational(314,100))),float((rho0-1).evalf(30)),
        bool(sp.Rational(5,1)/(6*sp.Rational(314,100)) < (rho0-1))))
print("  等价地只需 pi > 5/(6(rho_0-1)) = %.10f"%float((5/(6*(rho0-1))).evalf(30)))
print("  L_w <= 1+5/(3pi) = %.20f ;  锥族上确界 (13+sqrt5)/10 = %.20f ; 差 %.6e"
      %(float((1+2*Negb).evalf(40)),float(((13+sqrt5)/10).evalf(40)),
        float((1+2*Negb).evalf(40))-float(((13+sqrt5)/10).evalf(40))))
print("  STW Phi = %.20f ; 本界改进 = %.6e"%(float(phi.evalf(40))*0+float(((3+sqrt5)/4).evalf(40)),
      float(((3+sqrt5)/4).evalf(40))-float((1+Negb).evalf(40))))
print("  冠军 P = (3+sqrt5)/5 = %.20f ; 上界 10/(3pi) = %.20f ; 松弛比 %.6f"
      %(float(((3+sqrt5)/5).evalf(40)),float(Pb.evalf(40)),float(Pb.evalf(40))/float(((3+sqrt5)/5).evalf(40))))

print()
print("="*100); print("C3  sup P 的数值探索（P 只依赖 (rho,phi)，约束 sum rho=2, sum rho e^{2i phi}=0）"); print("="*100)
def Pval(rho,ang):
    c=np.cos(ang[:,None]-ang[None,:])
    return float(np.sum((rho[:,None]*rho[None,:])*np.where(c<0,c**2,0.0)))
def negpen(x,N,mu):
    rho=x[:N]; ang=x[N:]
    if np.any(rho<0): return 1e6+abs(rho.min())*1e6
    v=(rho.sum()-2.0)**2+abs(np.sum(rho*np.exp(2j*ang)))**2
    return -Pval(rho,ang)+mu*v
res={}
for N in [3,4,5,6,7,8,9,10,12,14,16,20,24]:
    b=-1;bx=None
    nstart=260 if N<=10 else 120
    for t in range(nstart):
        x0=np.concatenate([2*rng.random(N)/N, rng.random(N)*2*np.pi])
        r=minimize(negpen,x0,args=(N,2e6),method='Nelder-Mead',
                   options={'maxiter':20000,'maxfev':20000,'xatol':1e-11,'fatol':1e-14})
        x=r.x;rho=x[:N];ang=x[N:]
        if np.any(rho<0): continue
        viol=abs(rho.sum()-2)+abs(np.sum(rho*np.exp(2j*ang)))
        if viol<1e-7:
            val=Pval(rho,ang)
            if val>b: b=val;bx=(rho.copy(),ang.copy(),viol)
    res[N]=b
    print("  N=%2d  sup P >= %.12f  (违例 %.1e)"%(N,b,bx[2] if bx else -1))
best=max(res.values()); Nb=max(res,key=lambda k:res[k])
print("  最大值 = %.12f 在 N=%d"%(best,Nb))
print("  冠军环 P = (3+sqrt5)/5 = %.12f"%float(((3+sqrt5)/5).evalf(40)))
print("  上界 10/(3pi) = %.12f"%float(Pb.evalf(40)))
print("  上界 - 数值sup = %.6e  (%.3f%%)"%(float(Pb.evalf(40))-best,100*(float(Pb.evalf(40))-best)/best))
print("  上界 - 冠军    = %.6e  (%.3f%%)"%(float(Pb.evalf(40))-float(((3+sqrt5)/5).evalf(40)),
      100*(float(Pb.evalf(40))-float(((3+sqrt5)/5).evalf(40)))/float(((3+sqrt5)/5).evalf(40))))
print("ALL C-SECTION CHECKS DONE")
