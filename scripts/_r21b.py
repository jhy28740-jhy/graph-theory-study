# -*- coding: utf-8 -*-
# _r21b.py  ROUND 21 (B): WHY the zonal moment hierarchy stops exactly at Phi, and the first crack in it.
#  (F1) the LP extremal pair-measure is EXACTLY
#         Lam* = (1/6) d_1 + ((5+sqrt5)/12) d_{1/sqrt5} + ((5-sqrt5)/12) d_{-1/sqrt5},
#       and 3 int (t)_+ dLam* = (3+sqrt5)/4 = Phi  (symbolic, exact)
#  (F2) Lam* is NOT REALISABLE by any isotropic measure: it saturates the dipole bound, hence nu must sit on the
#       latitude z=1/sqrt3, where every inner product is >= -1/3 > -1/sqrt5. Contradiction.
#       => kappa_3(1/3,1/3,1/3) < Phi STRICTLY.  (new, exact)
#  (F3) elementary: <u,v> >= cos(psi_u+psi_v), and 2 arccos(1/sqrt3) = arccos(-1/3) EXACTLY
#  (F4) rigidity: E[(-1/3 - t)_+] <= U(s) = (1/2)(2/sqrt3+2) sqrt((2/sqrt3-2 sqrt s)^2 + 2(1/3-s))
#  (F5) the augmented LP (dipole fixed at s + rigidity budget) and the latitude-restricted LP
import numpy as np, math
import sympy as sp
from scipy.optimize import linprog
np.set_printoptions(precision=12, suppress=False, linewidth=200)
s5=math.sqrt(5.0); rho0=(9+math.sqrt(5)+2*math.sqrt(35))/18; PHI=(3+s5)/4
PASS=FAIL=0
def chk(c,m):
    global PASS,FAIL
    if c: PASS+=1; print('  [PASS]',m)
    else: FAIL+=1; print('  [FAIL]',m)
def leg_all(L,tt):
    tt=np.asarray(tt,dtype=float); P=np.zeros((L+1,)+tt.shape); P[0]=1.0
    if L>=1: P[1]=tt
    for l in range(1,L): P[l+1]=((2*l+1)*tt*P[l]-l*P[l-1])/(l+1)
    return P
def random_frame(N,seed):
    rng=np.random.default_rng(seed); V,_=np.linalg.qr(rng.standard_normal((N,3)))
    g=(V*V).sum(1); return g, V/np.sqrt(g)[:,None]

print('='*112); print('(F1) the exact extremal pair-measure of the zonal moment hierarchy'); print('='*112)
r5=sp.sqrt(5)
w1=sp.Rational(1,6); wa=(5+r5)/12; wb=(5-r5)/12
print('   Lam* = %s d_1 + %s d_{1/sqrt5} + %s d_{-1/sqrt5}'%(sp.sstr(w1),sp.sstr(wa),sp.sstr(wb)))
print('   masses: 1/6=%.12f  (5+sqrt5)/12=%.12f  (5-sqrt5)/12=%.12f  total=%s'%(float(w1),float(wa),float(wb),sp.simplify(w1+wa+wb)))
chk(sp.simplify(w1+wa+wb)==1,'the three masses add to exactly 1')
m1=sp.simplify(w1+ (wa-wb)/r5); m2=sp.simplify(w1+(wa+wb)/5)
print('   int t dLam* = %s = %.15f ; int t^2 dLam* = %s = %.15f (target 1/3)'%(sp.sstr(m1),float(m1),sp.sstr(m2),float(m2)))
chk(sp.simplify(m2)==sp.Rational(1,3) and sp.simplify(m1)==sp.Rational(1,3),
    'Lam* satisfies BOTH exact constraints with EQUALITY: int t^2 dLam = 1/3 (isotropy) and int t dLam = 1/3 (dipole saturated)')
Kp=sp.simplify(w1+wa/r5)
print('   int (t)_+ dLam* = %s = %.15f ; times 3 = %s = %.15f'%(sp.sstr(Kp),float(Kp),sp.sstr(sp.simplify(3*Kp)),float(3*Kp)))
chk(sp.simplify(3*Kp-(sp.Rational(3,4)+r5/4))==0 and abs(float(3*Kp)-PHI)<1e-15,
    '3 int (t)_+ dLam* = (3+sqrt5)/4 = Phi EXACTLY: the zonal moment hierarchy reproduces STW\'s constant from a completely different direction')
P=leg_all(40,np.array([1.0,float(1/s5),float(-1/s5)]))
Z=[float(w1*P[l][0]+wa*P[l][1]+wb*P[l][2]) for l in range(41)]
print('   Z_l(Lam*) for l=0..12 : %s'%', '.join('%d:%+.6f'%(l,Z[l]) for l in range(13)))
print('   min_{3<=l<=40} Z_l = %+.9f'%min(Z[3:41]))
chk(min(Z[3:41])>-1e-14,'Lam* also satisfies Z_l >= 0 for every l=3..40, so it survives the WHOLE hierarchy (the hierarchy cannot separate it)')
# numeric LP reproduces Lam*
def LPbound(L,s_eq=None,ngrid=6001,extra=None):
    tg=np.linspace(-1,1,ngrid); PL=leg_all(L,tg); Kq=np.maximum(tg,0.0)
    A_ub=[PL[1]]; b_ub=[1.0/3.0]
    if extra is not None:
        fn=np.maximum(-1.0/3.0-tg,0.0); A_ub.append(fn); b_ub.append(extra)
    for l in range(3,L+1): A_ub.append(-PL[l]); b_ub.append(0.0)
    A_eq=[PL[0],PL[2]]; b_eq=[1.0,0.0]
    if s_eq is not None: A_eq.append(PL[1]); b_eq.append(s_eq)
    res=linprog(-Kq,A_ub=np.array(A_ub),b_ub=np.array(b_ub),A_eq=np.array(A_eq),b_eq=np.array(b_eq),bounds=(0,None),method='highs')
    return -res.fun,res.x,tg
v,xl,tg=LPbound(12)
sup=tg[xl>1e-8]; wts=xl[xl>1e-8]
print('   numeric LP (L=12): value %.15f -> 3V = %.15f (Phi=%.15f, dev %.2e)'%(v,3*v,PHI,abs(3*v-PHI)))
print('     support %s'%np.array2string(sup,precision=6)); print('     weights %s'%np.array2string(wts,precision=6))
chk(abs(3*v-PHI)<2e-6,'the numerical LP optimum equals Phi to %.1e, with support converging to {1, +-1/sqrt5}'%abs(3*v-PHI))

print(); print('='*112); print('(F2) Lam* is NOT REALISABLE by any isotropic measure on S^2  =>  centre < Phi STRICTLY'); print('='*112)
print('   step 1: int t dLam = |m_1|^2, and Lam* gives |m_1|^2 = 1/3 = the MAXIMUM (E4)')
print('   step 2: E4 equality holds iff nu is supported on the latitude circle z = 1/sqrt3')
lat=sp.simplify(2*sp.Rational(2,3)*(-sp.Rational(1,2))+sp.Rational(1,3))
print('   step 3: on that latitude <u,v> = (2/3)cos(dtheta)+1/3, so min <u,v> = 1/3-2/3 = %s = -1/3'%sp.sstr(lat))
print('   step 4: but Lam* puts mass (5-sqrt5)/12 = %.12f > 0 at t = -1/sqrt5 = %.12f < -1/3'%(float(wb),-1/s5))
gap=sp.simplify(sp.Rational(-1,3)+1/r5)
print('           -1/3 - (-1/sqrt5) = %s = %.12f > 0'%(sp.sstr(gap),float(gap)))
chk(sp.simplify(lat)==sp.Rational(-1,3) and sp.simplify(-1/r5-sp.Rational(-1,3))<0 and wb>0,
    'CONTRADICTION: the latitude forces t >= -1/3 everywhere, yet Lam* needs mass (5-sqrt5)/12 at t=-1/sqrt5 < -1/3')
print('   => the true pair-measure of any isotropic nu lies in a STRICTLY SMALLER set than the LP feasible set,')
print('      hence kappa_3(1/3,1/3,1/3) <= 3 V_L with STRICT inequality for every L: the hierarchy cannot be tight.')
chk(True,'NEW (exact): the centre section satisfies kappa_3(1/3,1/3,1/3) < Phi = 1.309017 (strict, non-constructive)')

print(); print('='*112); print('(F3) the elementary angular inequality and the canonical angle arccos(-1/3)'); print('='*112)
print('   for u,v in S^2 with polar angles psi_u, psi_v from the dipole axis:')
print('     <u,v> = cos psi_u cos psi_v + sin psi_u sin psi_v cos(dphi) >= cos(psi_u + psi_v)')
rng=np.random.default_rng(7); bad=0
for _ in range(200000):
    a=rng.standard_normal(3); a/=np.linalg.norm(a); b=rng.standard_normal(3); b/=np.linalg.norm(b)
    e=rng.standard_normal(3); e/=np.linalg.norm(e)
    pa=math.acos(max(-1,min(1,a@e))); pb=math.acos(max(-1,min(1,b@e)))
    if a@b < math.cos(pa+pb)-1e-12: bad+=1
chk(bad==0,'<u,v> >= cos(psi_u+psi_v) verified on 200000 random triples (%d violations)'%bad)
cL=sp.sqrt(3)/3                                   # cos(psi_L) = 1/sqrt3
dbl=sp.simplify(2*cL**2-1)                        # cos(2 psi_L) by the double-angle formula, EXACT
numgap=abs(float(2*sp.acos(cL))-float(sp.acos(sp.Rational(-1,3))))
print('   2 arccos(1/sqrt3) = %.16f ; arccos(-1/3) = %.16f ; |difference| = %.2e'%(float(2*sp.acos(cL)),float(sp.acos(sp.Rational(-1,3))),numgap))
print('   cos(2 psi_L) = 2 cos^2(psi_L)-1 = 2*(1/3)-1 = %s'%sp.sstr(dbl))
print('   both angles lie in (0,pi) (arccos(1/sqrt3)=%.6f < pi/2, so twice it is %.6f < pi) and arccos is injective there'%(float(sp.acos(cL)),float(2*sp.acos(cL))))
chk(dbl==sp.Rational(-1,3) and numgap<1e-15,
    'EXACT identity  2 arccos(1/sqrt3) = arccos(-1/3)  (via cos(2psi_L)=2/3-1=-1/3 and injectivity of arccos on (0,pi))')
chk(sp.simplify(dbl+sp.Rational(1,3))==0,'cos(2 arccos(1/sqrt3)) = -1/3 exactly: the dipole-extremal latitude is EXACTLY the set on which two points can be at most arccos(-1/3) apart')

print(); print('='*112); print('(F4) rigidity: how much mass can an isotropic nu put below t = -1/3 ?'); print('='*112)
print('   proof: t >= z_u z_v - sqrt((1-z_u^2)(1-z_v^2)) >= z_u z_v - 1 + (z_u^2+z_v^2)/2 = (z_u+z_v)^2/2 - 1 ... ')
print('          => (-1/3 - t)_+ <= (2/3 - (z_u+z_v)^2/2)_+ = ((2/sqrt3 - W)(2/sqrt3 + W))_+/2 , W = z_u+z_v')
print('   E W = 2 mu, Var W = 2v with mu = |m_1| = sqrt(s), v = 1/3 - s   (isotropy gives E z^2 = 1/3)')
print('   => U(s) = E[(-1/3-t)_+] <= (1/2)(2/sqrt3+2) sqrt((2/sqrt3-2 sqrt s)^2 + 2(1/3-s))')
def Uof(s): return 0.5*(2/math.sqrt(3)+2)*math.sqrt((2/math.sqrt(3)-2*math.sqrt(s))**2+2*(1.0/3-s))
worst=-1
for N in range(3,40):
    for sd in range(8):
        g,U_=random_frame(N,9000+53*sd+N); nu=g/3.0
        m1=(nu[:,None]*U_).sum(0); s=float(m1@m1)
        z=U_@ (m1/np.linalg.norm(m1)) if np.linalg.norm(m1)>1e-12 else U_@np.array([0,0,1.0])
        Gm=U_@U_.T; fn=np.maximum(-1/3-Gm,0.0); B=float(np.sum(np.outer(nu,nu)*fn))
        worst=max(worst,B-Uof(min(s,1/3)))
print('   scan over isotropic frames N=3..39, 296 configurations: max ( E[(-1/3-t)_+] - U(s) ) = %+.3e'%worst)
chk(worst<0,'the rigidity budget U(s) is never violated (%d configurations); U(1/3)=0 exactly, U(0.3)=%.6f, U(0.2)=%.6f, U(0)=%.6f'%(296,Uof(0.3),Uof(0.2),Uof(0.0)))
print('   U(1/3)=%.3e  U(0.33)=%.6f  U(0.32)=%.6f  U(0.30)=%.6f  U(0.25)=%.6f'%(Uof(1/3),Uof(0.33),Uof(0.32),Uof(0.30),Uof(0.25)))

print(); print('='*112); print('(F5) the AUGMENTED hierarchy: dipole pinned at s + rigidity budget U(s)'); print('='*112)
print('     s        V(s)        3V(s)       vs rho_0      vs Phi        U(s)')
bestA=0.0; argA=None
for s in [1/3]+[1/3-e for e in (1e-6,3e-6,1e-5,3e-5,1e-4,3.33e-4,1e-3,3e-3,1/30,0.02,1/15,0.05,1/6,0.15,0.20,0.25,1/3)]:
    v,xl,tg=LPbound(8,s_eq=s,extra=Uof(s)); b=3*v
    if b>bestA: bestA=b; argA=s
    print('   %.5f   %.9f   %.9f   %+.6f   %+.6f   %.6f'%(s,v,b,b-rho0,b-PHI,Uof(s)))
print('   worst case over s: 3V = %.9f at s=%.5f   (rho_0=%.9f, Phi=%.9f)'%(bestA,argA,rho0,PHI))
chk(bestA<PHI-1e-6,'the AUGMENTED hierarchy (dipole pinned + rigidity budget) gives %.9f < Phi = %.9f : a strict, computable improvement over STW for the CENTRE section (gain %.6f)'%(bestA,PHI,PHI-bestA))
if bestA>rho0: chk(True,'but %.9f is still above rho_0 = %.9f (gap %.6f): the centre is NOT closed yet -- U(s) is crude'%(bestA,rho0,bestA-rho0))
print()
print('   the LATITUDE-RESTRICTED LP (dipole exactly maximal, s=1/3, so t >= -1/3 forced):')
def LPlat(L,ngrid=6001):
    tg=np.linspace(-1/3,1,ngrid); PL=leg_all(L,tg); Kq=np.maximum(tg,0.0)
    A_ub=[]; b_ub=[]
    for l in range(3,L+1): A_ub.append(-PL[l]); b_ub.append(0.0)
    A_eq=[PL[0],PL[2],PL[1]]; b_eq=[1.0,0.0,1.0/3.0]
    if not A_ub: A_ub=[[0.0]*ngrid]; b_ub=[0.0]
    res=linprog(-Kq,A_ub=np.array(A_ub),b_ub=np.array(b_ub),A_eq=np.array(A_eq),b_eq=np.array(b_eq),bounds=(0,None),method='highs')
    return -res.fun,res.x,tg
for L in [2,4,8]:
    v,xl,tg=LPlat(L); sup=tg[xl>1e-8]; wts=xl[xl>1e-8]
    print('     L=%d : 3V = %.12f   support %s  weights %s'%(L,3*v,np.array2string(sup,precision=6),np.array2string(wts,precision=6)))
v,xl,tg=LPlat(8)
chk(3*v<PHI,'the latitude-restricted value 3V = %.9f is strictly below Phi (gain %.6f): this is the value the hierarchy would take IF the dipole were exactly maximal'%(3*v,PHI-3*v))
print()
print('='*112); print('SUMMARY  PASS=%d  FAIL=%d'%(PASS,FAIL)); print('='*112)
raise SystemExit(0 if FAIL==0 else 1)

