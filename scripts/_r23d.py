# -*- coding: utf-8 -*-
# _r23d.py  Round 23 part E: the equi-latitude circle (A*^2 = 3) solved
import numpy as np, math, itertools, sys
import sympy as sp
import mpmath as mp
from scipy.optimize import minimize
mp.mp.dps=50
PASS=FAIL=0
def chk(n,c,e=""):
    global PASS,FAIL
    if c: PASS+=1; print("   [PASS] "+n+("   "+e if e else ""))
    else:    FAIL+=1; print("   [FAIL] "+n+("   "+e if e else ""))
def hd(t): print("\n"+"="*112); print(t); print("="*112)
S5=mp.sqrt(5); rho0=(mp.mpf(9)+S5+2*mp.sqrt(35))/18
c1=(20+12*S5)/225; c0=6*rho0-(51*S5-21)/16

hd("(E0) exact constants of the coupled inequality (corrected threshold)")
print("   c1 = 1/(-9beta) = (20+12sqrt5)/225 =",mp.nstr(c1,18))
print("   c0 = 6rho_0-(51sqrt5-21)/16        =",mp.nstr(c0,18))
print("   threshold: A*^2 <= c0 + 1/(5 c1) =",mp.nstr(c0+1/(5*c1),18),"  -> the floor m4>=1/5 alone closes it")
print("   open window: A*^2 in (c0+1/(5c1), 3]  i.e. V=6-2A*^2 in [0, %.10f)"%float(6-2*(c0+1/(5*c1))))
print("   at A*^2=3: need m4 >=",mp.nstr(c1*(3-c0),18),";  exact equi-latitude floor 19/81 =",mp.nstr(mp.mpf(19)/81,18),
      ";  margin =",mp.nstr(mp.mpf(19)/81-c1*(3-c0),10))

hd("(E1) EXACT lemma: on the equi-latitude circle  m4 = 19/81 + (8|nu3|^2+2|nu4|^2)/729")
t0=1/math.sqrt(3); r0=math.sqrt(1-t0*t0)
def lat(angles,masses):
    U=np.stack([r0*np.cos(angles),r0*np.sin(angles),np.full(len(angles),t0)],1)
    return masses,U
def nhat(angles,masses,k): return complex(np.sum(masses*np.exp(-1j*k*angles)))
rng=np.random.default_rng(3); worst=0.0; worstL=0.0
def cons_mat(angles):
    n=len(angles)
    rows=[np.ones(n),np.cos(angles),np.sin(angles),np.cos(2*angles),np.sin(2*angles)]
    return np.array(rows)          # 5 x n ;  m must satisfy Cm = (3,0,0,0,0)
for trial in range(400):
    n=int(rng.integers(6,13)); ang=np.sort(rng.uniform(0,2*math.pi,n))
    C=cons_mat(ang); b=np.array([3.0,0,0,0,0])
    # affine solution set
    sol,res_,rk,sv=np.linalg.lstsq(C,b,rcond=None)
    if rk<5: continue
    Ns=np.linalg.svd(C)[2][rk:]                 # null space basis (n-rk) x n
    if Ns.shape[0]==0: continue
    for _ in range(6):
        coef=rng.normal(size=Ns.shape[0])*0.3
        m=sol+coef@Ns
        if m.min()<=1e-6: continue
        if abs(C@m-b).max()>1e-9: continue
        U=lat(ang,m)[1]; Cm=U@U.T
        m4=float(np.einsum('i,j,ij->',m,m,Cm**4))/9.0
        n3=nhat(ang,m,3); n4=nhat(ang,m,4)
        pred=19/81+(8*abs(n3)**2+2*abs(n4)**2)/729
        worst=max(worst,abs(m4-pred))
        L=float(np.einsum('i,j,ij->',m,m,np.abs(Cm)))
        worstL=max(worstL,L)
        break
chk("m4 = 19/81 + (8|nu_3|^2 + 2|nu_4|^2)/729 on random equi-latitude measures",worst<1e-11,"max dev %.2e"%worst)
print("   => m4 >= 19/81 EXACTLY on the equi-latitude circle, equality iff nu_3=nu_4=0")
chk("19/81 exceeds the coupled requirement at A*^2=3", mp.mpf(19)/81 > c1*(3-c0),
    "margin %.6e"%float(mp.mpf(19)/81-c1*(3-c0)))

hd("(E2) vertex reduction: L(m)=m^T A m is convex in m => the max is at <=5 atoms")
def Lmax_k(k,nstart=400,seed=0):
    r2=np.random.default_rng(seed); best=(-1,None,None)
    starts=[]
    starts.append(np.arange(k)*2*math.pi/k)                      # regular k-gon
    for _ in range(nstart): starts.append(np.sort(r2.uniform(0,2*math.pi,k)))
    b=np.array([3.0,0,0,0,0])
    for th0 in starts:
        def negL(th):
            C=cons_mat(th)
            if k<5:
                m,res_,rk,sv=np.linalg.lstsq(C,b,rcond=None)
                if rk<k or np.abs(C@m-b).max()>1e-8: return 1e6
            else:
                try: m=np.linalg.solve(C,b)
                except np.linalg.LinAlgError: return 1e6
            if m.min()<=0: return 1e6
            U=lat(th,m)[1]; return -float(np.einsum('i,j,ij->',m,m,np.abs(U@U.T)))
        if k<5:
            rr=minimize(negL,th0,method='Nelder-Mead',
                        options={'maxiter':8000,'maxfev':8000,'xatol':1e-12,'fatol':1e-14})
        else:
            rr=minimize(negL,th0,method='Nelder-Mead',
                        options={'maxiter':20000,'maxfev':20000,'xatol':1e-12,'fatol':1e-14})
        if -rr.fun>best[0]:
            th=rr.x; C=cons_mat(th)
            m=np.linalg.lstsq(C,b,rcond=None)[0] if k<5 else np.linalg.solve(C,b)
            best=(-rr.fun,th.copy(),m.copy())
    return best
tot=(-1,None,None,None)
for k in (3,4,5):
    L,th,m=Lmax_k(k,nstart=(300 if k<5 else 500),seed=17*k)
    U=lat(th,m)[1]; Cm=U@U.T
    print("   k=%d : L_max = %.12f   masses=%s   angles/2pi=%s"%(k,L,np.array2string(m,precision=6),
          np.array2string((th/(2*math.pi))%(1.0),precision=6)))
    print("          m4=%.12f  |nu3|=%.3e |nu4|=%.3e  isotropy=%.1e"%(
        float(np.einsum('i,j,ij->',m,m,Cm**4))/9.0,abs(nhat(th,m,3)),abs(nhat(th,m,4)),
        abs((m[:,None,None]*(U[:,:,None]*U[:,None,:])).sum(0)-np.eye(3)).max()))
    if L>tot[0]: tot=(L,k,th,m)
print("   BEST equi-latitude L = %.12f  (k=%d)"%(tot[0],tot[1]))
print("   regular square (4 atoms, mass 3/4, angles 0,90,180,270):")
thsq=np.array([0,math.pi/2,math.pi,3*math.pi/2]); msq=np.full(4,0.75)
Usq=lat(thsq,msq)[1]; Csq=Usq@Usq.T
Lsq=float(np.einsum('i,j,ij->',msq,msq,np.abs(Csq)))
m4sq=float(np.einsum('i,j,ij->',msq,msq,Csq**4))/9.0
print("      L = %.15f  (=9/2 ? %.2e)   m4 = %.15f  (=7/27 ? %.2e)"%(Lsq,abs(Lsq-4.5),m4sq,abs(m4sq-7/27)))
print("      isotropy defect = %.2e"%abs((msq[:,None,None]*(Usq[:,:,None]*Usq[:,None,:])).sum(0)-np.eye(3)).max())
chk("the regular square attains L = 9/2 exactly on the equi-latitude circle",abs(Lsq-4.5)<1e-13,"L=%.15f"%Lsq)
chk("equi-latitude L_max <= 6rho_0-3 = %.12f (margin %.6f)"%(float(6*rho0-3),float(6*rho0-3)-max(tot[0],Lsq)),
    max(tot[0],Lsq)<float(6*rho0)-3,"margin %.6f"%(float(6*rho0)-3-max(tot[0],Lsq)))
chk("equi-latitude L_max <= 9/2 (numerical, k<=5 vertex search)",max(tot[0],Lsq)<4.5+1e-9,
    "L_max=%.12f"%max(tot[0],Lsq))
print("\n   CONSEQUENCE: at A*^2=3,  J = L + A*^2 <= 9/2+3 = 15/2 = 7.5 < 6rho_0 = %.9f (margin %.9f)"%(
      float(6*rho0),float(6*rho0)-7.5))
sys.stdout.flush(); print("\nPART E SUMMARY PASS=%d FAIL=%d"%(PASS,FAIL))
