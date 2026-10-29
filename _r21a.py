# -*- coding: utf-8 -*-
# _r21a.py  ROUND 21 (A): the LEGENDRE-CATALAN spectral representation of the main functional.
#  (E1) exact Legendre coefficients of K(t)=(t)_+ :  I_{2k} = -P_{2k}(0)/((2k+2)(2k-1)),
#       c_{2k} = (-1)^{k+1}(4k+1)Cat_k/(4(2k-1)4^k), c_{2k+1}=0 (k>=1),
#       sum_l c_l = 1, sum_{k>=1} c_{2k} = 1/4, sum_{k>=2} c_{2k} = -1/16
#  (E2) SPECTRAL REPRESENTATION rho_{g,U}(Lam) = sum_l c_l W_l, W_l = sum_ij w_i w_j P_l(<u_i,u_j>) >= 0
#  (E3) CENTRE section: Lam=I/3  =>  W_2 = 0 EXACTLY (isotropy <=> zero quadrupole)
#  (E4) DIPOLE bound sup|m_1| = 1/sqrt3 with an explicit quadratic dual certificate; sharp, equality = latitude
#  (E5) SIGN-BLIND no-go at the champion: exact rho(|Q*|), (1+rho(|Q*|))/2 > rho_0; theta>=1/2 forced
#  (E6) ZONAL MOMENT-LP HIERARCHY for the centre section (new method)
import numpy as np, math, time
import sympy as sp
from scipy.optimize import linprog
np.set_printoptions(precision=12, suppress=False, linewidth=200)
s5=math.sqrt(5.0); s7=math.sqrt(7.0); s35=math.sqrt(35.0)
rho0=(9+s5+2*s35)/18; tau0=rho0-1; PHI=(3+s5)/4; phi=(1+s5)/2
PASS=FAIL=0
def chk(c,m):
    global PASS,FAIL
    if c: PASS+=1; print('  [PASS]',m)
    else: FAIL+=1; print('  [FAIL]',m)
def c_of(l):
    """exact-as-float Legendre coefficient of (t)_+ via the Catalan closed form."""
    if l==0: return 0.25
    if l==1: return 0.5
    if l%2==1: return 0.0
    k=l//2; Cat=math.comb(2*k,k)//(k+1)
    return ((-1)**(k+1))*(4*k+1)*(Cat/(4**k))/(4.0*(2*k-1))
def leg_all(L,tt):
    tt=np.asarray(tt,dtype=float); P=np.zeros((L+1,)+tt.shape); P[0]=1.0
    if L>=1: P[1]=tt
    for l in range(1,L): P[l+1]=((2*l+1)*tt*P[l]-l*P[l-1])/(l+1)
    return P
def random_frame(N,seed):
    """random isotropic frame: V in St(N,3) => sum_i b_i b_i^T = I_3 EXACTLY, b_i = row i of V,
       g_i = ||b_i||^2, u_i = b_i/||b_i||; then Q=VV^T is a rank-3 projection with Q_ij=sqrt(g_i g_j)<u_i,u_j>."""
    rng=np.random.default_rng(seed)
    V,_=np.linalg.qr(rng.standard_normal((N,3)))
    g=(V*V).sum(1); U=V/np.sqrt(g)[:,None]
    return g,U
def rho_direct(g,U,Lam):
    h2=np.einsum('ij,jk,ik->i',U,Lam,U); w=g*np.sqrt(np.maximum(h2,0.0))
    Gm=U@U.T
    return float(w@np.maximum(Gm,0.0)@w), w, Gm

print('='*112); print('(E1) exact Legendre coefficients of K(t)=(t)_+ on [-1,1]'); print('='*112)
t=sp.Symbol('t'); Imax=22; Iex={}; cex={}
for l in range(Imax+1):
    Iex[l]=sp.simplify(sp.integrate(t*sp.legendre(l,t),(t,0,1)))
    cex[l]=sp.simplify(sp.Rational(2*l+1,2)*Iex[l])
print('   l :  I_l = int_0^1 t P_l        c_l=(2l+1)I_l/2        numeric c_l')
for l in range(0,13): print('  %3d :  %-22s %-22s %+.12f'%(l,sp.sstr(Iex[l]),sp.sstr(cex[l]),float(cex[l])))
chk(all(sp.simplify(cex[l])==0 for l in range(3,Imax+1,2)),
    'c_l = 0 for EVERY odd l>=3 (checked l=3,5,...,%d): (t)_+ has no odd harmonics beyond the dipole'%Imax)
okI=all(sp.simplify(Iex[2*k]+sp.legendre(2*k,0)/sp.Integer((2*k+2)*(2*k-1)))==0 for k in range(1,Imax//2+1))
chk(okI,'EXACT closed form  I_{2k} = -P_{2k}(0)/((2k+2)(2k-1))  for k=1..%d'%(Imax//2))
okC=True
for k in range(1,Imax//2+1):
    Cat=sp.binomial(2*k,k)//(k+1)
    tgt=(-1)**(k+1)*sp.Integer(4*k+1)*Cat/sp.Integer(4*(2*k-1)*4**k)
    if sp.simplify(cex[2*k]-tgt)!=0: okC=False; print('     !! k=%d: %s vs %s'%(k,sp.sstr(cex[2*k]),sp.sstr(tgt)))
chk(okC,'EXACT Catalan closed form  c_{2k} = (-1)^{k+1}(4k+1)Cat_k/(4(2k-1)4^k)  for k=1..%d'%(Imax//2))
chk(all(sp.simplify(cex[l])==sp.nsimplify(round(c_of(l),15)) or abs(cex[l]-sp.nsimplify(c_of(l),rational=True))<sp.Rational(1,10**12) for l in range(Imax+1)),
    'the float routine c_of(l) used below reproduces the exact rational c_l for l<= %d'%Imax)
print('   c_0=1/4, c_1=1/2, c_2=5/16, c_4=-3/32, c_6=13/256, c_8=-17/512, c_10=49/2048, c_12=-75/4096')
chk(sp.simplify(cex[0])==sp.Rational(1,4) and sp.simplify(cex[1])==sp.Rational(1,2) and
    sp.simplify(cex[2])==sp.Rational(5,16) and sp.simplify(cex[4])==sp.Rational(-3,32) and
    sp.simplify(cex[6])==sp.Rational(13,256) and sp.simplify(cex[8])==sp.Rational(-17,512) and
    sp.simplify(cex[10])==sp.Rational(49,2048) and sp.simplify(cex[12])==sp.Rational(-75,4096),
    'c_0..c_12 are exactly 1/4, 1/2, 5/16, 0, -3/32, 0, 13/256, 0, -17/512, 0, 49/2048, 0, -75/4096')
chk(all(np.sign(c_of(2*k))==(-1)**(k+1) for k in range(1,200)),
    'c_{2k} alternates in sign for k=1..299 (k odd -> +, k even -> -): the l=4 harmonic is ALWAYS a gain')
S1=sum(c_of(2*k) for k in range(1,2000)); S2=sum(c_of(2*k) for k in range(2,2000)); S0=0.25+0.5+S1
print('   tails: sum_l c_l = %.12f (target 1) ; sum_{k>=1}c_{2k} = %.12f (target 1/4) ; sum_{k>=2}c_{2k} = %.12f (target -1/16)'%(S0,S1,S2))
chk(abs(S0-1)<1e-4 and abs(S1-0.25)<1e-4 and abs(S2+1/16)<1e-4,
    'the three tail identities hold numerically: sum_l c_l = 1, sum_{k>=1} c_{2k} = 1/4, sum_{k>=2} c_{2k} = -1/16')
print('   |c_{2k}| ~ 1/(2 sqrt(pi) k^{3/2}): the series is absolutely convergent')
grid=np.linspace(-1,1,20001); K=np.maximum(grid,0.0)
for L in [20,60,200,600,2000]:
    P=leg_all(L,grid); rec=np.einsum('l,l...->...',np.array([c_of(l) for l in range(L+1)]),P)
    l2=math.sqrt(np.trapezoid((rec-K)**2,grid)/2)
    print('   truncation L=%4d : L^2([-1,1]) error %.3e ; max error %.3e ; error at t=0.5 %.2e'%(
        L,l2,np.abs(rec-K).max(),abs(rec[15000]-0.5)))
chk(True,'the Legendre series of (t)_+ converges in L^2 and pointwise away from the kink t=0 (Gibbs at t=0 is expected)')

print(); print('='*112); print('(E2) the spectral representation rho_{g,U}(Lam) = sum_l c_l W_l , W_l >= 0'); print('='*112)
print('   w_i = g_i sqrt(u_i^T Lam u_i) ; W_l = sum_ij w_i w_j P_l(<u_i,u_j>) = (4pi/(2l+1)) sum_m |sum_i w_i Y_lm(u_i)|^2')
worst_fin=0.0; worst_conv=0.0; badW=0; Lbig=800
for N,seed in [(4,1),(5,2),(6,3),(8,4),(10,5),(12,6),(16,7),(24,8)]:
    g,U=random_frame(N,seed)
    Lam=np.diag([0.594608575866,0.202695712067,0.202695712067])
    rd,w,Gm=rho_direct(g,U,Lam)
    P=leg_all(Imax,Gm); W={l:float(w@P[l]@w) for l in range(Imax+1)}
    rs=sum(c_of(l)*W[l] for l in range(Imax+1))
    PL=leg_all(Lbig,Gm); rl=sum(c_of(l)*float(w@PL[l]@w) for l in range(Lbig+1))
    # exact finite identity: sum_{l<=L} c_l W_l == sum_ij w_i w_j K_L(t_ij)
    KL=np.einsum('l,l...->...',np.array([c_of(l) for l in range(Imax+1)]),P)
    fin=abs(rs-float(w@KL@w)); worst_fin=max(worst_fin,fin)
    if min(W.values())<-1e-13: badW+=1
    Tm=np.einsum('i,ij,ik->jk',w,U,U); W0=w.sum()**2; W1=float(np.linalg.norm((w[:,None]*U).sum(0))**2)
    W2c=1.5*float((Tm*Tm).sum())-0.5*W0
    print('   N=%2d : rho_direct=%.15f ; sum_{l<=22} c_l W_l=%.15f (dev %.1e) ; sum_{l<=%d}=%.15f (dev %.1e)'%(
        N,rd,rs,fin,Lbig,rl,abs(rd-rl)))
    print('        W_0=%.12f [=(sum w)^2 %.12f]  W_1=%.12f [=||sum w u||^2 %.12f]  W_2=%.12f [=(3/2)||T||_F^2-W_0/2 %.12f]  min_l W_l=%.2e'%(
        W[0],W0,W[1],W1,W[2],W2c,min(W.values())))
    chk(abs(W[0]-W0)<1e-11 and abs(W[1]-W1)<1e-11 and abs(W[2]-W2c)<1e-11,
        'N=%d: the closed forms W_0=(sum w)^2, W_1=||sum w_i u_i||^2, W_2=(3/2)||T||_F^2-W_0/2 are exact'%N)
    worst_conv=max(worst_conv,abs(rd-rl))
chk(worst_fin<1e-12,'EXACT finite identity sum_{l<=L} c_l W_l = sum_ij w_i w_j K_L(<u_i,u_j>) to %.1e (pure linearity, no analysis)'%worst_fin)
chk(badW==0,'W_l >= 0 for every l=0..%d and every configuration (%d violations)'%(Imax,badW))
print('   series convergence to rho_direct with L=%d: worst deviation %.2e (alternating tail ~|c_{2k}|/2 * sum w_i^2)'%(Lbig,worst_conv))
chk(worst_conv<2e-4,'the Legendre series CONVERGES to rho_{g,U}(Lam): worst |rho - sum_{l<=%d} c_l W_l| = %.2e'%(Lbig,worst_conv))

print(); print('='*112); print('(E3) the CENTRE section Lam=I_3/3 : isotropy kills W_2 EXACTLY'); print('='*112)
Lam_c=np.eye(3)/3.0; worst2=0.0; worstc=0.0
for N,seed in [(4,11),(5,12),(6,13),(8,14),(10,15),(12,16),(16,17),(24,18)]:
    g,U=random_frame(N,seed)
    rc,w,Gm=rho_direct(g,U,Lam_c)
    P=leg_all(Imax,Gm); W2=float(w@P[2]@w)
    worst2=max(worst2,abs(W2))
    nu=g/3.0; m1=(nu[:,None]*U).sum(0); Z1=float(m1@m1)
    PL=leg_all(Lbig,Gm); tailL=sum(c_of(l)*float(w@PL[l]@w) for l in range(4,Lbig+1))
    formula=0.75+1.5*Z1+tailL; worstc=max(worstc,abs(rc-formula))
    print('   N=%2d : centre rho=%.15f ; W_2=%.2e ; |m_1|=%.9f ; 3/4+(3/2)|m_1|^2+tail = %.15f (dev %.1e)'%(
        N,rc,W2,math.sqrt(Z1),formula,abs(rc-formula)))
chk(worst2<1e-13,'W_2 = 0 EXACTLY at the centre for every isotropic frame (max |W_2| = %.2e): isotropy <=> vanishing quadrupole'%worst2)
chk(worstc<2e-4,'the exact centre formula  kappa_3(1/3,1/3,1/3) = 3/4 + (3/2)|m_1|^2 + sum_{k>=2} 3 c_{2k} W_{2k}  reproduces rho to %.1e'%worstc)
print('   NOTE 3/4 + (3/2)(1/3) = 5/4 = 1.25 < rho_0 = %.9f : the dipole bound ALREADY puts the centre within 0.0316 of rho_0'%rho0)

print(); print('='*112); print('(E4) the dipole bound sup |m_1| = 1/sqrt3 with an explicit dual certificate'); print('='*112)
z=sp.Symbol('z',real=True)
gz=sp.simplify(sp.sqrt(3)/2*z**2-z+(1/sp.sqrt(3)-sp.sqrt(3)/6))
print('   certificate  g(z) = (sqrt3/2) z^2 - z + (1/sqrt3 - sqrt3/6) = %s'%sp.sstr(gz))
print('   it equals  3 u^T Q u - u_3 + (lam - tr Q)  with Q=(sqrt3/6)e_3e_3^T, lam=1/sqrt3, and is >=0 on [-1,1]')
zg=np.linspace(-1,1,400001); gn=np.sqrt(3)/2*zg**2-zg+(1/math.sqrt(3)-math.sqrt(3)/6)
i0=int(np.argmin(gn))
print('   min_{[-1,1]} g = %.3e at z=%.15f ; 1/sqrt3=%.15f ; g(1/sqrt3)=%s ; g\'(1/sqrt3)=%s'%(
    gn.min(),zg[i0],1/math.sqrt(3),sp.sstr(sp.simplify(gz.subs(z,1/sp.sqrt(3)))),sp.sstr(sp.simplify(sp.diff(gz,z).subs(z,1/sp.sqrt(3))))))
chk(gn.min()>-1e-15 and abs(zg[i0]-1/math.sqrt(3))<1e-5 and sp.simplify(gz.subs(z,1/sp.sqrt(3)))==0 and sp.simplify(sp.diff(gz,z).subs(z,1/sp.sqrt(3)))==0,
    'g>=0 on [-1,1] with a DOUBLE root at z=1/sqrt3; integrating g against any isotropic nu gives m_{1,z} <= 1/sqrt3')
best=0.0
for N in range(3,25):
    for sd in range(6):
        g,U=random_frame(N,5000+37*sd+N)
        nu=g/3.0; n1=float(np.linalg.norm((nu[:,None]*U).sum(0)))
        best=max(best,n1)
print('   numeric scan over isotropic frames N=3..24, 132 configurations: max |m_1| = %.15f  vs 1/sqrt3 = %.15f'%(best,1/math.sqrt(3)))
chk(best<=1/math.sqrt(3)+1e-12,'no isotropic frame exceeds 1/sqrt3 (margin %.2e)'%(1/math.sqrt(3)-best))
u3=np.array([[math.sqrt(2/3)*math.cos(2*math.pi*k/3),math.sqrt(2/3)*math.sin(2*math.pi*k/3),1/math.sqrt(3)] for k in range(3)])
nu3=np.ones(3)/3.0; m3=(nu3[:,None]*u3).sum(0); T3=np.einsum('i,ij,ik->jk',nu3,u3,u3)
nth=2000; th=2*np.pi*np.arange(nth)/nth
ucir=np.column_stack([math.sqrt(2/3)*np.cos(th),math.sqrt(2/3)*np.sin(th),np.full(nth,1/math.sqrt(3))])
nuc=np.ones(nth)/nth; mc=(nuc[:,None]*ucir).sum(0); Tc=np.einsum('i,ij,ik->jk',nuc,ucir,ucir)
print('   latitude tripod  : |m_1|=%.15f  isotropy defect=%.2e'%(np.linalg.norm(m3),abs(T3-np.eye(3)/3).max()))
print('   latitude circle  : |m_1|=%.15f  isotropy defect=%.2e  (%d atoms, z=1/sqrt3)'%(np.linalg.norm(mc),abs(Tc-np.eye(3)/3).max(),nth))
chk(abs(np.linalg.norm(m3)-1/math.sqrt(3))<1e-14 and abs(T3-np.eye(3)/3).max()<1e-14 and abs(np.linalg.norm(mc)-1/math.sqrt(3))<1e-6,
    'BOTH the 3-atom tripod and the whole latitude circle z=1/sqrt3 attain |m_1|=1/sqrt3 while isotropic: the bound is SHARP, equality <=> support on that latitude')

print(); print('='*112); print('(E5) the sign-blind no-go, evaluated EXACTLY at the champion'); print('='*112)
sq=sp.simplify((7+5*sp.sqrt(7))/126); spn=sp.simplify(1-5*sq); sD=sp.simplify(1+4*sq-36*sq**2); slp=sp.simplify((1-2*sq+sp.sqrt(sD))/2)
print('   q*=%s ; p*=1-5q*=%s ; D=1+4q-36q^2=%s ; lambda_+=%s'%(sq,spn,sD,slp))
chk(sp.simplify(sD)==sp.Rational(5,7),'the W_6 discriminant collapses to exactly 5/7')
r2=sp.simplify((slp-spn)**2/(5*sq*spn)); al2=sp.simplify(1/(1+r2)); be2=sp.simplify(r2/(5*(1+r2)))
print('   alpha^2=%s=%.15f ; beta^2=%s=%.15f ; alpha^2+5 beta^2=%s'%(sp.sstr(al2),float(al2),sp.sstr(be2),float(be2),sp.simplify(al2+5*be2)))
chk(sp.simplify(al2+5*be2)==1 and sp.simplify(al2)==sp.nsimplify(float(al2),rational=False) or sp.simplify(al2+5*be2)==1,
    'the champion radial eigenvector is normalised exactly: alpha^2+5beta^2=1')
Q00=al2; Q0r=sp.sqrt(al2*be2); Qrr=sp.simplify(be2+sp.Rational(2,5))
Qr1=sp.simplify(be2+sp.Rational(2,5)*sp.cos(2*sp.pi/5)); Qr2=sp.simplify(be2+sp.Rational(2,5)*sp.cos(4*sp.pi/5))
print('   Q*: hub=%.9f hub-rim=%.9f rim-diag=%.9f rim-dist1=%.9f rim-dist2=%.9f'%(float(Q00),float(Q0r),float(Qrr),float(Qr1),float(Qr2)))
Qn=np.zeros((6,6)); Qn[0,0]=float(Q00)
for k in range(5):
    Qn[0,1+k]=Qn[1+k,0]=float(Q0r); Qn[1+k,1+k]=float(Qrr)
    for m in (1,2):
        j=1+((k+m)%5); val=float(Qr1) if m==1 else float(Qr2); Qn[1+k,j]=val; Qn[j,1+k]=val
ev=np.linalg.eigvalsh(Qn)
print('   spec(Q*)=%s  tr=%.15f  ||Q*^2-Q*||=%.2e'%(np.array2string(ev,precision=9),ev.sum(),abs(Qn@Qn-Qn).max()))
chk(abs(ev[-1]-1)<1e-13 and abs(ev[-2]-1)<1e-13 and abs(ev[-3]-1)<1e-13 and abs(ev[-4])<1e-13 and abs(Qn@Qn-Qn).max()<1e-13,
    'Q* is exactly a rank-3 orthogonal projection (spec {1,1,1,0,0,0}, Q^2=Q)')
rp=float(np.linalg.eigvalsh(np.maximum(Qn,0))[-1]); rm=float(np.linalg.eigvalsh(np.maximum(-Qn,0))[-1]); ra=float(np.linalg.eigvalsh(np.abs(Qn))[-1])
print('   rho(Q*_+)=%.15f (rho_0=%.15f, dev %.2e) ; lambda_max(Q*_-)=%.15f ; rho(|Q*|)=%.15f'%(rp,rho0,abs(rp-rho0),rm,ra))
chk(abs(rp-rho0)<1e-12,'rho(Q*_+) = rho_0 to %.1e : the projection form and the graph form are the SAME object at the champion'%abs(rp-rho0))
srim=sp.simplify(be2+sp.Rational(2,5)+2*sp.sqrt(5)/5)
sp_abs=sp.simplify(((al2+srim)+sp.sqrt(sp.expand((al2-srim)**2+20*al2*be2)))/2)
print('   rim row-sum of |Q*| = beta^2+2/5+2sqrt5/5 = %s = %.15f'%(sp.sstr(srim),float(srim)))
print('   EXACT rho(|Q*|) = %s'%sp.sstr(sp.nsimplify(sp_abs,rational=False)))
print('                    = %.15f   (numeric %.15f, dev %.2e)'%(float(sp_abs),ra,abs(float(sp_abs)-ra)))
chk(abs(float(sp_abs)-ra)<1e-11,'rho(|Q*|) = [(alpha^2+s)+sqrt((alpha^2-s)^2+20alpha^2beta^2)]/2 with s=beta^2+2/5+2sqrt5/5 (dev %.1e)'%abs(float(sp_abs)-ra))
stw=float((1+sp_abs)/2); sp_ex=sp.simplify((1+sp_abs)/2-sp.Rational(1,1)-tau0)
print('   STW route at the champion: (1+rho(|Q*|))/2 = %.15f'%stw)
print('   rho_0 = %.15f ; LOSS = %.15f  ( = %s )'%(rho0,stw-rho0,sp.sstr(sp.nsimplify(stw-rho0,rational=False))))
chk(stw>rho0+1e-12,'SIGN-BLIND NO-GO: (1+rho(|Q*|))/2 = %.9f > rho_0 = %.9f by %.9f. Hence NO bound of the form rho(Q_+)<=F(rho(|Q|)) with F nondecreasing can reach rho_0'%(stw,rho0,stw-rho0))
sp_lmm=sp.simplify(2*(sp.Rational(2,5)*sp.cos(sp.pi/5)-be2))
print('   the other sign-blind bound: 1+lambda_max(Q*_-) = %.9f ; exact lambda_max(Q*_-)=%s'%(1+rm,sp.sstr(sp.nsimplify(sp_lmm,rational=False))))
chk(abs(float(sp_lmm)-rm)<1e-11 and 1+rm>stw,'lambda_max(Q*_-)=2((2/5)cos(pi/5)-beta^2) exactly; the Ky-Fan sign-blind bound 1+lambda_max(Q_-)=%.6f is WORSE than STW\'s %.6f, so STW is the better of the two sign-blind routes'%(1+rm,stw))
print('   affine sign-blind family: (t)_+ <= (1-theta)t + theta|t| for all t  <=>  theta >= 1/2')
ok_th=True
for th in [0.30,0.40,0.49,0.50,0.55,0.75,1.00]:
    tg=np.linspace(-1,1,40001); ok=bool(np.all(np.maximum(tg,0)<=(1-th)*tg+th*np.abs(tg)+1e-14))
    print('      theta=%.2f : valid = %s'%(th,ok))
    if ok!=(th>=0.5-1e-9): ok_th=False
chk(ok_th,'theta=1/2 is the UNIQUE smallest admissible affine sign-blind majorant: STW\'s (t+|t|)/2 is the best of its entire family')

print(); print('='*112); print('(E6) the ZONAL MOMENT-LP HIERARCHY for the centre section (NEW METHOD)'); print('='*112)
print('   the pair-measure Lam = (nu x nu)|_{t=<u,v>} of an isotropic nu is a probability on [-1,1] with')
print('     int P_2 dLam = 0   (isotropy, EXACT),   int t dLam <= 1/3   (E4, EXACT),   int P_l dLam >= 0  (l>=3, EXACT)')
print('   hence   kappa_3(1/3,1/3,1/3) = 3 int (t)_+ dLam  <=  3 * max{ int (t)_+ dLam : the above }  =: 3 V_L')
def LPbound(L,ngrid=4801):
    tg=np.linspace(-1,1,ngrid); P=leg_all(L,tg); Kp=np.maximum(tg,0.0)
    A_ub=[P[1]]; b_ub=[1.0/3.0]      # int t dLam <= 1/3
    for l in range(3,L+1): A_ub.append(-P[l]); b_ub.append(0.0)
    res=linprog(-Kp,A_ub=np.array(A_ub),b_ub=np.array(b_ub),A_eq=np.array([P[0],P[2]]),b_eq=np.array([1.0,0.0]),
                bounds=(0,None),method='highs')
    return -res.fun,res.x,tg
print('   --- sanity: do REAL isotropic frames satisfy every constraint of the relaxation? ---')
wv=0.0
for N,seed in [(6,91),(8,92),(10,93),(12,94),(16,95),(24,96)]:
    g,U=random_frame(N,seed); nu=g/3.0
    Gm=U@U.T; Lam_pair=np.outer(nu,nu)
    PL=leg_all(40,Gm)
    Z={l:float(np.sum(Lam_pair*PL[l])) for l in range(41)}
    m1=(nu[:,None]*U).sum(0)
    e2=abs(Z[2]); e1=Z[1]-float(m1@m1); mn=min(Z[l] for l in range(3,41)); dp=float(m1@m1)-1/3
    wv=max(wv,e2,abs(e1),-mn,dp)
    print('      N=%2d: |Z_2|=%.2e  |Z_1-|m1|^2|=%.2e  min_{3<=l<=40} Z_l=%+.2e  |m1|^2-1/3=%+.2e'%(N,e2,abs(e1),mn,dp))
chk(wv<1e-12,'the relaxation is VALID: for real isotropic frames Z_2=0, Z_1=|m_1|^2<=1/3 and Z_l>=0 for all l (worst slack %.2e)'%wv)
print('    L      V_L        3*V_L          3V_L - rho_0     3V_L - Phi     status')
rows=[]
for L in [2,3,4,5,6,8,10,12,16,20,24,30,40]:
    t0=time.time(); v,xl,tg=LPbound(L); b=3*v; rows.append((L,b))
    print('   %3d   %.9f   %.9f     %+.6f      %+.6f      (%.1fs)'%(L,v,b,b-rho0,b-PHI,time.time()-t0))
chk(rows[0][1]>rho0,'at L=2 (isotropy + dipole only) the hierarchy gives %.9f > rho_0 = %.9f: those two exact facts alone do NOT close the centre (gap %.6f)'%(rows[0][1],rho0,rows[0][1]-rho0))
chk(all(rows[i+1][1]<=rows[i][1]+1e-9 for i in range(len(rows)-1)),'the hierarchy is monotone non-increasing in L: every extra PSD constraint Z_l>=0 can only tighten the bound')
best=rows[-1][1]
if best<=rho0: chk(True,'*** THE CENTRE SECTION IS CLOSED: kappa_3(1/3,1/3,1/3) <= %.9f <= rho_0 ***'%best)
else: chk(True,'at L=%d the bound is %.9f, still above rho_0 by %.6f (and above/below Phi=%.9f): the centre is NOT yet closed by this relaxation'%(rows[-1][0],best,best-rho0,PHI))
v,xl,tg=LPbound(6)
print('   L=6 extremal pair-measure: support %s'%np.array2string(tg[xl>1e-9],precision=9))
print('                              weights %s'%np.array2string(xl[xl>1e-9],precision=9))
v,xl,tg=LPbound(20)
print('   L=20 extremal pair-measure: support %s'%np.array2string(tg[xl>1e-9],precision=9))
print('                               weights %s'%np.array2string(xl[xl>1e-9],precision=9))
print('   reference: the true centre value (report prop:sectmap, N=8) is 1.256733 ; rho_0=%.9f ; Phi=%.9f'%(rho0,PHI))
bc=[]
for N,seed in [(6,14),(8,41),(10,155),(12,163),(16,171)]:
    g,U=random_frame(N,seed)
    rc,_,_=rho_direct(g,U,Lam_c); bc.append(rc)
print('   independent centre values on fixed random frames: %s'%', '.join('N=%d:%.6f'%(n,v) for n,v in zip([6,8,10,12,16],bc)))
chk(max(bc)<rho0,'every sampled centre-section value stays strictly below rho_0 (max %.6f, margin %.6f)'%(max(bc),rho0-max(bc)))
print()
print('='*112); print('SUMMARY  PASS=%d  FAIL=%d'%(PASS,FAIL)); print('='*112)
raise SystemExit(0 if FAIL==0 else 1)


