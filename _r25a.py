# -*- coding: utf-8 -*-
"""R25-a: 单轴截面的谐波分解恒等式 + 极值剖面探测."""
import numpy as np, itertools, sys
from numpy.linalg import eigh, norm

S5 = np.sqrt(5.0); S35 = np.sqrt(35.0); S3 = np.sqrt(3.0)
rho0 = (9+S5+2*S35)/18
gam = (1+S5)/16; alp = (7*S5-5)/8; bet = (25-15*S5)/16
cW = (2*S5-1)/14; cS = (10+S5)/14; cZ = (5-3*S5)/14
K1 = 9*gam+3*alp; B = -9*bet
del0 = 3*(1+S5)/2+3-6*rho0
print("rho0=%.15f  gam=%.10f alp=%.10f bet=%.10f" % (rho0,gam,alp,bet))
print("cW=%.10f cS=%.10f cZ=%.10f  K1=%.10f B=%.10f delta0=%.12f" % (cW,cS,cZ,K1,B,del0))
print("check cW=gam-3bet/35: %.3e ; cS=alp+6bet/7: %.3e ; cZ=8bet/35: %.3e"
      % (cW-(gam-3*bet/35), cS-(alp+6*bet/7), cZ-8*bet/35))
print("check phi+1-2rho0 = delta0/3: %.3e" % (( (1+S5)/2 +1-2*rho0) - del0/3))

def P2(t): return (3*t*t-1)/2
def P4(t): return (35*t**4-30*t*t+3)/8

def signmax(vecs, w):
    """A_w = max over sign patterns eps of || sum eps_i w_i u_i || ; returns (A, eps, z)."""
    n = len(w); best = (-1,None,None)
    for eps in itertools.product([1,-1], repeat=n-1):
        e = np.array(eps+(1,))
        s = (e*w)@vecs; ns = norm(s)
        if ns > best[0]: best = (ns, e, s/ns if ns>0 else np.array([0,0,1.]))
    return best

def analyse(g, U, rhoa, e=None, verbose=False):
    """uniaxial Lambda = b I + rho e e^T, b=(1-rho)/3.  returns dict."""
    if e is None: e = np.array([0,0,1.0])
    b = (1-rhoa)/3
    t = U@e
    h = np.sqrt(b+rhoa*t*t)
    w = g*h
    W = w.sum()
    Sw = (w[:,None,None]*(U[:,:,None]*U[:,None,:])).sum(0)
    nrmSF = (Sw*Sw).sum()
    n = len(g)
    G = U@U.T
    M4 = (np.outer(w,w)*G**4).sum()
    Z4 = (np.outer(w,w)*P4(G)).sum()
    Z2 = (np.outer(w,w)*P2(G)).sum()
    Lw = (np.outer(w,w)*abs(G)).sum()
    A, eps, z = signmax(U, w)
    aw = abs(U@z)
    Vw = (np.outer(w,w)*(aw[:,None]-aw[None,:])**2).sum()
    zet = z@Sw@z
    qw = gam+alp*G**2+bet*G**4
    Dw = (np.outer(w,w)*(qw-abs(G))).sum()
    Jw = Lw+A*A
    lhs = cW*W*W+cS*nrmSF+W*zet
    rhs = -cZ*Z4+Dw+0.5*Vw
    return dict(W=W,Sw=Sw,nrmSF=nrmSF,lam=eigh(Sw)[0],zet=zet,M4=M4,Z4=Z4,Z2=Z2,
                Lw=Lw,A=A,Aw2=A*A,Vw=Vw,Dw=Dw,Jw=Jw,half=Jw/2,
                id1=M4-(-3*W*W/35+6*nrmSF/7+8*Z4/35),
                id2=Sw-nrmSF, id3=Lw-(gam*W*W+alp*nrmSF+bet*M4-Dw),
                id4=A*A-(W*zet-0.5*Vw),
                id5=Jw-(lhs-rhs), id6=(2*rho0-Jw)-(2*rho0-lhs-rhs),
                t=t,h=h,w=w,eps=eps,z=z)

# ---------- champion: pole + regular 5-gon ----------
astar = (3*S35+5*np.sqrt(7)-7-7*S5)/14
zstar = 0.5004400497682123079777569143187991903444
rho_star = (3*astar-1)/2
bstar = (1-rho_star)/3
Gcon = 2/(5*(1-zstar**2)); gh = 1-5*Gcon*zstar**2
print("\nchampion: a*=%.15f rho*=%.15f b=%.15f g_hub=%.12f g_con=%.12f"
      % (astar,rho_star,bstar,gh,Gcon))
U6 = [np.array([0,0,1.0])]
for k in range(5):
    ph = 2*np.pi*k/5
    U6.append(np.array([np.sqrt(1-zstar**2)*np.cos(ph), np.sqrt(1-zstar**2)*np.sin(ph), zstar]))
U6 = np.array(U6); g6 = np.array([gh]+[Gcon]*5)
print("isotropy defect: %.3e ; sum g-3 = %.3e" % (np.abs((g6[:,None,None]*(U6[:,:,None]*U6[:,None,:])).sum(0)-np.eye(3)).max(), g6.sum()-3))
d = analyse(g6,U6,rho_star)
print("champion J_w/2 = %.15f   (rho0=%.15f)  diff=%.3e" % (d['half'],rho0,d['half']-rho0))
for k in ['W','nrmSF','zet','M4','Z4','Lw','Aw2','Vw','Dw']:
    print("   %-7s = %.12f" % (k, d[k]))
print("   lam(S_w) =", np.round(d['lam'],9), " sqrt(b)=%.9f sqrt(a)=%.9f"%(np.sqrt(bstar),np.sqrt(astar)))
print("   identity residuals: M4decomp=%.2e  Lwdecomp=%.2e  Aw2decomp=%.2e  Jwdecomp=%.2e"
      % (d['id1'],d['id3'],d['id4'],d['id5']))
print("   main-term lhs=%.12f   corrections rhs=%.12f   2rho0-Jw=%.12f" % (d['zet']*d['W']+cW*d['W']**2+cS*d['nrmSF'], -cZ*d['Z4']+d['Dw']+0.5*d['Vw'], 2*rho0-d['Jw']))
