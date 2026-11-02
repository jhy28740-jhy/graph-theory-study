import numpy as np, sympy as sp
from numpy.linalg import eigh, norm
sq5=np.sqrt(5); sq3=np.sqrt(3)
R=lambda a,b: a/b
# ---------- exact constants ----------
d0=-R(25624284539519,53655078400000)+sq5/2-41*sq3/120
d1=d0+R(103,9000000); d=d1+R(1,3000000)
c=(3+sq5)/4-d/2
alpha=(2*sq5-1)/14; beta=(sq5+10)/14; rho=(15*sq5-25)/16
a_=(1+sq5)/6; b_=(7*sq5-5)/8
a3=R(12928,41503); a4=R(10240,41503)
c0=R(8352,41503); c1=R(43806,41503); c2=R(54016,41503)
be=R(16386,18865); bV=R(64608,290521); bJ=R(316672,290521)
dq=(3+sq5)/2-R(544886,207515)
g4c=R(9,4)-511*sq3/420; g3c=R(37,20)-11*sq3/12
eta=0.1; th=0.24; ell=R(937,1400)
E=R(5183683,5659500); Vc=R(145858703,581042000); Jc=R(7989533,14526050)
Bc=R(42565069,653672250); Cc=R(88944217,232416800)
BIG=R(725029049,278853120)
Qe=R(148831,181104000); QV=R(42847433,5229378000000); QJ=R(65082103,20917512000)
kap=R(3,125000)
I3=np.eye(3)
d_ij=np.eye(3)
OMEGA=(np.einsum('ij,kl->ijkl',d_ij,d_ij)+np.einsum('ik,jl->ijkl',d_ij,d_ij)+np.einsum('il,jk->ijkl',d_ij,d_ij))/15.0
def L3(y): return (np.einsum('i,jk->ijk',y,d_ij)+np.einsum('j,ik->ijk',y,d_ij)+np.einsum('k,ij->ijk',y,d_ij))/5.0
def L4(Z): return (np.einsum('ij,kl->ijkl',d_ij,Z)+np.einsum('ik,jl->ijkl',d_ij,Z)+np.einsum('il,jk->ijkl',d_ij,Z)
                  +np.einsum('jk,il->ijkl',d_ij,Z)+np.einsum('jl,ik->ijkl',d_ij,Z)+np.einsum('kl,ij->ijkl',d_ij,Z))/7.0
def Psi3(u): return np.einsum('i,j,k->ijk',u,u,u)-L3(u)
def Psi4(u): return np.einsum('i,j,k,l->ijkl',u,u,u,u)-OMEGA-L4(np.outer(u,u)-I3/3)
def f_(w): return c0+c1*w+c2*w**2-a3*w**3-a4*w**4
def sig(w): return f_(w)-np.abs(w)-w
def h0(z): return (1-z)*(sq5*z-1)**2*(sq5*z+sq5+2)/(2*(sq5+1)**2)
def g3f(t): return t**3-3*t/5
def g4f(t): return t**4-6*t**2/7+R(3,35)
def zstar(t): return t**4-R(85,64)*t**2+R(15,128)*t+R(27,128)
def gfun(w): return (1-9*w**2)/25

# ---------- sanity of harmonic kernels ----------
rng=np.random.default_rng(7)
def rnd_u():
    x=rng.normal(size=3); return x/norm(x)
err=0
for _ in range(400):
    u,v=rnd_u(),rnd_u(); t=u@v
    err=max(err,abs(np.sum(Psi3(u)*Psi3(v))-g3f(t)), abs(np.sum(Psi3(v)**2)-0.4),
               abs(np.sum(Psi4(u)*Psi4(v))-g4f(t)), abs(np.sum(Psi4(v)**2)-R(8,35)))
print("harmonic kernels g3,g4 max err =", f"{err:.3e}")

def analyse(p,Vm):
    """p masses, Vm rows = atoms v_i. returns dict of all quantities."""
    o={}
    C2=(p[:,None]*Vm*Vm).sum(0).sum()   # tr E XX^T
    Ri=norm(Vm,axis=1); Ui=Vm/np.where(Ri>0,Ri,1)[:,None]
    w=Ui@Ui.T                            # cos matrix
    A=(p[:,None]*p[None,:]*np.abs(Vm@Vm.T)).sum()
    m=(p[:,None]*Vm).sum(0); B=m@m; e=1-B
    C=(p*Ri).sum(); V=3-C**2
    D=(3+sq5)/2-A-B
    Mmat=np.einsum('i,ij,ik->jk',p,Vm,Ui)   # E[R U U^T] = sum p v u^T
    Z=Mmat-(C/3)*I3
    J=1-np.sum(Mmat*Mmat)
    S=(p[:,None]*p[None,:]*Ri[:,None]*Ri[None,:]*h0(np.abs(w))).sum()
    K2=(p[:,None]*p[None,:]*Ri[:,None]*Ri[None,:]*(a_+b_*(w**2-1/3)-rho*(w**4-1/5)-np.abs(w))).sum()
    P=(p[:,None]*p[None,:]*Ri[:,None]*Ri[None,:]*sig(w)).sum()
    T3=np.einsum('i,imno->mno',p*Ri,np.stack([Psi3(u) for u in Ui]))
    T4=np.einsum('i,imnop->mnop',p*Ri,np.stack([Psi4(u) for u in Ui]))
    H=np.sqrt(np.sum(T3**2)); K=np.sqrt(np.sum(T4**2))
    o.update(A=A,B=B,e=e,C=C,V=V,D=D,J=J,S=S,K2=K2,P=P,H=H,K=K,Mmat=Mmat,Z=Z,
             Ri=Ri,Ui=Ui,w=w,p=p,Vm=Vm,m=m,T3=T3,T4=T4)
    # regression
    if B>1e-12:
        x=np.sqrt(B); vh=m/x; q0=1/np.sqrt(3); q=C/(3*x)
        t=Ui@vh; Y=Ri*t
        L=(p*Ri*(t-q0)**2).sum()
        xi=Ui-(Vm@Mmat.T)
        w0=(p[:,None]*xi).sum(0); Jv=(p*(xi@vh)**2).sum()
        F0=2*V/9+J/3+q**2*e+2*q*np.sqrt(max(e*J,0))
        G0=V/3+3*q**2*e+2*q*np.sqrt(max(e*V,0))
        Lbar=(7*q**2/4)*e+(5/18)*V+(1/4)*J+(2*q/3)*np.sqrt(max(e*V,0))+(3*q/2)*np.sqrt(max(e*J,0))
        Lq=(p*Ri*(t-q)**2).sum()
        p3=(p*Ri*g3f(t)).sum(); p4=(p*Ri*g4f(t)).sum()
        r_reg=Lbar-L; r4=-p4-g4c+e/15+19*V/1400+53*L/28; r3=-p3-g3c-9*e/40+V/8+2*L
        R3=a3*np.sum((T3+(th/(2*a3))*Psi3(vh))**2); R4=a4*np.sum((T4+(eta/(2*a4))*Psi4(vh))**2)
        R0=R3+R4+eta*r4+th*r3
        o.update(x=x,vh=vh,q=q,L=L,F0=F0,G0=G0,Lbar=Lbar,Lq=Lq,p3=p3,p4=p4,
                 r_reg=r_reg,r4=r4,r3=r3,R0=R0,R3=R3,R4=R4,w0=w0,Jv=Jv,
                 p3b=np.sum(T3*Psi3(vh)),p4b=np.sum(T4*Psi4(vh)),
                 Lq_minus_L=(1/np.sqrt(3))*((C/np.sqrt(3)-x)**2)*((C/np.sqrt(3))+2*x)/x**2)
    # radial remainders
    RR=Ri[:,None]*Ri[None,:]
    zz=3*(Ri[:,None]-Ri[None,:])**2/(4*RR)
    r_rad=P-(3/50)*(V-9*J)
    r_tr=(p[:,None]*p[None,:]*RR*((8/9-zz)*sig(w)+zz*(sig(w)-gfun(w)))).sum()
    Zs=(p[:,None]*p[None,:]*RR*zstar(w)).sum()
    Zs_f=K**2+R(7,320)-R(15,128)*e-R(561,4480)*V+R(211,448)*J
    Cc_=Zs/125-3*C**2/125000
    r_ct=P-Cc_
    Cc_f=Zs_f/125-3*C**2/125000
    o.update(zz=zz,r_rad=r_rad,r_tr=r_tr,Zs=Zs,Zs_f=Zs_f,Ccal=Cc_,Ccalf=Cc_f,r_ct=r_ct,
             ratio=Ri.max()/Ri[Ri>0].min())
    if B>1e-12:
        q=o['q']; u_=np.sqrt(e); z_=np.sqrt(V); j_=np.sqrt(max(J,0))
        Q=(13/200)*(z_-(937/273)*q*u_)**2+(19/50)*(j_-(2811/2128)*q*u_)**2+(183/200-BIG*q**2)*e+Qe*e+QV*V+QJ*J
        Q7=(13/200)*(z_-(937/273)*q*u_)**2+(19/50)*(j_-(2811/2128)*q*u_)**2+(183/200-BIG*q**2)*e+R(10481,11319000)*e+R(305491,2614689000)*V+R(625833,232416800)*J
        o.update(Q=Q,Q7=Q7)
    return o

def resid(o):
    r={}
    r['budget4: D-(e+aV+bJ+rhoK2+S)'] = o['D']-(o['e']+alpha*o['V']+beta*o['J']+rho*o['K2']+o['S'])
    r['budget7: D-(dq+be e+bV V+bJ J+a3H2+a4K2+P)'] = o['D']-(dq+be*o['e']+bV*o['V']+bJ*o['J']+a3*o['H']**2+a4*o['K']**2+o['P'])
    r['Zstar: Zs - formula'] = o['Zs']-o['Zs_f']
    r['payment: r_rad>=0 check  P-(3/50)(V-9J)-r_rad'] = o['P']-(3/50)*(o['V']-9*o['J'])-o['r_rad']
    r['reserve: r_rad-(P/9+r_tr)'] = o['r_rad']-(o['P']/9+o['r_tr'])
    r['Ccal exact vs expanded'] = o['Ccal']-o['Ccalf']
    if 'Q' in o:
        r['p3 = <T3,Psi3(v)>'] = o['p3b']-o['p3']
        r['p4 = <T4,Psi4(v)>'] = o['p4b']-o['p4']
        r['a3H2 = -th p3 - th^2/(10a3)+R3'] = a3*o['H']**2-(-th*o['p3']-th**2/(10*a3)+o['R3'])
        r['a4K2 = -eta p4 - 2eta^2/(35a4)+R4'] = a4*o['K']**2-(-eta*o['p4']-2*eta**2/(35*a4)+o['R4'])
        r['Lq-L formula'] = o['Lq']-o['L']-o['Lq_minus_L']
        r['Lbar-(3/4F0+1/3G0)'] = o['Lbar']-(0.75*o['F0']+o['G0']/3)
        r['squares: Ee+VcV+JcJ-ell*Lbar - Q7'] = (E*o['e']+Vc*o['V']+Jc*o['J']-ell*o['Lbar'])-o['Q7']
        r['cert-start residual'] = o['D']-(d0+E*o['e']+Vc*o['V']+Jc*o['J']-ell*o['Lbar']+o['R0']+o['r_rad']+ell*o['r_reg'])
        r['cert8 residual'] = o['D']-(d1+o['Q']+o['R0']+ell*o['r_reg']+o['K']**2/1125+o['r_tr']+o['r_ct']/9)
        r['Q - Q7 - (paid corr)'] = o['Q']-o['Q7']-(-o['e']/9600+(-R(187,1680000)+R(1,375000))*o['V']+R(211,504000)*o['J'])
    return r

def rand_law(k,seed):
    g=np.random.default_rng(seed)
    Y=g.normal(size=(k,3)); w0=g.random(k)+0.2; w0/=w0.sum()
    S=np.einsum("i,ij,ik->jk",w0,Y,Y); ev,U=eigh(S); Yw=Y@(U@np.diag(1/np.sqrt(ev))@U.T)
    return w0,Yw

print("\n=== identity residuals on 12 random isotropic laws (k=5..12) ===")
worst={}
mins={'r_reg':1e9,'r4':1e9,'r3':1e9,'R0':1e9,'r_rad':1e9,'r_tr':1e9,'r_ct':1e9,'Q':1e9,'Q7':1e9,'D-d':1e9}
for s in range(12):
    k=5+(s%8); p,Vm=rand_law(k,100+s)
    o=analyse(p,Vm); r=resid(o)
    for kk,vv in r.items(): worst[kk]=max(worst.get(kk,0),abs(vv))
    for kk in mins:
        if kk=='D-d': mins[kk]=min(mins[kk],o['D']-d)
        elif kk in o: mins[kk]=min(mins[kk],o[kk])
    mins['ratio_check']=min(mins.get('ratio_check',1e9),o['ratio'])
for kk,vv in worst.items(): print(f"  max |residual| {kk:52s} = {vv:.3e}")
print("\n  minima over sample: r_reg=%.3e r4=%.3e r3=%.3e R0=%.3e"%(mins['r_reg'],mins['r4'],mins['r3'],mins['R0']))
print("                     r_rad=%.3e r_tr=%.3e r_ct=%.3e Q=%.3e"%(mins['r_rad'],mins['r_tr'],mins['r_ct'],mins['Q']))
print("  min(D-d) = %.6f   (d=%.10f)   max radius ratio seen=%.3f"%(mins['D-d'],d,mins['ratio_check']))
print("  c =",f"{c:.15f}","  d0=",f"{d0:.12f}"," d1=",f"{d1:.12f}"," d=",f"{d:.12f}")



