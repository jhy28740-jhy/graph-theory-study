# -*- coding: utf-8 -*-
# _r12.py  R19: (A) correct FD check of the Danskin gradient ; (B) NEW SECTION DECOMPOSITION
#   kappa_3 = sup_{Lam} sup_{frames} rho_{g,U}(Lam) = sup_{mu in Bloch simplex} kappa_3(mu),
#   rho(V;Lam) = sum_ij (<v_i,v_j>)_+ sqrt(v_i'Lam v_i) sqrt(v_j'Lam v_j)   [x_i = sqrt(v_i'Lam v_i)]
import numpy as np, math
np.set_printoptions(precision=12, suppress=False, linewidth=200)
exec(open('_o1.py').read().split('print("="*104)')[0])
s5=math.sqrt(5.0); rho0=(9+s5+2*math.sqrt(35.0))/18
PASS=FAIL=0
def chk(c,m):
    global PASS,FAIL
    if c: PASS+=1; print('  [PASS]',m)
    else: FAIL+=1; print('  [FAIL]',m)
def rhoLam(U,g,L):
    h=np.sqrt(np.maximum(np.einsum('ij,jk,ik->i',U,L,U),0.0)); gh=g*h
    m=(U*gh[:,None]).sum(0); return float(m@m)+float((np.outer(gh,gh)*np.maximum(-(U@U.T),0.0)).sum())
def gradformula(U,g,L):
    Lin=g*np.einsum('ij,jk,ik->i',U,L,U); xx=np.sqrt(np.maximum(Lin,1e-300))
    Qp=np.maximum(np.sqrt(np.outer(g,g))*(U@U.T),0.0); y=(Qp@xx)/xx
    return np.einsum('i,ij,ik->jk',y,U,U)
print('='*104); print('(A) Danskin closed form  grad_Lam rho = sum_i ((Q+ x)_i/x_i) g_i u_i u_i^T   (central FD on PSD-interior Lam)'); print('='*104)
rng=np.random.default_rng(4242); w=0.0
for t in range(40):
    B=rng.normal(size=(3,3)); L=B@B.T; L=L/np.trace(L)
    ev,EV=np.linalg.eigh(L); ev=np.maximum(ev,0.06); L=(EV*ev)@EV.T; L/=np.trace(L)
    Gf=gradformula(Uc,gc,L); e=1e-6; Gn=np.zeros((3,3))
    for a in range(3):
        for b in range(a,3):
            E=np.zeros((3,3)); E[a,b]=1.0
            if a!=b: E[b,a]=1.0
            d=(rhoLam(Uc,gc,L+e*E)-rhoLam(Uc,gc,L-e*E))/(2*e)
            Gn[a,b]=d; Gn[b,a]=d
    dev=np.abs(Gf-Gn).max()/max(1.0,np.abs(Gf).max()); w=max(w,dev)
print('  champion: max rel dev (formula vs FD) over 40 Lam = %.3e'%w); chk(w<1e-5,'Danskin gradient formula exact (rel dev %.2e)'%w)
mn=min(np.linalg.eigvalsh(gradformula(Uc,gc,(lambda A:A@A.T/np.trace(A@A.T))(rng.normal(size=(3,3))))).min() for _ in range(80))
print('  min eig grad rho(Lam) over 80 random Lam = %+.6e  (>=0)'%mn); chk(mn>-1e-12,'grad rho(Lam) is PSD for every Lam (min eig %+.2e)'%mn)
print('  NOTE: grad rho(Lam) = S(T*(Lam)) with T*(Lam)_ij = x_j(Lam)/x_i(Lam); Danskin applies because rho=min_T <Lam,S(T)>.')
print()
print('='*104); print('(B) NEW: section decomposition  kappa_3 = sup_{mu} kappa_3(mu),  kappa_3(mu)=sup_{V in St(N,3)} rho(V;diag mu)'); print('='*104)
def rhoVL(V,L):
    Q=V@V.T; s=np.sqrt(np.maximum(np.einsum('ij,jk,ik->i',V,L,V),0.0))
    return float((np.maximum(Q,0.0)*np.outer(s,s)).sum())
def gradVL(V,L):
    Q=V@V.T; M=(Q>0).astype(float); s=np.sqrt(np.maximum(np.einsum('ij,jk,ik->i',V,L,V),1e-300))
    Lv=V@L.T; A=(M*s)@V; B=(M*(Q*np.outer(s,s)))@Lv/s[:,None]
    return 2.0*(s[:,None]*A+B)
def ascend(V,L,iters=2500,step=0.05):
    best=(-1e18,None)
    for it in range(iters):
        G=gradVL(V,L); G=G-V@((V.T@G+G.T@V)/2.0); V=V+step*G
        Qr,_=np.linalg.qr(V); V=Qr[:,:3]
        if it%25==0:
            val=rhoVL(V,L)
            if val>best[0]: best=(val,V.copy())
    val=rhoVL(V,L)
    if val>best[0]: best=(val,V.copy())
    return best
def sec(mu,N,seeds):
    L=np.diag(mu); b=(-1e18,None)
    for sd in seeds:
        r=np.random.default_rng(sd); Q0,_=np.linalg.qr(r.normal(size=(N,3))); V=Q0[:,:3]
        val,VV=ascend(V,L); 
        if val>b[0]: b=(val,VV)
    return b
mus=[('center (1/3,1/3,1/3)',[1/3,1/3,1/3]),
     ('uniaxial a=0.45',[0.45,0.275,0.275]),
     ('uniaxial a=0.50',[0.50,0.25,0.25]),
     ('uniaxial a=0.5946 (=lam_max(Lam~) of champion)',[0.594608575866,0.202695712067,0.202695712067]),
     ('uniaxial a=0.70',[0.70,0.15,0.15]),
     ('uniaxial a=0.85',[0.85,0.075,0.075]),
     ('rank-2 (1/2,1/2,0)',[0.5,0.5,0.0]),
     ('rank-2 (0.7,0.3,0)',[0.7,0.3,0.0]),
     ('rank-1 (1,0,0)',[1.0,0.0,0.0]),
     ('asym (0.60,0.30,0.10)',[0.60,0.30,0.10]),
     ('asym (0.50,0.35,0.15)',[0.50,0.35,0.15])]
print('   section mu                        N=6        N=8       N=10       N=14      (kappa_3(mu), numerical)')
res={}
for tag,mu in mus:
    row=[]
    for N in (6,8,10,14):
        v,_=sec(mu,N,[1000+7*N+k for k in range(4)])
        row.append(v)
    res[tag]=max(row)
    print('   %-42s %s   max=%.12f  (%+.2e vs rho0)'%(tag,' '.join('%.9f'%z for z in row),max(row),max(row)-rho0))
print()
print('   champion check: rho(V_champ; Lam~) = %.15f  (rho0 = %.15f)'%(rhoVL(Uc*np.sqrt(gc)[:,None],np.diag([0.594608575866,0.202695712067,0.202695712067]) if False else None) if False else 0.0,rho0))
U,g=Uc,gc; Kf=U@U.T; ac=Kf>1e-13; Vc=U*np.sqrt(g)[:,None]
Qp=np.maximum(Vc@Vc.T,0.0); wv,ev=np.linalg.eigh(Qp); k=int(np.argmax(wv)); xx=ev[:,k].copy()
if xx.sum()<0: xx=-xx
h=xx/np.sqrt(g); W=np.array([((g*h)[ac[i]][:,None]*U[ac[i]]).sum(0) if ac[i].any() else np.zeros(3) for i in range(len(g))])
Lt=(1.0/wv[k])*np.einsum('i,ij,ik->jk',g*h,W,U)
print('   rho(V_champ; Lam~) = %.15f ;  spec(Lam~) = %s'%(rhoVL(Vc,Lt),np.array2string(np.linalg.eigvalsh(Lt)[::-1],precision=12)))
print('   => the champion realises kappa_3(mu*) >= rho0 at the UNIAXIAL section mu* = spec(Lam~).')
print()
print('第十二段汇总：PASS=%d FAIL=%d'%(PASS,FAIL))
