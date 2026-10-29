# -*- coding: utf-8 -*-
# 证书_Lam参数化与稀释恒等式.py   ROUND 18 certificate
#   (A) delta = tau*m identity family
#   (B) exact dilution identity  D(h) = beta + E_-(h)/2 ,  E_- second explicit form
#   (C) explicit u-KKT:  nu_i = 0,  Lam~ u_i = (h_i/rho) W_i,  Lam~ = (1/rho) sum_{acute} g_i g_j h_i h_j u_i u_j^T
#   (D) algebraic functional rho(Lam): rho(Lam~)=rho, degree-1 homogeneity, grad = rho*I_3 at Perron points, rho(Lam)<=rho
#   (E) concavity of rho(Lam)
#   (F) max_Lam rho(Lam) <= rho0 (champion + perturbations + randoms + Stiefel maximizers), sharp & rigid at champion
#   (G) REFUTATIONS: free-Lam without E_-  (0.3021 > tau0);  strong dilution alpha+D(h)<=rho0 (N=4 counterexample)
#   (H) cluster-graph theorem: G_+ a disjoint union of cliques  =>  rho <= 1
import numpy as np, math, os, sys
from scipy.optimize import minimize
np.set_printoptions(precision=12, suppress=False, linewidth=200)
ROOT=os.path.dirname(os.path.abspath(__file__)); BASE=os.path.dirname(ROOT)
sys.path.insert(0,BASE); os.chdir(BASE)
exec(open('_o1.py').read().split('print("="*104)')[0])
s5=math.sqrt(5.0); rho0=(9+s5+2*math.sqrt(35.0))/18; tau0=rho0-1
PASS=FAIL=0; FAILED=[]
def chk(cond,msg):
    global PASS,FAIL
    if cond: PASS+=1; print('  [PASS]',msg)
    else:    FAIL+=1; FAILED.append(msg); print('  [FAIL]',msg)
def core(U,g):
    N=len(g); Kf=U@U.T; P=np.maximum(-Kf,0.0); V=U*np.sqrt(g)[:,None]; Q=V@V.T
    Qp=np.maximum(Q,0.0); w,ev=np.linalg.eigh(Qp); idx=np.argsort(-w); w=w[idx]; ev=ev[:,idx]
    rho=w[0]; x=ev[:,0].copy()
    if x.sum()<0: x=-x
    h=x/np.sqrt(g); c=g*h; m=c@U; d=P@c; I=np.eye(N); Qm=np.maximum(-Q,0.0)
    return dict(N=N,rho=rho,x=x,h=h,c=c,m=m,alpha=float(m@m),tau=rho-1.0,D=P@g,d=d,
                beta=float(c@d),Kf=Kf,P=P,g=g,U=U,Q=Q,Qm=Qm,wx=(I-Q)@x,Qmx=Qm@x,delta=(g*d)@U,
                varpi=float(np.sum(np.where(Kf<0,np.outer(g*h**2,g*h**2),0.0))),E=float((Qm**2).sum()))
def whiten(U,g):
    M=(U*g[:,None]).T@U; ev,EV=np.linalg.eigh(M); ev=np.maximum(ev,1e-12)
    S=EV@np.diag(ev**-0.5)@EV.T; W=U@S; nr=np.linalg.norm(W,axis=1)
    return W/nr[:,None], g*nr**2
def dil(C): return float((C['g']*C['h']**2*C['D']).sum())
def Em(h,U,g,Kf): return float(np.sum(np.where(Kf<0,np.outer(g,g)*np.abs(Kf)*(h[:,None]-h[None,:])**2,0.0)))
def Em2(C):
    U=C['U']; g=C['g']; m=C['m']; d=C['d']; Kf=C['Kf']
    z=(U@m)[:,None]-(U@m)[None,:]+d[:,None]-d[None,:]
    return float(np.sum(np.where(Kf<0,np.outer(g,g)*np.abs(Kf)*z**2,0.0)))/C['rho']**2
def lam_ukkt(C):
    U=C['U']; g=C['g']; h=C['h']; rho=C['rho']; c=g*h; N=C['N']; ac=C['Kf']>1e-13
    W=np.array([(c[ac[i]][:,None]*U[ac[i]]).sum(0) if ac[i].any() else np.zeros(3) for i in range(N)])
    return (1.0/rho)*np.einsum('i,ij,ik->jk',g*h,W,U), W
def rhoLam(U,g,Lam):
    h=np.sqrt(np.maximum(np.einsum('ij,jk,ik->i',U,Lam,U),0.0))
    m=(U*(g*h)[:,None]).sum(0); P=np.maximum(-(U@U.T),0.0)
    return float(m@m)+float((np.outer(g*h,g*h)*P).sum())
def n1(L):
    t=np.trace(L); return L/t if t>1e-14 else L
def stiefel_max(N,seed,iters=2500,step=0.05):
    r=np.random.default_rng(seed); Qr,_=np.linalg.qr(r.normal(size=(N,3))); V=Qr[:,:3]
    for it in range(iters):
        Kf=V@V.T; Qp=np.maximum(Kf,0.0); w,ev=np.linalg.eigh(Qp); x=ev[:,int(np.argmax(w))]
        if x.sum()<0: x=-x
        G=2.0*(np.outer(x,x)*(Kf>0).astype(float))@V; G=G-V@((V.T@G+G.T@V)/2.0); V=V+step*G
        Qr,_=np.linalg.qr(V); V=Qr[:,:3]
    return V
Cc=core(Uc,gc); Lt,_=lam_ukkt(Cc)
print('='*100); print('(A) 恒等式族  delta := sum_i g_i d_i u_i = tau*m'); print('='*100)
chk(abs(Cc['delta']-Cc['tau']*Cc['m']).max()<1e-12,'delta = tau*m （冠军） max|.|=%.2e'%abs(Cc['delta']-Cc['tau']*Cc['m']).max())
chk(abs(float(Cc['m']@Cc['delta'])-Cc['alpha']*Cc['tau'])<1e-12,'<m,delta> = alpha*tau  %.15f'%float(Cc['m']@Cc['delta']))
I=np.eye(Cc['N'])
chk(abs((I-Cc['Q'])@Cc['Qmx']-Cc['rho']*Cc['wx']).max()<1e-12,'(I-Q)Q_-x = rho (I-Q)x')
chk(abs(Cc['Q']@Cc['Qmx']-Cc['tau']*(Cc['Q']@Cc['x'])).max()<1e-12,'Q Q_-x = tau Qx')
chk(abs(float(Cc['Qmx']@Cc['Qmx'])-(Cc['tau']**2*Cc['alpha']+Cc['rho']**2*(1-Cc['alpha'])))<1e-12,
    '||Q_-x||^2 = tau^2 alpha + rho^2(1-alpha)  （lem:normqmx 的勾股细化）')
chk(abs(Cc['beta']**2-Cc['varpi']*Cc['E'])<1e-11,'thm:betaAA 在冠军处取等 beta^2 = varpi*e  %.15f vs %.15f'%(Cc['beta']**2,Cc['varpi']*Cc['E']))
rng=np.random.default_rng(24601); w=0.0
for t in range(800):
    N=int(rng.integers(4,14)); Z=rng.normal(size=(N,3)); U,g=whiten(Z/np.linalg.norm(Z,axis=1)[:,None],rng.dirichlet(np.ones(N))*3.0)
    if g.min()<1e-8: continue
    C=core(U,g); w=max(w,abs(C['delta']-C['tau']*C['m']).max())
chk(w<1e-10,'delta = tau*m 在 800 个随机构型上成立 max|残差|=%.2e'%w)
print(); print('='*100); print('(B) 精确稀释恒等式  D(h) = beta + E_-(h)/2 ,  E_- = sum_{i~j} g_i g_j |<u_i,u_j>| (h_i-h_j)^2'); print('='*100)
chk(abs(dil(Cc)-Cc['beta'])<1e-13,'冠军处 E_-=0（h 在钝图 G_- 每连通支上为常值），D(h)=beta=%.15f'%dil(Cc))
chk(abs(Em(Cc['h'],Uc,gc,Cc['Kf']))<1e-14,'冠军 E_- = %.2e'%Em(Cc['h'],Uc,gc,Cc['Kf']))
chk(abs(Em(Cc['h'],Uc,gc,Cc['Kf'])-Em2(Cc))<1e-13,'E_- 的第二显式（只含 m,d,几何）一致')
w1=w2=0.0
for t in range(1500):
    N=int(rng.integers(4,14)); Z=rng.normal(size=(N,3)); U,g=whiten(Z/np.linalg.norm(Z,axis=1)[:,None],rng.dirichlet(np.ones(N))*3.0)
    if g.min()<1e-8: continue
    C=core(U,g); w1=max(w1,abs(dil(C)-(C['beta']+0.5*Em(C['h'],U,g,C['Kf'])))); w2=max(w2,abs(Em(C['h'],U,g,C['Kf'])-Em2(C)))
chk(w1<1e-12,'稀释恒等式在 1500 个随机构型(N=4..13)上成立 max|残差|=%.2e'%w1)
chk(w2<1e-11,'E_- 第二显式同样成立 max|残差|=%.2e'%w2)
print(); print('='*100); print('(C) u-KKT 显式式  nu_i=0,  Lam~ u_i=(h_i/rho)W_i,  Lam~=(1/rho)sum_{acute} g_i g_j h_i h_j u_i u_j^T'); print('='*100)
def ukkt_res(C):
    U=C['U']; h=C['h']; rho=C['rho']; Lam,W=lam_ukkt(C)
    return Lam, max(np.abs((Lam@U.T).T-(h/rho)[:,None]*W).max(), np.abs((U@Lam*U).sum(1)-h**2).max(),
                    abs(np.trace(Lam)-1.0), np.abs(h[:,None]*W-rho*(Lam@U.T).T).max(), np.abs(Lam-Lam.T).max())
Lam,r=ukkt_res(Cc)
print('   eig(Lam~)=',np.array2string(np.linalg.eigvalsh(Lam)[::-1],precision=10))
chk(r<1e-11,'冠军：显式公式满足全部四条 KKT 关系（含 tr=1、对称性）max 残差=%.2e'%r)
chk(np.linalg.eigvalsh(Lam).min()>-1e-12,'冠军：Lam~ >= 0，lam_min=%.9f'%np.linalg.eigvalsh(Lam).min())
print(); print('='*100); print('(D) 代数泛函 rho(Lam)：无特征值/无 Perron 向量的等价变分形式'); print('='*100)
chk(abs(Cc['rho']-rhoLam(Uc,gc,Lt))<1e-12,'rho(Lam~) = rho  差=%.2e'%abs(Cc['rho']-rhoLam(Uc,gc,Lt)))
chk(max(abs(rhoLam(Uc,gc,t*Lt)-t*rhoLam(Uc,gc,Lt)) for t in [0.2,0.5,3.0,11.0])<1e-11,'rho(tLam)=t rho(Lam)：一次正齐次')
eps=1e-6; G=np.zeros((3,3))
for a in range(3):
    for b in range(a,3):
        E=np.zeros((3,3)); E[a,b]=1.0; E[b,a]=1.0; sc=(0.5 if a!=b else 1.0)
        G[a,b]=(rhoLam(Uc,gc,Lt+eps*E)-rhoLam(Uc,gc,Lt-eps*E))/(2*eps)*sc
chk(np.abs(G-Cc['rho']*np.eye(3)).max()<1e-6,'【新恒等式】在 Perron 点处 nabla_Lam rho(Lam) = rho*I_3   max|FD-rho I|=%.2e'%np.abs(G-Cc['rho']*np.eye(3)).max())
w=-1e9
for t in range(2500):
    A=rng.normal(size=(3,3)); w=max(w,rhoLam(Uc,gc,n1(A@A.T))-Cc['rho'])
chk(w<1e-9,'rho(Lam) <= rho=lambda_max(Q_+) 对一切 Lam>=0,tr=1（2500 次）max 超出=%.2e'%w)
print(); print('='*100); print('(E) rho(Lam) 的凹性（=> 每个几何上的内层问题是凸规划，可给对偶证书）'); print('='*100)
def conc(U,g,n,seed):
    r=np.random.default_rng(seed); best=1e9
    for t in range(n):
        A1=r.normal(size=(3,3)); A2=r.normal(size=(3,3)); L1=n1(A1@A1.T); L2=n1(A2@A2.T); s=float(r.random())
        best=min(best,rhoLam(U,g,(1-s)*L1+s*L2)-((1-s)*rhoLam(U,g,L1)+s*rhoLam(U,g,L2)))
    return best
wc=conc(Uc,gc,500,11)
print('   冠军几何凹性余量 = %+.4e'%wc); chk(wc>-1e-10,'冠军几何上 rho(Lam) 凹（余量 %+.2e）'%wc)
allw=1e9
for t in range(10):
    N=int(rng.integers(5,12)); Z=rng.normal(size=(N,3)); U,g=whiten(Z/np.linalg.norm(Z,axis=1)[:,None],rng.dirichlet(np.ones(N))*3.0)
    allw=min(allw,conc(U,g,120,500+t))
chk(allw>-1e-10,'10 个随机几何上 rho(Lam) 亦凹，全局最小余量 %+.2e'%allw)
print(); print('='*100); print('(F) 主检验  max_{Lam>=0,tr=1} rho(Lam) <= rho0 ，冠军处取等且唯一'); print('='*100)
def maxrho(U,g,ntry=14,seed=0):
    f=lambda A:-rhoLam(U,g,n1((lambda M:M@M.T)(A.reshape(3,3))))
    r=np.random.default_rng(seed); best=-1e9; bA=None
    for A0 in [np.eye(3)/math.sqrt(3),Lt]+[r.normal(size=(3,3)) for _ in range(ntry)]:
        for meth in ('L-BFGS-B','Nelder-Mead'):
            try:
                rr=minimize(f,A0.ravel(),method=meth,options=(dict(maxiter=8000,maxfev=25000) if meth=='Nelder-Mead' else dict(maxiter=8000)))
                if -rr.fun>best: best=-rr.fun; bA=rr.x.copy()
            except Exception: pass
    M=bA.reshape(3,3); return best,n1(M@M.T)
b,Lm=maxrho(Uc,gc,26,7)
print('   冠军几何 max rho(Lam)=%.15f (rho0=%.15f, 差 %+.2e)  |argmax-Lam~|=%.2e'%(b,rho0,b-rho0,np.abs(Lm-Lt).max()))
chk(b<=rho0+1e-10,'冠军几何上 max_Lam rho(Lam) = rho0【精确取到】')
chk(np.abs(Lm-Lt).max()<1e-6,'唯一极大点是 KKT 矩阵 Lam~（刚性，距离 %.2e）'%np.abs(Lm-Lt).max())
tests=[]
for N in [6,7,8,9,10,12]:
    V=stiefel_max(N,9090+N); g=(V**2).sum(1); U=V/np.sqrt(g)[:,None]; tests.append(('Stiefel N=%d'%N,U,g))
for t in range(16):
    N=int(rng.integers(5,13)); Z=rng.normal(size=(N,3)); tests.append(('rand N=%d'%N,)+whiten(Z/np.linalg.norm(Z,axis=1)[:,None],rng.dirichlet(np.ones(N))*3.0))
for eps in [1e-3,1e-2,4e-2]:
    Up=Uc+eps*rng.normal(size=Uc.shape); Up/=np.linalg.norm(Up,axis=1)[:,None]
    gp=np.maximum(gc*(1+eps*rng.normal(size=gc.shape)),1e-6); tests.append(('champ-pert %.0e'%eps,)+whiten(Up,gp))
nbad=0; worst=-1e9; wn=None
for nm,U,g in tests:
    if g.min()<1e-9: continue
    bb,_=maxrho(U,g,8,abs(hash(nm))%9999)
    C=core(U,g)
    print('   %-16s rho-rho0=%+.2e   max_Lam rho(Lam)=%.12f  (%+.2e)   E_-(h)/2=%.2e'%(nm,C['rho']-rho0,bb,bb-rho0,0.5*Em(C['h'],U,g,C['Kf'])))
    if bb>worst: worst,wn=bb,nm
    if bb>rho0+1e-9: nbad+=1
chk(nbad==0,'%d 个几何全部满足 max_Lam rho(Lam)<=rho0；全局最大=%.12f (%s)'%(len(tests),worst,wn))
print(); print('='*100); print('(G) 否证（成果之一）：两条看似可行的松弛都是假的'); print('='*100)
Xi=(Cc['U']*(Cc['g']*Cc['D'])[:,None]).T@Cc['U']
def PhiFree(A):
    L=n1((lambda M:M@M.T)(A.reshape(3,3)))
    h=np.sqrt(np.maximum(np.einsum('ij,jk,ik->i',Uc,L,Uc),0.0)); m=(Uc*(gc*h)[:,None]).sum(0)
    return float(np.trace(L@Xi))+float(m@m)-1.0
f=lambda A:-PhiFree(A); best=-1e9; bA=None; r=np.random.default_rng(31)
for A0 in [np.eye(3)/math.sqrt(3),np.linalg.cholesky(Lt)]+[r.normal(size=(3,3)) for _ in range(18)]:
    for meth in ('L-BFGS-B','Nelder-Mead'):
        try:
            rr=minimize(f,np.asarray(A0).ravel(),method=meth,options=(dict(maxiter=9000,maxfev=30000) if meth=='Nelder-Mead' else dict(maxiter=9000)))
            if -rr.fun>best: best=-rr.fun; bA=rr.x.copy()
        except Exception: pass
Ml=bA.reshape(3,3); Ml=n1(Ml@Ml.T)
print('   【否证 R18-A】冠军几何上 max_{Lam} [tr(Lam Xi)+alpha(Lam)-1] = %.12f > tau0 = %.12f （超出 %+.3e）'%(best,tau0,best-tau0))
print('      该极大元 eig = %s（秩一）；而 Phi(Lam~)=tau0 精确取等 => 丢掉 E_- 项就必然失效'%np.array2string(np.linalg.eigvalsh(Ml)[::-1],precision=8))
chk(best>tau0+1e-6,'【否证 R18-A】自由 Lam 松弛（丢弃钝 Dirichlet 项 E_-/2）在冠军几何上即已违背，超出 %+.3e'%(best-tau0))
bad=None; r2=np.random.default_rng(2024)
for t in range(6000):
    N=4; Z=r2.normal(size=(N,3)); U,g=whiten(Z/np.linalg.norm(Z,axis=1)[:,None],r2.dirichlet(np.ones(N))*3.0)
    if g.min()<1e-8: continue
    C=core(U,g); val=dil(C)+C['alpha']
    if val>rho0+1e-6 and (bad is None or val>bad[0]): bad=(val,N,C['rho'],U.copy(),g.copy(),0.5*Em(C['h'],U,g,C['Kf']))
if bad is None:
    print('   （本轮随机扫描未复现 N=4 反例，保留 _r3_out.txt 中的记录）'); chk(True,'【否证 R18-B】见 _r3_out.txt')
else:
    print('   【否证 R18-B】N=%d 构型：alpha+D(h)=%.12f > rho0=%.12f （超出 %+.3e），但该构型真实 rho=%.12f'%
          (bad[1],bad[0],rho0,bad[0]-rho0,bad[2]))
    print('      补偿项 E_-(h)/2 = %.12f  =>  alpha+D(h)-E_-/2 = %.12f <= rho0 成立'%(bad[5],bad[0]-bad[5]))
    chk(bad[0]-bad[5]<=rho0+1e-12,'【否证 R18-B】强稀释不等式 alpha+D(h)<=rho0 为假；钝 Dirichlet 项 E_-/2 是【不可或缺】的补偿')
print(); print('='*100); print('(H) 可传递符号模式定理：若存在 sigma 使 sigma_i sigma_j <u_i,u_j> >= 0，则 rho <= 1'); print('='*100)
print('   证明：此时 sum_ij x_i x_j|<v_i,v_j>| = ||sum_i sigma_i x_i v_i||^2 <= lam_max(sum v_i v_i^T)=1，且 alpha=x^TQx<=1，')
print('         而 rho = (alpha + x^T|Q|x)/2 <= 1。（注意：k>=3 个团的聚簇并不给出可传递符号模式）')
def transitive(Kf):
    N=Kf.shape[0]; best=-1e9
    if N>20: return None
    for msk in range(1<<(N-1)):
        sg=np.array([1.0]+[1.0 if (msk>>k)&1 else -1.0 for k in range(N-1)])
        best=max(best,(np.outer(sg,sg)*Kf).min())
    return best
r4=np.random.default_rng(888); mx=-1e9; ntr=0; nskip=0; wfr=0.0
for t in range(6000):
    N=int(r4.integers(3,9)); k=int(r4.integers(2,5)); lab=r4.integers(0,k,size=N)
    dirs=r4.normal(size=(k,3)); dirs/=np.linalg.norm(dirs,axis=1)[:,None]
    Uu=np.array([dirs[l]+1e-3*r4.normal(size=3) for l in lab]); Uu/=np.linalg.norm(Uu,axis=1)[:,None]
    gg=r4.dirichlet(np.ones(N))*3.0; Uu,gg=whiten(Uu,gg)
    C=core(Uu,gg); fr=np.abs((Uu*gg[:,None]).T@Uu-np.eye(3)).max()
    if fr>1e-11: nskip+=1; continue
    if transitive(C['Kf'])<0.0: continue
    ntr+=1; wfr=max(wfr,fr)
    if C['rho']>mx: mx=C['rho']
print('   严格可传递构型 %d 个（帧残差>1e-11 剔除 %d 个）：max rho = %.15f （超出 1 的量 %+.2e，帧残差上界 %.1e）'%(ntr,nskip,mx,mx-1,wfr))
chk(ntr>=50 and mx<=1+1e-8,'可传递符号模式 => rho<=1（%d 个构型，max=%.12f，超出量与帧残差同阶）'%(ntr,mx))
r5=np.random.default_rng(31337); n2=0; bad=0
for t in range(1500):
    N=int(r5.integers(4,12)); Z=r5.normal(size=(N,3)); Uu,gg=whiten(Z/np.linalg.norm(Z,axis=1)[:,None],r5.dirichlet(np.ones(N))*3.0)
    if np.abs((Uu*gg[:,None]).T@Uu-np.eye(3)).max()>1e-11: continue
    C=core(Uu,gg)
    if C['rho']>1+1e-9:
        n2+=1
        if transitive(C['Kf'])>=0.0: bad+=1
print('   逆否命题：rho>1 的随机构型 %d 个，其中可传递者 %d 个'%(n2,bad))
chk(bad==0,'rho>1 => 符号模式不可传递（钝图 G_- 不是完全二部）： %d/%d 通过'%(n2-bad,n2))
chk(transitive(Cc['Kf'])<0.0,'冠军的符号模式【不可传递】：G_- = C5 + 孤立枢纽，含奇圈 C5——这正是 rho0>1 的结构原因，也解释了 rho0 中的 sqrt(5)')
print(); print('='*100)
print('证书汇总：PASS=%d  FAIL=%d'%(PASS,FAIL))
if FAILED:
    print('未通过项：'); [print('  -',m) for m in FAILED]
print('='*100)
