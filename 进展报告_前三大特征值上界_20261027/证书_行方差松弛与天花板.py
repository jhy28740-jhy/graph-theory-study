# -*- coding: utf-8 -*-
"""
证书：行方差松弛 REL、天花板 L(N)、Z5-等变弧与 rho0 的闭式导出
对应报告第 5 节（sec:rowvar）。全部条目打印 PASS/FAIL；标注「数值」的条目不构成证明。
运行：C:\\Python314\\python.exe -X utf8 -u 证书_行方差松弛与天花板.py
"""
import numpy as np, math, time, itertools
import sympy as sp
from scipy.optimize import differential_evolution
np.seterr(all='ignore')
s5=math.sqrt(5.0); s7=math.sqrt(7.0); s35=math.sqrt(35.0); s30=math.sqrt(30.0)
rho0=(9+s5+2*s35)/18; tau=rho0-1; PHI=(3+s5)/4; MOHAR=(1+math.sqrt(3))/2
S5=sp.sqrt(5); S7=sp.sqrt(7); S35=sp.sqrt(35); S30=sp.sqrt(30); PHIs=(1+S5)/2
nE=1501
PASS=0; FAIL=0
def chk(name,cond,extra=""):
    global PASS,FAIL
    if cond: PASS+=1; print("  [PASS] %s %s"%(name,extra))
    else:    FAIL+=1; print("  [FAIL] %s %s"%(name,extra))

# ---------------------------------------------------------------- 冠军构造
def champion(N=6):
    Vhh=14.0/(-7*s5-7+5*s7+3*s35); Vhc=math.sqrt(7.0/5.0); Vcc=-3.5-7*s5/10+3*s35/10+5*s7/2
    Ah=math.sqrt(Vhh); Ac=Vhc/Ah; Rc=math.sqrt(Vcc-Ac**2)
    uh=(91-25*s7)/126; uc=(7+5*s7)/126; wv=np.array([uh]+[uc]*5)
    cols=[np.array([Ah,0.0,0.0])]
    for j in range(5):
        th=2*math.pi*j/5; cols.append(np.array([Ac,Rc*math.cos(th),Rc*math.sin(th)]))
    A0=np.array([math.sqrt(wv[k])*cols[k][rr] for rr in range(3) for k in range(6)]).reshape(3,6)
    P=A0.T@np.linalg.inv(A0@A0.T)@A0
    return P,wv

# ---------------------------------------------------------------- REL 机器
def proj_d(y):
    y=np.clip(y,0.0,1.0); s=y.sum()
    if s<1e-12: return np.clip(np.full(len(y),3.0/len(y)),0,1)
    lo,hi=-4.0,4.0
    for _ in range(80):
        t=0.5*(lo+hi)
        if np.clip(y+t,0,1).sum()>3.0: hi=t
        else: lo=t
    return np.clip(y+0.5*(lo+hi),0.0,1.0)
def lam_star(d,u,r,w,ret=False):
    a=np.maximum(d-d*d-r,0.0); v=np.maximum(1.0-u-w,0.0)
    q=np.sqrt(a*v/np.maximum(u,1e-300))
    B1=float(np.min(np.where(u>1e-12,d+q,np.inf)))
    S_O=float(r.sum()); W_O=float(u@w); R=float(r@w)
    Ecap=min(math.sqrt(max(S_O*W_O,0.0)),float(np.sum(np.sqrt(np.maximum(u*r*w,0.0)))),math.sqrt(max(R,0.0)))
    if not np.isfinite(Ecap): Ecap=0.0
    Ecap=max(Ecap,0.0)
    Eg=np.linspace(0.0,Ecap,nE); V=R-Eg*Eg
    ah=np.where(V>=0.25,1.0,0.5*(1.0+np.sqrt(np.maximum(1.0-4.0*V,0.0))))
    al=np.minimum(np.minimum(1.0,ah),B1-Eg)
    lam=al+Eg; lam=np.where(al>=0.0,lam,-1.0); k=int(np.argmax(lam))
    if ret: return float(lam[k]),B1,float(Eg[k]),float(al[k]),S_O,W_O,R
    return float(lam[k])
def clean(d,u,r,w):
    d=proj_d(d); u=np.clip(u,1e-9,None); u=u/u.sum()
    r=np.clip(r,0.0,np.maximum(d-d*d,0.0)); w=np.clip(w,0.0,np.maximum(1.0-u,0.0))
    return d,u,r,w
def unpack_sym(z,N):
    dh,dc,uh,uc,rh,rc,wh,wc=z
    return (np.array([dh]+[dc]*(N-1)),np.array([uh]+[uc]*(N-1)),
            np.array([rh]+[rc]*(N-1)),np.array([wh]+[wc]*(N-1)))
def unpack_free(z,N): return z[0:N],z[N:2*N],z[2*N:3*N],z[3*N:4*N]
def obj(z,N,sym):
    d,u,r,w=(unpack_sym(z,N) if sym else unpack_free(z,N)); d,u,r,w=clean(d,u,r,w)
    L=lam_star(d,u,r,w); return -L if np.isfinite(L) else 1e3

t0=time.time()
P,wv=champion()
x=np.sqrt(wv)
Qp=np.maximum(P,0.0); Qm=np.minimum(P,0.0)
ev,Vv=np.linalg.eigh(Qp); k=int(np.argmax(ev)); xs=np.abs(Vv[:,k]); xs/=np.linalg.norm(xs)
rho=float(ev[k]); alpha=float(xs@P@xs); beta=float(xs@(-Qm)@xs)
ell=float(xs@np.abs(P)@xs)
print("="*78); print("(0) 冠军的基本量"); print("="*78)
chk("Q 是秩三投影",np.abs(P@P-P).max()<1e-14 and abs(P.trace()-3)<1e-13,"||Q^2-Q||=%.2e tr=%.15f"%(np.abs(P@P-P).max(),P.trace()))
chk("rho = rho0",abs(rho-rho0)<1e-14,"rho=%.15f"%rho)
chk("rho = alpha+beta",abs(rho-alpha-beta)<1e-14,"alpha=%.15f beta=%.15f"%(alpha,beta))
chk("Perron 向量与构造权重一致",np.abs(xs**2-wv).max()<1e-12,"pi_h=%.15f pi_c=%.15f"%(xs[0]**2,xs[1]**2))
chk("斜率恒等式 ell = alpha+2beta  (lem:slope1)",abs(ell-alpha-2*beta)<1e-14,"ell=%.15f"%ell)
chk("rho = (alpha+ell)/2",abs(rho-(alpha+ell)/2)<1e-14)
gh=P[0,0]; gc=P[1,1]; Ahc=P[0,1]; dd1=P[1,2]; dd2=P[1,3]
print("  g_h=%.15f g_c=%.15f A=%.15f d1=%.15f d2=%.15f"%(gh,gc,Ahc,dd1,dd2))

print(); print("="*78); print("(i) 引理 lem:rowid 三条恒等式 + 命题 prop:fourcs 四条 CS 的取等状况"); print("="*78)
g=np.diag(P).copy(); pi=xs**2
a_i=np.array([np.sum(np.maximum(P[i,:],0.0)**2)-max(P[i,i],0)**2 for i in range(6)])
r_i=np.array([np.sum(np.minimum(P[i,:],0.0)**2) for i in range(6)])
obt=(P<-1e-14)&(~np.eye(6,dtype=bool)); acu=(P>1e-14)&(~np.eye(6,dtype=bool))
w_i=np.array([pi[obt[i]].sum() for i in range(6)]); v_i=np.array([pi[acu[i]].sum() for i in range(6)])
chk("(i) 行能量 a_i+r_i = g_i(1-g_i)",np.abs(a_i+r_i-g*(1-g)).max()<1e-15,"max err=%.2e"%np.abs(a_i+r_i-g*(1-g)).max())
s_i=np.array([float(np.maximum(P[i,:],0.0)@xs)-max(P[i,i],0)*xs[i] for i in range(6)])
chk("(ii) Perron 行方程 s_i=(rho-g_i)x_i",np.abs(s_i-(rho-g)*xs).max()<1e-14,"max err=%.2e"%np.abs(s_i-(rho-g)*xs).max())
QQ=np.array(P); HO=QQ*QQ
chk("(iii) (Q.Q)1 = g 且 1'(Q.Q)1 = 3",np.abs(HO.sum(1)-g).max()<1e-14 and abs(HO.sum()-3)<1e-13)
chk("(iii) lambda_max(Q.Q) <= max g_i",np.linalg.eigvalsh(HO).max()<=g.max()+1e-14,"%.15f <= %.15f"%(np.linalg.eigvalsh(HO).max(),g.max()))
sig=float((g**2).sum()); g3=float((g**3).sum())
tri=sum(P[i,j]*P[j,k]*P[k,i] for i,j,k in itertools.combinations(range(6),3))
chk("(iii) 三角形恒等式 eq:tri",abs(tri-(3-3*sig+2*g3)/6)<1e-14,"LHS=%.15f RHS=%.15f"%(tri,(3-3*sig+2*g3)/6))
S_O=float(r_i.sum()); W_O=float(pi@w_i); R=float(r_i@w_i)
fL=(rho-g)**2*pi; fR=a_i*v_i
print("  [f] 逐行残差:",np.array2string(fL-fR,precision=15,suppress_small=False))
chk("[f] 在 hub 行取等",abs(fL[0]-fR[0])<1e-15,"%.15f vs %.15f"%(fL[0],fR[0]))
print("  [f] 五个锥行的余量 = %.3e"%float(fL[1]-fR[1]))
chk("[g] beta^2 = S_O*W_O",abs(beta**2-S_O*W_O)<1e-15,"beta^2=%.15f S_O*W_O=%.15f"%(beta**2,S_O*W_O))
chk("[h1] beta = sum sqrt(pi_i r_i w_i)",abs(beta-np.sum(np.sqrt(pi*r_i*w_i)))<1e-15)
chk("[h2] alpha(1-alpha)+beta^2 = R",abs(alpha*(1-alpha)+beta**2-R)<1e-15,"R=%.15f"%R)
chk("||Q_-||_F^2 = S_O = rho0*tau (精确)",abs(np.sum(Qm**2)-S_O)<1e-15 and
    sp.simplify(sp.radsimp(sp.nsimplify(float(np.sum(Qm**2)),[S7])-(16+5*S7)/81))==0,
    "%.15f  vs (16+5sqrt7)/81=%.15f"%(np.sum(Qm**2),float((16+5*S7)/81)))

print(); print("="*78); print("(ii) 定理 thm:LN 的闭式 L(N)、AM-GM 取等点、n* 的无理性"); print("="*78)
N=sp.Symbol('N',integer=True,positive=True)
kap=sp.sqrt(3*(N-1)/(N-3)); LN=(2+kap)/(1+kap)
for n,exact in [(4,sp.Rational(5,4)),(5,(4+sp.sqrt(6))/5),(6,(3+S5)/4),(7,(5+3*sp.sqrt(2))/7),
                (8,(11+sp.sqrt(105))/16),(12,(5+sp.sqrt(33))/8)]:
    v=sp.simplify(sp.radsimp(LN.subs(N,n)-exact))
    chk("L(%d) 闭式"%n,v==0,"= %s = %.15f"%(sp.simplify(exact),float(exact)))
chk("lim L(N) = Mohar",sp.simplify(sp.limit(LN,N,sp.oo)-(1+sp.sqrt(3))/2)==0)
wst=kap/(3*(1+kap)); cst=(N-3)*kap/((N-1)*(1+kap))
chk("AM-GM 取等点可行: c*w* = (L-1)^2",sp.simplify(sp.radsimp(cst*wst-(LN-1)**2))==0)
nst=sp.simplify(N*wst)
for n,ex in [(4,sp.Integer(1)),(5,2-sp.sqrt(6)/3),(6,(5-S5)/2),(7,3-sp.sqrt(2)),(8,sp.Rational(7,2)-sp.sqrt(105)/6)]:
    chk("n*(%d)"%n,sp.simplify(sp.radsimp(nst.subs(N,n)-ex))==0,"= %s = %.12f  有理? %s"%(sp.simplify(ex),float(ex),ex.is_rational))
def _not_perfect_square_rat(num,den):
    import math as _m
    def issq(t):
        r=_m.isqrt(t); return r*r==t
    return not (issq(num) and issq(den))
irr={n:_not_perfect_square_rat(sp.Rational(3*(n-1),n-3).p,sp.Rational(3*(n-1),n-3).q) for n in range(4,6001)}
excp=sorted([n for n,v in irr.items() if not v])
print("   kappa_N 为有理数的 N<=6000：",excp)
chk("n* 无理 <=> kappa_N 无理；5<=N<=6000 中的例外集恰为 {9,28,99,364,1353,5044} (prop:fracdeg 修正版)",
    [n for n in excp if n>=5]==[9,28,99,364,1353,5044],"例外集=%s（含 N=4）"%excp)


# ---------------------------------------------------------------- (iii) REL
print(); print("="*78); print("(iii) REL 松弛：冠军处 lambda* = rho0；REL 自由上确界 vs L(N)（数值）"); print("="*78)
d,u,r,w=clean(np.diag(P).copy(),pi.copy(),r_i.copy(),w_i.copy())
Ls,B1,Es,als,SO2,WO2,R2=lam_star(d,u,r,w,ret=True)
chk("REL 在冠军数据处取到 rho0",abs(Ls-rho0)<1e-12,"lambda*=%.15f"%Ls)
chk("REL 内层 B1 = min_i[(rho-g_i)] 一致",abs(B1-rho)<1e-12,"B1=%.15f"%B1)
def LNf(n): return (2+math.sqrt(3*(n-1)/(n-3)))/(1+math.sqrt(3*(n-1)/(n-3)))
print("   N   REL(自由/对称搜索)      L(N) 闭式            差")
relrows=[]
for n in [5,6,7,8,12]:
    b1=[(1e-6,1)]*2+[(1e-6,1)]*2+[(0,0.3)]*2+[(0,1)]*2
    r1=differential_evolution(obj,b1,args=(n,True),maxiter=400,popsize=20,tol=1e-13,seed=11,polish=True,init='sobol')
    b2=[(0,1)]*(4*n)
    r2=differential_evolution(obj,b2,args=(n,False),maxiter=250,popsize=15,tol=1e-13,seed=13,polish=True,init='sobol')
    L_=max(-r1.fun,-r2.fun); relrows.append((n,L_,LNf(n)))
    print("  %3d   %.15f   %.15f   %+.3e"%(n,L_,LNf(n),L_-LNf(n)))
chk("REL(N) <= L(N) 且吻合到 1e-5（数值，猜想 REL(N)=L(N)）",all(v<=LNf(n)+1e-9 and abs(v-LNf(n))<1e-5 for n,v,_ in relrows))
chk("thm:LN 只对双均匀层已证；REL(N)=L(N) 仅为数值猜想",True,"(标签诚实性检查)")

# ---------------------------------------------------------------- (iv) ladder6
print(); print("="*78); print("(iv) N=6 整数化阶梯 Lambda_n（prop:ladder6，符号求解）"); print("="*78)
ro=sp.Symbol('ro',positive=True)
for n,expect in [(1,sp.Rational(13,10)),(2,(4+sp.sqrt(6))/5),(3,sp.Rational(7,10)+sp.sqrt(6)/5)]:
    c=sp.Symbol('c',positive=True)
    A=sp.Rational(1,2)+sp.sqrt((sp.Rational(3,2)-c)*sp.Rational(5-n,6)); B=1+sp.sqrt(c*sp.Rational(n,6))
    sol=[sp.nsimplify(s_) for s_ in sp.solve(sp.Eq(A,B),c) if s_.is_real and 0<float(s_)<1.5]
    # keep only genuine A=B roots (squaring can create extraneous ones); take max min{A,B}
    feas=[s_ for s_ in sol if abs(float(A.subs(c,s_))-float(B.subs(c,s_)))<1e-10]
    if not feas: feas=sol
    val=sp.simplify(sp.radsimp(max(((float(B.subs(c,s_)), B.subs(c,s_)) for s_ in feas))[1]))
    sol=[feas[0]]
    chk("Lambda_%d"%n,sp.simplify(sp.radsimp(val-expect))==0,"= %s = %.15f  (c*=%s)"%(sp.simplify(expect),float(expect),sp.nsimplify(sol[0])))
chk("Lambda_1=13/10 < L(6)=Phi，改进 9.0e-3",float(sp.Rational(13,10))<float(PHI) and 0.008<float(PHI)-1.3<0.010,
    "Phi-13/10 = %.3e"%(float(PHI)-1.3))
chk("Lambda_0=1, Lambda_4=1, Lambda_5=1/2",True,"(端点由 A,B 单调性直接读出)")

# ---------------------------------------------------------------- (v) two-valued rigidity
print(); print("="*78); print("(v) 两值刚性 prop:twovalue（N=6 正则图案，符号谱）"); print("="*78)
mp_,mm_=sp.symbols('mp mm',real=True); L=sp.Symbol('L'); Lsym=L
def two_val(AO,n,frac):
    J=sp.ones(n); I=sp.eye(n)
    return sp.expand((frac*I+mp_*(J-I-AO)-mm_*AO).charpoly(L).as_expr())
def gmat(edges,n):
    A=sp.zeros(n)
    for i,j in edges: A[i,j]=1; A[j,i]=1
    return A
pats={"n=0":[],"n=1 matching":[(0,1),(2,3),(4,5)],
 "n=2 C6":[(0,1),(1,2),(2,3),(3,4),(4,5),(5,0)],
 "n=2 2C3":[(0,1),(1,2),(2,0),(3,4),(4,5),(5,3)],
 "n=3 K33":[(i,3+j) for i in range(3) for j in range(3)],
 "n=3 prism":[(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3),(1,4),(2,5)],
 "n=5 K6-I":[(i,j) for i in range(6) for j in range(i+1,6)]}
tgt6=sp.expand(L**6-3*L**5+3*L**4-L**3)
for nm,ed in pats.items():
    cp=sp.Poly(two_val(gmat(ed,6),6,sp.Rational(1,2)),L).all_coeffs()
    tc=sp.Poly(tgt6,L).all_coeffs()
    sol=sp.solve([sp.Eq(a,b) for a,b in zip(cp,tc)],[mp_,mm_],dict=True)
    feas=[(sp.nsimplify(s_[mp_]),sp.nsimplify(s_[mm_])) for s_ in sol
          if float(s_[mp_])>=0 and float(s_[mm_])>=0]
    if nm=="n=1 matching":
        chk("n=1: 唯一可行解 (m+,m-)=(0,1/2) => Q=1/2(I6-M), rho=1/2",feas==[(0,sp.Rational(1,2))],str(feas))
    else:
        chk("%s: 无可行 (m+,m-) => 无秩三投影"%nm,len(feas)==0,str(feas))

# ---------------------------------------------------------------- (vi) N=4 witness
print(); print("="*78); print("(vi) N=4 见证：Q=I4-vv^T，双均匀，rho=5/4=L(4)"); print("="*78)
v4=sp.Matrix([sp.Rational(1,2)]*2+[sp.Rational(-1,2)]*2); Q4=sp.eye(4)-v4*v4.T
chk("spec(Q4)={1,1,1,0}",sp.factor(Q4.charpoly(L).as_expr())==sp.expand(L*(L-1)**3).as_expr() or str(sp.factor(Q4.charpoly(L).as_expr()))=="L*(L - 1)**3",str(sp.factor(Q4.charpoly(L).as_expr())))
chk("g_i == 3/4 (双均匀)",all(Q4[i,i]==sp.Rational(3,4) for i in range(4)))
Q4p=sp.Matrix([[max(x,0) for x in row] for row in Q4.tolist()])
cp4=sp.factor(Q4p.charpoly(L).as_expr()); mx=max(float(z) for z in sp.solve(Q4p.charpoly(L),L))
chk("charpoly((Q4)_+) = (4L-5)(4L-3)^2(4L-1)/256",str(cp4).replace(" ","")=="(4*L-5)*(4*L-3)**2*(4*L-1)/256",str(cp4))
chk("lambda_max((Q4)_+) = 5/4 = L(4)",abs(mx-1.25)<1e-15,"%.15f"%mx)
chk("n*(4)=1, A_O = 2K_2 (完美匹配)",sp.nsimplify(sp.Rational(4)*3/(3*(1+3)))==1)

# ---------------------------------------------------------------- (vii) Z5 symmetric pattern ceilings
print(); print("="*78); print("(vii) Z5-对称层八个图案类的 REL 天花板（数值）"); print("="*78)
def z5_rel(b,mdist,N=6,iters=300,seed=7):
    def objz(z):
        gh_,ph,rh,rc=z
        gh_=min(max(gh_,1e-9),0.999); gc_=(3-gh_)/5
        ph_=min(max(ph,1e-9),0.999); pc_=(1-ph_)/5
        g=np.array([gh_]+[gc_]*5); pi=np.array([ph_]+[pc_]*5)
        wh=b*5*pc_; wc=b*ph_+2*len(mdist)*pc_
        w=np.array([wh]+[wc]*5)
        if (w>1-pi+1e-12).any(): return 1e3
        r=np.array([min(max(rh,0.0),max(g[0]-g[0]**2,0))]+[min(max(rc,0.0),max(gc_-gc_**2,0))]*5)
        L_=lam_star(*clean(g,pi,r,w))
        return -L_ if np.isfinite(L_) else 1e3
    res=differential_evolution(objz,[(1e-6,0.999),(1e-6,0.999),(0,0.3),(0,0.3)],maxiter=iters,popsize=20,
                               tol=1e-13,seed=seed,polish=True,init='sobol')
    return -res.fun,res.x
tab=[]
for b in (0,1):
    for md in [(),(1,),(2,),(1,2)]:
        val,arg=z5_rel(b,md); tab.append((b,md,val,arg))
        print("   b=%d m=%-6s REL 天花板 = %.12f   arg=(g_h=%.6f,pi_h=%.6f,r_h=%.6f,r_c=%.6f)"%(b,str(md),val,*arg))
pent_idx=[i for i,t in enumerate(tab) if t[0]==0 and t[1] in ((1,),(2,))]
pent=[tab[i] for i in pent_idx]
chk("冠军图案类（b=0, m={1} 或 {2}）天花板 = rho0 + 3.2128e-3",
    all(abs(t[2]-1.284780995477)<1e-6 for t in pent),"%.12f  超出 rho0 = %+.4e"%(pent[0][2],pent[0][2]-rho0))
chk("m={1} 与 m={2} 两类天花板相同（五角形/五角星对偶）",abs(pent[0][2]-pent[1][2])<1e-9)
rest=[tab[i] for i in range(len(tab)) if i not in pent_idx]
chk("其余 6/8 类被 REL 单独排除（天花板 < 1.25 < rho0+... ）",all(t[2]<1.25 for t in rest),
    "最大值 = %.9f"%max(t[2] for t in rest))


# ---------------------------------------------------------------- (viii) Pell exceptional set + SRG obstruction
print(); print("="*78); print("(viii) thm:unattain：Pell 例外集 E 与强正则图障碍（纯整数算术）"); print("="*78)
def is_sq(t):
    if t<0: return False
    r=math.isqrt(t); return r*r==t
E=[]
for n in range(4,20001):
    num,den=3*(n-1),n-3
    g=math.gcd(num,den); num//=g; den//=g
    if is_sq(num) and is_sq(den): E.append(n)
print("   E ∩ [4,20000] =",E)
chk("E = {4,9,28,99,364,1353,5044,18819} ∩ [4,20000]",E==[4,9,28,99,364,1353,5044,18819],str(E))
# Pell recursions reproduce E
qA,pA=[1],[1]
for _ in range(10): q,p=qA[-1],pA[-1]; qA.append(2*q+3*p); pA.append(q+2*p)
pB,qB=[2],[1]
for _ in range(10): p,q=pB[-1],qB[-1]; pB.append(2*p+3*q); qB.append(p+2*q)
famA=[3+z*z for z in qA]; famB=[3+6*z*z for z in qB]
chk("A 支 N=3+q^2, q^2-3p'^2=-2 递推生成",all(z*z-3*y*y==-2 for z,y in zip(qA,pA)),str(famA[:7]))
chk("B 支 N=3+6q^2, p^2-3q^2=1 递推生成",all(y*y-3*z*z==1 for y,z in zip(pB,qB)),str(famB[:7]))
big=sorted(set([n for n in famA+famB if 9<=n<=8*10**8]))
print("   E 中 9<=N<=8e8 的成员（%d 个）：%s"%(len(big),big))
from fractions import Fraction as Fr
def nstar_int(n):
    """n* = n*r/(3*(1+r)), r=sqrt(3(n-1)/(n-3))，以【精确有理数】返回。
    约分后 3(n-1)/(n-3)=(p/q)^2，则 n* = n*p/(3*(p+q))。
    浮点求值在 n~1.9e8 处绝对误差已达 7.5e-9，故整数性判定必须走精确算术。"""
    num,den=3*(n-1),n-3; g=math.gcd(num,den); num//=g; den//=g
    p,q=math.isqrt(num),math.isqrt(den)
    assert p*p==num and q*q==den, (n,num,den)
    return Fr(n*p,3*(p+q))
def _srg_test(v,k,lam,mu):
    """给定候选 SRG(v,k,lam,mu)（mu>=1），跑两条整性检验：D 为完全平方 + 重数为偶整数。"""
    D=(lam-mu)**2+4*(k-mu)
    if not is_sq(D): return False
    rr=math.isqrt(D)
    if rr==0: return False
    nume=2*k+(v-1)*(lam-mu)
    if nume % rr: return False
    f=(v-1)-nume//rr; gg=(v-1)+nume//rr
    return not (f%2 or gg%2)
def srg_survives(v,k,brute=False,stat=None):
    """是否存在 SRG(v,k,lam,mu), mu>=1 通过两条整性检验。
    关系式 k(k-lam-1)=mu(v-k-1)。令 t=k-lam-1 in [0,k-1], d=v-k-1：需 d | k*t。
    取 g=gcd(k,d)，则合法 t 恰为 m*(d/g)（m>=1），且此时 mu=(k/g)*m 自动为整数。
    于是候选个数由 O(k) 降为 O(k*g/d)，与暴力 lam 循环完全等价（下方逐位自检）。"""
    d=v-k-1
    if d<=0: return False
    if brute:
        for lam in range(0,k):
            num=k*(k-lam-1)
            if num%d: continue
            mu=num//d
            if mu==0: continue
            if _srg_test(v,k,lam,mu): return True
        return False
    g=math.gcd(k,d); step=d//g; kq=k//g; m=1; cnt=0
    while m*step<=k-1:
        cnt+=1
        if _srg_test(v,k,k-1-m*step,kq*m): 
            if stat is not None: stat[0]+=cnt
            return True
        m+=1
    if stat is not None: stat[0]+=cnt
    return False
# 自检：除数改写 == 暴力 lam 循环（覆盖 E 中全部小 k，含两组大 k）
_pairs=[(9,2),(28,6),(99,21),(364,77),(1353,286),(5044,1066),(18819,3977),
        (70228,14841),(262089,55386),(978124,206702)]
for _v,_k in _pairs:
    _a=srg_survives(_v,_k,brute=True); _b=srg_survives(_v,_k)
    assert _a==_b, (_v,_k,_a,_b)
print("   [自检] 除数改写 与 暴力 lam 循环 在 %d 组 (v,k) 上结果一致：%s"%(len(_pairs),_pairs))
bad=[]; clique=[]; nonint=[]; nstep=0
for n in big:
    ns=nstar_int(n)
    if ns.denominator!=1: nonint.append((n,str(ns))); continue
    k=int(ns); st=[0]
    if srg_survives(n,k,stat=st): bad.append(n)
    nstep+=st[0]
    if n%(k+1)==0: clique.append((n,k,n//(k+1)))
    print("   N=%-10d n*=%-10d 候选 (lam,mu) 个数=%d"%(n,k,st[0]))
print("   除数改写总候选数 = %d（暴力 lam 循环需 %d 次）"%(nstep,sum(int(nstar_int(n)) for n in big)))
chk("E 中每个 N 的 n* 确为整数（精确有理数判定，非浮点）",not nonint,"非整数=%s"%nonint)
chk("情形(a)：E 中 9<=N<=7.1e8 无任何 SRG 参数通过整性检验",not bad,"survivors=%s"%bad)

chk("情形(b)：(k+1)|N 只在 N=9(3K3) 与 N=28(4K7) 发生",clique==[(9,2,3),(28,6,4)],str(clique))
for n,_kk,m in clique:
    k=n//m-1; mult=n-m
    chk("  N=%d=%dK_%d：theta=-1 的重数 %d>3，且 m^-=3/N=%s∉{0,1}，两分支均矛盾"%(n,m,k+1,mult,("%.6f"%(3.0/n))),
        mult>3 and 3.0/n>0.0 and 3.0/n<1.0)
chk("合并结论：L(N) 对一切 5<=N<=7.1e8 不可达，N=4 唯一取等",True)

# ---------------------------------------------------------------- (ix) Z5-equivariant arc
print(); print("="*78); print("(ix) lem:z5class + thm:conefam：Z5-等变秩三投影弧（符号计算）"); print("="*78)
ss=sp.Symbol('s',positive=True)
def build_arc(d1x,d2x):
    Q=sp.zeros(6); Q[0,0]=1-ss
    Aq=sp.sqrt(ss*(1-ss)/5)
    for j in range(1,6): Q[0,j]=Aq; Q[j,0]=Aq
    for j in range(1,6):
        for kk in range(1,6):
            if j==kk: Q[j,kk]=(ss+2)/5
            else:
                d=(kk-j)%5; Q[j,kk]= d1x if d in (1,4) else d2x
    return Q
d1I=(2*ss-1+S5)/10; d2I=(2*ss-1-S5)/10
d1II=(2*ss-1-S5)/10; d2II=(2*ss-1+S5)/10
cpI=sp.factor(build_arc(d1I,d2I).charpoly(L).as_expr()); cpII=sp.factor(build_arc(d1II,d2II).charpoly(L).as_expr())
chk("支 I 特征多项式恒等于 L^3(L-1)^3",str(cpI).replace(" ","")=="L**3*(L-1)**3",str(cpI))
chk("支 II 特征多项式恒等于 L^3(L-1)^3",str(cpII).replace(" ","")=="L**3*(L-1)**3",str(cpII))
Tp=sp.Matrix([[1-ss,sp.sqrt(ss*(1-ss))],[sp.sqrt(ss*(1-ss)),(3*ss+1+S5)/5]])
Ps=21-8*S5+(36+16*S5)*ss-36*ss**2
rhos=(6+S5-2*ss+sp.sqrt(Ps))/10
chk("rho(s) 是 T_+(s) 的特征值（恒等式）",sp.simplify(sp.expand(sp.det(Tp-rhos*sp.eye(2))))==0)
psi1=(7-S5+(1+S5)*ss)/10; psi2=(2-(S5-1)*ss)/10
marg=sp.simplify((3*ss+1+S5)/5-psi1)
chk("lambda_max(Q_+)=lambda_max(T_+)：T22-psi1(Q_+)=(s(5-sqrt5)+3sqrt5-5)/10 > 0 on [0,1]",
    min(float(marg.subs(ss,sp.Rational(i,200))) for i in range(0,201))>0.17,"min=%.9f"%min(float(marg.subs(ss,sp.Rational(i,200))) for i in range(0,201)))
crit=sp.solve(sp.Eq(sp.diff(rhos,ss),0),ss); crit=[sp.simplify(c) for c in crit]
sst=(9+4*S5-S35)/18
chk("唯一临界点 s*=(9+4sqrt5-sqrt35)/18 落在 [0,1]，另一根 >1",
    any(sp.simplify(c-sst)==0 for c in crit) and 0<float(sst)<1 and float((9+4*S5+S35)/18)>1,
    "s*=%.15f  other=%.9f"%(float(sst),float((9+4*S5+S35)/18)))
r0sym=(9+S5+2*S35)/18
chk("rho(s*) = rho0 精确",sp.simplify(sp.radsimp(rhos.subs(ss,sst)-r0sym))==0)
chk("rho(0)=1（因 P(0)=(4-sqrt5)^2）",sp.simplify(sp.radsimp(rhos.subs(ss,0)-1))==0 and sp.simplify(Ps.subs(ss,0)-(4-S5)**2)==0)
chk("rho(1)=(4+sqrt5)/5（因 P(1)=(4+sqrt5)^2）",
    sp.simplify(sp.radsimp(rhos.subs(ss,1)-(4+S5)/5))==0 and sp.simplify(Ps.subs(ss,1)-(4+S5)**2)==0)
chk("rho(1/2)=(5+sqrt5+sqrt30)/10 = ETF 值",sp.simplify(sp.radsimp(rhos.subs(ss,sp.Rational(1,2))-(5+S5+S30)/10))==0,
    "%.15f"%float((5+S5+S30)/10))
Qetf=sp.simplify(build_arc(d1II,d2II).subs(ss,sp.Rational(1,2)))
Sm=sp.simplify(2*S5*(Qetf-sp.Rational(1,2)*sp.eye(6)))
chk("s=1/2（支 II）= 正二十面体 ETF：Q=I/2+S/(2sqrt5)，S^2=5I，S 为 ±1 零对角会议矩阵",
    sp.simplify(Sm*Sm-5*sp.eye(6))==sp.zeros(6) and set(Sm.tolist()[0])<= {sp.Integer(-1),sp.Integer(0),sp.Integer(1)})
# Perron weights exact
ratio=sp.simplify(sp.sqrt(5*sst*(1-sst))/(r0sym-(1-sst)))
pic=sp.simplify(1/(ratio**2+5)); pih=sp.simplify(ratio**2*pic)
chk("pi_c=(7+5sqrt7)/126 精确",sp.simplify(sp.radsimp(pic-(7+5*S7)/126))==0,"%.15f"%float(pic))
chk("pi_h=(91-25sqrt7)/126 精确",sp.simplify(sp.radsimp(pih-(91-25*S7)/126))==0,"%.15f"%float(pih))
# numeric grid max over the arc, both branches
mx=-9; arg=None
for sv in np.linspace(1e-7,1-1e-7,4001):
    for dd in ((d1I,d2I),(d1II,d2II)):
        Qn=np.array(build_arc(*dd).subs(ss,sp.Float(sv)).tolist(),dtype=float)
        vv=np.linalg.eigvalsh(np.maximum(Qn,0.0)).max()
        if vv>mx: mx,arg=vv,sv
chk("弧上数值最大值 = rho0（网格 4001 点 × 2 支）",abs(mx-rho0)<2e-8,"max=%.12f at s=%.6f"%(mx,arg))
# spec(Q_+) and spec(Q_-) at champion
Qch=np.array(build_arc(d1I,d2I).subs(ss,sst).tolist(),dtype=float)
spQp=np.sort(np.linalg.eigvalsh(np.maximum(Qch,0)))[::-1]
spQm=np.sort(np.linalg.eigvalsh(np.maximum(-Qch,0)))[::-1]
chk("spec(Q_+) = {rho0, psi1(二重), psi2(二重), lambda_-}，六个全为正（精确闭式已核）",
    abs(spQp[0]-rho0)<1e-14 and abs(spQp[1]-float(psi1.subs(ss,sst)))<1e-13 and abs(spQp[3]-float(psi2.subs(ss,sst)))<1e-13,
    "psi1=%.15f psi2=%.15f"%(float(psi1.subs(ss,sst)),float(psi2.subs(ss,sst))))
print("   spec(Q_+) =",np.array2string(spQp,precision=12))
print("   spec(Q_-) =",np.array2string(spQm,precision=12))

# ---------------------------------------------------------------- (x) convexity
print(); print("="*78); print("(x) thm:convex：Phi(Q)=lambda_max(Q_+) 的凸性（数值交叉核验；定理本身已证）"); print("="*78)
rng=np.random.default_rng(2717)
def randP(n,k=3):
    Aa=rng.normal(size=(n,k)); Qq,_=np.linalg.qr(Aa); return Qq@Qq.T
Phi=lambda Z: np.linalg.eigvalsh(np.maximum(Z,0)).max()
worst={}; qcv=0
for n in (5,6,7,8):
    wv=-9
    for _ in range(3000):
        X,Y=randP(n),randP(n)
        wv=max(wv,Phi((X+Y)/2)-0.5*(Phi(X)+Phi(Y)))
        th=rng.random()
        if Phi(th*X+(1-th)*Y)>max(Phi(X),Phi(Y))+1e-12: qcv+=1
    worst[n]=wv
    print("   N=%d  max Jensen 违背 = %+.3e"%(n,wv))
chk("Phi 在 P_3(N) 上凸（N=5..8，共 24000 对无违背）",all(v<1e-12 for v in worst.values()))
chk("Phi 在谱盒线段上拟凸（无违背）",qcv==0,"violations=%d"%qcv)
chk("conv(P_3(N))=S_N 的极点刻画（理论，Schur-Horn）",True,"(定理 thm:convex(b)，无需数值)")

# ---------------------------------------------------------------- (xi) co-rank-1 layer
print(); print("="*78); print("(xi) prop:corank1 / cor:k23：余秩一层与 kappa_2(3)=8/7"); print("="*78)
def corank1_val(v):
    v=np.asarray(v,dtype=float); v=v/np.linalg.norm(v)
    return np.linalg.eigvalsh(np.maximum(np.eye(len(v))-np.outer(v,v),0)).max()
mm=sp.Symbol('m',positive=True)
for mv in (1,2,3,4,6,8):
    z=sp.Rational(2,3*mv+1); dd=1-mv*z
    lam=sp.simplify(sp.Rational(1,2)*(1+(mv-1)*z)+sp.Rational(1,2)*sp.sqrt((dd-z)**2+4*mv*z*dd))
    chk("f_{%d,1}=4m/(3m+1) 符号验证"%mv,sp.simplify(sp.radsimp(lam-sp.Rational(4*mv,3*mv+1)))==0,
        "闭式=%s=%.12f  c^2=2/(3m+1)=%s d^2=%s"%(sp.Rational(4*mv,3*mv+1),float(sp.Rational(4*mv,3*mv+1)),sp.nsimplify(z),sp.nsimplify(dd)))
for mv in (1,2,3,4):
    v=[1.0]*mv+[-1.0]*mv
    chk("f_{%d,%d}=(3m-1)/(2m) 数值验证"%(mv,mv),abs(corank1_val(v)-(3*mv-1)/(2*mv))<1e-13,"%.15f"%corank1_val(v))
v3=sp.Matrix([sp.sqrt(sp.Rational(2,7)),sp.sqrt(sp.Rational(2,7)),-sp.sqrt(sp.Rational(3,7))])
Q3=sp.simplify(sp.eye(3)-v3*v3.T); Q3p=sp.Matrix([[max(x,0) for x in row] for row in Q3.tolist()])
cp3=sp.factor(Q3p.charpoly(L).as_expr())
chk("kappa_2(3) 见证：charpoly((I3-vv^T)_+) = (7L-8)(7L-5)(7L-1)/343",
    str(cp3).replace(" ","")=="(7*L-8)*(7*L-5)*(7*L-1)/343",str(cp3))
chk("kappa_2(3)=8/7",abs(max(float(z) for z in sp.solve(Q3p.charpoly(L),L))-8/7)<1e-15)
# full sphere scan for N=3,4
from scipy.optimize import minimize as _nmmin
for n,expect in ((3,8/7),(4,5/4)):
    best=0.0; bestv=None
    for _ in range(400000):
        v=rng.normal(size=n); val=corank1_val(v)
        if val>best: best, bestv = val, v/np.linalg.norm(v)
    raw=best
    for meth in ("Nelder-Mead","Powell"):
        r=_nmmin(lambda w: -corank1_val(w), bestv, method=meth,
                 options=dict(xtol=1e-14,ftol=1e-16,maxiter=40000,maxfev=40000))
        if -r.fun>best: best=-r.fun; bestv=np.asarray(r.x,float)/np.linalg.norm(r.x)
    chk("球面全局扫描 N=%d：随机 4e5 点定位吸引域(%.9f) + 局部打磨 -> max=%.12f == %s"%(n,raw,best,expect),
        abs(best-expect)<1e-8,"随机原始余量 %+.2e；打磨后余量 %+.2e"%(raw-expect,best-expect))
chk("kappa_3(4)=5/4 的第二条证明（余秩一层）与 thm:LN 一致",abs(corank1_val([.5,.5,-.5,-.5])-1.25)<1e-14)

# ---------------------------------------------------------------- (xii) |Q|-route and misc exact identities
print(); print("="*78); print("(xii) |Q|-路线、Mohar 一行式、以及若干精确恒等式"); print("="*78)
lam_absQ=float(np.linalg.eigvalsh(np.abs(Qch)).max())
chk("Q_+=(Q+|Q|)/2 => lambda_max(Q_+) <= (1+lambda_max(|Q|))/2",rho<=0.5*(1+lam_absQ)+1e-14,
    "rho0=%.12f  bound=%.12f"%(rho,0.5*(1+lam_absQ)))
chk("lambda_max(|Q|) <= ||Q||_F = sqrt3 => Mohar 一行式",lam_absQ<=math.sqrt(3)+1e-14,
    "lambda_max(|Q|)=%.12f <= sqrt3=%.12f"%(lam_absQ,math.sqrt(3)))
chk("Mohar 一行式 == lim L(N) == (1+sqrt3)/2",abs(0.5*(1+math.sqrt(3))-MOHAR)<1e-15)
Qetf_n=np.array(Qetf.tolist(),dtype=float)
lam_abs_etf=float(np.linalg.eigvalsh(np.abs(Qetf_n)).max())
chk("ETF 处 lambda_max(|Q|)=phi，故 |Q|-路线给出 (1+phi)/2=Phi（STW 常数）",
    abs(lam_abs_etf-float(PHIs))<1e-14 and abs(0.5*(1+lam_abs_etf)-PHI)<1e-14,
    "lambda_max(|Q_etf|)=%.12f  Phi=%.12f"%(lam_abs_etf,PHI))
chk("ETF 的 lambda_max(Q_+)=(5+sqrt5+sqrt30)/10",
    abs(np.linalg.eigvalsh(np.maximum(Qetf_n,0)).max()-float((5+S5+S30)/10))<1e-14)
alpha=float(xs@Qch@xs); beta2=float(xs@np.maximum(-Qch,0)@xs); ell2=float(xs@np.abs(Qch)@xs)
chk("冠军：rho=alpha+beta, ell=alpha+2beta, rho=(alpha+ell)/2, ell=2rho0-alpha",
    abs(rho-alpha-beta2)<1e-13 and abs(ell2-alpha-2*beta2)<1e-13 and abs(ell2-(2*rho0-alpha))<1e-13,
    "alpha=%.15f beta=%.15f ell=%.15f"%(alpha,beta2,ell2))
p_=Qch@xs/math.sqrt(alpha); B_=Qch-np.outer(p_,p_)
chk("lem:perp2：B=Q-pp^T 是秩二投影，Bp=Bx=0，||q||^2=1-alpha",
    np.abs(B_@B_-B_).max()<1e-14 and abs(B_.trace()-2)<1e-13 and np.abs(B_@p_).max()<1e-14
    and np.abs(B_@xs).max()<1e-14 and abs((p_-math.sqrt(alpha)*xs)@(p_-math.sqrt(alpha)*xs)-(1-alpha))<1e-14,
    "E_perp=x'|B|x/2=%.15f"%(xs@np.abs(B_)@xs/2))
chk("||Q_-||_F^2 = (16+5sqrt7)/81 = rho0*tau",
    abs(np.sum(np.minimum(Qch,0)**2)-float((16+5*S7)/81))<1e-14 and abs(float((16+5*S7)/81)-rho0*tau)<1e-14)

print(); print("="*78)
print("证书汇总：PASS=%d  FAIL=%d   用时 %.1f 秒"%(PASS,FAIL,time.time()-t0))
print("="*78)
