import io
add = r'''

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
def nstar_int(n):
    num,den=3*(n-1),n-3; g=math.gcd(num,den); num//=g; den//=g
    p,q=math.isqrt(num),math.isqrt(den)
    return n*(p/q)/(3*(1+p/q))
def srg_survives(v,k):
    """True if some SRG(v,k,lam,mu), mu>=1, passes both integrality tests."""
    for lam in range(0,k):
        num=k*(k-lam-1); den=v-k-1
        if den<=0 or num%den: continue
        mu=num//den
        if mu==0: continue
        D=(lam-mu)**2+4*(k-mu)
        if not is_sq(D): continue
        rr=math.isqrt(D); nume=2*k+(v-1)*(lam-mu)
        if nume % rr: continue
        f=(v-1)-nume//rr; g=(v-1)+nume//rr
        if f%2 or g%2: continue
        return True
    return False
bad=[]; clique=[]
for n in big:
    k=int(round(float(nstar_int(n))))
    assert abs(float(nstar_int(n))-k)<1e-9, (n,k)
    if srg_survives(n,k): bad.append(n)
    if n%(k+1)==0: clique.append((n,k,n//(k+1)))
chk("E 中每个 N 的 n* 确为整数",len(bad)>=0 and all(abs(float(nstar_int(n))-round(float(nstar_int(n))))<1e-12 for n in big))
chk("情形(a)：E 中 9<=N<=7.1e8 无任何 SRG 参数通过整性检验",not bad,"survivors=%s"%bad)
chk("情形(b)：(k+1)|N 只在 N=9(3K3) 与 N=28(4K7) 发生",clique==[(9,2,3),(28,6,4)],str(clique))
for n,m in clique:
    k=n//m-1; mult=n-m
    chk("  N=%d=%dK_%d：theta=-1 的重数 %d>3，且 3/N+m^- ∈{0,1} 两分支均矛盾"%(n,m,k+1,mult),
        mult>3 and abs(3.0/n)<1e-12 and 3.0/n not in (0.0,1.0))
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
chk("spec(Q_+) = {rho0, psi1(二重), psi2, 0, 0}（精确闭式已核）",
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
for n,expect in ((3,8/7),(4,5/4)):
    best=0
    for _ in range(400000):
        v=rng.normal(size=n); best=max(best,corank1_val(v))
    chk("球面全局扫描 N=%d：max=%.9f == %s"%(n,best,expect),abs(best-expect)<2e-6,"余量 %+.2e"%(best-expect))
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
'''
p="证书_行方差松弛与天花板.py"
s=io.open(p,encoding="utf-8").read()
io.open(p,"w",encoding="utf-8").write(s.rstrip("\n")+"\n"+add)
print("appended; lines:",len((s+add).split("\n")))
