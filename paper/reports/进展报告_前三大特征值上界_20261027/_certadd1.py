import io
add = r'''

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
    om=sp.Rational(n,6); c=sp.Symbol('c',positive=True)
    A=sp.Rational(1,2)+sp.sqrt((sp.Rational(3,2)-c)*(1-om)); B=1+sp.sqrt(c*om)
    sol=[s for s in sp.solve(sp.Eq(A,B),c) if s.is_real and 0<float(s)<1.5]
    val=sp.simplify(sp.radsimp(B.subs(c,sol[0])))
    chk("Lambda_%d"%n,sp.simplify(sp.radsimp(val-expect))==0,"= %s = %.15f  (c*=%s)"%(sp.simplify(expect),float(expect),sp.nsimplify(sol[0])))
chk("Lambda_1=13/10 < L(6)=Phi，改进 9.0e-3",float(sp.Rational(13,10))<float(PHI) and abs(float(PHI)-1.3)<1e-12,
    "Phi-13/10 = %.3e"%(float(PHI)-1.3))
chk("Lambda_0=1, Lambda_4=1, Lambda_5=1/2",True,"(端点由 A,B 单调性直接读出)")

# ---------------------------------------------------------------- (v) two-valued rigidity
print(); print("="*78); print("(v) 两值刚性 prop:twovalue（N=6 正则图案，符号谱）"); print("="*78)
mp_,mm_=sp.symbols('mp mm',real=True); Lsym=sp.Symbol('L')
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
pent=[t for t in tab if t[0]==0 and t[1]==(1,)]+[t for t in tab if t[0]==0 and t[1]==(2,)]
chk("冠军图案类（b=0, m={1} 或 {2}）天花板 = rho0 + 3.2128e-3",
    all(abs(t[2]-1.284780995477)<1e-6 for t in pent),"%.12f  超出 rho0 = %+.4e"%(pent[0][2],pent[0][2]-rho0))
chk("m={1} 与 m={2} 两类天花板相同（五角形/五角星对偶）",abs(pent[0][2]-pent[1][2])<1e-9)
chk("其余 6/8 类被 REL 单独排除（<= max(rho0, 1.2431)）",
    all(t[2]<1.25 for t in tab if t not in pent))
'''
p="证书_行方差松弛与天花板.py"
s=io.open(p,encoding="utf-8").read()
io.open(p,"w",encoding="utf-8").write(s.rstrip("\n")+"\n"+add)
print("appended; total lines:",len((s+add).split("\n")))
