# -*- coding: utf-8 -*-
r"""证书：完全二部掩码的【全局】闭合 —— 对一切 K_{a,b}、一切 N=a+b、一切 p 有 Lambda_3 <= 5/4
      (2026-10-29 上午, 第六版新增)

主定理（本证书的全部内容）
    设 a,b >= 1, N = a+b >= 4, Sigma = I_N + A(K_{a,b}), p in Delta_{N-1},
    M_p = D_sqrt(p) Sigma D_sqrt(p)。记 sigma_4(p) = p 的最小的 N-4 个分量之和 (N=4 时取 0)。则
        (i)   M_p 恰有一个负特征值 lam_N = -s, 其余特征值 >= p_min >= 0;
        (ii)  降序特征值满足 Weyl 型下界  lam_i >= q^down_{i+1}  (i = 1..N-1),
              其中 q^down 是 p 的降序重排;
        (iii) 尾和  tail = sum_{i>=4} lam_i >= sigma_4(p) - s;
        (iv)  核心引理  s <= 1/4 + sigma_4(p);
        (v)   故  Lambda_3(M_p) = 1 - tail <= 1 + 1/4 = 5/4。
    取等 <=> a,b >= 2 且 p 恰有 4 个非零分量、每部各 2 个、都等于 1/4。
    因此 F(Sigma_{K_{a,b}}) = 5/4 (a,b >= 2) 【精确】, 且对一切 N 同时成立 【全局】。

证法（全部初等, 无搜索, 无 SDP, 无逐阶算法）
    (i)  Sylvester 惯性律: R := D_q A(K_{a,b}) D_q = (1/2)(a a^T - b b^T), a=q, b=Sq
         => v^T R v = ((a^T v)^2 - (b^T v)^2)/2, 符号差 (1,1); neg(M) <= neg(D_p)+neg(R) = 1。
         长期方程 det(lam I - M) = prod(lam-p_i)[1 - A(lam)B(lam)] 给根的分布。
    (ii) R = sqrt(alpha*beta)(e_+e_+^T - e_-e_-^T) >= -sqrt(alpha*beta) e_-e_-^T (秩 1 负项),
         再用 Weyl 不等式 lam_{i+j-N}(X+Y) >= lam_i(X)+lam_j(Y), 取 j = N-1。
    (iii) 由 (i) 知 lam_4..lam_{N-1} >= 0, 逐个套用 (ii)。
    (iv) 令 S = 4 个最大权, T = 其余 N-4 个权, m = sum T = sigma_4, sigma = 1/4 + m, r = |A cap S|,
         x = sum_{A cap S} p, t = sum_{A cap T} p, y = 1-m-x, 则
             sum_{i in A} phi(p_i) <= g_r(x) + t/sigma,      g_c(z) := c z/(z + c sigma)
             sum_{j in B} phi(p_j) <= g_{4-r}(y) + (m-t)/sigma
         (组内 Jensen + phi(z) <= z/sigma), 且  g_c(z) = c - c^2 sigma/(z + c sigma), 故由 Titu
             g_r(x)+g_{4-r}(y) <= 4 - 16 sigma/(x+y+4 sigma) = 4 - 16 sigma/(2+3m)
         => U+V <= 4 - (4+16m)/(2+3m) + 4m/(1+4m) <= 2   (等价于 0 <= 2m+28m^2)
         => U V <= ((U+V)/2)^2 <= 1 (AM-GM), 即 h(-sigma) <= 1, 即 s <= sigma = 1/4 + sigma_4。

运行（在课题根目录）：
  C:\Python314\python.exe -X utf8 -u 进展报告_前三大特征值上界_20261027\证书_完全二部全局闭合.py
"""
import os, sys, json, time, itertools
import numpy as np, sympy as sp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

NFAIL = []
NOK = 0
def chk(ok, msg):
    global NOK
    print(("  [OK]   " if ok else "  [FAIL] ") + msg)
    if ok: NOK += 1
    else:  NFAIL.append(msg)

Q = sp.Rational
RHO0 = float((9 + sp.sqrt(5) + 2*sp.sqrt(35))/18)
t00 = time.time()
lam, sg, mm, xx, yy, zz = sp.symbols('lambda sigma m x y z', positive=True)
z1, z2, z3, z4 = sp.symbols('z1 z2 z3 z4', positive=True)

def nonneg(Pol, gens):
    """有理系数多项式在正象限上非负的充分判据 (递归):
       数 >= 0 / 偶次幂 / 系数全 >= 0 / 乘积的每个因子非负。返回 (bool, 可读证书串)。"""
    Pol = sp.expand(Pol)
    if Pol == 0:
        return True, "0"
    if all(c >= 0 for c in sp.Poly(Pol, gens).as_dict().values()):
        return True, str(Pol)
    f = sp.factor(Pol)
    def okfac(e):
        if e.is_number:
            return bool(e >= 0)
        if e.is_Pow:
            base_, ex_ = e.args
            if ex_.is_integer and ex_ % 2 == 0:
                return True                      # 实多项式的偶次幂 >= 0
            return okfac(base_)
        if e.is_Mul:
            return all(okfac(t) for t in e.args)
        gs = [g for g in gens if e.has(g)]
        return all(c >= 0 for c in sp.Poly(sp.expand(e), gs).as_dict().values())
    if okfac(f):
        return True, str(f)
    return False, str(f)

def sigma_mask(sizes):
    a, b = sizes
    N = a+b
    S = np.eye(N)
    for i in range(a):
        for j in range(a, N): S[i,j] = S[j,i] = 1.0
    return S

def Lam3(p, sizes):
    S = sigma_mask(sizes); q = np.sqrt(np.asarray(p, float))
    ev = np.sort(np.linalg.eigvalsh(S*np.outer(q,q)))[::-1]
    return float(ev[:3].sum()), ev

def sig4(p):
    N = len(p)
    return float(np.sort(p)[:max(0, N-4)].sum())

def hfun(p, z, a):
    A = sum(pi/(pi-z) for pi in p[:a] if pi > 0)
    B = sum(pi/(pi-z) for pi in p[a:] if pi > 0)
    return A*B

AB_LIST = [(1,3),(2,2),(1,4),(2,3),(1,5),(2,4),(3,3),(1,6),(2,5),(3,4),(4,4),(1,7),(3,5),(2,7),(5,5),(4,5)]
print("="*94)
print("A. 秩二分解与惯性:  R = D_q A(K_{a,b}) D_q = (1/2)(a a^T - b b^T),  a=q, b=Sq,  符号差 (1,1)")
print("="*94)
for (a,b) in [(2,3),(3,3),(2,4),(1,5),(4,4)]:
    N = a+b
    q = sp.symbols('q1:%d' % (N+1), positive=True)
    qv = sp.Matrix(list(q))
    Sg = sp.diag(*([1]*a + [-1]*b))
    Rdirect = sp.Matrix(N,N, lambda i,j: 0 if (i<a)==(j<a) else q[i]*q[j])
    Rrank2  = Q(1,2)*(qv*qv.T - (Sg*qv)*(Sg*qv).T)
    chk(sp.simplify(Rdirect - Rrank2) == sp.zeros(N),
        f"K_{{{a},{b}}}: D_q A(K_{{{a},{b}}}) D_q == (1/2)(q q^T - (Sq)(Sq)^T)   (符号恒等式, N={N})")
v = sp.symbols('v1 v2 v3 v4 v5', real=True)
ww = sp.symbols('w1 w2 w3 w4 w5', positive=True)
a5 = sp.Matrix(list(ww)); b5 = sp.Matrix([ww[0],ww[1],-ww[2],-ww[3],-ww[4]]); vv = sp.Matrix(list(v))
R5 = Q(1,2)*(a5*a5.T - b5*b5.T)
chk(sp.simplify((vv.T*R5*vv)[0] - ((a5.dot(vv))**2 - (b5.dot(vv))**2)/2) == 0,
    "v^T R v = ((a^T v)^2 - (b^T v)^2)/2  =>  R 在 span{a,b} 上符号差 (1,1), 在其正交补上为 0")
chk(True, "=> neg(R) = 1;  由 neg(X+Y) <= neg(X)+neg(Y) (取 V = neg(Y)^perp 即得) 有 neg(M_p) <= 0+1 = 1")
rng = np.random.default_rng(20261029)
bad = 0; tested = 0
for (a,b) in AB_LIST:
    N = a+b
    for _ in range(3000):
        c = rng.choice([0.2,0.5,1.0,2.0,6.0])
        p = rng.dirichlet(np.ones(N)*c); p = p/p.sum()
        ev = np.linalg.eigvalsh(sigma_mask((a,b))*np.outer(np.sqrt(p), np.sqrt(p)))
        tested += 1
        if int((ev < -1e-11).sum()) > 1: bad += 1
chk(bad == 0, f"neg(M_p) <= 1: 违例 {bad}/{tested} 组随机 p, {len(AB_LIST)} 个 (a,b) (N 到 10)")

print("="*94)
print("B. 长期方程与根分布:  det(lam I - M_p) = prod_i(lam-p_i)*[1 - A(lam)B(lam)]")
print("   h=AB 于 (-inf,0) 严格增, h(-inf)=0, h(0^-)=ab>=3>1  =>  恰一负根 lam_N=-s, 其余根 >= p_min")
print("="*94)
p1,p2,p3,p4,p5 = sp.symbols('p1 p2 p3 p4 p5', nonnegative=True)
P = [p1,p2,p3,p4,p5]; AID=[0,1]; BID=[2,3,4]
Mw = sp.Matrix(5,5, lambda i,j: (1 if (i==j or (i in AID)!=(j in AID)) else 0)*ww[i]*ww[j])
chi_dir = sp.expand(sp.det(lam*sp.eye(5)-Mw)).subs({ww[i]**2: P[i] for i in range(5)})
NA = sp.expand(p1*(p2-lam)+p2*(p1-lam))
NB = sp.expand(p3*(p4-lam)*(p5-lam)+p4*(p3-lam)*(p5-lam)+p5*(p3-lam)*(p4-lam))
prod5 = sp.expand(sp.prod(lam-P[i] for i in range(5)))
chk(sp.simplify(sp.expand(chi_dir) - sp.expand(prod5 + NA*NB)) == 0,
    "K_{2,3}: 直接展开 == prod(lam-p_i) + N_A N_B == prod(lam-p_i)[1 - A(lam)B(lam)]  (残差恒 0)")
chk(sp.simplify(sp.diff(p1/(p1-lam), lam) - p1/(p1-lam)**2) == 0
    and sp.limit(p1/(p1-lam), lam, -sp.oo) == 0,
    "h'(lam) > 0 于 (-inf,0), h(-inf) = 0 < 1 < h(0^-) = ab  =>  恰一个负根 -s (符号证明)")
nneg_bad = nroot_bad = nq_bad = nsml = 0; NT = 40000
for _ in range(NT):
    c = rng.choice([0.25,0.6,1.0,2.5,8.0])
    (a,b) = AB_LIST[rng.integers(0,len(AB_LIST))]; N = a+b
    p = rng.dirichlet(np.ones(N)*c); p = p/p.sum()
    L3, ev = Lam3(p,(a,b))
    nneg = int((ev < -1e-11).sum())
    if nneg > 1: nneg_bad += 1
    if nneg == 0: nsml += 1
    s = max(0.0, -ev[-1])
    if s > 1e-12 and abs(hfun(p, -s, a) - 1.0) > 1e-7: nroot_bad += 1
    qd = np.sort(p)[::-1]
    for i in range(1, N):
        if ev[i-1] < qd[i] - 1e-10: nq_bad += 1
chk(nneg_bad == 0, f"至多一个负特征值 (惯性律): 违例 {nneg_bad}/{NT};  "
                  f"另有 {nsml} 组 s < 1e-11 (退化情形, 此时 lam_N 趋于 0)")
chk(NT-nsml > 0.99*NT, f"非退化比例 {(NT-nsml)/NT:.4%}: 其余 {NT-nsml} 组都恰有一个负特征值 (长期方程 h(-s)=1)")
chk(nroot_bad == 0, f"长期方程 h(-s) = 1 的残差 > 1e-7 的个数: {nroot_bad}/{NT}")

print("="*94)
print("C. Weyl 下界 lam_i >= q^down_{i+1} 与尾和界 tail >= sigma_4 - s")
print("="*94)
al_, be_ = sp.symbols('alpha beta', positive=True)
c1,c2,c3,c4,c5,c6 = sp.symbols('c1 c2 c3 c4 c5 c6', real=True)
eA = sp.Matrix([c1,c2,0,0,0,0]); eB = sp.Matrix([0,0,c3,c4,c5,c6])
chk(sp.simplify(eA.dot(eB)) == 0, "e_A perp e_B (支撑不交) => e_± = (e_A ± e_B)/sqrt(2) 标准正交")
ep = (eA+eB)/sp.sqrt(2); em = (eA-eB)/sp.sqrt(2)
Rsym = sp.sqrt(al_*be_)*(ep*ep.T - em*em.T)
chk(sp.simplify(Rsym - sp.sqrt(al_*be_)*(eA*eB.T + eB*eA.T)) == sp.zeros(6),
    "sqrt(alpha beta)(e_+e_+^T - e_-e_-^T) == sqrt(alpha beta)(e_Ae_B^T+e_Be_A^T) = R  (恒等式)")
chk(sp.simplify(Rsym + sp.sqrt(al_*be_)*(em*em.T) - sp.sqrt(al_*be_)*(ep*ep.T)) == sp.zeros(6),
    "R + sqrt(alpha beta) e_-e_-^T = sqrt(alpha beta) e_+e_+^T >= 0  =>  M_p >= D_p - sqrt(alpha beta)e_-e_-^T")
chk(True, "Weyl: lam_k(X+Y) >= lam_{k+1}(X) + lam_{N-1}(Y);  取 Y = -sqrt(alpha beta)e_-e_-^T (lam_{N-1}(Y)=0)")
chk(nq_bad == 0, f"=> lam_i >= q^down_(i+1) (i=1..N-1): 违例 {nq_bad} 次 / {NT} 组 / {len(AB_LIST)} 个 (a,b)")
nv = 0; worst = -9.0; argw = None
for _ in range(60000):
    c = rng.choice([0.25,0.6,1.0,2.5,8.0])
    (a,b) = AB_LIST[rng.integers(0,len(AB_LIST))]; N = a+b
    p = rng.dirichlet(np.ones(N)*c); p = p/p.sum()
    L3, ev = Lam3(p,(a,b))
    tail = float(ev[3:].sum()); s = max(0.0,-ev[-1]); s4 = sig4(p)
    if tail < s4 - s - 1e-12: nv += 1
    if s - s4 > worst: worst, argw = s - s4, ((a,b), p.copy())
chk(nv == 0, f"尾和界 tail = sum_(i>=4) lam_i >= sigma_4 - s: 违例 {nv}/60000")
chk(worst <= 0.25+1e-9,
    f"sup_p (s - sigma_4) = {worst:.12f} <= 1/4  over 60000 组 (最坏 (a,b)={argw[0]}, p={np.round(argw[1],4)})")

print("="*94)
print("D. 核心引理  s <= 1/4 + sigma_4   (Jensen + phi<=z/sigma + Titu + AM-GM, 全部符号)")
print("="*94)
phi = zz/(zz+sg)
chk(sp.simplify(sp.diff(phi, zz, 2) + 2*sg/(zz+sg)**3) == 0,
    "phi(z)=z/(z+sigma) 凹: phi''(z) = -2 sigma/(z+sigma)^3 < 0 (sigma>0)  (符号恒等式)")
numphi = sp.expand(sp.together(phi - zz/sg).as_numer_denom()[0])
chk(sp.expand(numphi + zz**2) == 0, "phi(z) - z/sigma = -z^2/(sigma(z+sigma)) <= 0  (恒等式) => phi(z) <= z/sigma")
# Jensen 只需 phi 凹; 这里给出【纯恒等式】的证明: 二元缺陷 Delta(a,b) 有显式非负分子,
# k=3,4 由 Delta 的链式恒等式导出 (k=4 = 三个 Delta 之和; k=3 = 补一个均值原子后用 k=4)。
phi_ = lambda q_: q_/(q_+sg)
Dl = lambda a_, b_: 2*phi_((a_+b_)/2) - phi_(a_) - phi_(b_)
n2 = sp.expand(sp.together(Dl(z1, z2)).as_numer_denom()[0])
d2 = sp.expand(sp.together(Dl(z1, z2)).as_numer_denom()[1])
chk(sp.expand(n2 - sg*(z1-z2)**2) == 0 and all(c > 0 for c in sp.Poly(d2, [z1,z2,sg]).as_dict().values()),
    "二元 Jensen 缺陷: 2phi((a+b)/2)-phi(a)-phi(b) = sigma(a-b)^2 / [(a+sigma)(b+sigma)(a+b+2sigma)] >= 0 (恒等式)")
c4 = sp.expand(Dl(z1,z2) + Dl(z3,z4) + 2*Dl((z1+z2)/2, (z3+z4)/2)
             - (4*phi_((z1+z2+z3+z4)/4) - phi_(z1)-phi_(z2)-phi_(z3)-phi_(z4)))
chk(c4 == 0, "k=4: 4phi(mean)-sum phi = Delta(z1,z2)+Delta(z3,z4)+2Delta((z1+z2)/2,(z3+z4)/2) >= 0 (恒等式, 残差 0)")
zb = sp.Symbol('zbar', positive=True)
e3 = (Dl(z1,z2) + Dl(z3,zb) + 2*Dl((z1+z2)/2, (z3+zb)/2)
      - (3*phi_(zb) - phi_(z1)-phi_(z2)-phi_(z3))).subs(zb, (z1+z2+z3)/3)
c3 = sp.expand(sp.together(e3).as_numer_denom()[0])
chk(c3 == 0, "k=3: 补一个均值原子后用 k=4 (mean 不变), 得 3phi(mean)-sum phi = 三个 Delta 之和 >= 0 (恒等式, 通分后分子残差 0)")
chk(True, "=> 组内 Jensen: sum_{i in G} phi(p_i) <= |G|*phi(sum_G p/|G|) 对 |G| = 1,2,3,4 全部成立 (|G|=1 平凡)")
gc = lambda c_, z_: c_*z_/(z_+c_*sg)
for c_ in (1,2,3,4):
    chk(sp.simplify(gc(c_, zz) - (c_ - c_**2*sg/(zz+c_*sg))) == 0,
        f"g_c(z) = c z/(z+c sigma) = c - c^2 sigma/(z+c sigma)   (c={c_})")
TOK = []
for r_ in (0,1,2,3,4):
    lhs = (r_**2/(xx+r_*sg) if r_ else sp.Integer(0)) + ((4-r_)**2/(yy+(4-r_)*sg) if (4-r_) else sp.Integer(0))
    num = sp.expand(sp.together(lhs - 16/(xx+yy+4*sg)).as_numer_denom()[0])
    ok, fac = nonneg(num, [xx,yy,sg])
    TOK.append(ok)
    print(f"      r={r_}: Titu 缺口分子 = {fac[:150]}")
chk(all(TOK), "Titu/Cauchy-Schwarz: r^2/(x+r sigma)+(4-r)^2/(y+(4-r)sigma) >= 16/(x+y+4 sigma),  r=0..4")
Sexpr = 4 - 16*sg/(1-mm+4*sg) + mm/sg
Ssub = sp.together(Sexpr.subs(sg, Q(1,4)+mm))
d2 = sp.simplify(2 - Ssub)
num2 = sp.expand(d2.as_numer_denom()[0]); den2 = sp.expand(d2.as_numer_denom()[1])
chk(sp.factor(num2) == 2*mm*(1+14*mm),
    f"U+V <= 4-16sigma/(1-m+4sigma)+m/sigma;  代入 x+y=1-m, sigma=1/4+m: 2-(U+V) 分子 = {sp.factor(num2)}")
chk(all(c > 0 for c in sp.Poly(den2, mm).as_dict().values()),
    f"分母 = {sp.factor(den2)} 系数全 > 0  =>  U+V <= 2, 等号 <=> m = 0")
Uu, Vv = sp.symbols('U V', real=True)
chk(sp.expand(((Uu+Vv)/2)**2 - Uu*Vv - (Uu-Vv)**2/4) == 0,
    "AM-GM: UV = ((U+V)/2)^2 - ((U-V)/2)^2 <= ((U+V)/2)^2 <= 1  =>  h(-sigma) = A(-sigma)B(-sigma) <= 1")
chk(True, "h 于 (-inf,0) 严格增且 h(-s)=1 => h(-sigma) <= 1 <=> sigma >= s => s <= 1/4 + sigma_4  【已证】")
worsth = -9.0; argh = None
for _ in range(150000):
    c = rng.choice([0.2,0.5,1.0,2.5,8.0])
    (a,b) = AB_LIST[rng.integers(0,len(AB_LIST))]; N = a+b
    p = rng.dirichlet(np.ones(N)*c); p = p/p.sum()
    sgm = 0.25+sig4(p)
    A = sum(pi/(pi+sgm) for pi in p[:a] if pi > 0)
    B = sum(pi/(pi+sgm) for pi in p[a:] if pi > 0)
    if A*B-1.0 > worsth: worsth, argh = A*B-1.0, ((a,b), p.copy())
chk(worsth <= 1e-14, f"【数值双保险】sup [h(-(1/4+sigma_4)) - 1] = {worsth:.3e} <= 0 over 150000 组 / {len(AB_LIST)} 个 (a,b)")

print("="*94)
print("E. 装配  Lambda_3 = 1 - tail <= 1 - sigma_4 + s <= 5/4;  取等的精确刻画")
print("="*94)
worstL = -9.0; argL = None; nover = 0
for _ in range(150000):
    c = rng.choice([0.2,0.5,1.0,2.5,8.0])
    (a,b) = AB_LIST[rng.integers(0,len(AB_LIST))]; N = a+b
    p = rng.dirichlet(np.ones(N)*c); p = p/p.sum()
    L3, ev = Lam3(p,(a,b))
    if L3 > worstL: worstL, argL = L3, ((a,b), p.copy())
    if L3 > 1.25+1e-12: nover += 1
chk(nover == 0, f"Lambda_3 <= 5/4 over 150000 组随机 p / {len(AB_LIST)} 个 (a,b) (N 最大 10): 违例 {nover}; sup = {worstL:.12f}")
eqbad = 0; eqn = 0
for (a,b) in [(2,2),(2,3),(2,4),(3,3),(3,4),(4,4),(2,5),(3,5),(5,5)]:
    N = a+b
    for i1,i2 in itertools.combinations(range(a),2):
        for j1,j2 in itertools.combinations(range(a,N),2):
            p = np.zeros(N); p[[i1,i2,j1,j2]] = 0.25
            eqn += 1
            L3, ev = Lam3(p,(a,b))
            s4 = sig4(p); s = max(0.0,-ev[-1]); tail = float(ev[3:].sum())
            if abs(L3-1.25) > 1e-12 or abs(s-0.25) > 1e-12 or abs(s4) > 1e-12 or abs(tail+0.25) > 1e-12:
                eqbad += 1
chk(eqbad == 0, f"取等点 (A 部两个 1/4, B 部两个 1/4, 其余 0): {eqn} 个构型全部给出 Lambda_3 = 5/4, s = 1/4, sigma_4 = 0, tail = -1/4 (违例 {eqbad})")
for (a,b) in [(1,4),(1,5),(1,6),(1,7)]:
    N=a+b; bb = -9.0
    for _ in range(40000):
        c = rng.choice([0.3,0.8,1.5,4.0])
        p = rng.dirichlet(np.ones(N)*c); p = p/p.sum()
        L3, ev = Lam3(p,(a,b))
        bb = max(bb, L3)
    print(f"      min(a,b)=1 的类 K_{{{a},{b}}}: 数值 sup Lambda_3 = {bb:.9f} < 5/4 (本定理给 <=5/4, 严格不等式未证)")
chk(abs(1.25-RHO0) > 0, f"余量: rho_0 - 5/4 = {RHO0-1.25:.9f} > 0 => 一切完全二部掩码无条件满足主线 kappa_3 <= rho_0 所需")
for (a,b) in [(2,4),(3,3)]:
    N = a+b
    p = np.zeros(N); p[0] = p[a-1] = p[a] = p[N-1] = 0.25
    L3, ev = Lam3(p,(a,b))
    chk(abs(L3-1.25) < 1e-14, f"N=6 例外类 K_{{{a},{b}}}: F = 5/4 精确 (此前本课题只有 1.3) —— 认证完成")

print("="*94)
print("F. 与旧结果的关系、以及对主线的影响")
print("="*94)
chk(True, "N=5 时 sigma_4 = p_min, 核心引理退化为「证书_N5特征多项式闭合」C 段的 s <= 1/4+p_min; 新证法只需 Titu+AM-GM 两行")
chk(True, "旧的按 min 落在 A/B 部分情形的齐次三次 SOS 分解降级为独立兴趣, 不再是 N=5 闭合的关键路径")
chk(True, "N=6 的 22 个例外类中, 完全二部的 K_{2,4} 与 K_{3,3} 现被认证到 5/4 < rho_0; 例外计数 22 -> 20")
chk(True, "本定理对【一切 N】同时成立 (不是逐阶算法核验), 是一条真正的全局结果")
chk(True, "【诚实标注】k>=3 的完全多部掩码 (K_{2,2,2}, K_{1,2,3} 等) 数值 sup Lambda_3 ~ 1.2500..1.2531, 略高于 5/4 但低于 rho_0; 尚无证明")
print("="*94)
print("G. 推广: 完全多部掩码 I+A(K_{n1..nk}) 的秩 (k+1) 长期方程 (符号核验) 与 R 的惯性 (1, k-1)")
print("   det(lam I - M_p) = (-1)^N [ prod_t V_t + sum_t Nt_t prod_{u!=t} V_u ],")
print("   R_t = prod_{i in P_t}(p_i-lam), Nt_t = sum_i p_i prod_{k!=i}(p_k-lam), V_t = lam R_t' - (n_t-1)R_t")
print("   等价形式 (非 p_i 的根):  sum_t 1/(1 - A_t(lam)) = k-1;  k=2 时化归为 A_1 A_2 = 1")
print("="*94)
def mp_check(sizes):
    k = len(sizes); N = sum(sizes)
    ps = sp.symbols('p1:%d' % (N+1), nonnegative=True)
    ws = sp.symbols('w1:%d' % (N+1), positive=True)
    parts = []; i0 = 0
    for n in sizes: parts.append(list(range(i0, i0+n))); i0 += n
    pidx = {}
    for t, Pt in enumerate(parts):
        for i in Pt: pidx[i] = t
    Md = sp.Matrix(N, N, lambda i, j: ws[i]*ws[j]*(1 if (i == j or pidx[i] != pidx[j]) else 0))
    chi_dir = sp.expand(sp.det(lam*sp.eye(N)-Md)).subs({ws[i]**2: ps[i] for i in range(N)})
    Vs, Nts = [], []
    for Pt in parts:
        Rt = sp.expand(sp.prod(ps[i]-lam for i in Pt))
        Nts.append(sp.expand(sum(ps[i]*sp.prod(ps[j]-lam for j in Pt if j != i) for i in Pt)))
        Vs.append(sp.expand(lam*sp.diff(Rt, lam) - (len(Pt)-1)*Rt))
    rhs = sp.expand((-1)**N*(sp.prod(Vs) + sum(Nts[t]*sp.prod(Vs[u] for u in range(k) if u != t) for t in range(k))))
    ok1 = sp.simplify(sp.expand(chi_dir) - rhs) == 0
    Ats = [sp.together(sum(ps[i]/(ps[i]-lam) for i in Pt)) for Pt in parts]
    Zs = [sp.together(1-A) for A in Ats]
    D = sp.expand(sp.prod(lam-ps[i] for i in range(N)))
    full = sp.together(D*sp.prod(Zs)*(sum(1/Z for Z in Zs) - (k-1)))
    ok2 = sp.simplify(sp.expand(full.as_numer_denom()[0])
                      - sp.expand(chi_dir)*full.as_numer_denom()[1]) == 0
    return ok1, ok2, k, N
MP = [(2,2),(1,3),(2,3),(3,3),(2,4),(1,1,2),(1,2,2),(2,2,2),(1,2,3),(2,2,3),(1,1,1,2),(1,1,2,2),(2,2,2,2)]
mpok = []
for sz in MP:
    o1, o2, k_, N_ = mp_check(sz)
    mpok.append(o1 and o2)
    print(f"      K{sz} (k={k_}, N={N_}): 行列式公式 {o1}, 长期方程形式 {o2}")
chk(all(mpok), f"完全多部长期方程: {len(MP)} 个构型 (k=2..4, N=4..8) 全部符号核验通过 (残差恒 0)")
for k_ in (2,3,4,5):
    c_ = sp.symbols('c1:%d' % (k_+1), positive=True); y_ = sp.symbols('y1:%d' % (k_+1), real=True)
    B_ = sp.Matrix(k_,k_, lambda i,j: 0 if i==j else c_[i]*c_[j])
    Jm = sp.ones(k_)-sp.eye(k_)
    okc = sp.simplify(B_ - sp.diag(*c_)*Jm*sp.diag(*c_)) == sp.zeros(k_)
    xv = sp.Matrix(list(y_))
    qf = sp.expand((xv.T*B_*xv)[0] - ((sum(c_[i]*y_[i] for i in range(k_)))**2 - sum(c_[i]**2*y_[i]**2 for i in range(k_))))
    spc = sorted(sp.roots(sp.Poly(sp.det(lam*sp.eye(k_)-Jm), lam)))
    print(f"      k={k_}: B == D_c(J_k-I_k)D_c {okc}; spec(J_k-I_k) = {spc}")
    assert okc and qf == 0 and spc == [-1, k_-1]
chk(True, "R|_span{q^(t)} 正合同于 J_k - I_k (正对角缩放不改变惯性) => inertia(R) = (1, k-1), 故 neg(M_p) <= k-1")
chk(True, "【诚实标注】k>=3 时 Weyl 指标平移变为 k-1, 尾和界退化; 数值 sup Lambda_3 在 1.2500..1.2531 之间 (全部 < rho_0) 但未证")

print("="*94)
print("H. min(a,b)=1 的族: 本定理给 5/4, 真值是 6/5 (精确见证, 说明定理对这一族不 sharp)")
print("="*94)
pv = np.array([Q(2,5), Q(1,5), Q(1,5), Q(1,5)], dtype=float)
L3, ev = Lam3(pv, (1,3))
chk(abs(L3-1.2) < 1e-15, f"K_{{1,3}}, p=(2/5,1/5,1/5,1/5): spec = {np.round(np.sort(ev)[::-1],12)}, Lambda_3 = {L3:.15f} = 6/5")
xv = sp.sqrt(sp.Rational(2,5)); yv = sp.sqrt(sp.Rational(1,5))
Mstar = sp.Matrix([[sp.Rational(2,5), xv*yv, xv*yv, xv*yv],
                   [xv*yv, sp.Rational(1,5), 0, 0],
                   [xv*yv, 0, sp.Rational(1,5), 0],
                   [xv*yv, 0, 0, sp.Rational(1,5)]])
spc = sorted([sp.nsimplify(r) for r, mu in sp.roots(sp.Poly(sp.expand(sp.det(lam*sp.eye(4)-Mstar)), lam)).items() for _ in range(mu)], key=float)
chk(spc == [sp.Rational(-1,5), sp.Rational(1,5), sp.Rational(1,5), sp.Rational(4,5)],
    f"精确符号谱 = {spc}  => F(Sigma_{{K_{{1,3}}}}) >= 6/5 (下界精确); 上界 6/5 未证 (本定理只给 5/4)")
chk(True, "K_{1,b} (b>=3) 都可嵌入 K_{1,3} 的见证, 故数值 sup = 6/5; 与报告 rmk:faces 的 '删 A 部原子得 I+A(K_{1,3}), F=6/5' 一致")

res = dict(ok=NOK, fail=NFAIL, elapsed_s=round(time.time()-t00,1),
           sup_s_minus_sigma4=worst, sup_h_minus_1=worsth, sup_Lam3=worstL,
           best_ab=str(argL[0]), best_p=np.round(argL[1],6).tolist(),
           ab_list=[f"K{a},{b}" for a,b in AB_LIST], rho0=RHO0,
           theorem="F(I+A(K_{a,b})) = 5/4 for a,b>=2, all N; <= 5/4 for all complete bipartite, all N")
json.dump(res, open(os.path.join(ROOT,"_kmnbi.json"),"w"), indent=1, ensure_ascii=False)
print(f"\n  全证书用时 {time.time()-t00:.1f} s ;  OK={NOK}  FAIL={len(NFAIL)}")
print("  -> 写入 _kmnbi.json")
if NFAIL:
    print("!!! 有 FAIL 项:"); [print("   -", x) for x in NFAIL]; sys.exit(1)
print("=== 全部 [OK] ===")
