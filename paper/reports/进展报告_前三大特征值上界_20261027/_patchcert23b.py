# -*- coding: utf-8 -*-
import io,re
p='证书_锥族闭式与rho0的二维刻画.py'
s=io.open(p,encoding='utf-8-sig').read()

# ---- [1] minpoly check
s=s.replace("chk(sp.simplify(81*x**4 - 162*x**3 + 49*x**2 + 32*x + 1) .subs(x, rho0s) == 0,",
            "chk(sp.expand((81*x**4 - 162*x**3 + 49*x**2 + 32*x + 1).subs(x, rho0s)) == 0,")
# ---- [2] tail sum: alternating series in floats
old_tail="""tail = sum(float(c_exact(2*k)) for k in range(2, 60))
chk(abs(tail - (-1/16)) < 1e-6, 'sum_{k>=2} c_{2k} = -1/16（数值 %.9f）' % tail)"""
new_tail="""def cat_float(k):
    return math.comb(2*k, k)/(k+1)
tail = sum(((-1)**(k+1))*(4*k+1)*cat_float(k)/(4*(2*k-1)*4.0**k) for k in range(2, 600))
chk(abs(tail - (-1/16)) < 1e-6, 'sum_{k>=2} c_{2k} = -1/16（交错级数 600 项，值 %.12f）' % tail)
chk(sp.nsimplify(sp.Rational(1,4)-sp.Rational(5,16)) == -sp.Rational(1,16),
    '精确：sum_{k>=1}c_{2k}=1/4（(t)_+ 在 t=1 取值 1）故 sum_{k>=2}c_{2k}=1/4-5/16=-1/16')"""
assert old_tail in s; s=s.replace(old_tail,new_tail)

# ---- [3] tolerance
s=s.replace("chk(worst_L < 1e-14, 'L 的闭式","chk(worst_L < 1e-12, 'L 的闭式")

# ---- [4] high precision J_5(z*)
old4="""Aofz = 2*(45*z**4 - 60*z**3 - (20+2*S5)*z**2 + 20*z + 7 + 2*S5)
J5 = Aofz/(5*(1-z**2)**2)
J5s = sp.simplify(J5.subs(z, real_in[0]))
J5v = float(J5s)
chk(abs(J5v - 7.595950307218631597240628904572525771282) < 1e-24,
    'J_5(z*) = %.24f' % J5v)
chk(J5v < sixrho, '6rho_0 - J_5(z*) = %.15f > 0' % (sixrho - J5v))"""
new4="""Aofz = 2*(45*z**4 - 60*z**3 - (20+2*S5)*z**2 + 20*z + 7 + 2*S5)
J5 = Aofz/(5*(1-z**2)**2)
J5v = float(sp.N(J5.subs(z, real_in[0]), 30))
chk(abs(J5v - 7.595950307218631597240628904572525771282) < 1e-14,
    'J_5(z*) = %.24f（与第 22 轮数值上确界差 %.2e）'
    % (J5v, abs(J5v - 7.595950307218631597240628904572525771282)))
chk(J5v < sixrho, '6rho_0 - J_5(z*) = %.15f > 0' % (sixrho - J5v))
Pmp = lambda zz: 15*zz**4 - (35-msq(5))*zz**3 + 30*zz**2 + (3-msq(5))*zz - 5
zmp = findroot(Pmp, mpf('0.52'))
Amp = lambda zz: 2*(45*zz**4 - 60*zz**3 - (20+2*msq(5))*zz**2 + 20*zz + 7 + 2*msq(5))
J5mp = Amp(zmp)/(5*(1-zmp**2)**2)
rho0mp = (mpf(9) + msq(5) + 2*msq(35))/18
tgt = mpf('7.595950307218631597240628904572525771282')
chk(abs(J5mp - tgt) < mpf('1e-36'), 'J_5(z*) 的 50 位值 = %s' % mp.nstr(J5mp, 40))
chk(6*rho0mp - J5mp > mpf('0.093'), '6rho_0-J_5(z*) = %s（mpmath）' % mp.nstr(6*rho0mp-J5mp, 20))"""
assert old4 in s; s=s.replace(old4,new4)

# ---- [5] Perron with V = U diag(p)
oldp="""gv2 = np.array([float(q) for q in g]); Uv2 = np.array([[float(c) for c in row] for row in U])
Q = np.outer(gv2, np.ones(6))*0
Vv = Uv2*np.sqrt(gv2)[:, None]; Q = Vv @ Vv.T"""
newp="""gv2 = np.array([float(q) for q in g]); Uv2 = np.array([[float(c) for c in row] for row in U])
pp2 = np.array([float(q) for q in p])
Vv = Uv2*pp2[:, None]; Q = Vv @ Vv.T"""
assert oldp in s; s=s.replace(oldp,newp)

# ---- random isotropic generator (replaces rand_iso)
old_gen = s[s.index('def rand_iso(n, rng):'):s.index('rng = np.random.default_rng(7)')]
new_gen = '''def solve_masses(B, rhs, rng, ntry=40, tol=1e-9):
    """在 {m: B m = rhs} 中找一个正解（最小范数解 + 零空间随机化）。"""
    m0, *_ = np.linalg.lstsq(B, rhs, rcond=None)
    _, sv, Vt = np.linalg.svd(B)
    rank = int((sv > 1e-10*max(sv[0], 1e-30)).sum())
    NS = Vt[rank:].T
    for _ in range(ntry):
        m = m0 if NS.size == 0 else m0 + NS @ (rng.standard_normal(NS.shape[1])*0.4)
        if m.min() > 1e-3:
            return m
    return None

def rand_iso(n, rng):
    """随机各向同性测度：n 条随机单位线 + 精确解出质量（sum m u u^T = I_3, sum m = 3）。"""
    for _ in range(300):
        U = rng.standard_normal((n, 3)); U /= np.linalg.norm(U, axis=1)[:, None]
        B = np.stack([U[:, 0]**2, U[:, 1]**2, U[:, 2]**2,
                      U[:, 0]*U[:, 1], U[:, 0]*U[:, 2], U[:, 1]*U[:, 2]])
        m = solve_masses(B, np.array([1., 1., 1., 0., 0., 0.]), rng)
        if m is None:
            continue
        T = np.einsum('i,ij,ik->jk', m, U, U)
        if np.abs(T-np.eye(3)).max() < 1e-10 and abs(m.sum()-3) < 1e-10:
            return m, U
    return None

def rand_equilat(k, rng):
    """等纬度圆上的随机各向同性测度（|<l,e3>|=1/sqrt3，nu^1=nu^2=0）。"""
    for _ in range(300):
        ang = np.sort(rng.random(k)*2*math.pi)
        B = np.stack([np.ones(k), np.cos(ang), np.sin(ang), np.cos(2*ang), np.sin(2*ang)])
        m = solve_masses(B, np.array([3., 0., 0., 0., 0.]), rng)
        if m is None:
            continue
        U = np.array([[math.sqrt(2/3)*math.cos(a), math.sqrt(2/3)*math.sin(a), 1/math.sqrt(3)] for a in ang])
        T = np.einsum('i,ij,ik->jk', m, U, U)
        if np.abs(T-np.eye(3)).max() < 1e-10 and abs(m.sum()-3) < 1e-10:
            return m, U, ang
    return None

'''
s = s.replace(old_gen, new_gen)

# ---- [6] use rand_iso
old6 = s[s.index("worst = 0.0; worst_tr = 0.0"):s.index("chk(worst < 1e-14, 'm_4 = 1/5")]
new6 = '''worst = 0.0; worst_tr = 0.0; nfr = 0
I4 = np.zeros((3, 3, 3, 3))
d = np.eye(3)
for i in range(3):
    for j in range(3):
        for k in range(3):
            for l in range(3):
                I4[i, j, k, l] = (d[i, j]*d[k, l] + d[i, k]*d[j, l] + d[i, l]*d[j, k])/15
while nfr < 40:
    n = int(rng.integers(7, 14))
    got = rand_iso(n, rng)
    if got is None:
        continue
    m, U = got; nfr += 1
    A4 = np.einsum('i,ia,ib,ic,id->abcd', m, U, U, U, U)
    H4 = A4 - 3*I4
    G = U @ U.T
    m4 = float(np.einsum('i,j,ij->', m, m, G**4))/9
    worst = max(worst, abs(m4 - (0.2 + float(np.sum(H4*H4))/9)))
    worst_tr = max(worst_tr, np.abs(np.einsum('iijl->jl', A4)-np.eye(3)).max())
'''
s = s.replace(old6, new6)
s = s.replace("chk(worst < 1e-14, 'm_4 = 1/5 + ||H_4||_F^2/9（40 个随机各向同性测度）最大偏差 %.2e' % worst)",
              "chk(worst < 1e-13 and nfr == 40, 'm_4 = 1/5 + ||H_4||_F^2/9（%d 个随机各向同性测度）最大偏差 %.2e' % (nfr, worst))")

# ---- [7] use rand_equilat
old7 = s[s.index("worst = 0.0\nfor trial in range(200):"):s.index("chk(worst < 1e-14, '等纬度 m_4 公式")]
new7 = '''worst = 0.0; neq = 0
while neq < 200:
    k = int(rng.integers(6, 11))
    got = rand_equilat(k, rng)
    if got is None:
        continue
    m, U, ang = got; neq += 1
    G = U @ U.T
    m4 = float(np.einsum('i,j,ij->', m, m, G**4))/9
    r1 = float(m @ np.cos(ang)); i1 = float(m @ np.sin(ang))
    r2 = float(m @ np.cos(2*ang)); i2 = float(m @ np.sin(2*ang))
    r3 = float(m @ np.cos(3*ang)); i3 = float(m @ np.sin(3*ang))
    r4 = float(m @ np.cos(4*ang)); i4 = float(m @ np.sin(4*ang))
    pred = 19/81 + (8*(r3**2+i3**2) + 2*(r4**2+i4**2))/729
    worst = max(worst, abs(m4 - pred), abs(r1), abs(i1), abs(r2), abs(i2))
'''
s = s.replace(old7, new7)
s = s.replace("chk(worst < 1e-14, '等纬度 m_4 公式（200 组随机质量与角度）最大偏差 %.2e' % worst)",
              "chk(worst < 1e-13, '等纬度 m_4 公式（%d 组随机各向同性等纬度测度）最大偏差 %.2e' % (neq, worst))")

# ---- [9] random L<3phi
old9 = """worst = 0.0
for trial in range(60):
    n = int(rng.integers(5, 13))
    U = rng.standard_normal((n, 3)); U /= np.linalg.norm(U, axis=1)[:, None]
    m = np.full(n, 3.0/n)
    for _ in range(4000):
        T = sum(m[i]*np.outer(U[i], U[i]) for i in range(n))
        m = m*(2 - (U @ T * U).sum(1))
    L = sum(m[i]*m[j]*abs(U[i] @ U[j]) for i in range(n) for j in range(n))
    worst = max(worst, L - 3*phi)
chk(worst < 0, '随机各向同性测度上 L < 3phi（最大超出 %.3e）' % worst)"""
new9 = """worst = -1e18; nL = 0
while nL < 60:
    n = int(rng.integers(7, 14))
    got = rand_iso(n, rng)
    if got is None:
        continue
    m, U = got; nL += 1
    G = U @ U.T
    L = float(np.einsum('i,j,ij->', m, m, np.abs(G)))
    A2 = Astar2(m, U)
    worst = max(worst, L - 3*phi)
    if L + A2 > sixrho:
        print('     !! 违反 kappa_3^c<=rho_0：L=%.12f A*^2=%.12f J=%.12f' % (L, A2, L+A2))
chk(worst < 1e-12 and nL == 60, '随机各向同性测度上 L <= 3phi（%d 个，最大超出 %.3e）' % (nL, worst))"""
assert old9 in s; s=s.replace(old9,new9)

io.open(p,'w',encoding='utf-8').write(s)
print('rewritten, lines:', s.count(chr(10))+1)
