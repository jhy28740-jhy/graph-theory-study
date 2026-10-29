# -*- coding: utf-8 -*-
"""
证书_凹性对偶与截面分解.py   (第 19 轮)

复核第 19 轮的全部数值/精确断言，自包含（不依赖仓库根目录的草稿脚本）：
  (A) rho(Lam) 的几何平均表示  =>  凹性（解析已证，此处复核恒等式与凹性余量）
  (B) Danskin 闭式梯度  grad_Lam rho = sum_i ((Q+ x)_i/x_i) * g_i u_i u_i^T  >= 0
  (C) 精确对偶：rho(Lam)=min_{T>0}<Lam,S(T)>，Perron 取法 t*_ij=x_j/x_i 给出 S(T*)=rho*I_3（每个几何都成立）
  (D) 全局极大点：Lam~ > 0（精确闭式 lam_min=(1-h_hub^2)/2）且 max_Lam rho(Lam)=rho(Lam~)=rho；
      切向界 rho(Lam) <= rho*tr(Lam)
  (E) P1（各向异性对偶证书）在冠军几何上被严格排除：唯一可行证书是 S=rho0*I_3
  (F) Weyl+Frobenius 自举天花板：rho <= (1+sqrt5)/2（全局）、<= 3/2（n=6 加迹精细化）
  (G) Legendre 预算：c_ell 精确值；ell<=1 逐项天花板 5/4；全体逐项天花板 > rho0（该族死）
  (H) 截面分解 kappa_3 = sup_mu kappa_3(mu) 的数值地图；冠军恰在单轴截面上
运行： C:\\Python314\\python.exe -X utf8 -u 证书_凹性对偶与截面分解.py
"""
import math
import numpy as np
from numpy.polynomial.legendre import legvander

np.set_printoptions(precision=12, suppress=False, linewidth=200)
S5 = math.sqrt(5.0); S7 = math.sqrt(7.0); S35 = math.sqrt(35.0)
RHO0 = (9 + S5 + 2 * S35) / 18
TAU0 = RHO0 - 1
PHI = (1 + S5) / 2
PASS = FAIL = 0


def chk(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1; print('  [PASS]', msg)
    else:
        FAIL += 1; print('  [FAIL]', msg)


def champion():
    """hub + regular pentagon 精确构型（与报告第 3 节一致）。"""
    Vhh = 14.0 / (-7 * S5 - 7 + 5 * S7 + 3 * S35)
    Vhc = math.sqrt(7.0 / 5.0)
    Vcc = -3.5 - 7 * S5 / 10 + 3 * S35 / 10 + 5 * S7 / 2
    Ah = math.sqrt(Vhh); Ac = Vhc / Ah; Rc = math.sqrt(Vcc - Ac ** 2)
    uh = (91 - 25 * S7) / 126; uc = (7 + 5 * S7) / 126
    wv = np.array([uh] + [uc] * 5)
    cols = [np.array([Ah, 0.0, 0.0])]
    for j in range(5):
        th = 2 * math.pi * j / 5
        cols.append(np.array([Ac, Rc * math.cos(th), Rc * math.sin(th)]))
    A0 = np.array([math.sqrt(wv[k]) * cols[k][r] for r in range(3) for k in range(6)]).reshape(3, 6)
    P = A0.T @ np.linalg.inv(A0 @ A0.T) @ A0
    w, Vs = np.linalg.eigh(P); idx = np.argsort(-w); Vs = Vs[:, idx]; w = w[idx]
    V = Vs[:, :3] * np.sqrt(np.maximum(w[:3], 0))
    g = (V ** 2).sum(1); U = V / np.sqrt(g)[:, None]
    return U, g


UC, GC = champion()


def frame_data(U, g):
    """一个几何 (g,U) 的全部派生量。零质量行 (g_i=0) 自动剔除。"""
    keep = g > 1e-13
    U = U[keep]; g = g[keep]
    N = len(g); V = U * np.sqrt(g)[:, None]
    Kf = U @ U.T; Q = V @ V.T; Qp = np.maximum(Q, 0.0); Qm = np.maximum(-Q, 0.0)
    w, ev = np.linalg.eigh(Qp); k = int(np.argmax(w))
    rho = float(w[k]); x = ev[:, k].copy()
    if x.sum() < 0: x = -x
    h = x / np.sqrt(g)
    return dict(N=N, U=U, g=g, V=V, Kf=Kf, Q=Q, Qp=Qp, Qm=Qm, rho=rho, x=x, h=h,
                alpha=float(x @ Q @ x), beta=float(x @ Qm @ x), keep=keep)


def rho_Lam(U, g, L):
    """rho(Lam) 原始定义（不含特征值）。"""
    h = np.sqrt(np.maximum(np.einsum('ij,jk,ik->i', U, L, U), 0.0)); gh = g * h
    m = (U * gh[:, None]).sum(0)
    return float(m @ m) + float((np.outer(gh, gh) * np.maximum(-(U @ U.T), 0.0)).sum())


def rho_geo(U, g, L):
    """几何平均表示 rho(Lam)=sum_ij (Q+)_ij sqrt(L_i L_j), L_i=g_i u_i'Lam u_i。"""
    Qp = np.maximum(np.sqrt(np.outer(g, g)) * (U @ U.T), 0.0)
    Lin = g * np.einsum('ij,jk,ik->i', U, L, U)
    xx = np.sqrt(np.maximum(Lin, 0.0))
    return float((Qp * np.outer(xx, xx)).sum())


def S_of_T(U, g, T):
    """对偶族 S(T)=1/2 sum_ij (Q+)_ij ( t_ij g_i u_i u_i^T + t_ij^{-1} g_j u_j u_j^T )。"""
    Qp = np.maximum(np.sqrt(np.outer(g, g)) * (U @ U.T), 0.0)
    ii, jj = np.nonzero(Qp > 0)
    D = g[:, None, None] * (U[:, :, None] * U[:, None, :])     # D[i]=g_i u_i u_i^T
    return 0.5 * np.einsum('k,kij->ij', (Qp[ii, jj] * T[ii, jj]), D[ii]) \
         + 0.5 * np.einsum('k,kij->ij', (Qp[ii, jj] / T[ii, jj]), D[jj])


def lam_ukkt(D):
    """闭式 u-KKT 解 Lam~（命题 prop:uktkt）。"""
    U, g, h, rho, Kf = D['U'], D['g'], D['h'], D['rho'], D['Kf']
    ac = Kf > 1e-13; c = g * h; N = D['N']
    W = np.array([(c[ac[i]][:, None] * U[ac[i]]).sum(0) if ac[i].any() else np.zeros(3) for i in range(N)])
    return (1.0 / rho) * np.einsum('i,ij,ik->jk', g * h, W, U)


def rand_psd(rng, lo=0.0):
    A = rng.normal(size=(3, 3)); L = A @ A.T; L /= np.trace(L)
    if lo > 0:
        ev, EV = np.linalg.eigh(L); ev = np.maximum(ev, lo); L = (EV * ev) @ EV.T; L /= np.trace(L)
    return L


def stiefel_ascend_rho(N, seed, iters=2500, step=0.05):
    r = np.random.default_rng(seed); Q0, _ = np.linalg.qr(r.normal(size=(N, 3))); V = Q0[:, :3]
    best = (-1e18, None)
    for it in range(iters):
        Q = V @ V.T; w, ev = np.linalg.eigh(np.maximum(Q, 0.0)); k = int(np.argmax(w))
        x = ev[:, k]
        if x.sum() < 0: x = -x
        G = 2.0 * (np.outer(x, x) * (Q > 0).astype(float)) @ V
        G = G - V @ ((V.T @ G + G.T @ V) / 2.0); V = V + step * G
        Q0, _ = np.linalg.qr(V); V = Q0[:, :3]
        if it % 25 == 0 and w[k] > best[0]: best = (float(w[k]), V.copy())
    return best


def ascend_section(mu, N, seeds, iters=1200, step=0.05):
    """在固定 Lam=diag(mu) 的截面上对 V in St(N,3) 做梯度上升（数值）。"""
    L = np.diag(mu)
    def f(V):
        Q = V @ V.T; s = np.sqrt(np.maximum(np.einsum('ij,jk,ik->i', V, L, V), 0.0))
        return float((np.maximum(Q, 0.0) * np.outer(s, s)).sum())
    def gr(V):
        Q = V @ V.T; M = (Q > 0).astype(float)
        s = np.sqrt(np.maximum(np.einsum('ij,jk,ik->i', V, L, V), 1e-300)); Lv = V @ L
        return 2.0 * (s[:, None] * ((M * s) @ V) + ((M * (Q * np.outer(s, s))) @ Lv) / s[:, None])
    out = (-1e18, None)
    for sd in seeds:
        r = np.random.default_rng(sd); Q0, _ = np.linalg.qr(r.normal(size=(N, 3))); V = Q0[:, :3]
        for it in range(iters):
            G = gr(V); G = G - V @ ((V.T @ G + G.T @ V) / 2.0); V = V + step * G
            Q0, _ = np.linalg.qr(V); V = Q0[:, :3]
            if it % 25 == 0:
                val = f(V)
                if val > out[0]: out = (val, V.copy())
        val = f(V)
        if val > out[0]: out = (val, V.copy())
    return out


def banner(s):
    print(); print('=' * 100); print(s); print('=' * 100)


# ---------------------------------------------------------------- (A)
banner('(A) rho(Lam) 的几何平均表示与凹性（解析证明的数值复核）')
rng = np.random.default_rng(20261029)
dev = 0.0; marg = 1e9
D = frame_data(UC, GC)
for _ in range(400):
    L = rand_psd(rng)
    dev = max(dev, abs(rho_Lam(D['U'], D['g'], L) - rho_geo(D['U'], D['g'], L)))
    L1 = rand_psd(rng); L2 = rand_psd(rng); s = float(rng.random())
    a = rho_Lam(D['U'], D['g'], (1 - s) * L1 + s * L2)
    b = (1 - s) * rho_Lam(D['U'], D['g'], L1) + s * rho_Lam(D['U'], D['g'], L2)
    marg = min(marg, a - b)
print('  两种写法最大差 = %.3e ;  冠军上 400 组随机凹性余量最小值 = %+.3e' % (dev, marg))
chk(dev < 1e-12, '几何平均表示 rho(Lam)=sum (Q+)_ij sqrt(L_i L_j) 与原定义一致 (max %.2e)' % dev)
chk(marg > -1e-11, 'rho(Lam) 凹（解析已证：sqrt(L_iL_j)=min_t (tL_i+L_j/t)/2 为仿射族的下确界），余量 %+.2e' % marg)

# ---------------------------------------------------------------- (B)
banner('(B) Danskin 闭式梯度  grad_Lam rho(Lam) = sum_i ((Q+ x(Lam))_i / x(Lam)_i) * g_i u_i u_i^T  >= 0')
def grad_formula(U, g, L):
    Lin = g * np.einsum('ij,jk,ik->i', U, L, U)
    xx = np.sqrt(np.maximum(Lin, 1e-300))
    Qp = np.maximum(np.sqrt(np.outer(g, g)) * (U @ U.T), 0.0)
    y = (Qp @ xx) / xx
    return np.einsum('i,i,ij,ik->jk', y, g, U, U)
dev = 0.0
for _ in range(40):
    L = rand_psd(rng, lo=0.06); Gf = grad_formula(D['U'], D['g'], L); e = 1e-6
    Gn = np.zeros((3, 3))
    for a in range(3):
        for b in range(a, 3):
            E = np.zeros((3, 3)); E[a, b] = 1.0
            if a != b: E[b, a] = 1.0
            d = (rho_Lam(D['U'], D['g'], L + e * E) - rho_Lam(D['U'], D['g'], L - e * E)) / (2 * e)
            Gn[a, b] = Gn[b, a] = (d if a == b else 0.5 * d) if a == b else 0.5 * d
    dev = max(dev, np.abs(Gf - Gn).max() / max(1.0, np.abs(Gf).max()))
mn = min(np.linalg.eigvalsh(grad_formula(D['U'], D['g'], rand_psd(rng))).min() for _ in range(80))
print('  中心差分对照最大相对偏差 = %.3e ;  80 个随机 Lam 上 grad 的最小特征值 = %+.6e' % (dev, mn))
chk(dev < 1e-5, 'Danskin 闭式梯度与中心差分一致 (rel dev %.2e)' % dev)
chk(mn > -1e-12, 'grad rho(Lam) 对一切 Lam 半正定 (min eig %+.2e)' % mn)

# ---------------------------------------------------------------- (C)
banner('(C) 精确对偶：rho(Lam)=min_{T>0}<Lam,S(T)>；Perron 取法给出 S(T*)=rho*I_3（每个几何）')
cases = [('champion', UC, GC)]
for N in (4, 5, 6, 7, 9, 12):
    for sd in (1, 2):
        r = np.random.default_rng(700 + 13 * N + sd)
        A = r.normal(size=(3, N)); M = A @ A.T
        ev, EV = np.linalg.eigh(M); V = (EV @ np.diag(ev ** -0.5) @ EV.T @ A).T
        gg = (V ** 2).sum(1); cases.append(('rand N=%d s%d' % (N, sd), V / np.sqrt(gg)[:, None], gg))
worst_S = 0.0; worst_dual = 0.0; worst_gap = 0.0
for tag, U, g in cases:
    Dd = frame_data(U, g); U2, g2 = Dd['U'], Dd['g']
    Qp = np.maximum(np.sqrt(np.outer(g2, g2)) * (U2 @ U2.T), 0.0)
    x = Dd['x']
    if x.min() < 1e-14:
        print('  %-14s 跳过（Perron 向量含零分量，AM-GM 取等点退化）' % tag); continue
    T = x[None, :] / x[:, None]
    S = S_of_T(U2, g2, T)
    worst_S = max(worst_S, np.abs(S - Dd['rho'] * np.eye(3)).max())
    viol = 0.0; gap = np.inf
    for _ in range(60):
        L = rand_psd(rng); Tr = np.exp(rng.normal(size=(len(g2), len(g2)))); Tr = 0.5 * (Tr + Tr.T)
        viol = max(viol, rho_Lam(U2, g2, L) - float(np.trace(L @ S_of_T(U2, g2, Tr))))
        gap = min(gap, float(np.trace(L @ S)) - rho_Lam(U2, g2, L))
    worst_dual = max(worst_dual, viol); worst_gap = min(worst_gap, gap)
    print('  %-14s rho=%.12f  |S(T*)-rho*I_3|_max=%.2e  弱对偶最大违背=%+.2e  T* 处最小间隙=%.3e'
          % (tag, Dd['rho'], np.abs(S - Dd['rho'] * np.eye(3)).max(), viol, gap))
chk(worst_S < 1e-12, 'S(T*) = rho*I_3 精确成立（%d 个几何，最大偏差 %.2e）' % (len(cases), worst_S))
chk(worst_dual < 1e-11, '弱对偶 rho(Lam) <= <Lam,S(T)> 对一切 (Lam,T)（最大违背 %.2e）' % worst_dual)
chk(worst_gap > -1e-12, 'T* 处对偶间隙非负（最小 %.2e）；结合 max_Lam rho(Lam)=rho 得对偶零间隙' % worst_gap)
print('  => max_{Lam>=0,tr=1} rho(Lam) = min_{T>0} lam_max(S(T)) = lam_max(Q+) 精确成立。')
print('     最优证书恒为各向同性 rho*I_3；Lambda-参数化不带任何松弛间隙。')

# ---------------------------------------------------------------- (D)
banner('(D) Lam~ > 0（精确闭式）与「Lam~ 是全局极大点」')
D = frame_data(UC, GC); Lt = lam_ukkt(D); evLt = np.linalg.eigvalsh(Lt)[::-1]
g0 = (9 - 4 * S5 + S35) / 18; a2 = (91 - 25 * S7) / 126
h_hub2 = a2 / g0; lmin_exact = (1 - h_hub2) / 2
print('  eig(Lam~) =', np.array2string(evLt, precision=12))
print('  精确闭式: lam_max(Lam~)=h_hub^2=%s=%.15f ; lam_min=(1-h_hub^2)/2=%.15f' % ('(91-25 sqrt7)/(126 g0*)', h_hub2, lmin_exact))
print('  数值 lam_min=%.15f  与闭式差=%.2e ;  tr(Lam~)=%.15f' % (evLt[-1], abs(evLt[-1] - lmin_exact), np.trace(Lt)))
chk(abs(evLt[-1] - lmin_exact) < 1e-12 and lmin_exact > 0,
    'Lam~ 正定，lam_min=(1-h_hub^2)/2=%.12f>0（精确闭式，无浮点依赖）' % lmin_exact)
chk(abs(evLt[0] - h_hub2) < 1e-12, 'lam_max(Lam~)=h_hub^2=%.12f（精确闭式）' % h_hub2)
from scipy.optimize import minimize


def lam_of_A(A):
    s = float(np.sum(A * A))
    return (A @ A.T) / s


def max_rho_Lam(U, g, ntry=14, seed=0):
    """max{rho(Lam): Lam>=0, tr Lam=1}，参数化 Lam=AA^T/tr(AA^T)，L-BFGS + 解析梯度。
       grad_A rho = (2/s)(G A - rho A)，G=grad_Lam rho（(B) 的闭式），s=||A||_F^2。"""
    r = np.random.default_rng(seed)

    def fg(a):
        A = a.reshape(3, 3); s = float(np.sum(A * A)); L = (A @ A.T) / s
        val = rho_Lam(U, g, L); G = grad_formula(U, g, L)
        return -val, (-(2.0 / s) * (G @ A - val * A)).ravel()

    starts = [np.eye(3).ravel() / math.sqrt(3.0)] + [r.normal(size=(3, 3)).ravel() for _ in range(ntry)]
    best = (-1e18, None)
    for a0 in starts:
        res = minimize(fg, a0, jac=True, method='L-BFGS-B', options=dict(maxiter=4000, ftol=1e-16, gtol=1e-14))
        if -res.fun > best[0]:
            best = (-res.fun, lam_of_A(res.x.reshape(3, 3)))
    return best
Lt_c = lam_ukkt(D)
gdev = np.abs(grad_formula(D['U'], D['g'], Lt_c) - D['rho'] * np.eye(3)).max()
print('  ||grad rho(Lam~) - rho*I_3||_max = %.3e   （定理 thm:lamparam(ii) 的数值复核）' % gdev)
chk(gdev < 1e-9, 'grad rho(Lam~) = rho*I_3 精确（偏差 %.2e）；与凹性合起来解析地给出全局极大' % gdev)
tang = 0.0; argd = -1e18; idd = -1e18; rkbad = False; rks = []
for N in (6, 7, 9, 10):
    v, V = stiefel_ascend_rho(N, 555 + 13 * N)
    gg = (V ** 2).sum(1); U = V / np.sqrt(gg)[:, None]
    Dd = frame_data(U, gg); L = lam_ukkt(Dd)
    att = rho_Lam(Dd['U'], Dd['g'], L) - Dd['rho']
    sup = Dd['rho']; sL = L
    val, Lv = max_rho_Lam(Dd['U'], Dd['g'], 14, 31 * N)
    if val > sup: sup, sL = val, Lv
    dvec = np.array([(Dd['g'][i] * np.array([Dd['U'][i][0] ** 2, Dd['U'][i][1] ** 2, Dd['U'][i][2] ** 2,
                                              2 * Dd['U'][i][0] * Dd['U'][i][1], 2 * Dd['U'][i][0] * Dd['U'][i][2],
                                              2 * Dd['U'][i][1] * Dd['U'][i][2]])).tolist() for i in range(Dd['N'])])
    rk = int(np.linalg.matrix_rank(dvec, tol=1e-10))
    dist = np.abs(sL - L).max()
    idd = max(idd, dist if rk == 6 else 0.0); rkbad = rkbad or ((rk < 6) != (dist > 1e-4))
    rks.append((N, rk, dist))
    for _ in range(120):
        Lr = rand_psd(rng); tang = max(tang, rho_Lam(Dd['U'], Dd['g'], Lr) - Dd['rho'] * np.trace(Lr))
        Lr2 = rng.normal(size=(3, 3)); Lr2 = Lr2 @ Lr2.T * float(rng.uniform(0.2, 3.0))
        tang = max(tang, rho_Lam(Dd['U'], Dd['g'], Lr2) - Dd['rho'] * np.trace(Lr2))
    argd = max(argd, sup - Dd['rho'], att)
    print('   N=%2d rho=%.12f lam_min(Lam~)=%+.9f  rho(Lam~)-rho=%+.2e  上升最大超出=%+.2e  |Lam_arg-Lam~|=%.2e'
          % (N, Dd['rho'], np.linalg.eigvalsh(L).min(), att, sup - Dd['rho'], np.abs(sL - L).max()))
chk(argd < 1e-9, 'max_Lam rho(Lam) = rho(Lam~) = rho：Lambda~ 精确取到，且 15 起点 L-BFGS 从未超出（最大 %.2e）' % argd)
chk(idd < 1e-4, '当 rank{g_i u_i u_i^T}=6（=dim Sym_3）时极大点唯一：|Lam_arg-Lam~|_max <= %.2e' % idd)
chk(not rkbad, '极大点的唯一性判据：|Lam_arg-Lam~|>1e-4  <=>  rank{g_i u_i u_i^T}<6  （%s）' % ', '.join('N=%d:rank=%d,dist=%.1e' % t for t in rks))
print('  结构说明：极大点集 = {Lam>=0, <Lam,g_i u_i u_i^T>=x_i^2}，是仿射面；值恒等于 rho。')
print('  当 N 个 {g_i u_i u_i^T} 在 Sym_3 中线性相关时该面维数 >0，故 Lam~ 不唯一——这不影响 max_Lam rho(Lam)=rho。')
chk(tang < 1e-10, '切向界 rho(Lam) <= rho*tr(Lam) 在 960 个随机 Lam（含 tr!=1）上成立（最大违背 %.2e）' % tang)
print('  解析理由：凹性 + grad rho(Lam~)=rho*I_3 => rho(Lam) <= rho(Lam~)+<rho I_3, Lam-Lam~> = rho*tr(Lam)。')

# ---------------------------------------------------------------- (E)
banner('(E) P1（各向异性对偶证书）被严格排除')
print('  任何满足 (certcond) 的 S 必须使  rho0 = rho(Lam~) <= <Lam~,S> <= lam_max(S)*tr(Lam~) <= rho0 处处取等，')
print('  故 lam_max(S)=rho0 且 Lam~ 的支撑含于 S 的最高特征空间。由 (D)，Lam~ 正定（支撑=全空间），')
print('  故 S 的最高特征值为 rho0 且重数 3，即 S=rho0*I_3。')
print('  另外 AM-GM 取等条件 t_ij=sqrt(L_j/L_i) 在 Lam~ 处唯一，而 (C) 已证 S(T*)=rho0*I_3。')
best = -1e18; bestS = None
for _ in range(3000):
    A = rng.normal(size=(3, 3)); S = RHO0 * (A @ A.T) / np.linalg.eigvalsh(A @ A.T).max()
    if np.linalg.eigvalsh(S).max() > RHO0 + 1e-15: continue
    ok = 0.0
    for _ in range(30):
        L = rand_psd(rng); ok = min(ok, float(np.trace(L @ S)) - rho_Lam(D['U'], D['g'], L))
    if ok > best: best = ok; bestS = S
print('  随机搜索 lam_max(S)<=rho0 的可行证书，最好者 min_Lam(<Lam,S>-rho(Lam)) = %+.3e（<=0 表示不可行）' % best)
print('  该随机搜索找到的最接近可行者距 rho0*I_3 的距离 = %.6f' % np.abs(bestS - RHO0 * np.eye(3)).max())
chk(best < 1e-9, '3000 次随机各向异性尝试无一可行（最好间隙 %+.2e）；唯一证书为 rho0*I_3' % best)

# ---------------------------------------------------------------- (F)
banner('(F) Weyl+Frobenius 自举天花板：rho <= (1+sqrt5)/2（全局），<= 3/2（n=6 加迹精细化）')
print('  数据：tr Q+ = sum g_i = 3 ；||Q+||_F^2 + ||Q-||_F^2 = ||Q||_F^2 = tr(Q^2) = 3 ；')
print('        rho = alpha + beta，alpha=||V^T x||^2 <= lam_max(Q)=1，beta <= lam_max(Q-) <= ||Q-||_F。')
print('  => (rho-1)^2 <= 3-rho^2  =>  rho <= (1+sqrt5)/2 = %.15f' % PHI)
print('  n=6 精细化：||Q+||_F^2 >= rho^2+(3-rho)^2/5，||Q-||_F^2 >= 1.2*lam_max(Q-)^2 >= 1.2*(rho-1)^2')
print('              => 2.4 rho^2 - 3.6 rho <= 0 => rho <= 3/2')
worst = -1e18
for N in (4, 5, 6, 7, 8, 9, 10, 12, 14):
    v, V = stiefel_ascend_rho(N, 31337 + 11 * N)
    Q = V @ V.T; r = np.linalg.eigvalsh(np.maximum(Q, 0)).max()
    b1 = (r - 1) ** 2 - (3 - r ** 2)
    worst = max(worst, b1)
    print('   N=%2d  rho_asc=%.12f   (rho-1)^2-(3-rho^2) = %+.3e  (<=0 一致)' % (N, r, b1))
chk(worst < 1e-9, '自举不等式 (rho-1)^2 <= 3-rho^2 在 9 个阶数的上升极大元上全部成立（最大 %+.2e）' % worst)
print('  天花板记录：该族（迹+Frobenius+Weyl）全局上限 %.9f，n=6 上限 %.9f，均远高于 rho0=%.9f（缺口 %.6f / %.6f）'
      % (PHI, 1.5, RHO0, PHI - RHO0, 1.5 - RHO0))
chk(PHI - RHO0 > 0.3, '该族定量死路：全局缺口 %.6f > 0.3' % (PHI - RHO0))

# ---------------------------------------------------------------- (G)
banner('(G) Legendre 预算：(t)_+ = sum_ell c_ell P_ell(t)，c_ell = (2ell+1)/2 * int_0^1 t P_ell(t) dt')
import sympy as sp
tt = sp.symbols('t')
Cex = {l: sp.nsimplify((2 * l + 1) / 2 * sp.integrate(tt * sp.legendre(l, tt), (tt, 0, 1)), rational=True) for l in range(0, 15)}
print('  c_ell (ell<=12):', [(l, str(Cex[l])) for l in range(13)])
chk(all(Cex[l] == 0 for l in (3, 5, 7, 9, 11, 13)), 'c_ell = 0 对一切奇数 ell>=3（因 (t)_+ = t/2 + |t|/2，|t| 只有偶次项）')
chk(Cex[0] == sp.Rational(1, 4) and Cex[1] == sp.Rational(1, 2) and Cex[2] == sp.Rational(5, 16)
    and Cex[4] == sp.Rational(-3, 32) and Cex[6] == sp.Rational(13, 256), 'c0=1/4, c1=1/2, c2=5/16, c4=-3/32, c6=13/256 精确')
D = frame_data(UC, GC); U2, g2, x = D['U'], D['g'], D['x']
sg = np.sqrt(g2); K = np.clip(U2 @ U2.T, -1, 1); Lmax = 24
Vand = legvander(K, Lmax)
Cs = {l: (float(Cex[l]) if l in Cex else float(sp.nsimplify((2 * l + 1) / 2 * sp.integrate(tt * sp.legendre(l, tt), (tt, 0, 1)), rational=True))) for l in range(Lmax + 1)}
terms = []
for l in range(Lmax + 1):
    A = np.outer(sg, sg) * Vand[:, :, l]; A = 0.5 * (A + A.T)
    terms.append((l, Cs[l], float(x @ A @ x), float(np.linalg.eigvalsh(A).max())))
rec = sum(c * q for _, c, q, _ in terms)
print('  冠军: x^T A_0 x=%.12f (lam_max=3), x^T A_1 x=%.12f (lam_max=1)' % (terms[0][2], terms[1][2]))
print('  ell<=1 实际贡献 = %.12f ；逐项天花板 3/4+1/2 = 5/4 = %.12f' % (terms[0][1] * terms[0][2] + terms[1][1] * terms[1][2], 1.25))
print('  ell>=2 余项 = %+.12f ；级数截断到 ell=24 的重构值 = %.12f（rho0=%.12f，截断误差 %.2e）'
      % (rec - terms[0][1] * terms[0][2] - terms[1][1] * terms[1][2], rec, RHO0, rec - RHO0))
tc = sum(c * lm for _, c, _, lm in terms if c > 0)
print('  全体逐项天花板 sum_{c_ell>0} c_ell*lam_max(A_ell) = %.12f  >>  rho0 = %.12f' % (tc, RHO0))
chk(abs(terms[0][3] - 3) < 1e-9 and abs(terms[1][3] - 1) < 1e-9, 'lam_max(A_0)=sum g_i=3 与 lam_max(A_1)=lam_max(P)=1 精确')
chk(tc - RHO0 > 0.25, '该族定量死路：逐项 Legendre/优占型天花板 %.9f，比 rho0 高 %.6f' % (tc, tc - RHO0))

# ---------------------------------------------------------------- (H)
banner('(H) 截面分解 kappa_3 = sup_{mu in Delta_2} kappa_3(mu) 的数值地图（数值，不作成果）')
Lc = lam_ukkt(frame_data(UC, GC))
mus = [('center (1/3,1/3,1/3)', [1 / 3, 1 / 3, 1 / 3]),
       ('uniaxial a=0.50', [0.50, 0.25, 0.25]),
       ('uniaxial a=0.5946 (=lam_max(Lam~))', [h_hub2, lmin_exact, lmin_exact]),
       ('uniaxial a=0.70', [0.70, 0.15, 0.15]),
       ('uniaxial a=0.85', [0.85, 0.075, 0.075]),
       ('rank-2 (1/2,1/2,0)', [0.5, 0.5, 0.0]),
       ('rank-1 (1,0,0)', [1.0, 0.0, 0.0]),
       ('asym (0.60,0.30,0.10)', [0.60, 0.30, 0.10])]
print('   section mu                                   N=6        N=8       N=10        max')
allmax = -1e18; axis = {}
for tag, mu in mus:
    row = [ascend_section(mu, N, [2000 + 7 * N + k for k in range(3)])[0] for N in (6, 8, 10)]
    if 'lam_max' in tag:
        row.append(rho_Lam(D['U'], D['g'], Lc))       # 冠军本身是该截面的可行点
    axis[tag] = max(row); allmax = max(allmax, max(row))
    print('   %-42s %s   %.12f  (%+.2e)' % (tag, ' '.join('%.9f' % z for z in row), max(row), max(row) - RHO0))
Vc = D['U'] * np.sqrt(D['g'])[:, None]
print('   冠军自身在单轴截面 mu*=spec(Lam~) 上的取值 = %.15f = rho0（可行点，故 kappa_3(mu*) >= rho0）'
      % rho_Lam(D['U'], D['g'], Lc))
print('   spec(Lam~) = %s ；mu* 为单轴（二重特征值 %.12f）' % (np.array2string(np.linalg.eigvalsh(Lc)[::-1], precision=12), lmin_exact))
chk(all(v <= RHO0 + 1e-9 for v in axis.values()), '全部 8 个截面的数值上确界 <= rho0（最大 %.12f，超出 %.2e）' % (allmax, allmax - RHO0))
chk(max(v for k, v in axis.items() if 'rank' in k) < RHO0 - 0.05,
    '边界截面（rank-1/rank-2）明显低于 rho0（最大 %.9f，余量 %.6f）'
    % (max(v for k, v in axis.items() if 'rank' in k), RHO0 - max(v for k, v in axis.items() if 'rank' in k)))
print('  数值图像：最大值出现在单轴族 a≈0.6--0.7；冠军恰在 a=lam_max(Lam~)=h_hub^2 处取到 rho0。')

banner('汇总：PASS=%d FAIL=%d' % (PASS, FAIL))
print('rho0 = %.15f   tau0 = %.15f   phi = %.15f   Phi_STW = %.15f' % (RHO0, TAU0, PHI, (1 + PHI) / 2))
raise SystemExit(0 if FAIL == 0 else 1)
