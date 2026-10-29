# -*- coding: utf-8 -*-
# _r24c.py  ROUND 24 (C): the general-section identity  rho(Lambda) = (L_w + A_w^2)/2,  w_i = g_i h_i.
#   This generalises thm:r22id (centre section) to EVERY section mu, and it isolates exactly what
#   breaks when one tries to carry the R24 block machinery to the uniaxial section:
#       objective weights  w_i = g_i h_i          (h-tilted)
#       isotropy weights   g_i                     (bare),  sum g_i u_i u_i^T = I_3, sum g_i = 3
#   In the centre section h_i == 1/sqrt(3) is constant, so the two systems coincide up to scale.
import numpy as np, math, itertools
np.set_printoptions(precision=12, suppress=False, linewidth=200)
PASS = FAIL = 0
def chk(c, m):
    global PASS, FAIL
    if c: PASS += 1; print('  [PASS]', m)
    else: FAIL += 1; print('  [FAIL]', m)

def solve_masses(Bm, rhs, rng, ntry=80):
    m0, *_ = np.linalg.lstsq(Bm, rhs, rcond=None)
    _, sv, VT = np.linalg.svd(Bm); rank = int((sv > 1e-10*max(sv[0], 1e-30)).sum()); NS = VT[rank:].T
    for _ in range(ntry):
        m = m0 if NS.size == 0 else m0 + NS @ (rng.standard_normal(NS.shape[1])*0.4)
        if m.min() > 1e-3: return m
    return None

def rand_frame(n, rng):
    """g_i>0, sum g_i u_i u_i^T = I_3 (hence sum g_i = 3)."""
    for _ in range(400):
        U = rng.standard_normal((n, 3)); U /= np.linalg.norm(U, axis=1)[:, None]
        Bm = np.stack([U[:,0]**2, U[:,1]**2, U[:,2]**2, U[:,0]*U[:,1], U[:,0]*U[:,2], U[:,1]*U[:,2]])
        g = solve_masses(Bm, np.array([1.,1.,1.,0.,0.,0.]), rng)
        if g is None: continue
        T = np.einsum('i,ij,ik->jk', g, U, U)
        if np.abs(T-np.eye(3)).max() < 1e-11: return g, U
    return None

def rand_Lam(rng, kind='generic'):
    if kind == 'uniaxial':
        a = rng.uniform(0.3, 1.0); return np.diag([a, (1-a)/2, (1-a)/2])
    if kind == 'centre':
        return np.eye(3)/3
    if kind == 'rank1':
        w = rng.standard_normal(3); w /= np.linalg.norm(w); return np.outer(w, w)
    M = rng.standard_normal((3, 3)); S = M@M.T
    ev, Q = np.linalg.eigh(S); ev = np.abs(ev)+1e-3; ev = ev/ev.sum()
    return Q@np.diag(ev)@Q.T

def Amax(w, U):
    """A_w = max_{||z||=1} sum_i w_i |<u_i,z>|  ==  max_sigma ||sum sigma_i w_i u_i||."""
    n = len(w); best = -1e18; bz = None
    for sg in itertools.product([1,-1], repeat=n-1):
        s = np.array((1,)+sg); v = (s*w) @ U; t = float(v@v)
        if t > best: best = t; bz = v
    return math.sqrt(best), (bz/math.sqrt(best) if best > 1e-300 else np.array([0.,0.,1.]))

def bestsign(w, U):
    """max over sign flips of ||sum sigma_i w_i u_i||^2  (= A_w^2), and the argmax signs."""
    n = len(w); best = -1e18; bs = None
    for sg in itertools.product([1,-1], repeat=n-1):
        s = np.array((1,)+sg); v = (s*w) @ U; t = float(v@v)
        if t > best: best = t; bs = s
    return best, bs

def both(g, U, Lam):
    h = np.sqrt(np.maximum(np.einsum('ij,jk,ik->i', U, Lam, U), 0))
    V = U*np.sqrt(g)[:, None]; Q = V@V.T; Qp = (Q+np.abs(Q))/2
    x = np.sqrt(g)*h
    rho_dir = float(x @ Qp @ x)
    w = g*h
    Lw = float(w @ np.abs(U@U.T) @ w)
    aw = w @ U                             # sum_i w_i u_i for the GIVEN representatives
    alw = float(aw@aw)
    rho_fix = 0.5*(Lw+alw)                 # exact per-frame identity
    Aw2, bs = bestsign(w, U)
    rho_id = 0.5*(Lw+Aw2)                  # sign-optimised => the section value
    Us = U*bs[:, None]
    hs = np.sqrt(np.maximum(np.einsum('ij,jk,ik->i', Us, Lam, Us), 0))
    Vs = Us*np.sqrt(g)[:, None]; Qs = Vs@Vs.T; Qsp = (Qs+np.abs(Qs))/2
    rho_flip = float(x @ Qsp @ x)          # rho after applying the optimal sign flip
    lam_max = float(np.linalg.eigvalsh(Qp)[-1])
    return rho_dir, rho_id, lam_max, float(x@x), h, w, rho_fix, rho_flip, alw, Aw2

print('='*116)
print('(H) 任意截面的精确恒等式：rho_{g,U}(Lam) = (L_w + A_w^2)/2,  w_i = g_i h_i,  h_i = sqrt(u_i^T Lam u_i)')
print('='*116)
rng = np.random.default_rng(777)
worst = {'fix':0., 'flip':0., 'nrm':0., 'le':0.}
for kind in ('centre', 'uniaxial', 'rank1', 'generic'):
    wf = 0.; wp = 0.; n = 0
    for _ in range(120):
        nn = int(rng.integers(5, 12)); fr = rand_frame(nn, rng)
        if fr is None: continue
        g, U = fr; Lam = rand_Lam(rng, kind)
        rd, ri, lm, nx, h, w, rfix, rflip, alw, Aw2 = both(g, U, Lam)
        wf = max(wf, abs(rd-rfix)/max(rd, 1e-12))
        wp = max(wp, abs(rflip-ri)/max(ri, 1e-12))
        worst['nrm'] = max(worst['nrm'], abs(nx-np.trace(Lam)))
        worst['le'] = max(worst['le'], rd-lm)
        n += 1
    worst['fix'] = max(worst['fix'], wf); worst['flip'] = max(worst['flip'], wp)
    print('     %-9s : %d 个随机帧 x 随机 Lam；固定代表元恒等式偏差 %.3e，符号最优恒等式偏差 %.3e'%(kind, n, wf, wp))
chk(worst['fix'] < 1e-13, '逐帧精确恒等式 rho_{g,U}(Lam) = (L_w + ||sum w_i u_i||^2)/2（最大相对偏差 %.2e）'%worst['fix'])
chk(worst['flip'] < 1e-13, '符号最优形式 max_{sigma} rho_{g,sigma*U}(Lam) = (L_w + A_w^2)/2（最大相对偏差 %.2e）'%worst['flip'])
chk(worst['nrm'] < 1e-12, 'sum_i w_i^2/... 归一化核对：||x(Lam)||^2 = tr Lam（最大偏差 %.2e）'%worst['nrm'])
chk(worst['le'] < 1e-12, 'rho(Lam) <= lambda_max(Q_+)（最大超出 %.2e）'%worst['le'])

print()
print('='*116)
print('(I) 两个权重系统的分离：目标用 w_i=g_i h_i，各向同性用 g_i（这正是推广的障碍）')
print('='*116)
# centre: h constant => w proportional to g, and sum g = 3 reproduces eq:r22iso with m_i = g_i
fr = rand_frame(8, rng); g, U = fr
rd, ri, lm, nx, h, w, rfix, rflip, alw, Aw2 = both(g, U, np.eye(3)/3)
chk(abs(h-h[0]).max() < 1e-15 and abs(np.linalg.norm(w-g/np.sqrt(3))) < 1e-14,
    '中心截面：h_i == 1/sqrt(3) 常数，w = g/sqrt(3)，故 m_i := g_i 满足 sum m = 3 与 sum m l l^T = I_3')
print('     中心截面：rho=%.15f  (L_w+alpha_w^2)/2=%.15f  符号最优 (L_w+A_w^2)/2=%.15f'%(rd, rfix, ri))
# uniaxial: h varies => w is NOT proportional to g; sum_i (w_i/h_i) u_i u_i^T = I but sum w_i u_i u_i^T != I
a = 0.594608575866
Lam = np.diag([a, (1-a)/2, (1-a)/2])
rd, ri, lm, nx, h, w, rfix, rflip, alw, Aw2 = both(g, U, Lam)
Tg = np.einsum('i,ij,ik->jk', g, U, U)
Tw = np.einsum('i,ij,ik->jk', w, U, U)
print('     单轴截面 a=h_hub^2=%.12f：||T_g - I||_max = %.2e（各向同性对 g 成立）'%(a, np.abs(Tg-np.eye(3)).max()))
print('                                ||T_w - I||_max = %.2e（对 w 不成立 => 中心截面的矩计算不可直接搬运）'%np.abs(Tw-np.eye(3)).max())
chk(np.abs(Tg-np.eye(3)).max() < 1e-11 and np.abs(Tw-np.eye(3)).max() > 1e-3,
    '单轴截面下 w 权重不满足各向同性：R24 的分块机器需要新的归一化')
# the two exact constraints that DO survive for a general section
c1 = float(w@h); c2 = np.abs(np.einsum('i,ij,ik->jk', w/h, U, U)-np.eye(3)).max()
print('     仍然精确成立的两个约束：sum_i w_i h_i = tr Lam = %.15f ；||sum_i (w_i/h_i) u_i u_i^T - I||_max = %.2e'%(c1, c2))
chk(abs(c1-np.trace(Lam)) < 1e-13 and c2 < 1e-12, '一般截面的精确约束组：sum w_i h_i = tr Lam 且 sum (w_i/h_i)u_iu_i^T = I_3')

print()
print('='*116)
print('SUMMARY  PASS=%d  FAIL=%d'%(PASS, FAIL))
print('='*116)
raise SystemExit(0 if FAIL == 0 else 1)
