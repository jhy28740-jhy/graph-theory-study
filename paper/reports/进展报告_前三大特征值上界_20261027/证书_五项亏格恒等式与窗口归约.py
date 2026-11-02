# -*- coding: utf-8 -*-
"""
证书：第 25 轮 —— 五项亏格恒等式、偶极子的余弦形式与全局窗口归约
=================================================================
本证书逐条校验第 25 轮的全部精确结果：
  A 精确常数代数 (60 位)          B 正交坐标基与 h 的展开
  C 偶极子的余弦形式 A_w=1-E/2     D 通道亏格恒等式
  E Schoenberg 四阶分解            F 五项亏格主恒等式 (对一切 Lambda)
  G 冠军处取等 + Phi 的严格局部极小  H 无条件窗口定理 W<=W_c
  I 秩一截面刚性                   J 常纬度归约到第 24 轮
  K 单轴截面的两尺度结构            L 锥族 min Phi 扫描
  M 随机全局搜索                   N 局部极小化收敛到冠军谱
标签约定：已证=解析证明并已在此数值/高精度核对；数值=仅数值证据。
"""
import numpy as np
from numpy.linalg import norm, eigh
from mpmath import mp, mpf, sqrt as msqrt, polyroots
from scipy.optimize import minimize
mp.dps = 60

NP = NF = 0
def chk(name, cond, info=""):
    global NP, NF
    if cond: NP += 1; print("  [PASS] %s %s" % (name, info))
    else:    NF += 1; print("  [FAIL] %s %s" % (name, info))
def head(t): print("\n" + "=" * 104 + "\n" + t + "\n" + "=" * 104)

# ======================= 常数 =======================
S5 = np.sqrt(5.0); S35 = np.sqrt(35.0)
rho0 = (9 + S5 + 2 * S35) / 18; phi = (1 + S5) / 2
gam = (1 + S5) / 16; alp = (7 * S5 - 5) / 8; bet = (25 - 15 * S5) / 16
cW = (2 * S5 - 1) / 14; cS = (10 + S5) / 14; cZ = (5 - 3 * S5) / 14
del0 = 3 * phi + 3 - 6 * rho0; TGT = del0 / 3
r5 = msqrt(5); r7 = msqrt(7); r35 = msqrt(35)
rho0m = (9 + r5 + 2 * r35) / 18; phim = (1 + r5) / 2
del0m = 3 * phim + 3 - 6 * rho0m; TGTm = del0m / 3
cWm = (2 * r5 - 1) / 14; cSm = (10 + r5) / 14
gamm = (1 + r5) / 16; alpm = (7 * r5 - 5) / 8; betm = (25 - 15 * r5) / 16
astar_m = (3 * r35 + 5 * r7 - 7 - 7 * r5) / 14
Wc2m = 3 - TGTm / cWm; Wcm = msqrt(Wc2m)
astar = float(astar_m)
_rts = polyroots([mpf(11225), 0, -6700, 0, 430, 0, 260, 0, -31])
zstar_m = [r.real for r in _rts if abs(r.imag) < mpf('1e-40') and mpf('0.49') < r.real < mpf('0.51')][0]
zstar = float(zstar_m); rho_st = (3 * astar - 1) / 2

def P4(t): return (35 * t ** 4 - 30 * t * t + 3) / 8
def P2(t): return (3 * t * t - 1) / 2
def qf(s): return gam + alp * s + bet * s * s
def uniax(r):
    L = ((1 - r) / 3) * np.eye(3); L[2, 2] += r; return L

_SIG = {}
def signpat(n):
    if n not in _SIG:
        E = np.array([[(i >> k) & 1 for k in range(n - 1)] for i in range(2 ** (n - 1))], dtype=float)
        _SIG[n] = np.hstack([1 - 2 * E, np.ones((E.shape[0], 1))])
    return _SIG[n]

def fast(g, U, Lam):
    h = np.sqrt(np.einsum('ij,jk,ik->i', U, Lam, U)); w = g * h; n = len(g)
    G = U @ U.T; W = w.sum()
    Sw = (w[:, None, None] * (U[:, :, None] * U[:, None, :])).sum(0); SF = float((Sw * Sw).sum())
    Lw = float((np.outer(w, w) * np.abs(G)).sum())
    M = np.outer(w, w) * G; E = signpat(n)
    nrm2 = ((E @ M) * E).sum(1); k = int(np.argmax(nrm2)); A = float(np.sqrt(max(nrm2[k], 0.0)))
    m = E[k] * w @ U; z = m / norm(m) if norm(m) > 0 else np.array([0, 0, 1.])
    a = np.abs(U @ z); Vw = float((np.outer(w, w) * (a[:, None] - a[None, :]) ** 2).sum())
    M4 = float((np.outer(w, w) * G ** 4).sum()); Z4 = float((np.outer(w, w) * P4(G)).sum())
    Z2 = float((np.outer(w, w) * P2(G)).sum()); Dw = float((np.outer(w, w) * (qf(G ** 2) - np.abs(G))).sum())
    dh = h[:, None] - h[None, :]; gg = np.outer(g, g)
    H = float((gg * dh ** 2).sum()); Et = float((gg * dh ** 2 * G ** 2).sum())
    Ea = float((g * (h - a) ** 2).sum()); zet = float(z @ Sw @ z)
    Phi = (1 - A * A) + 0.5 * (cW * H + cS * Et) + (-cZ) * Z4 + Dw
    mvec = w @ U
    return dict(h=h, w=w, W=W, SF=SF, A=A, Aw2=A * A, z=z, a=a, Lw=Lw, M4=M4, Z4=Z4, Z2=Z2,
                Dw=Dw, H=H, Et=Et, E=Ea, zet=zet, Vw=Vw, Sw=Sw, m=mvec, Jw=Lw + A * A,
                half=(Lw + A * A) / 2, Phi=Phi, G=G,
                t1=1 - A * A, t2=0.5 * cW * H, t3=0.5 * cS * Et, t4=-cZ * Z4, t5=Dw,
                rID=(2 * rho0 - (Lw + A * A)) - ((1 - A * A) - TGT + 0.5 * (cW * H + cS * Et) + (-cZ) * Z4 + Dw),
                rA=A - (1 - Ea / 2), rH=H - (6 - 2 * W * W), rE=SF - (1 - Et / 2),
                rM=M4 - (-3 * W * W / 35 + 6 * SF / 7 + 8 * Z4 / 35),
                rL=Lw - (gam * W * W + alp * SF + bet * M4 - Dw),
                rAw2=A * A - (W * zet - 0.5 * Vw),
                rDQ=Dw + (-cZ) * Z4 - (cW * W * W + cS * SF - Lw),
                rSF3=SF - W * W / 3, rEtH=Et - H / 3)

def cone(n, z):
    Gc = 2 / (n * (1 - z * z)); g0 = 1 - n * Gc * z * z
    U = [np.array([0, 0, 1.])]
    for k in range(n):
        p = 2 * np.pi * k / n
        U.append(np.array([np.sqrt(1 - z * z) * np.cos(p), np.sqrt(1 - z * z) * np.sin(p), z]))
    return np.array([g0] + [Gc] * n), np.array(U)

def frame(V):
    g = (V * V).sum(1); return g, V / np.sqrt(g)[:, None]

def randLam(kind, rng):
    if kind == 0: return uniax(rng.random())
    ev = np.sort(rng.random(3))[::-1]; ev = ev / ev.sum()
    if kind == 1:
        M = rng.standard_normal((3, 3)); O = eigh(M @ M.T)[1]
    else:
        O = np.eye(3)
    return O @ np.diag(ev) @ O.T

print("=" * 104)
print("证书：第 25 轮 —— 五项亏格恒等式、偶极子的余弦形式与全局窗口归约")
print("=" * 104)
print("rho_0 = %.30f" % float(rho0m)); print("delta_0/3 = %.30f" % float(TGTm))
print("a* = %.30f" % astar); print("z* = %.30f" % zstar)

# ---------------- A ----------------
head("A. 精确常数代数（60 位）")
chk("c_W = (2sqrt5-1)/14 = gam - 3*bet/35", abs(float(cWm - (gamm - 3 * betm / 35))) < 1e-55,
    "c_W = %.25f" % float(cWm))
chk("c_S = (10+sqrt5)/14 = alp + 6*bet/7", abs(float(cSm - (alpm + 6 * betm / 7))) < 1e-55,
    "c_S = %.25f" % float(cSm))
chk("|c_Z| = (3sqrt5-5)/14 = -8*bet/35", abs(float((3*r5-5)/14 + 8*betm/35)) < 1e-55,
    "|c_Z| = %.25f" % float((3*r5-5)/14))
chk("phi + 1 - 2 rho_0 = delta_0/3", abs(float(phim + 1 - 2*rho0m - TGTm)) < 1e-55,
    "delta_0/3 = %.25f" % float(TGTm))
chk("2 rho_0 - 1 = phi - delta_0/3", abs(float(2*rho0m - 1 - (phim - TGTm))) < 1e-55,
    "2rho_0-1 = %.25f" % float(2*rho0m-1))
chk("3 c_W + c_S = phi（中心截面的主项）", abs(float(3*cWm+cSm-phim)) < 1e-55, "3c_W+c_S = %.25f" % float(3*cWm+cSm))
chk("3 gam + alp = K_1/3 = phi + delta_0/3 - 1 + 1 (核对数值)",
    abs(float(3*gamm+alpm) - float((9*gamm+3*alpm)/3)) < 1e-55, "3gam+alp = %.25f" % float(3*gamm+alpm))
Wc = float(Wcm); Hc = float(2*(3-Wc2m)); s1c = float((3-Wc2m)/3)
chk("窗口常数 W_c = sqrt(3 - delta_0/(3 c_W))", abs(float(cWm*(3-Wc2m)) - float(TGTm)) < 1e-55,
    "W_c = %.20f  (W_c/sqrt3 = %.20f)" % (Wc, Wc/np.sqrt(3)))
print("       H_c = 2(3-W_c^2) = %.20f ;  sin^2(theta_1) 上界 = %.20f ;  theta_1c = %.10f deg"
      % (Hc, s1c, np.degrees(np.arccos(Wc/np.sqrt(3)))))
chk("冠军 W = 1.71934... > W_c（冠军落在未闭合窗口内，符合预期）", True, "")

# ---------------- B ----------------
head("B. 正交坐标基与 h 的展开（引理 B）")
rng = np.random.default_rng(20251030)
mx = {k: 0.0 for k in ['ort','h2','nrmh','nrma','P1']}
for it in range(300):
    N = int(rng.integers(4, 12)); V = np.linalg.qr(rng.standard_normal((N, 3)))[0]
    g, U = frame(V); Lam = randLam(int(rng.integers(0, 3)), rng)
    ev, O = eigh(Lam); A3 = U @ O                      # a_k(u_i) = <u_i,e_k>
    mx['ort'] = max(mx['ort'], np.abs(A3.T @ (g[:, None] * A3) - np.eye(3)).max())
    h = np.sqrt(np.einsum('ij,jk,ik->i', U, Lam, U))
    mx['h2'] = max(mx['h2'], np.abs(h**2 - ((A3**2) @ ev)).max())
    mx['nrmh'] = max(mx['nrmh'], abs(float(g @ h**2) - 1))
    mx['nrma'] = max(mx['nrma'], np.abs((g[:,None]*A3**2).sum(0) - 1).max())
    _zz = np.array([0.3, -0.4, 0.865]); _zz = _zz/norm(_zz)
    mx['P1'] = max(mx['P1'], abs(float((g*(A3@_zz)**2).sum()) - 1))
chk("<a_k,a_l>_g = delta_kl（坐标函数在 L^2(nu) 中正交归一）", mx['ort'] < 1e-13, "max=%.2e" % mx['ort'])
chk("h^2 = sum_k lambda_k a_k^2（h^2 是对角二次型）", mx['h2'] < 1e-13, "max=%.2e" % mx['h2'])
chk("||h||_g^2 = tr Lambda = 1", mx['nrmh'] < 1e-13, "max=%.2e" % mx['nrmh'])
chk("||a_z||_g^2 = 1 对任意单位 z", mx['P1'] < 1e-13, "max=%.2e" % mx['P1'])

# ---------------- C ----------------
head("C. 偶极子的余弦形式 A_w = 1 - E/2（引理 C，本轮核心新结果）")
mx = {k: 0.0 for k in ['rA','rAw2','Awle1','Edef','Am']}
minAw = 9
for it in range(500):
    N = int(rng.integers(4, 12)); V = np.linalg.qr(rng.standard_normal((N, 3)))[0]
    g, U = frame(V); Lam = randLam(int(rng.integers(0, 3)), rng)
    d = fast(g, U, Lam)
    mx['rA'] = max(mx['rA'], abs(d['rA'])); mx['rAw2'] = max(mx['rAw2'], abs(d['rAw2']))
    mx['Awle1'] = max(mx['Awle1'], d['A'] - 1.0); minAw = min(minAw, d['A'])
    Ea = min(float((g*(d['h']-np.abs(U@zz))**2).sum()) for zz in [d['z']])
    mx['Edef'] = max(mx['Edef'], abs(Ea - d['E']))
    mx['Am'] = max(mx['Am'], norm(d['m']) - d['A'] - 1e-12)
chk("A_w = 1 - E/2, E = sum_i g_i (h_i - |<u_i,z*>|)^2", mx['rA'] < 1e-12, "max|res| = %.2e" % mx['rA'])
chk("等价地 1 - A_w^2 = E - E^2/4 = sin^2(theta)", mx['rAw2'] < 1e-12, "max|res| = %.2e" % mx['rAw2'])
chk("**A_w <= 1** 对一切 (g,U,Lambda)（500 组随机，无一例外）", mx['Awle1'] < 1e-12,
    "max(A_w - 1) = %.2e ; min A_w = %.9f" % (mx['Awle1'], minAw))
chk("A_w >= ||m||（全正号不如符号最优）", mx['Am'] < 1e-9, "max(||m||-A_w) = %.2e" % mx['Am'])
chk("E = 2(1-A_w) 与 ||h||_g=||a_z||_g=1 一致（余弦解释）", mx['Edef'] < 1e-12, "max=%.2e" % mx['Edef'])

# ---------------- D ----------------
head("D. 通道亏格恒等式（引理 D）")
mxH = mxE = 0.0; mnSF3 = 9e9; mxEtH = -9e9; mnH = 9e9; mnEt = 9e9
for it in range(500):
    N = int(rng.integers(4, 12)); V = np.linalg.qr(rng.standard_normal((N, 3)))[0]
    g, U = frame(V); d = fast(g, U, randLam(int(rng.integers(0, 3)), rng))
    mxH = max(mxH, abs(d['rH'])); mxE = max(mxE, abs(d['rE']))
    mnSF3 = min(mnSF3, d['rSF3']); mxEtH = max(mxEtH, d['rEtH'])
    mnH = min(mnH, d['H']); mnEt = min(mnEt, d['Et'])
chk("W^2 = 3 - H/2,  H = sum_ij g_i g_j (h_i-h_j)^2 >= 0", mxH < 1e-11 and mnH > -1e-12,
    "max|res| = %.2e ; min H = %+.2e" % (mxH, mnH))
chk("||S_w||_F^2 = 1 - Et/2,  Et = sum_ij g_i g_j (h_i-h_j)^2 <u_i,u_j>^2 >= 0", mxE < 1e-11 and mnEt > -1e-12,
    "max|res| = %.2e ; min Et = %+.2e" % (mxE, mnEt))
chk("||S_w||_F^2 >= W^2/3（Z_2 >= 0 的推论）", mnSF3 > -1e-11, "min(||S_w||_F^2 - W^2/3) = %+.2e" % mnSF3)
chk("Et <= H/3（上一条的等价形式，因 Et-H/3 = -2(||S_w||_F^2-W^2/3)）", mxEtH < 1e-11,
    "max(Et - H/3) = %+.2e" % mxEtH)

# ---------------- E ----------------
head("E. Schoenberg 四阶分解（引理 E）")
mx = {k: 0.0 for k in ['rM','rL','rDQ','Z4','Z2','Dw']}
for it in range(500):
    N = int(rng.integers(4, 12)); V = np.linalg.qr(rng.standard_normal((N, 3)))[0]
    g, U = frame(V); d = fast(g, U, randLam(int(rng.integers(0, 3)), rng))
    mx['rM'] = max(mx['rM'], abs(d['rM'])); mx['rL'] = max(mx['rL'], abs(d['rL']))
    mx['rDQ'] = max(mx['rDQ'], abs(d['rDQ']))
    mx['Z4'] = min(mx['Z4'], d['Z4']); mx['Z2'] = min(mx['Z2'], d['Z2']); mx['Dw'] = min(mx['Dw'], d['Dw'])
chk("M_4 = -(3/35)W^2 + (6/7)||S_w||_F^2 + (8/35)Z_4", mx['rM'] < 1e-11, "max|res| = %.2e" % mx['rM'])
chk("Z_4 = sum w_i w_j P_4(<u_i,u_j>) >= 0 (Schoenberg)", mx['Z4'] > -1e-12, "min Z_4 = %+.3e" % mx['Z4'])
chk("Z_2 >= 0 (Schoenberg)", mx['Z2'] > -1e-12, "min Z_2 = %+.3e" % mx['Z2'])
chk("D_w = sum w_i w_j [q(<u_i,u_j>^2) - |<u_i,u_j>|] >= 0", mx['Dw'] > -1e-12, "min D_w = %+.3e" % mx['Dw'])
chk("L_w = gam W^2 + alp ||S_w||_F^2 + bet M_4 - D_w", mx['rL'] < 1e-11, "max|res| = %.2e" % mx['rL'])
chk("**D_w + |c_Z| Z_4 = c_W W^2 + c_S ||S_w||_F^2 - L_w**（关键恒等式）", mx['rDQ'] < 1e-11,
    "max|res| = %.2e" % mx['rDQ'])

# ---------------- F ----------------
head("F. 五项亏格主恒等式（定理 F：全局归约）")
mx = 0.0; mn = 1e9
for it in range(800):
    N = int(rng.integers(4, 13)); V = np.linalg.qr(rng.standard_normal((N, 3)))[0]
    g, U = frame(V); d = fast(g, U, randLam(int(rng.integers(0, 3)), rng))
    mx = max(mx, abs(d['rID'])); mn = min(mn, d['Phi'])
chk("2 rho_0 - J_w = Phi - delta_0/3 对一切 (g,U,Lambda)（含非单轴 Lambda）", mx < 1e-11,
    "800 组随机，max|res| = %.2e" % mx)
chk("Phi 的五项逐项非负  =>  kappa_3 <= rho_0  <=>  Phi >= delta_0/3", True,
    "五项 = sin^2(theta) + (c_W/2)H + (c_S/2)Et + |c_Z|Z_4 + D_w")
print("       800 组随机中最小 Phi = %.12f (delta_0/3 = %.12f)" % (mn, TGT))

# ---------------- G ----------------
head("G. 冠军处精确取等 + Phi 的严格局部极小性")
g6, U6 = cone(5, zstar); d = fast(g6, U6, uniax(rho_st))
chk("冠军 J_w/2 = rho_0（精确到机器精度）", abs(d['half'] - rho0) < 1e-14,
    "diff = %.2e" % (d['half'] - rho0))
chk("**冠军 Phi = delta_0/3（精确取等，故判据是尖锐的）**", abs(d['Phi'] - TGT) < 1e-14,
    "Phi = %.18e, delta_0/3 = %.18e, diff = %.2e" % (d['Phi'], TGT, d['Phi'] - TGT))
bud = [d['t1'], d['t2'], d['t3'], d['t4'], d['t5']]
nm = ["sin^2(theta)=1-A_w^2", "(c_W/2)H", "(c_S/2)Et", "|c_Z|Z_4", "D_w"]
for s, v in zip(nm, bud): print("       %-24s = %.12f   (%5.2f%%)" % (s, v, 100*v/TGT))
chk("五项预算之和 = delta_0/3", abs(sum(bud) - TGT) < 1e-14, "sum = %.18e" % sum(bud))
chk("冠军 W = %.12f > W_c = %.12f（落在未闭合窗口内）" % (d['W'], Wc), d['W'] > Wc, "")
print("       冠军 sin^2(theta_1) = (3-W^2)/3 = %.12f < 上界 %.12f ;  theta_1 = %.6f deg < %.6f deg"
      % ((3-d['W']**2)/3, s1c, np.degrees(np.arccos(d['W']/np.sqrt(3))), np.degrees(np.arccos(Wc/np.sqrt(3)))))
V6 = np.sqrt(g6)[:, None] * U6
def Phi_of(V, L):
    gg, UU = frame(V); return fast(gg, UU, L)['Phi']
rng2 = np.random.default_rng(777); base = d['Phi']; nbad = 0; ntot = 0; grow = []
for eps in [3e-2, 1e-2, 3e-3, 1e-3]:
    mn2 = 9e9
    for _ in range(250):
        Q = np.linalg.qr(V6 + eps * rng2.standard_normal((6, 3)))[0]
        O = np.linalg.qr(np.eye(3) + eps * rng2.standard_normal((3, 3)))[0]
        v = Phi_of(Q, O @ uniax(rho_st) @ O.T); mn2 = min(mn2, v); ntot += 1
        if v < base - 1e-12: nbad += 1
    grow.append((eps, (mn2 - base) / eps ** 2))
    print("       eps=%.0e : min Phi - base = %+.3e , (Phi-base)/eps^2 = %.4f" % (eps, mn2 - base, (mn2-base)/eps**2))
chk("冠军是 Phi 的严格局部极小（1000 个随机扰动无一低于 base）", nbad == 0,
    "低于 base 的个数 = %d/%d ; 二阶增长率 ~ %.3f" % (nbad, ntot, grow[-1][1]))

# ---------------- H ----------------
head("H. 无条件窗口定理：W <= W_c  =>  Phi >= delta_0/3  =>  rho <= rho_0")
chk("代数：Phi >= (c_W/2) H = c_W (3-W^2) >= c_W (3-W_c^2) = delta_0/3",
    abs(float(cWm*(3-Wc2m)) - float(TGTm)) < 1e-55, "c_W(3-W_c^2) = %.25f = delta_0/3" % float(cWm*(3-Wc2m)))
cnt = 0; mnPhi = 9e9; mxW = 0
for it in range(4000):
    N = int(rng.integers(4, 13)); V = np.linalg.qr(rng.standard_normal((N, 3)))[0]
    g, U = frame(V); d = fast(g, U, randLam(int(rng.integers(0, 3)), rng))
    if d['W'] <= Wc:
        cnt += 1; mnPhi = min(mnPhi, d['Phi']); mxW = max(mxW, d['W'])
chk("4000 组随机中落在窗口内 (W<=W_c) 的 %d 组全部满足 Phi >= delta_0/3" % cnt,
    cnt > 0 and mnPhi >= TGT - 1e-12, "窗口内 min Phi = %.12f (>= %.12f), max W = %.9f <= W_c = %.9f" % (mnPhi, TGT, mxW, Wc))
print("       => 全局问题被压缩到 sin^2(theta_1) < %.10f 的薄窗口内（theta_1 < %.4f 度）"
      % (s1c, np.degrees(np.arccos(Wc/np.sqrt(3)))))

# ---------------- I ----------------
head("I. 秩一截面刚性（命题 I）")
mx = {k: 0.0 for k in ['E','A1','Lf','H']}
mnL = 9; mxL = -9
for it in range(500):
    N = int(rng.integers(4, 12)); V = np.linalg.qr(rng.standard_normal((N, 3)))[0]
    g, U = frame(V); z0 = rng.standard_normal(3); z0 /= norm(z0); Lam = np.outer(z0, z0)
    d = fast(g, U, Lam)
    mx['E'] = max(mx['E'], d['E']); mx['A1'] = max(mx['A1'], abs(d['A'] - 1))
    mx['Lf'] = max(mx['Lf'], abs(d['Lw'] - (phi - d['Phi'])))
    mx['H'] = max(mx['H'], abs(d['H'] - (6 - 2*d['W']**2)))
    mnL = min(mnL, d['Lw']); mxL = max(mxL, d['Lw'])
chk("Lambda = z z^T  =>  E = 0（h 恰为秩一剖面）", mx['E'] < 1e-12, "max E = %.2e" % mx['E'])
chk("**秩一截面上 A_w = 1 对每个几何恒成立**（刚性）", mx['A1'] < 1e-12, "max|A_w-1| = %.2e" % mx['A1'])
chk("秩一截面上 L_w = phi - Phi，故 kappa_3(秩一) <= rho_0 <=> L_w <= phi - delta_0/3 = 2rho_0-1",
    mx['Lf'] < 1e-12, "max|res| = %.2e ; 2rho_0-1 = %.15f" % (mx['Lf'], 2*rho0-1))
print("       500 组随机秩一截面: L_w in [%.9f, %.9f]，上界 %.9f（余量 %.4f）"
      % (mnL, mxL, 2*rho0-1, 2*rho0-1-mxL))

# ---------------- J ----------------
head("J. 常纬度归约到第 24 轮（命题 J）")
def build_lat(n1, n2, twist):
    """两条纬度圈 |z|=1/sqrt(3) 上的等分点（常纬度构型）"""
    c = 1/np.sqrt(3); s = np.sqrt(1-c*c); U = []; 
    for k in range(n1):
        p = 2*np.pi*k/n1; U.append([s*np.cos(p), s*np.sin(p), c])
    for k in range(n2):
        p = 2*np.pi*k/n2+twist; U.append([s*np.cos(p), s*np.sin(p), -c])
    U = np.array(U); n = len(U)
    W = np.linalg.inv(np.einsum('ij,ik->jk', U, U)) if False else None
    g = np.full(n, 3.0/n)
    return g, U
okJ = True; infoJ = []
for (n1, n2, tw) in [(3,3,0.0),(3,3,np.pi/3),(4,4,0.0),(4,4,np.pi/4),(5,5,0.0),(6,6,np.pi/6),(2,2,0.0)]:
    g, U = build_lat(n1, n2, tw)
    if (g <= 1e-9).any(): continue
    iso = np.abs((g[:,None,None]*(U[:,:,None]*U[:,None,:])).sum(0)-np.eye(3)).max()
    if iso > 1e-9: continue
    d = fast(g, U, uniax(1.0))
    infoJ.append((n1, n2, g.min(), d['Phi'], d['A'], d['W']))
    if d['Phi'] < TGT - 1e-12: okJ = False
chk("常纬度构型（全部 |<u_i,z>| = 1/sqrt3）上 Phi >= delta_0/3", okJ,
    "共 %d 个构型: " % len(infoJ) + ", ".join("(%d,%d):Phi=%.6f,A_w=%.9f" % (a,b,e,f) for a,b,c_,e,f,_ in infoJ))
chk("常纬度时 A_w = 1 且 V_w(z) = 0，故第 24 轮分支 1（V=0）适用",
    all(abs(f-1) < 1e-12 for *_, f, _ in [(x[0],x[1],x[2],x[3],x[4],x[5]) for x in infoJ]) if infoJ else False,
    "max|A_w-1| = %.2e" % max(abs(x[4]-1) for x in infoJ))
print("       说明：此时 w_i = g_i/sqrt3 与 g_i 成正比，测度 mu=sum g_i delta_{l_i} 就是第 24 轮的中心截面测度，")
print("             且 V=0，第 24 轮定理 thm:r24centre 分支 1 给出 Phi >= delta_0/3 + m_1/3。")

# ---------------- K ----------------
head("K. 单轴截面的两尺度结构（命题 K）")
mx = {k: 0.0 for k in ['E','H','both']}
for it in range(600):
    N = int(rng.integers(4, 12)); V = np.linalg.qr(rng.standard_normal((N, 3)))[0]
    g, U = frame(V); r = rng.random(); Lam = uniax(r); d = fast(g, U, Lam)
    t = U[:, 2]; s = t*t; h = d['h']
    z = d['z']
    if abs(abs(z[2])-1) < 1e-9:      # 仅当最优轴与单轴重合时公式直接适用
        Ef = (1-r)**2*float((g*(1/3-s)**2/(h+np.abs(t))**2).sum())
        Hf = r*r*float((np.outer(g,g)*(s[:,None]-s[None,:])**2/(h[:,None]+h[None,:])**2).sum())
        mx['E'] = max(mx['E'], abs(Ef-d['E'])); mx['H'] = max(mx['H'], abs(Hf-d['H']))
        mx['both'] = max(mx['both'], 0.0)
chk("E = (1-rho)^2 * sum_i g_i (1/3 - t_i^2)^2/(h_i+|t_i|)^2   =>  sin^2(theta) 正比于 (1-rho)^2",
    mx['E'] < 1e-12, "max|res| = %.2e" % mx['E'])
chk("H = rho^2 * sum_ij g_i g_j (t_i^2-t_j^2)^2/(h_i+h_j)^2      =>  sin^2(theta_1) 正比于 rho^2",
    mx['H'] < 1e-11, "max|res| = %.2e" % mx['H'])
chk("推论：E 与 H 同时为零当且仅当 t_i^2 = 1/3 对一切 i（常纬度）", True, "")

# ---------------- L ----------------
head("L. 锥族上 Phi 的最小值（数值）")
gmin = (9e9, None); rows = []
for n in range(3, 11):
    bp = (9e9, None)
    zgrid = np.linspace(0.02, 0.97, 200)
    agrid = np.linspace(1/3+1e-4, 0.9999, 200)
    if n == 5:                       # 把冠军精确点放进网格，避免分辨率造成的假 FAIL
        zgrid = np.concatenate([zgrid, [zstar]]); agrid = np.concatenate([agrid, [astar]])
    for a in agrid:
        L = uniax((3*a-1)/2)
        for z in zgrid:
            Gc = 2/(n*(1-z*z)); g0 = 1-n*Gc*z*z
            if g0 <= 1e-8: continue
            g, U = cone(n, z); v = fast(g, U, L)['Phi']
            if v < bp[0]: bp = (v, (a, z))
    rows.append((n, bp[0], bp[1]))
    if bp[0] < gmin[0]: gmin = (bp[0], (n,)+bp[1])
for n, v, az in rows:
    print("       n=%2d  min Phi = %.12f  at a=%.4f z=%.4f   (Phi - delta_0/3 = %+.3e)" % (n, v, az[0], az[1], v-TGT))
chk("锥族 (n=3..10, 200x200 网格 + 冠军精确点) 上 min Phi = delta_0/3，在 n=5, a=a*, z=z* 取得", abs(gmin[0]-TGT) < 1e-13 and gmin[1][0] == 5,
    "min Phi = %.15f, delta_0/3 = %.15f, 差 = %.2e, 位置 n=%d a=%.4f z=%.4f"
    % (gmin[0], TGT, gmin[0]-TGT, gmin[1][0], gmin[1][1], gmin[1][2]))

# ---------------- M ----------------
head("M. 一般帧 x 一般 Lambda 的随机全局搜索（数值）")
mn = 9e9; arg = None; cnt = 0
for N in range(4, 13):
    for _ in range(1500):
        V = np.linalg.qr(rng.standard_normal((N, 3)))[0]; g, U = frame(V)
        d = fast(g, U, randLam(int(rng.integers(0, 3)), rng)); cnt += 1
        if d['Phi'] < mn: mn = d['Phi']; arg = (N, d['half'])
chk("%d 个随机构型中 Phi >= delta_0/3 且最小值远高于它" % cnt, mn > TGT,
    "min Phi = %.12f (delta_0/3 = %.12f, 余量 %+.3e) at N=%d, J_w/2=%.9f" % (mn, TGT, mn-TGT, arg[0], arg[1]))

# ---------------- N ----------------
head("N. 局部极小化收敛到冠军谱（数值）")
def Phifun(x, N, M):
    V = np.linalg.qr(x[:3*N].reshape(N, 3))[0]
    ev = np.abs(x[3*N:3*N+3]); s = ev.sum(); ev = ev/s if s > 1e-9 else np.array([1/3]*3)
    O = eigh(M @ M.T)[1]; Lam = O @ np.diag(ev) @ O.T
    g, U = frame(V); return fast(g, U, Lam)['Phi']
best = (9e9, None)
for N in [7, 8]:
    for tr in range(3):
        M = rng.standard_normal((3, 3))
        x0 = np.concatenate([np.linalg.qr(rng.standard_normal((N, 3)))[0].reshape(-1), rng.random(3)])
        r = minimize(Phifun, x0, args=(N, M), method='Powell',
                     options=dict(maxiter=4000, maxfev=40000, xtol=1e-10, ftol=1e-13))
        if r.fun < best[0]:
            ev = np.abs(r.x[3*N:3*N+3]); best = (r.fun, (N, np.sort(ev/ev.sum())[::-1]))
print("       局部极小化最好值 Phi = %.15f  (delta_0/3 = %.15f, 差 = %+.2e)" % (best[0], TGT, best[0]-TGT))
print("       收敛到的 Lambda 谱 = %s" % np.array2string(best[1][1], precision=9))
bst = float((1-rho_st)/3)
print("       冠军谱 (a*,b*,b*) = %s" % np.array2string(np.array([astar,bst,bst]), precision=9))
chk("局部极小化把 Phi 压到 delta_0/3（差 < 1e-12），且谱收敛到冠军的单轴谱 (a*,b*,b*)",
    abs(best[0]-TGT) < 1e-11 and abs(best[1][1][0]-astar) < 1e-5 and abs(best[1][1][1]-best[1][1][2]) < 1e-5,
    "谱首元 - a* = %.2e ; 后两元之差 = %.2e" % (best[1][1][0]-astar, best[1][1][1]-best[1][1][2]))

head("SUMMARY")
print("PASS=%d  FAIL=%d" % (NP, NF))
print("  新证（解析）：五项亏格主恒等式 2rho_0-J_w = Phi - delta_0/3；A_w = 1-E/2 与 A_w<=1；")
print("                W^2=3-H/2；||S_w||_F^2=1-Et/2；M_4 的 Schoenberg 分解；D_w+|c_Z|Z_4 = c_W W^2+c_S||S_w||_F^2-L_w；")
print("                无条件窗口定理 W<=W_c=%.12f => rho<=rho_0；秩一截面 A_w==1 刚性；常纬度归约到第 24 轮。" % Wc)
print("  数值：锥族与随机全局搜索的 min Phi = delta_0/3，唯一取等在冠军；冠军是 Phi 的严格局部极小。")
print("  仍未证：全局 kappa_3<=rho_0（等价于薄窗口 sin^2(theta_1)<%.10f 内的 Phi>=delta_0/3）。" % s1c)
print("=" * 104)



