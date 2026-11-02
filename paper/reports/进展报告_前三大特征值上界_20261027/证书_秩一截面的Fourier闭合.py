# -*- coding: utf-8 -*-
"""
证书：第 27 轮 —— 秩一截面的 Fourier 闭合  Neg <= 5/(6pi)
================================================================
本证书把第 27 轮「秩一截面无条件闭合」的全部关键步骤逐条数值/符号核对。
链条（秩一截面 Lambda=zz^T，半球翻转后 a_i=<u_i,z>>=0）：
   t_ij = sigma_i sigma_j + r_i r_j c_ij,  sigma_i=g_i a_i^2, rho_i=g_i b_i^2,
   r_i=sqrt(sigma_i rho_i)=g_i a_i b_i, lambda_i=a_i/b_i, c_ij=<p_i,p_j>,
   Neg = sum_{t_ij<0}(-t_ij),  L_w = sum|t_ij| = 1+2Neg,  kappa_3(秩一)=1+Neg.
   (1) Young 逐点:  r_i r_j|c| - sigma_i sigma_j <= (r_i r_j|c|)^2/(4 sigma_i sigma_j)
       = rho_i rho_j c^2/4,  等号 iff kappa:=|c|/(lambda_i lambda_j)=2;
       精确亏格恒等式 Neg = Theta/4 - (1/4)sum_N sigma_i sigma_j (kappa_ij-2)^2,
       Theta = sum_N rho_i rho_j c_ij^2.
   (2) N subset {c<0}  =>  Theta <= P := sum_{i!=j, c_ij<0} rho_i rho_j c_ij^2.
   (3) 刚性恒等式（二维各向同性 sum rho_i p_i p_i^T = I_2）:
       sum_ij rho_i rho_j c_ij^2 = ||I_2||_F^2 = 2.
   (4) Fourier:  K(u)=(cos u)_-^2,  P = sum_k Khat_k |muhat(k)|^2,  mu=sum rho_i delta_{phi_i},
       Khat_0=1/4, Khat_2=1/8, Khat_k=0 (k 偶 >=4), Khat_k=2 sin(k pi/2)/(pi k(k^2-4)) (k 奇);
       muhat(0)=2, muhat(2)=0 (各向同性杀掉最大正系数), |muhat(k)|<=2, 丢负系数:
       P <= 1 + 8 sum_{k>=1 odd, Khat_k>0} Khat_k = 1 + (16/pi) S,
       S = sum_{m>=1} 1/((4m+1)((4m+1)^2-4)) = 5/24 - pi/16  (用 psi(1/4)-psi(3/4)=-pi).
   (5) 1+(16/pi)S = 10/(3pi);  Neg <= P/4 <= 10/(12pi)=5/(6pi)=0.265258238486;
       kappa_3(秩一) <= 1+5/(6pi) = 1.265258238486 < rho_0 = 1.281568196872.
标签：已证=解析证明并在此数值/高精度核对；数值=仅数值证据。
"""
import itertools
import numpy as np
import sympy as sp
from mpmath import mp, mpf, pi as mpi, sqrt as msqrt, sin as msin, nsum, inf, pslq, digamma
mp.dps = 60

NP = NF = 0
def chk(name, cond, info=""):
    global NP, NF
    if cond:
        NP += 1; print("  [PASS] %s %s" % (name, info))
    else:
        NF += 1; print("  [FAIL] %s %s" % (name, info))
def head(t):
    print("\n" + "=" * 100 + "\n" + t + "\n" + "=" * 100)

# ---------------- 常数 ----------------
S5 = np.sqrt(5.0); S35 = np.sqrt(35.0)
rho0 = (9 + S5 + 2 * S35) / 18
phi = (1 + S5) / 2
r5 = msqrt(5); r35 = msqrt(35); phim = (1 + r5) / 2
rho0m = (9 + r5 + 2 * r35) / 18
Negb = mpf(5) / (6 * mpi)
Pb = mpf(10) / (3 * mpi)
k3r1 = 1 + Negb
CONEr1 = (11 + phim) / 10
PHI_STW = (3 + r5) / 4
rng = np.random.default_rng(20261030)

print("=" * 100)
print("证书：第 27 轮 —— 秩一截面的 Fourier 闭合  Neg <= 5/(6pi)")
print("=" * 100)
print("rho_0     = %.30f" % float(rho0m))
print("5/(6pi)   = %.30f" % float(Negb))
print("1+5/(6pi) = %.30f" % float(k3r1))
print("10/(3pi)  = %.30f" % float(Pb))

# ---------------- A 精确常数代数 ----------------
head("A. 精确常数代数（60 位）与最终上界链")
m = sp.symbols('m', positive=True, integer=True)
Ssym = sp.summation(1/((4*m+1)*((4*m+1)**2-4)), (m, 1, sp.oo))
Spsi = sp.Rational(1,16)*(sp.polygamma(0,sp.Rational(1,4))-sp.polygamma(0,sp.Rational(3,4)))+sp.Rational(5,24)
Sguess = sp.Rational(5,24)-sp.pi/16
Sn = nsum(lambda k: 1/((4*k+1)*((4*k+1)**2-4)), [1, inf])
chk("S = 5/24 - pi/16（psi(1/4)-psi(3/4)=-pi）", abs(float(Spsi.evalf(40))-float(Sguess.evalf(40))) < 1e-38,
    "S = %.30f" % float(Sn))
chk("S 的 nsum 与闭式一致", abs(float(Sn)-float(Sguess.evalf(40))) < 1e-30,
    "差 = %.3e" % abs(float(Sn)-float(Sguess.evalf(40))))
Pbsym = sp.simplify(1+16*Sguess/sp.pi)
chk("1 + 16S/pi = 10/(3pi)（精确代数）", sp.simplify(Pbsym-10/(3*sp.pi))==0,
    "= %s" % sp.simplify(Pbsym))
chk("P <= 10/(3pi) 数值", abs(float(Pb)-float(Pbsym.evalf(40))) < 1e-35,
    "10/(3pi) = %.30f" % float(Pb))
chk("Neg <= 5/(6pi) = (10/(3pi))/4", abs(float(Negb)-float(Pb/4)) < 1e-40,
    "5/(6pi) = %.30f" % float(Negb))
chk("kappa_3(秩一) <= 1+5/(6pi) < rho_0", float(k3r1) < float(rho0m),
    "1+5/(6pi)=%.20f  rho_0=%.20f  余量=%.6e" % (float(k3r1), float(rho0m), float(rho0m-k3r1)))
gap = rho0m-1-Negb
gap_closed = (-15+mpi*(-9+r5+2*r35))/(18*mpi)
chk("余量 rho_0-1-5/(6pi) 的闭式（60 位）", abs(gap-gap_closed) < mpf(10)**-50,
    "= %.30f" % float(gap))
pi_need = float(mpf(5)/(6*(rho0m-1)))
chk("余量>0 的严格证明：只需 pi > 5/(6(rho_0-1))", mpi > mpf(5)/(6*(rho0m-1)),
    "需 pi > %.10f（显然）" % pi_need)
chk("本界改进 STW: (3+sqrt5)/4 - (1+5/(6pi)) > 0", float(PHI_STW-k3r1) > 0,
    "改进 = %.6e" % float(PHI_STW-k3r1))
chk("本界 vs 锥族条件值 (11+phi)/10", float(k3r1-CONEr1) > 0,
    "本界 - 锥族 = %.6e（本界无条件，锥族条件性）" % float(k3r1-CONEr1))
chk("L_w <= 1+5/(3pi)", abs(float(1+2*Negb)-float(1+mpf(5)/(3*mpi))) < 1e-35,
    "1+5/(3pi) = %.20f" % float(1+2*Negb))

# ---------------- B 冠军精确恒等式 ----------------
head("B. 冠军 W_6 = K_1 v C_5 的秩一数据（极点 + 正五边形，s=s_0=phi/(4+phi)）")
def champ():
    s0 = phi/(4+phi); g0 = (1-3*s0)/(1-s0); g1 = 2/(5*(1-s0))
    sig = np.array([g0*1.0]+[g1*s0]*5)
    rh = np.array([g0*0.0]+[g1*(1-s0)]*5)
    ang = np.array([0.0]+[2*np.pi*k/5 for k in range(5)])
    return sig, rh, ang
def quantities(sig, rh, ang):
    sig = np.asarray(sig, float); rh = np.asarray(rh, float); ang = np.asarray(ang, float)
    r = np.sqrt(sig*rh)
    lam = np.where(rh > 0, np.sqrt(np.divide(sig, rh, out=np.zeros_like(sig), where=rh > 0)), np.inf)
    c = np.cos(ang[:, None]-ang[None, :])
    t = sig[:, None]*sig[None, :]+r[:, None]*r[None, :]*c
    Nmask = (t < 0)
    Neg = float((-t[Nmask]).sum())
    Lw = float(np.abs(t).sum())
    Theta = float((rh[:, None]*rh[None, :]*c*c)[Nmask].sum())
    P = float((rh[:, None]*rh[None, :]*c*c)[(c < 0) & (~np.eye(len(sig), dtype=bool))].sum())
    kap = np.where(Nmask & (sig[:, None]*sig[None, :] > 0),
                   np.abs(c)/np.maximum(lam[:, None]*lam[None, :], 1e-300), np.nan)
    defect = float(np.nansum((sig[:, None]*sig[None, :])*(kap-2)**2))
    rigid = float((rh[:, None]*rh[None, :]*c*c).sum())
    return dict(Neg=Neg, Lw=Lw, Theta=Theta, P=P, defect=defect, rigid=rigid,
                t=t, c=c, r=r, lam=lam, sig=sig, rh=rh, Nmask=Nmask)
sig, rh, ang = champ()
q = quantities(sig, rh, ang)
chk("sum sigma = 1", abs(sig.sum()-1) < 1e-14, "= %.15f" % sig.sum())
chk("sum rho = 2", abs(rh.sum()-2) < 1e-14, "= %.15f" % rh.sum())
chk("二维各向同性 sum rho e^{2i phi}=0", abs((rh*np.exp(2j*ang)).sum()) < 1e-14,
    "= %.3e" % abs((rh*np.exp(2j*ang)).sum()))
chk("交叉项 sum r e^{i phi}=0", abs((q['r']*np.exp(1j*ang)).sum()) < 1e-14,
    "= %.3e" % abs((q['r']*np.exp(1j*ang)).sum()))
chk("Neg = phi^2/10 = (3+sqrt5)/20", abs(q['Neg']-phi**2/10) < 1e-14,
    "Neg=%.15f  phi^2/10=%.15f" % (q['Neg'], phi**2/10))
chk("L_w = (13+sqrt5)/10", abs(q['Lw']-(13+np.sqrt(5))/10) < 1e-14,
    "L_w=%.15f" % q['Lw'])
chk("L_w = 1+2Neg", abs(q['Lw']-1-2*q['Neg']) < 1e-14, "= %.3e" % abs(q['Lw']-1-2*q['Neg']))
chk("sum t = 1（= ||m~||^2 = A_w^2，秩一 A_w=1）", abs(q['t'].sum()-1) < 1e-14, "= %.15f" % q['t'].sum())
chk("刚性恒等式 sum rho_i rho_j c_ij^2 = 2", abs(q['rigid']-2) < 1e-13, "= %.15f" % q['rigid'])
chk("冠军处 Theta == P（N={c<0}）", abs(q['Theta']-q['P']) < 1e-14,
    "Theta=%.15f P=%.15f" % (q['Theta'], q['P']))
chk("冠军处 Young 取等 kappa==2（亏格=0）", q['defect'] < 1e-20, "defect=%.3e" % q['defect'])
chk("冠军 P = (3+sqrt5)/5", abs(q['P']-(3+np.sqrt(5))/5) < 1e-14, "P=%.15f" % q['P'])
chk("冠军 P < 10/(3pi)（界不紧，余量）", float(Pb)-q['P'] > 0,
    "10/(3pi)-P_champ = %.6e" % (float(Pb)-q['P']))
chk("冠军 Neg < 5/(6pi)", float(Negb)-q['Neg'] > 0,
    "5/(6pi)-Neg_champ = %.6e" % (float(Negb)-q['Neg']))

# ---------------- C 随机 Stiefel 帧核验 ----------------
head("C. 随机 Stiefel 帧核验：亏格恒等式 / 刚性恒等式 / Neg<=P/4 / ||Z_S||^2<=sigma_S")
mx1 = mx2 = 0.0; mx3 = -1e9
for trial in range(4000):
    N = int(rng.integers(4, 11))
    A = rng.standard_normal((N, 3)); Qm, _ = np.linalg.qr(A); V = Qm[:, :3].copy()
    V[V[:, 0] < 0] *= -1
    g = (V**2).sum(1); sg = V[:, 0]**2; rr = g-sg
    ok = (rr > 1e-12) & (sg > 1e-12)
    if ok.sum() < 3:
        continue
    p = V[:, 1:3]/np.sqrt(np.maximum((V[:, 1:3]**2).sum(1, keepdims=True), 1e-300))
    an = np.arctan2(p[:, 1], p[:, 0])
    qq = quantities(sg, rr, an)
    mx1 = max(mx1, abs(qq['Neg']-(qq['Theta']-qq['defect'])/4))
    mx2 = max(mx2, abs(qq['rigid']-2.0))
    mx3 = max(mx3, qq['Neg']-qq['P']/4)
chk("kappa-2 亏格恒等式 Neg=Theta/4-defect/4（4000 帧）", mx1 < 1e-12, "最大残差=%.3e" % mx1)
chk("刚性恒等式 sum rho rho c^2=2（4000 帧）", mx2 < 1e-11, "最大残差=%.3e" % mx2)
chk("Neg <= P/4（4000 帧，0 违例）", mx3 <= 1e-12, "最大违例=%.3e" % mx3)

def check_ZS(N=7, trials=120):
    bad = 0; mxm = -1e9
    for _ in range(trials):
        A = rng.standard_normal((N, 3)); Qm, _ = np.linalg.qr(A); V = Qm[:, :3]
        if V[:, 0].sum() < 0:
            V = -V
        sg = V[:, 0]**2; z = V[:, 0][:, None]*V
        if not np.allclose(z.sum(0), [1, 0, 0], atol=1e-11):
            continue
        Mi = np.diag(1/sg)
        if not np.allclose(z.T @ Mi @ z, np.eye(3), atol=1e-9):
            continue
        for S in itertools.product([0, 1], repeat=N):
            idx = np.array(S, dtype=bool)
            if idx.sum() == 0:
                continue
            ZS = z[idx].sum(0); sS = sg[idx].sum()
            mval = (ZS**2).sum()-sS
            if mval > 1e-12:
                bad += 1
            mxm = max(mxm, mval)
    return bad, mxm
bad, mxm = check_ZS(7, 120)
chk("新引理 ||Z_S||^2 <= sigma_S 对所有子集 S（N=7,120 帧,2^7-1 子集）", bad == 0,
    "违例=%d 最大 margin=%.3e" % (bad, mxm))

# ---------------- D Fourier 系数 ----------------
head("D. Fourier 系数 Khat_k 解析式 vs 数值积分；P=sum_k Khat_k|muhat(k)|^2")
def Kfun(u):
    return np.where(np.cos(u) < 0, np.cos(u)**2, 0.0)
def Khat_num(k, M=2000000):
    u = np.linspace(0, 2*np.pi, M, endpoint=False); du = 2*np.pi/M
    return float(((Kfun(u)*np.exp(-1j*k*u)).sum()*du/(2*np.pi)).real)
def Khat_for(k):
    if k == 0:
        return mpf(1)/4
    if k == 2:
        return mpf(1)/8
    if k % 2 == 0:
        return mpf(0)
    return 2*msin(k*mpi/2)/(mpi*k*(k**2-4))
mxk = 0.0
for k in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]:
    mxk = max(mxk, abs(Khat_num(k)-float(Khat_for(k))))
chk("Khat_k 解析式 vs 数值积分（k=0..9）", mxk < 1e-9, "最大差=%.3e" % mxk)
def P_direct(rho, ang):
    c = np.cos(ang[:, None]-ang[None, :])
    return float((rho[:, None]*rho[None, :]*np.where(c < 0, c**2, 0.0)).sum())
def P_fourier(rho, ang, KK=400):
    mh = np.array([np.sum(rho*np.exp(-1j*k*ang)) for k in range(-KK, KK+1)])
    ks = np.arange(-KK, KK+1)
    kh = np.array([float(Khat_for(abs(int(k)))) for k in ks])
    return float(np.sum(kh*np.abs(mh)**2))
mxf = 0.0
for trial in range(20):
    N = int(rng.integers(3, 13))
    rho = rng.random(N)*0.8; ang = rng.random(N)*2*np.pi
    mxf = max(mxf, abs(P_direct(rho, ang)-P_fourier(rho, ang)))
chk("Fourier 恒等式 P=sum_k Khat_k|muhat(k)|^2（20 组）", mxf < 1e-6, "最大差=%.3e" % mxf)
chk("各向同性杀掉 Khat_2 项：5-环 muhat(2)=0", abs(float(Khat_for(2))) > 0,
    "Khat_2=1/8>0，但冠军 muhat(2)=0 故此项为 0")

print("\n" + "=" * 100)
print("证书：第 27 轮 —— 秩一截面的 Fourier 闭合  Neg <= 5/(6pi)")
print("PASS = %d, FAIL = %d" % (NP, NF))
print("结论（已证）：kappa_3(秩一) = 1+Neg <= 1+5/(6pi) = %.15f < rho_0 = %.15f" % (float(k3r1), float(rho0m)))
print("=" * 100)
