# -*- coding: utf-8 -*-
"""_aud13.py : 修正 _aud9.py 中 budget4 残差的口径错误。
_aud9.py 把变量名 K2 赋成了 E[RR'(a+b(t^2-1/3)-rho(t^4-1/5)-|t|)]，
而这恰好等于 S=E[RR'h_0(|t|)]（上界函数恒等式），不是 K^2=||T_4||_F^2。
因此 "budget4" 那一行把 S 数了两次，残差 = rho*(K^2 - S)，与稿件无关。
本脚本用正确的 K^2 重算。"""
import numpy as np
from numpy.linalg import eigh, norm
sq5 = np.sqrt(5.0)
a_ = (1 + sq5) / 6.0
b_ = (7 * sq5 - 5) / 8.0
rho = (15 * sq5 - 25) / 16.0
alpha = (2 * sq5 - 1) / 14.0
beta = (sq5 + 10) / 14.0
I3 = np.eye(3)

def h0(z):
    return (1 - z) * (sq5 * z - 1) ** 2 * (sq5 * z + sq5 + 2) / (2 * (sq5 + 1) ** 2)

d_ij = np.eye(3)
OMEGA = (np.einsum('ij,kl->ijkl', d_ij, d_ij)
         + np.einsum('ik,jl->ijkl', d_ij, d_ij)
         + np.einsum('il,jk->ijkl', d_ij, d_ij)) / 15.0

def L4(Z):
    return (np.einsum('ij,kl->ijkl', d_ij, Z) + np.einsum('ik,jl->ijkl', d_ij, Z)
            + np.einsum('il,jk->ijkl', d_ij, Z) + np.einsum('kl,ij->ijkl', d_ij, Z)
            + np.einsum('jl,ik->ijkl', d_ij, Z) + np.einsum('jk,il->ijkl', d_ij, Z)) / 7.0

def Psi4(u):
    u4 = np.einsum('i,j,k,l->ijkl', u, u, u, u)
    uu = np.outer(u, u)
    return u4 - OMEGA - L4(uu - I3 / 3.0)

def analyse(p, Vm):
    Ri = norm(Vm, axis=1); Ui = Vm / np.where(Ri > 0, Ri, 1)[:, None]
    w = Ui @ Ui.T
    RR = Ri[:, None] * Ri[None, :]
    pp = p[:, None] * p[None, :]
    A = (pp * np.abs(Vm @ Vm.T)).sum()
    m = (p[:, None] * Vm).sum(0); B = m @ m; e = 1 - B
    C = (p * Ri).sum(); V = 3 - C ** 2
    D = (3 + sq5) / 2 - A - B
    M = np.einsum('i,ij,ik->jk', p, Vm, Ui)
    Z = M - (C / 3) * I3
    Znorm2 = float(np.sum(Z * Z))
    J = 1 - np.sum(M * M)
    S_direct = (pp * RR * h0(np.abs(w))).sum()
    S_majorant = (pp * RR * (a_ + b_ * (w ** 2 - 1 / 3) - rho * (w ** 4 - 1 / 5) - np.abs(w))).sum()
    T4 = np.einsum('i,imnop->mnop', p * Ri, np.stack([Psi4(u) for u in Ui]))
    K2_true = float(np.sum(T4 ** 2))
    # E[RR'(U.U')^4] = ||E[R U^{ox4}]||^2 = K^2 + (6/7)||Z||^2 + C^2/5
    RU4 = np.einsum('i,imnop->mnop', p * Ri, np.stack([np.einsum('i,j,k,l->ijkl', u, u, u, u) for u in Ui]))
    m4 = float(np.sum(RU4 ** 2))
    return dict(A=A, B=B, e=e, C=C, V=V, D=D, J=J, Znorm2=Znorm2,
                S_direct=S_direct, S_majorant=S_majorant, K2_true=K2_true,
                m4=m4, m4_pred=K2_true + 6 / 7 * Znorm2 + C ** 2 / 5,
                Jpred=V / 3 - Znorm2)

def champion():
    b = (7 + 5 * np.sqrt(7)) / 126; a = 1 - 5 * b
    wt = np.array([a] + [b] * 5)
    A6 = np.zeros((6, 6))
    for i in range(5): A6[0, i + 1] = A6[i + 1, 0] = 1
    for i in range(5):
        j = (i + 1) % 5; A6[i + 1, j + 1] = A6[j + 1, i + 1] = 1
    Mq = np.diag(np.sqrt(wt)) @ (A6 + np.eye(6)) @ np.diag(np.sqrt(wt))
    ev, Vv = eigh(Mq); o = np.argsort(-ev); Vv = Vv[:, o]
    F = Vv[:, :3] / np.sqrt(wt)[:, None]
    if F[:, 0].sum() < 0: F[:, 0] *= -1
    return wt, F


def _selftest():
    g = np.random.default_rng(1)
    Zs = g.normal(size=(3, 3)); Zs = Zs + Zs.T; Zs -= np.trace(Zs) / 3 * I3
    L = L4(Zs)
    ct = np.einsum('iikl->kl', L)
    print("  selftest ct L4(Z)-Z            = %.3e" % np.abs(ct - Zs).max())
    print("  selftest ||L4(Z)||^2-(6/7)|Z|^2 = %.3e" % abs(np.sum(L * L) - 6 / 7 * np.sum(Zs * Zs)))
    u = g.normal(size=3); u /= norm(u)
    print("  selftest |Psi3(u)|^2-2/5        = %.3e" % abs(np.sum(Psi3(u) ** 2) - 0.4))
    print("  selftest |Psi4(u)|^2-8/35       = %.3e" % abs(np.sum(Psi4(u) ** 2) - 8 / 35))
    v = g.normal(size=3); v /= norm(v)
    print("  selftest <Psi3(u),Psi3(v)>-g3   = %.3e" % abs(np.sum(Psi3(u) * Psi3(v)) - ((u@v)**3 - 3*(u@v)/5)))
    print("  selftest <Psi4(u),Psi4(v)>-g4   = %.3e" % abs(np.sum(Psi4(u) * Psi4(v)) - ((u@v)**4 - 6*(u@v)**2/7 + 3/35)))

def Psi3(u):
    return np.einsum('i,j,k->ijk', u, u, u) - (np.einsum('i,jk->ijk', u, d_ij)
            + np.einsum('j,ik->ijk', u, d_ij) + np.einsum('k,ij->ijk', u, d_ij)) / 5.0

print("=== tensor self-tests ===")
_selftest()

print("=== champion ===")
oc = analyse(*champion())
for k in ['A', 'B', 'e', 'C', 'V', 'D', 'J', 'Jpred', 'Znorm2', 'S_direct', 'S_majorant', 'K2_true', 'm4', 'm4_pred']:
    print(f"  {k:12s} = {oc[k]: .15f}")
print(f"  S_direct - S_majorant = {oc['S_direct']-oc['S_majorant']:+.3e}   (上界函数恒等式)")
print(f"  m4 - m4_pred          = {oc['m4']-oc['m4_pred']:+.3e}   (谐波分解)")
print(f"  J - Jpred             = {oc['J']-oc['Jpred']:+.3e}")
r_bad = oc['D'] - (oc['e'] + alpha * oc['V'] + beta * oc['J'] + rho * oc['S_direct'] + oc['S_direct'])
r_ok  = oc['D'] - (oc['e'] + alpha * oc['V'] + beta * oc['J'] + rho * oc['K2_true'] + oc['S_direct'])
print(f"  budget4 残差（_aud9.py 的错误口径, 用 S 冒充 K^2）= {r_bad:+.3e}")
print(f"  budget4 残差（正确口径, 用真 K^2）               = {r_ok:+.3e}")
print(f"  rho*(K^2-S) = {rho*(oc['K2_true']-oc['S_direct']):+.3e}   <- 应与错误口径一致")

print("\n=== 12 random isotropic laws ===")
g = np.random.default_rng(2026)
worst_ok = 0.0; worst_h0 = 0.0; worst_m4 = 0.0
for t in range(12):
    k = int(g.integers(4, 10))
    Y = g.normal(size=(k, 3)); w0 = g.random(k) + 0.2; w0 /= w0.sum()
    S2 = np.einsum("i,ij,ik->jk", w0, Y, Y); ev, U = eigh(S2)
    V_ = Y @ (U @ np.diag(1 / np.sqrt(ev)) @ U.T)
    o = analyse(w0, V_)
    r = o['D'] - (o['e'] + alpha * o['V'] + beta * o['J'] + rho * o['K2_true'] + o['S_direct'])
    worst_ok = max(worst_ok, abs(r))
    worst_h0 = max(worst_h0, abs(o['S_direct'] - o['S_majorant']))
    worst_m4 = max(worst_m4, abs(o['m4'] - o['m4_pred']))
print(f"  max |budget4 residual| (correct K^2)      = {worst_ok:.3e}")
print(f"  max |S_direct - S_majorant|               = {worst_h0:.3e}")
print(f"  max |E[RR'w^4] - (K^2+6||Z||^2/7+C^2/5)| = {worst_m4:.3e}")
