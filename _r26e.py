# -*- coding: utf-8 -*-
# _r26e.py  第 26 轮补充复核（向量化）：
#   (1) n=3..24 的 max_s L_w(s,n) 独立重算，验证片段所列 n>=7 数值表；
#   (2) 四个精确极大值 (n=3,4,5,6) 的 60 位符号/高精度复核；
#   (3) 极限锥族函数 L_inf 的高精度驻点、导数残差与端点值；
#   (4) lim T_n/n 的闭式 (2*sqrt3+pi/3)/pi 与收敛速率；
#   (5) 关键常数与两处余量。
import numpy as np
from mpmath import mp, mpf, sqrt, acos, pi, diff, findroot
mp.dps = 60

sqrt5 = np.sqrt(5.0)
rho0  = (9+sqrt5+2*np.sqrt(35))/18

# ---------- (1) ----------
def maxL(n, M=200001):
    d = np.arange(1, n)
    c = np.cos(2*np.pi*d/n)
    def L(sv):
        sv = np.atleast_1d(np.asarray(sv, float))
        A = np.abs((1-sv)[:, None]*c[None, :] + sv[:, None]).sum(axis=1)
        Rn = A
        return ((1-3*sv)*(1+sv) + 4*sv*(1+Rn)/n)/(1-sv)**2
    xs = np.linspace(1e-12, 1/3, M)
    v = L(xs); i = int(np.argmax(v))
    lo, hi = xs[max(i-2,0)], xs[min(i+2,M-1)]
    for _ in range(120):
        m1 = lo+(hi-lo)/3; m2 = hi-(hi-lo)/3
        if L(m1)[0] < L(m2)[0]: lo = m1
        else: hi = m2
    s = (lo+hi)/2
    return s, float(L(s)[0])

print("=== (1) n=3..24 的 max_s L_w(s,n)（向量化网格 2e5 + 三分精修）===")
res = {}
for n in range(3, 25):
    s, v = maxL(n); res[n] = (s, v)
    print("  n=%2d  s*=%.12f  max=%.9f" % (n, s, v))
tail = [res[n][1] for n in range(7, 25)]
print("  n=7..24 的最大值 = %.6f (在 n=%d)" % (max(tail), 7+int(np.argmax(tail))))
print("  片段所列: 1.463854,1.485702,1.460603,1.456940,1.473843,1.464273,1.454849,")
print("            1.468392,1.464989,1.454110,1.465352,1.464826,1.456417,1.463446,")
print("            1.464391,1.459143,1.462151,1.463886")
mine = ",".join("%.6f" % res[n][1] for n in range(7, 25))
print("  本次重算:", mine)

print()
print("=== (2) 四个精确极大值（mpmath 60 位）===")
def Rn_mp(s, n):
    return sum(abs((1-s)*mp.cos(2*pi*d/n) + s) for d in range(1, n))
def Lw_mp(s, n):
    return ((1-3*s)*(1+s) + 4*s*(1+Rn_mp(s, n))/n)/(1-s)**2
exact = {3:(mpf(1)/5, mpf(4)/3), 4:(mpf(1)/3, mpf(3)/2),
         5:((1+2*sqrt(5))/19, (13+sqrt(5))/10), 6:(mpf(1)/4, mpf(13)/9)}
for n in (3,4,5,6):
    s, v = exact[n]; got = Lw_mp(s, n)
    print("  n=%d  s*=%-22s L=%-22s |L-v|=%s" % (n, mp.nstr(s,20), mp.nstr(got,20), mp.nstr(abs(got-v),2)))
print("  s0=(1+2*sqrt5)/19 = %s ; phi/(4+phi) = %s ; 差=%s"
      % (mp.nstr((1+2*sqrt(5))/19,20), mp.nstr(((1+sqrt(5))/2)/(4+(1+sqrt(5))/2),20),
         mp.nstr(abs((1+2*sqrt(5))/19 - ((1+sqrt(5))/2)/(4+(1+sqrt(5))/2)),2)))
print("  (13+sqrt5)/10 = %s ; (6+phi)/5 = %s"
      % (mp.nstr((13+sqrt(5))/10,20), mp.nstr((6+(1+sqrt(5))/2)/5,20)))

print()
print("=== (3) 极限函数 L_inf(s)=[(1-3s)(1+s)+4s r(s)]/(1-s)^2 ===")
def r_mp(s):
    return (2/pi)*(sqrt(1-2*s) + s*acos(-s/(1-s))) - s
def L_mp(s):
    return ((1-3*s)*(1+s) + 4*s*r_mp(s))/(1-s)**2
f = lambda s: diff(L_mp, s)
s_star = mp.re(findroot(f, mpf('0.2828')))
print("  s*        =", mp.nstr(s_star, 30))
print("  L_inf(s*) =", mp.nstr(L_mp(s_star), 30))
print("  导数残差 |L'(s*)| =", mp.nstr(abs(f(s_star)), 3))
e13 = L_mp(mpf(1)/3); cl = (2*sqrt(3)+pi/3)/pi
print("  L_inf(1/3)              =", mp.nstr(e13, 30))
print("  (2*sqrt3+pi/3)/pi       =", mp.nstr(cl, 30))
print("  |差| =", mp.nstr(abs(e13-cl), 3))
print("  L_inf(s*) < (13+sqrt5)/10 ?", L_mp(s_star) < (13+sqrt(5))/10,
      " 余量 =", mp.nstr((13+sqrt(5))/10 - L_mp(s_star), 6))
print("  L_inf(s*) < 1.485703 ?", L_mp(s_star) < mpf('1.485703'))

print()
print("=== (4) 环族 lim T_n/n ===")
def Tn(n):
    return float(np.sum(np.abs(2*np.cos(2*np.pi*np.arange(n)/n)+1)))
for n in [3,4,5,6,8,10**2,10**3,10**4,10**5,10**6]:
    print("  n=%7d  T_n/n = %.12f" % (n, Tn(n)/n))
print("  闭式 (2*sqrt3+pi/3)/pi = %.15f" % float(cl))
print("  T_5 = %.9f ; T_5/5 = %.9f ; 3n/2 - T_5 = %.6e" % (Tn(5), Tn(5)/5, 7.5-Tn(5)))
print("  T_n <= 3n/2 对 n=3..4000 的违例数 =",
      sum(1 for n in range(3,4001) if Tn(n) > 1.5*n + 1e-9),
      " 取等集 =", [n for n in range(3,4001) if abs(Tn(n)-1.5*n) < 1e-9])

print()
print("=== (5) 关键常数 ===")
for nm, v in [("rho0", rho0), ("2rho0-1", 2*rho0-1), ("(13+sqrt5)/10", (13+sqrt5)/10),
              ("(23+sqrt5)/20", (23+sqrt5)/20), ("STW Phi", (3+sqrt5)/4)]:
    print("  %-14s = %.15f" % (nm, v))
print("  余量 2rho0-1 - (13+sqrt5)/10 = %.6e" % (2*rho0-1-(13+sqrt5)/10))
print("  余量 rho0 - (23+sqrt5)/20    = %.6e" % (rho0-(23+sqrt5)/20))
print("  STW 差距 Phi - rho0          = %.12f" % ((3+sqrt5)/4-rho0))
print("ALL CHECKS DONE")
