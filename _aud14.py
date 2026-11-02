# -*- coding: utf-8 -*-
"""_aud14.py : 在冠军律上精确分解 (certificate-final) 的 8 个非负项，
以及 Q 的 6 个子项，量化 "常数缺口 c-rho_0 究竟藏在哪一项里"。
全部量沿用 _aud9.py 的定义；这里只读不写。"""
import numpy as np
from numpy.linalg import eigh, norm
from fractions import Fraction as F
sq5=np.sqrt(5.0)
d0=-25624284539519/53655078400000+sq5/2-41*np.sqrt(3)/120
d1=d0+103/9000000; d=d1+1/3000000
c=(3+sq5)/4-d/2; cw=(9+sq5+2*np.sqrt(35))/18
I3=np.eye(3)
a_=(1+sq5)/6; b_=(7*sq5-5)/8; rho=(15*sq5-25)/16
a3=F(12928,41503); a4=F(10240,41503)
ell=F(937,1400); eta=F(1,10); th=F(6,25)
BIG=F(725029049,278853120)
kappa=3/125000

d_ij=np.eye(3)
OMEGA=(np.einsum('ij,kl->ijkl',d_ij,d_ij)+np.einsum('ik,jl->ijkl',d_ij,d_ij)
       +np.einsum('il,jk->ijkl',d_ij,d_ij))/15.0
def L4(Z):
    return (np.einsum('ij,kl->ijkl',d_ij,Z)+np.einsum('ik,jl->ijkl',d_ij,Z)
            +np.einsum('il,jk->ijkl',d_ij,Z)+np.einsum('jk,il->ijkl',d_ij,Z)
            +np.einsum('jl,ik->ijkl',d_ij,Z)+np.einsum('kl,ij->ijkl',d_ij,Z))/7.0
def Psi4(u):
    return np.einsum('i,j,k,l->ijkl',u,u,u,u)-OMEGA-L4(np.outer(u,u)-I3/3.0)
def sig(w):
    num=8352+43806*w+54016*w**2-12928*w**3-10240*w**4
    return num/41503-np.abs(w)-w
def h0(z): return (1-z)*(sq5*z-1)**2*(sq5*z+sq5+2)/(2*(sq5+1)**2)
def zstar(t): return t**4-F(85,64)*t**2+F(15,128)*t+F(27,128)

b=(7+5*np.sqrt(7))/126; a=1-5*b
wt=np.array([a]+[b]*5)
A6=np.zeros((6,6))
for i in range(5): A6[0,i+1]=A6[i+1,0]=1
for i in range(5):
    j=(i+1)%5; A6[i+1,j+1]=A6[j+1,i+1]=1
Mq=np.diag(np.sqrt(wt))@(A6+np.eye(6))@np.diag(np.sqrt(wt))
ev,Vv=eigh(Mq); o=np.argsort(-ev); Vv=Vv[:,o]
Vm=Vv[:,:3]/np.sqrt(wt)[:,None]
if Vm[:,0].sum()<0: Vm[:,0]*=-1
p=wt

Ri=norm(Vm,axis=1); Ui=Vm/np.where(Ri>0,Ri,1)[:,None]
w=Ui@Ui.T; RR=Ri[:,None]*Ri[None,:]; pp=p[:,None]*p[None,:]
A=(pp*np.abs(Vm@Vm.T)).sum(); m=(p[:,None]*Vm).sum(0); B=m@m; e=1-B
C=(p*Ri).sum(); V=3-C**2; D=(3+sq5)/2-A-B
M=np.einsum('i,ij,ik->jk',p,Vm,Ui); Z=M-(C/3)*I3
J=1-np.sum(M*M)
S=(pp*RR*h0(np.abs(w))).sum()
T4=np.einsum('i,imnop->mnop',p*Ri,np.stack([Psi4(u) for u in Ui]))
K2=float(np.sum(T4**2)); K=np.sqrt(K2)
P=(pp*RR*sig(w)).sum()
V_=3-C**2; J_=J
r_rad=P-(3/50)*(V-9*J)
x=np.sqrt(B); vh=m/x; q=C/(3*x); q2=q*q
Zs=(pp*RR*np.vectorize(lambda t: float(zstar(t)))(w)).sum()
Ccal=Zs/125-3*C**2/125000
r_ct=P-Ccal
chi=(p**2*Ri**2).sum()
kchi=kappa*chi; r_off=r_ct-kchi
r_tr=r_rad-P/9

# ---- Q 的六个子项 ----
se=np.sqrt(e); sv=np.sqrt(V); sj=np.sqrt(J)
Q1=(13/200)*(sv-float(937/273)*q*se)**2
Q2=(19/50)*(sj-float(2811/2128)*q*se)**2
Q3=(183/200-float(BIG)*q2)*e
Q4=float(F(148831,181104000))*e
Q5=float(F(42847433,5229378000000))*V
Q6=float(F(65082103,20917512000))*J
Q=Q1+Q2+Q3+Q4+Q5+Q6

terms=[("d (常数项)",d),("Qhat = Q - V/9e6",Q-V/9000000),("R_0",0.0),
       ("ell*r_reg",0.0),("K^2/1125",K2/1125),("r_tr",r_tr),
       ("kappa/9*(chi-C^2/24)",kappa/9*(chi-C**2/24)),("r_off/9",r_off/9)]
slack=D-d
print(f"c      = {c:.15f}")
print(f"rho_0  = {cw:.15f}")
print(f"c-rho_0= {c-cw:.15f}   2*(c-rho_0) = {2*(c-cw):.15f}")
print(f"D(champion) = {D:.15f}   d = {d:.15f}   slack D-d = {slack:.15f}")
print(f"  (slack 应等于 2*(c-rho_0) = {2*(c-cw):.12f})  一致? {abs(slack-2*(c-cw))<1e-12}")
print()
print("=== certificate-final 各项（R_0 与 ell*r_reg 见下）===")
r_reg=None
# r_reg = Lbar - L ;  需要 L, Lbar
t=Ui@vh; Y=Ri*t
w0=None
# L = E[R(t-q0)^2], q0=1/sqrt3
q0=1/np.sqrt(3)
L=(p*Ri*(t-q0)**2).sum()
# F0,G0
Zv=Z@vh; Jv=(p*( (np.stack([ (Ui[i]@vh - (M@vh)@Ui[i]) for i in range(len(p)) ]) )**2)*1.0).sum() if False else None
# 用 xi = U - M X 直接算
Xi=Vm  # X 的取值就是 v_i (质量 p_i)
xi=Ui-(M@Xi.T).T
w0v=(p[:,None]*xi).sum(0)
Jv=(p*(xi@vh)**2).sum()
F0=float(np.sum(Zv**2))+q2*e+Jv-2*q*(vh@w0v)
ERY=(p*Ri*Y).sum()
G0=1-2*q*ERY+3*q2
Lq=(p*Ri*(t-q)**2).sum()
Lbar=0.75*F0+G0/3
r_reg=Lbar-L
R0=None
# R_0 来自 certificate-start: D = d0 + E e + Vc V + Jc J - ell Lbar + R0 + r_rad + ell r_reg
E_=float(F(5183683,5659500)); Vc=float(F(145858703,581042000)); Jc=float(F(7989533,14526050))
R0=D-(d0+E_*e+Vc*V+Jc*J-float(ell)*Lbar+r_rad+float(ell)*r_reg)
print(f"  r_reg = Lbar-L           = {r_reg:.12f}   ell*r_reg = {float(ell)*r_reg:.12f}")
print(f"  R_0 (由 certificate-start 反解) = {R0:.12f}")
terms[2]=("R_0",R0); terms[3]=("ell*r_reg",float(ell)*r_reg)
tot=0.0
for nm,val in terms:
    tot+=val
print()
print(f"{'项':32s} {'数值':>16s} {'占 slack 比例':>14s}")
for nm,val in terms[1:]:
    print(f"{nm:32s} {val:16.12f} {100*val/slack:13.2f}%")
print(f"{'合计(应=D)':32s} {tot:16.12f}   vs D = {D:.12f}  残差 {tot-D:+.3e}")
print()
print("=== Q 的六个子项 ===")
for nm,val in [("Q1 (13/200)(sqrtV-937/273 q sqrte)^2",Q1),("Q2 (19/50)(sqrtJ-2811/2128 q sqrte)^2",Q2),
               ("Q3 (183/200-BIG q^2) e",Q3),("Q4 148831/181104000 e",Q4),
               ("Q5 42847433/5229378000000 V",Q5),("Q6 65082103/20917512000 J",Q6)]:
    print(f"  {nm:44s} {val:.12f}   占Q {100*val/Q:6.2f}%   占slack {100*val/slack:6.2f}%")
print(f"  {'Q 合计':44s} {Q:.12f}")
print()
print("=== r_tr 的来源分解 ===")
print(f"  P = E[RR' sigma]            = {P:.12f}")
print(f"  (8/9)P                      = {8*P/9:.12f}")
print(f"  -E[RR' z_rad g]  (=r_tr-8P/9) = {r_tr-8*P/9:.12f}")
print(f"  r_tr                        = {r_tr:.12f}")
print(f"  V-9J                        = {V-9*J:.12f}   (<0 => 径向支付为负, 是 r_tr 的主要来源)")
print(f"  2V-9||Z||^2 = 9J-V          = {9*J-V:.12f}")
print()
print("=== Q3 系数的松紧 ===")
print(f"  183/200 - BIG*q^2 (冠军 q^2={q2:.9f}) = {183/200-float(BIG)*q2:.12f}")
print(f"  证书能用的下界 (q^2<=20/57)          = {float(F(53754367,19868284800)):.12f}")
print(f"  松弛倍数                              = {(183/200-float(BIG)*q2)/float(F(53754367,19868284800)):.2f}x")
print()
print("=== 小分支窗口 ===")
thr=((3+sq5)/2-1/20)/2
print(f"  D*<1/20  <=>  M > {thr:.15f}")
print(f"  rho_0 = {cw:.15f} < {thr:.9f} ?  {cw<thr}   => 若 M=rho_0, 小分支为空")
print(f"  c     = {c:.15f} > {thr:.9f} ?  {c>thr}")
print(f"  小分支窗口宽度 = c - {thr:.9f} = {c-thr:.12f}")
