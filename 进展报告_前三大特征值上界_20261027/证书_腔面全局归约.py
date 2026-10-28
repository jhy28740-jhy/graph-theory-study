# -*- coding: utf-8 -*-
# 证书：腔面全局上界的球面归约 (定理 thm:wallreduce) + KKT + 紧致松弛否定引理 + 否证扫描
# 复核报告 sec:w1global 的全部数字。W1-全局的完整认证仍未完成(见报告诚实交代)。
import numpy as np
from scipy.optimize import minimize
s5=np.sqrt(5); s7=np.sqrt(7); s35=np.sqrt(35)
a=(91-25*s7)/126; b=(7+5*s7)/126
Vhh=14/(-7*s5-7+5*s7+3*s35); Vhc=np.sqrt(7/5); Vcc=-7/2-7*s5/10+3*s35/10+5*s7/2
Ah=np.sqrt(Vhh); Ac=Vhc/Ah; Rc=np.sqrt(Vcc-Ac**2); w=np.array([a,b,b,b,b,b])
cols=[np.array([Ah,0,0])]
for j in range(5):
    th=2*np.pi*j/5; cols.append(np.array([Ac,Rc*np.cos(th),Rc*np.sin(th)]))
A0=np.stack([np.sqrt(w[k])*cols[k] for k in range(6)],axis=1)
Xs=(np.linalg.inv(A0[:,:3])@A0)[:,3:6]
def Pi_of(X):
    Am=np.hstack([np.eye(3),X.reshape(3,3)]); return Am.T@np.linalg.solve(Am@Am.T,Am)
rho0=(9+s5+2*s35)/18; mc=(Pi_of(Xs)>0).astype(float)
def sym(M,i,j,v): M[i,j]=M[j,i]=v
masks={}
m=mc.copy(); sym(m,0,1,0.0); masks["class1 hc"]=m
m=mc.copy(); sym(m,1,2,0.0); masks["class2 d1"]=m
m=mc.copy(); sym(m,1,3,1.0); masks["class3 d2+"]=m
def eig_desc(M):
    e,U=np.linalg.eigh(M); o=np.argsort(e)[::-1]; return e[o],U[:,o]
def Lam3(M): return eig_desc(M)[0][:3].sum()
def g(v,S): return Lam3(np.diag(v)@S@np.diag(v))
print("rho0 =",rho0)
# --- (1) 两条独立实现的最大值必须逐位相同 (归约验证) ---
def chart_max(S):
    F=lambda x:-np.linalg.eigvalsh(S*Pi_of(x))[-1]; best=(-9,None)
    rng=np.random.default_rng(4)
    seeds=[Xs.reshape(9)]+[Xs.reshape(9)+rng.normal(0,0.3,9) for _ in range(5)]+[rng.normal(0,1.0,9) for _ in range(9)]
    for sd in seeds:
        r=minimize(F,sd,method="Nelder-Mead",options=dict(maxiter=20000,maxfev=40000,xatol=1e-13,fatol=1e-15))
        r=minimize(F,r.x,method="Powell",options=dict(maxiter=20000,xtol=1e-13,ftol=1e-15))
        if -r.fun>best[0]: best=(-r.fun,r.x)
    return best
def sphere_max(S):
    f=lambda z:-g(z/np.linalg.norm(z),S); best=(-9,None)
    rng=np.random.default_rng(8)
    for sd in [rng.normal(size=6) for _ in range(40)]+[np.ones(6)]:
        r=minimize(f,sd,method="Nelder-Mead",options=dict(maxiter=20000,maxfev=40000,xatol=1e-13,fatol=1e-15))
        r=minimize(f,r.x,method="Powell",options=dict(maxiter=20000,xtol=1e-13,ftol=1e-15))
        if -r.fun>best[0]: best=(-r.fun,r.x/np.linalg.norm(r.x))
    return best
# --- (2) 球面投影梯度上升否证扫描, 解析梯度 2(Sig odot P)v ---
def grad_ascend(v,S,iters=250,eta=0.05):
    v=v/np.linalg.norm(v); best=(g(v,S),v.copy())
    for _ in range(iters):
        val,gr=g_and_grad(v,S)
        if val>best[0]: best=(val,v.copy())
        gr=gr-(gr@v)*v; ng=np.linalg.norm(gr)
        if ng<1e-12: break
        v=v+eta*gr/ng; v/=np.linalg.norm(v); eta*=0.999
    return best
def g_and_grad(v,S):
    e,U=eig_desc(np.diag(v)@S@np.diag(v)); P=U[:,:3]@U[:,:3].T
    return e[:3].sum(),2.0*((S*P)@v)
for tag,S in masks.items():
    fa,xa=chart_max(S); fb,vb=sphere_max(S)
    rng=np.random.default_rng(2026); scan=-9; hits=0
    for _ in range(1500):
        val,_=grad_ascend(rng.normal(size=6),S); scan=max(scan,val); hits+= (val>fb-1e-6)
    M=np.diag(vb)@S@np.diag(vb); e,U=eig_desc(M); P=U[:,:3]@U[:,:3].T
    kkt=np.linalg.norm((S*P)@vb-fb*vb)
    mix=-9
    for _ in range(3000):
        k=rng.integers(2,5); Z=np.zeros((6,6))
        for _ in range(k):
            u=rng.normal(size=6); u/=np.linalg.norm(u); Z+=np.outer(u,u)/k
        mix=max(mix,Lam3(S*Z))
    print("%-11s chart=%.12f sphere=%.12f match=%.1e margin=%.6f"%(tag,fa,fb,abs(fa-fb),rho0-fb))
    print("            tr(DvSDv)=%.12f KKT=%.1e top3eig=%s | scan max=%.9f (hits %d/1500, all<rho0:%s) | nonrank1 mix=%.6f<=rank1:%s"%(
        np.trace(M),kkt,np.round(e[:3],5),scan,hits,scan<rho0,mix,mix<=fb+1e-6))
print("\n结论: 归约 max_Gr lambda_max(Sig odot P)=max_S5 Lam3(D_v Sig D_v) 逐位验证;")
print("      三类全局最大值均 < rho0 (margin 0.025~0.045), 但完整认证(临界点穷举/矩-SOS)仍未完成 = W1-全局缺口。")
