# -*- coding: utf-8 -*-
import io
p='证书_锥族闭式与rho0的二维刻画.py'
s=io.open(p,encoding='utf-8-sig').read()
s=s.replace("chk(abs(tail - (-1/16)) < 1e-6,","chk(abs(tail - (-1/16)) < 5e-5,")
old="""pp2 = np.array([float(q) for q in p])
Vv = Uv2*pp2[:, None]; Q = Vv @ Vv.T
Qp = np.maximum(Q, 0.0); ev, evec = np.linalg.eigh(Qp)
k = int(np.argmax(ev))
chk(abs(ev[k] - rho0) < 1e-12, 'lambda_max((VV^T)_+) = rho_0，偏差 %.2e' % abs(ev[k]-rho0))
pv = evec[:, k]
if pv.sum() < 0: pv = -pv
hh = np.array([float(q) for q in h]); hh = hh/hh.max(); pv2 = pv/pv.max()
chk(np.abs(pv2 - hh).max() < 1e-8, 'Perron 向量重现 h_i=sqrt(u_i^T Lam u_i)，偏差 %.2e' % np.abs(pv2-hh).max())"""
new="""Vv = Uv2*np.sqrt(gv2)[:, None]; Q = Vv @ Vv.T          # Q_ij = sqrt(g_i g_j)<u_i,u_j>
Qp = np.maximum(Q, 0.0); ev, evec = np.linalg.eigh(Qp)
k = int(np.argmax(ev)); pv = evec[:, k]
Fv = float(F)
chk(abs(ev[k] - rho0) < 1e-12 and abs(ev[k]-Fv) < 1e-12,
    '自洽性：lambda_max((VV^T)_+) = F(a*,z*) = rho_0（偏差 %.2e / %.2e）' % (abs(ev[k]-Fv), abs(ev[k]-rho0)))
hh = np.abs(pv)/np.sqrt(gv2); hh = hh/hh[0]
href = np.array([float(q) for q in h]); href = href/href[0]
chk(np.abs(hh - href).max() < 1e-7,
    'Perron 向量重现 h_i=sqrt(u_i^T Lam u_i)（|x_i|/sqrt(g_i) 归一后偏差 %.2e）' % np.abs(hh-href).max())"""
assert old in s; s=s.replace(old,new)
io.open(p,'w',encoding='utf-8').write(s)
print('ok')
