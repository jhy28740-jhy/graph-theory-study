# -*- coding: utf-8 -*-
import io
p='证书_锥族闭式与rho0的二维刻画.py'
s=io.open(p,encoding='utf-8-sig').read()
old="""def cat_float(k):
    return math.comb(2*k, k)/(k+1)
tail = sum(((-1)**(k+1))*(4*k+1)*cat_float(k)/(4*(2*k-1)*4.0**k) for k in range(2, 600))"""
new="""def catalan_ratio_terms(K):
    \"\"\"r_k = Cat_k/4^k（稳定递推），返回 c_{2k}=(-1)^{k+1}(4k+1)r_k/(4(2k-1))。\"\"\"
    out = {}; r = 0.25
    for k in range(1, K+1):
        out[k] = ((-1)**(k+1))*(4*k+1)*r/(4*(2*k-1))
        r *= (2*k+1)/(2.0*(k+2))
    return out
_ct = catalan_ratio_terms(600)
tail = sum(_ct[k] for k in range(2, 601))"""
assert old in s; s=s.replace(old,new)
io.open(p,'w',encoding='utf-8').write(s)
print('ok')
