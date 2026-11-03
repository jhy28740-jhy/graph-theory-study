# -*- coding: utf-8 -*-
import io
p='证书_锥族闭式与rho0的二维刻画.py'
s=io.open(p,encoding='utf-8-sig').read()
o1="chk(sp.simplify(2*(S5+1)/2*3 - (3+2*S35+S5)) == 0 or True, 'placeholder')\n"
n1="chk(sp.simplify((51*S5-21)/16 + (225-135*S5)/80 - 3*(1+S5)/2) == 0,\n    '(51sqrt5-21)/16 + 9beta/5 = 3phi（符号精确，即 L<=3phi 的锐常数）')\n"
assert o1 in s; s=s.replace(o1,n1)
o2="chk(sp.simplify(sum(c_exact(2*k) for k in range(2, 12))) == sp.Rational(-1, 16) - sum(c_exact(2*k) for k in range(12, 40), sp.Rational(0)) or True, 'placeholder2')\n"
assert o2 in s; s=s.replace(o2,"")
io.open(p,'w',encoding='utf-8').write(s)
print('ok')
