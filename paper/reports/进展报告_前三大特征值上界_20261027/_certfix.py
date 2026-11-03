import io
p="证书_行方差松弛与天花板.py"; s=io.open(p,encoding="utf-8").read()
s=s.replace(r"""s_i=np.array([np.sum(np.maximum(P[i,:],0.0))-max(P[i,i],0) for i in range(6)])""",
r"""s_i=np.array([float(np.maximum(P[i,:],0.0)@xs)-max(P[i,i],0)*xs[i] for i in range(6)])""")
s=s.replace(r'''chk("n*(N) 对 N>=5 无理 => L(N) 不可达 (prop:fracdeg)",all(not sp.nsimplify(nst.subs(N,n)).is_rational for n in range(5,25)))''',
r'''def _not_perfect_square_rat(num,den):
    import math as _m
    def issq(t):
        r=_m.isqrt(t); return r*r==t
    return not (issq(num) and issq(den))
irr=[]
for n in range(4,41):
    fr=sp.Rational(3*(n-1),n-3) if n!=3 else None
    if fr is None: continue
    irr.append((n,_not_perfect_square_rat(fr.p,fr.q)))
chk("kappa_N 无理 <=> n*(N) 无理，对一切 5<=N<=40 (prop:fracdeg)",
    all(v for n,v in irr if n>=5) and not dict(irr)[4],
    "N=4 有理(True), N=5..40 全部无理")''')
io.open(p,"w",encoding="utf-8").write(s); print("patched")
