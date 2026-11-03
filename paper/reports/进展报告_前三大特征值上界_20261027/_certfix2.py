import io
p="证书_行方差松弛与天花板.py"; s=io.open(p,encoding="utf-8").read()
old=r'''irr=[]
for n in range(4,41):
    fr=sp.Rational(3*(n-1),n-3) if n!=3 else None
    if fr is None: continue
    irr.append((n,_not_perfect_square_rat(fr.p,fr.q)))
chk("kappa_N 无理 <=> n*(N) 无理，对一切 5<=N<=40 (prop:fracdeg)",
    all(v for n,v in irr if n>=5) and not dict(irr)[4],
    "N=4 有理(True), N=5..40 全部无理")'''
new=r'''irr={n:_not_perfect_square_rat(sp.Rational(3*(n-1),n-3).p,sp.Rational(3*(n-1),n-3).q) for n in range(4,3001)}
excp=sorted([n for n,v in irr.items() if not v])
print("   kappa_N 为有理数的 N<=3000：",excp)
chk("n* 无理 <=> kappa_N 无理；5<=N<=3000 中的例外集恰为 {9,28,99,364,1353,5044} (prop:fracdeg 修正版)",
    [n for n in excp if n>=5]==[9,28,99,364,1353,5044],"例外集=%s（含 N=4）"%excp)'''
assert old in s; s=s.replace(old,new)
old2=r'''pent=[t for t in tab if t[0]==0 and t[1]==(1,)]+[t for t in tab if t[0]==0 and t[1]==(2,)]
chk("冠军图案类（b=0, m={1} 或 {2}）天花板 = rho0 + 3.2128e-3",
    all(abs(t[2]-1.284780995477)<1e-6 for t in pent),"%.12f  超出 rho0 = %+.4e"%(pent[0][2],pent[0][2]-rho0))
chk("m={1} 与 m={2} 两类天花板相同（五角形/五角星对偶）",abs(pent[0][2]-pent[1][2])<1e-9)
chk("其余 6/8 类被 REL 单独排除（<= max(rho0, 1.2431)）",
    all(t[2]<1.25 for t in tab if t not in pent))'''
new2=r'''pent_idx=[i for i,t in enumerate(tab) if t[0]==0 and t[1] in ((1,),(2,))]
pent=[tab[i] for i in pent_idx]
chk("冠军图案类（b=0, m={1} 或 {2}）天花板 = rho0 + 3.2128e-3",
    all(abs(t[2]-1.284780995477)<1e-6 for t in pent),"%.12f  超出 rho0 = %+.4e"%(pent[0][2],pent[0][2]-rho0))
chk("m={1} 与 m={2} 两类天花板相同（五角形/五角星对偶）",abs(pent[0][2]-pent[1][2])<1e-9)
rest=[tab[i] for i in range(len(tab)) if i not in pent_idx]
chk("其余 6/8 类被 REL 单独排除（天花板 < 1.25 < rho0+... ）",all(t[2]<1.25 for t in rest),
    "最大值 = %.9f"%max(t[2] for t in rest))'''
assert old2 in s; s=s.replace(old2,new2)
io.open(p,"w",encoding="utf-8").write(s); print("patched")
