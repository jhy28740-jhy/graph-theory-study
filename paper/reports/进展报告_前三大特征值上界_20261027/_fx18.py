p='进展报告_前三大特征值之和的精确上界.tex'
s=open(p,encoding='utf-8').read()
assert s.count('conj:main')==1
s=s.replace('conj:main','conj:A1',1)
open(p,'w',encoding='utf-8').write(s)
print('ok')
