import re,io
p='进展报告_前三大特征值之和的精确上界.tex'
s=open(p,encoding='utf-8').read()
bad='两者都不再需要「整条矩核直线上 $Q_+\\succeq0」这种难以处理的断言。'
good='两者都不再需要「整条矩核直线上 $Q_+\\succeq0$」这种难以处理的断言。'
assert bad in s, 'not found'
s=s.replace(bad,good,1)
open(p,'w',encoding='utf-8').write(s)
print('fixed 1')
# scan whole file for lines with odd number of unescaped $ (excluding $$ display lines)
lines=s.split('\n')
for i,l in enumerate(lines,1):
    t=re.sub(r'(?<!\\)\\\$','',l)          # drop escaped \$
    t=re.sub(r'\$\$','',t)                 # drop display math delimiters
    n=t.count('$')
    if n%2==1:
        print('ODD',i,l[:160])
