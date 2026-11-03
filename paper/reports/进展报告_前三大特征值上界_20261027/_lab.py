import io,re
p='进展报告_前三大特征值之和的精确上界.tex'
s=io.open(p,encoding='utf-8').read()
lines=s.split('\n')
for i,l in enumerate(lines,1):
    if re.search(r'\\label\{(thm|lem|prop|cor|conj|def|rmk|eq|sec|app):',l):
        m=re.search(r'\\label\{([^}]+)\}',l)
        env=re.search(r'\\begin\{([^}]+)\}',l)
        print('%5d %-28s %s'%(i,m.group(1),l.strip()[:60]))
