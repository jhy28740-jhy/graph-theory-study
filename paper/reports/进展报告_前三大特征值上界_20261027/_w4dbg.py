import io
s=io.open("_frag_rowvar_4b.tex",encoding="utf-8").read().split("\n")
for i,l in enumerate(s):
    if "prop:corank1" in l or "lambda_{\\max}(Q_+)=" in l or "end{proposition}" in l or "end{enumerate}" in l or "eq:fm1" in l:
        print(i,repr(l))
