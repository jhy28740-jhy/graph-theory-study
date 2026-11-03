import io
p="_frag_rowvar_4.tex"; s=io.open(p,encoding="utf-8").read()
s=s.replace("第~\\ref{sec:champ}\\ 节给出了它的构造","定理~\\ref{thm:champ} 给出了它的构造")
s=s.replace("第~\\ref{sec:lower}\\ 节的定理~\\ref{thm:lower} 说","定理~\\ref{thm:lower} 说")
assert "sec:champ" not in s and "sec:lower" not in s
io.open(p,"w",encoding="utf-8").write(s); print("ok")
