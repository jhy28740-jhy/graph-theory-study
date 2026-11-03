import io
L=io.open("_frag_rowvar_4b.tex",encoding="utf-8").read().split("\n")
for i in range(101,119): print(i,repr(L[i]))
