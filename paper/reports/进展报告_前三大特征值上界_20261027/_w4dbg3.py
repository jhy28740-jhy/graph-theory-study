import io
L=io.open("_frag_rowvar_4b.tex",encoding="utf-8").read().split("\n")
for i in range(101,122): print(i,repr(L[i]))
print("....")
for i in range(138,152): print(i,repr(L[i]))
