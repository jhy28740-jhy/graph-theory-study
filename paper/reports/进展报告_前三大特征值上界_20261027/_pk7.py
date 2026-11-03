import io
s=io.open("_frag_r26.tex",encoding="utf-8").read().split("\n")
print("### R26 status section")
for i in range(489,len(s)): print(i+1,":",s[i])
