import io
s=io.open("_frag_r25.tex",encoding="utf-8").read().split("\n")
print("### R25 status section (from line 500)")
for i in range(499,len(s)): print(i+1,":",s[i])
