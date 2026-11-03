import io,re
s=io.open("进展报告_前三大特征值之和的精确上界.tex",encoding="utf-8").read().split("\n")
# locate eq:frame18
for i,l in enumerate(s,1):
    if "label{eq:frame18}" in l: print("frame18 at",i)
    if "label{eq:stiefelkkt}" in l: print("stiefelkkt at",i)
    if re.search(r"\\label\{(thm:r23uni|prop:r24gen|thm:r22id|rem:r26struct|rem:r26one|prop:r26conecl|eq:r26cl|lem:r26Rn|prop:r26rank1|lem:r26hemi)\}", l): print(l.strip()[:40],"->",i)
