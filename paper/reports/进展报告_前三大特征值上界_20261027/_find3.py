import io,re
L=io.open("进展报告_前三大特征值之和的精确上界.tex",encoding="utf-8").read().split("\n")
for i,l in enumerate(L[:215]):
    if re.search(r"^(第[一二三四五六七八九]|\\end\{abstract\}|.*第六，|.*第七，)", l.strip()):
        print(i+1, l[:170])
