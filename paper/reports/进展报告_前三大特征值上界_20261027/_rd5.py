import io,re
s=io.open("进展报告_前三大特征值之和的精确上界.tex",encoding="utf-8").read()
for m in re.finditer(r"L_w\s*:?=|L_w:=|\\label\{eq:Lw\}", s):
    ln=s[:m.start()].count("\n")+1
    print(ln,":",s.split("\n")[ln-1][:160])
