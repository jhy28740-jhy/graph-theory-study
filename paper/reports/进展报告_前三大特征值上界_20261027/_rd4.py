import io,re
s=io.open("进展报告_前三大特征值之和的精确上界.tex",encoding="utf-8").read().split("\n")
for i in range(5680,5760): print(i+1,":",s[i])
