import io
s=io.open("进展报告_前三大特征值之和的精确上界.tex",encoding="utf-8").read().split("\n")
for i in range(9295,9440): print(i+1,":",s[i])
