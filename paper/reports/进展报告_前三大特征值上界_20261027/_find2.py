import io
L=io.open("进展报告_前三大特征值之和的精确上界.tex",encoding="utf-8").read().split("\n")
for i in list(range(2993,3018))+list(range(3065,3080)):
    print(i+1, L[i][:190])
