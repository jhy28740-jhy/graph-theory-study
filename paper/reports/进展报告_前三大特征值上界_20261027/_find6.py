import io
L=io.open("进展报告_前三大特征值之和的精确上界.tex",encoding="utf-8").read().split("\n")
for i in range(5031,5050): print(i+1, L[i][:180])
print("---- obstacles 4440-4470 ----")
for i in range(4439,4470): print(i+1, L[i][:180])
