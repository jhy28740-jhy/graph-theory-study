import io,re
L=io.open("进展报告_前三大特征值之和的精确上界.tex",encoding="utf-8").read().split("\n")
st=None
for i,l in enumerate(L):
    if "label{sec:honest}" in l: st=i
    if "label{sec:obstacles}" in l: ob=i
    if "subsection{一句话总结}" in l: yj=i
    if "label{app:repro}" in l: ar=i
print("honest",st+1,"obstacles",ob+1,"yijuhua",yj+1,"apprepro",ar+1)
for i in range(yj-14,yj+40): print(i+1, L[i][:170])
