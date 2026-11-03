import io,re
main="进展报告_前三大特征值之和的精确上界.tex"
L=io.open(main,encoding="utf-8").read().split("\n")
for i,l in enumerate(L):
    if "rmk:tab19read" in l or "eq:chamber6" in l or "label{sec:pq6}" in l or "rmk:deadroutes" in l or "label{sec:obstacles}" in l or "label{sec:honest}" in l or "一句话总结" in l:
        print(i+1, l[:150])
