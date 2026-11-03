import io
main="进展报告_前三大特征值之和的精确上界.tex"
s=io.open(main,encoding="utf-8").read()
old=r"\subsection{分数度障碍：$L(N)$ 对 $N\ge5$ 不可达；整数化阶梯}\label{sec:ladder}"
new=r"\subsection{分数度障碍：$L(N)$ 对 $N\ge5$ 不可达；整数化阶梯}\label{sec:fracdeg}"
assert s.count(old)==1
s=s.replace(old,new)
io.open(main,"w",encoding="utf-8").write(s); print("renamed sec:ladder -> sec:fracdeg")
