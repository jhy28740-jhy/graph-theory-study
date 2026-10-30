# -*- coding: utf-8 -*-
tex="进展报告_前三大特征值之和的精确上界.tex"
s=open(tex,encoding="utf-8").read()
# 1) version bump in \date
old1="第十二版，在第十一版基础上再增一节（\\S\\ref{sec:r28}，第 28 轮）。"
new1="第十三版，在第十二版基础上再增一节（\\S\\ref{sec:r29}，第 29 轮）。"
assert s.count(old1)==1, ("v1",s.count(old1))
s=s.replace(old1,new1,1)
# 2) insert R29 recap right before the "第十版原有的" tail of \date
tail="第十版原有的 \\S\\ref{sec:r25}"
assert s.count(tail)==1, ("tail",s.count(tail))
recap=("\\S\\ref{sec:r29}（第 29 轮）：导出对\\emph{任意}秩三投影成立的 \\textbf{Perron--钝角对齐恒等式} "
"$\\|Q_-x\\|^2=\\rho(1-\\rho)+\\mathrm{Neg}(2\\rho-1)=\\mathrm{Neg}^2+\\mathcal A(1-\\mathcal A)$（定理~\\ref{thm:r29align}），"
"给出近极值元必要条件 $\\mathrm{Neg}\\ge\\frac{\\rho(\\rho-1)}{2\\rho-1}$（推论~\\ref{cor:r29nec}）与秩一截面的共享特征向量刚性 "
"$Q_-x=\\mathrm{Neg}\\,x$（推论~\\ref{cor:r29rank1}）；\\textbf{精确否证}三条解耦型全局路线——Weyl 界 $1+\\lambda_{\\max}(Q_-)$"
"（冠军钝角图 $=C_5$、$\\lambda_{\\max}(Q_-)=2w=0.379920$，超 $\\rho_0$ 达 $0.098352$，命题~\\ref{prop:r29weyl}）、"
"绝对值 Gram $\\rho(|Q|)=1.602534$、对角分裂 $1.327442$（命题~\\ref{prop:r29other}）——证明全局闭合\\emph{必须}经由共享 Perron 耦合；"
"并\\textbf{精确定位} P1$'$ 的矩对偶缺口（点态支配迫使 $\\alpha\\ge\\frac12$、自相关松弛 $\\ge2$，命题~\\ref{prop:r29pointwise}--\\ref{prop:r29gap}），"
"把 P1$'$ 收窄到唯一的 Toeplitz--PSD/Fejér--Riesz 路线；另用投影梯度独立再确认 $\\kappa_3=\\rho_0$、极值帧为 $6$ 原子 $5$-锥（$N=4..9$）。")
s=s.replace(tail, recap+tail, 1)
open(tex,"w",encoding="utf-8").write(s)
print("DATE UPDATED. new len",len(s))
print("第十三版 count:",s.count("第十三版")," sec:r29 in date:", "第 29 轮" in s)