# -*- coding: utf-8 -*-
tex="进展报告_前三大特征值之和的精确上界.tex"
s=open(tex,encoding="utf-8").read()
dich=r"""\begin{corollary}[对齐二分，已证]\label{cor:r29dich}
设 $x$ 为 $Q_+$ 的 Perron 向量，则
\[
\mathrm{Neg}=\lambda_{\max}(Q_-)\ \Longleftrightarrow\ Q_-x=\mathrm{Neg}\cdot x
\quad\Longrightarrow\quad \mathcal A\in\{0,1\}.
\]
即 $x$ 同时为 $Q_-$ 的顶特征向量当且仅当偶极子饱和 $\mathcal A\in\{0,1\}$；等价地，凡 $0<\mathcal A<1$
（含冠军 $\mathcal A=0.9766$），$x$ 与 $Q_-$ 的顶特征方向\textbf{严格错位}、$\mathrm{Neg}<\lambda_{\max}(Q_-)$。
\end{corollary}
\begin{proof}
$\mathrm{Neg}=x^\top Q_-x\le\lambda_{\max}(Q_-)\|x\|^2=\lambda_{\max}(Q_-)$，取等（Rayleigh--Ritz）当且仅当
$Q_-x=\lambda_{\max}(Q_-)x=\mathrm{Neg}\,x$；此时 $\|Q_-x\|^2=\mathrm{Neg}^2$，代入定理~\ref{thm:r29align}~第二式得 $\mathcal A(1-\mathcal A)=0$。
\end{proof}
\begin{remark}
推论~\ref{cor:r29dich}~把冠军置于「$\mathcal A$ 接近 $1$ 但严格小于 $1$、$x$ 与 $Q_-$ 顶方向轻微错位」的状态：
错位量 $\lambda_{\max}(Q_-)-\mathrm{Neg}=0.379920-0.304973=0.074947$，正对应偶极子亏损 $\mathcal A(1-\mathcal A)=0.022857$。
$300$ 个随机秩三投影上，$0<\mathcal A<1$ 时 $\mathrm{Neg}<\lambda_{\max}(Q_-)$ 无一例外（最小间隔 $0.0185$，\texttt{\_r29f.py}）。
\end{remark}

"""
a1="\\subsection{解耦型全局路线的精确否证}\\label{sec:r29:weyl}"
assert s.count(a1)==1,("a1",s.count(a1))
s=s.replace(a1, dich+a1, 1)

pertable=r"""另对\textbf{单轴截面}做直接 P1 探针（\texttt{\_r29f.py}）：固定 $a$，锥族值 $F_{\rm cone}(a)=\max_{n,z}\rho(\Lam(a))$
由 $z$-网格（步长 $\sim7\times10^{-4}$、$n=3..8$）求得：
\begin{center}\small
\begin{tabular}{c|ccccccc}
\toprule
$a$ & $0.34$ & $0.40$ & $0.50$ & $0.5946$ & $0.70$ & $0.80$ & $0.90$\\
\midrule
$F_{\rm cone}(a)$ & $1.266773$ & $1.272985$ & $1.279617$ & $\mathbf{1.281568}$ & $1.279498$ & $1.274592$ & $1.268261$\\
最优锥 $n$ & $5$ & $5$ & $5$ & $5$ & $5$ & $5$ & $5$\\
\bottomrule
\end{tabular}
\end{center}
$F_{\rm cone}(a)$ 在 $a^*=0.5946$ 处取到 $\rho_0$、其余 $a$ 严格更小，且最优锥\textbf{恒为} $n=5$，
独立再确认定理~\ref{thm:r23uni}（单轴锥族 $\max F=\rho_0$）。
\emph{诚实标注}：同脚本对\textbf{一般帧}的多起点投影梯度（含 $\Lam(a)$-加权）因 $x(\Lam)$ 对 $V$ 的依赖采用了近似梯度而\textbf{欠收敛}，
所得 $F_{\rm gen}(a)\approx1.245$--$1.250$ 仅为\textbf{下界}、不足以佐证 P1；P1 的可靠数值证据是上文对
$\rho((VV^\top)_+)$（不含 $\Lam$-加权、梯度精确）的全局搜索（$N=4..9$ 均 $\le\rho_0$、极值为 $6$ 原子 $5$-锥）。

"""
a2="\\subsection{本轮状态、诚实标签与方法优先队列（第十六版）}\\label{sec:r29:status}"
assert s.count(a2)==1,("a2",s.count(a2))
s=s.replace(a2, pertable+a2, 1)

# add dichotomy to status table row + repro list mention
s=s.replace("秩一截面 $Q_-x=\\mathrm{Neg}\\,x$（共享特征向量刚性） & \\textbf{已证（精确）} & 推论~\\ref{cor:r29rank1}\\\\",
            "秩一截面 $Q_-x=\\mathrm{Neg}\\,x$（共享特征向量刚性） & \\textbf{已证（精确）} & 推论~\\ref{cor:r29rank1}\\\\\n对齐二分 $\\mathrm{Neg}=\\lambda_{\\max}(Q_-)\\Rightarrow\\mathcal A\\in\\{0,1\\}$ & \\textbf{已证（精确）} & 推论~\\ref{cor:r29dich}\\\\",1)
s=s.replace("\\texttt{\\_r29e.py}（$400$ 随机投影上恒等式与必要条件的普适核验）。",
            "\\texttt{\\_r29e.py}（$400$ 随机投影上恒等式与必要条件的普适核验）、\n\\texttt{\\_r29f.py}（单轴截面直接 P1 探针 $F_{\\rm cone}(a)$、对齐二分核验；输出 \\texttt{\\_r29f\\_out.txt}）。",1)
open(tex,"w",encoding="utf-8").write(s)
print("ADDENDUM OK. new len",len(s)," cor:r29dich count:",s.count("cor:r29dich"))