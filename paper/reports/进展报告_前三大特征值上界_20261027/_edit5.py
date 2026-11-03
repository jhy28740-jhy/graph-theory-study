import io
main="进展报告_前三大特征值之和的精确上界.tex"
s=io.open(main,encoding="utf-8").read()
def rep(old,new,cnt=1):
    global s
    assert s.count(old)==cnt, (s.count(old), old[:60])
    s=s.replace(old,new,cnt)

rep(r"""对这 $17$ 类，任何自洽见证都必须是正权的，故「$\notin\mathcal S_6^{\mathrm{int}}$」\emph{等价}于
「$\Sigma$ 完全不是腔掩码」，即引理~\ref{lem:chamber} 的纯几何陈述""",
r"""对这 $17$ 类，任何自洽见证都必须是正权的，故「$\Sigma$ 完全不是腔掩码」\emph{蕴含}「$\Sigma\notin\mathcal S_6^{\mathrm{int}}$」
（\textbf{反方向不成立}：$K_{2,2,2}$ 是腔掩码却与 $F$ 无关，见命题~\ref{prop:chamberdead}），即引理~\ref{lem:chamber} 的纯几何陈述""")

rep(r"""一个 $\mathrm{Gr}(3,6)$（$9$ 维）上的问题，与 $F$ 的极大化点无关。""",
r"""一个 $\mathrm{Gr}(3,6)$（$9$ 维）上的问题，与 $F$ 的极大化点无关。

\smallskip
\noindent\fbox{\parbox{0.95\linewidth}{\small\textbf{第 12 轮更新（否证）}：
上面这条路线已经\textbf{死掉}。命题~\ref{prop:chamberdead} 给出一个精确有理证书
$Q=\frac12(I_6-M)\in\mathcal P_3(6)$，使 $1\{Q\ge0\}=I_6+A(K_{2,2,2})$，
即 \eqref{eq:chamber6} 对 $\Sigma=I+A(K_{2,2,2})$ 为假；数值扫描进一步显示
$19$ 个未闭合类\textbf{全部}是腔掩码（符号违背量 $\le5\times10^{-19}$）。
因此 \eqref{eq:chamber6} 不能作为 \eqref{eq:s6int} 的后件使用，
必须换成第~\ref{sec:chamberdead} 小节末尾指出的\emph{耦合}条件
$\Sigma=1\{\Pi_3(D_p^{1/2}\Sigma D_p^{1/2})\ge0\}$。详见第~\ref{sec:rowvar} 节。}}""")

rep(r"""\subsection{更新后的优先队列}\label{sec:pq6}

按「期望收益 $\times$ 成功概率 $/$ 代价」排序，已同步写入 \texttt{研究方法优先级协议.md}。""",
r"""\subsection{更新后的优先队列（第六版，\textbf{已被第~\ref{sec:pq7} 小节取代}）}\label{sec:pq6}

\textbf{阅读提示（第 12 轮）}：本小节保留原文以记录思路演化，但\textbf{不再有效}。
下面的 P1 已被命题~\ref{prop:chamberdead} 否证（$19/19$ 类都是腔掩码），
P6 已被引理~\ref{lem:slope1} 吸收为与主猜想逐字等价。现行队列见第~\ref{sec:pq7} 小节。

按「期望收益 $\times$ 成功概率 $/$ 代价」排序，已同步写入 \texttt{研究方法优先级协议.md}。""")

rep(r"""\item \textbf{P1（队首）—— $N=6$ 的 $17$ 个无支配点类：证 \eqref{eq:chamber6}，即 $\Sigma\ne1\{Q\ge0\}$ 对一切 $Q\in\mathcal P_3(6)$。}""",
r"""\item \textbf{P1（队首，\textcolor{red}{第 12 轮已否证并弹出}）—— $N=6$ 的 $17$ 个无支配点类：证 \eqref{eq:chamber6}，即 $\Sigma\ne1\{Q\ge0\}$ 对一切 $Q\in\mathcal P_3(6)$。\ \textcolor{red}{【$19/19$ 类都是腔掩码，命题~\ref{prop:chamberdead}】}}""")
io.open(main,"w",encoding="utf-8").write(s)
print("edits applied; lines:",len(s.split("\n")))
