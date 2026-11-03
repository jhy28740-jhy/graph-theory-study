# -*- coding: utf-8 -*-
import io
main="进展报告_前三大特征值之和的精确上界.tex"
s=io.open(main,encoding="utf-8").read()
def rep(old,new,cnt=1):
    global s
    assert s.count(old)==cnt, ("COUNT %d"%s.count(old), old[:70]); s=s.replace(old,new,cnt)
# abstract fix
rep(r"""命题~\ref{prop:fracdeg} 解释了它为什么不能更好：取等要求的钝角度 $n^*=N\kappa_N/(3(1+\kappa_N))$
对一切 $N\ge5$ 都是\emph{无理数}（$n^*(6)=\frac{5-\sqrt5}2$），所以 $\Phi-\rho_0$ 这个缺口的本质是\textbf{度的整性}，
不是 Cauchy--Schwarz 太粗。""",
r"""命题~\ref{prop:fracdeg} 解释了它为什么不能更好：取等\emph{强制}钝角图正则、度为 $n^*=N\kappa_N/(3(1+\kappa_N))$，
并强制 $Q$ 两值，故 $n^*$ 必须是整数；$n^*\in\Z$ 当且仅当 $\frac{3(N-1)}{N-3}$ 是 $\Q$ 中的平方，
这只对一支\emph{无穷但极稀疏}的 Pell 例外集 $\mathcal E=\{4,9,28,99,364,1353,5044,18819,\dots\}$ 发生
（$n^*(6)=\frac{5-\sqrt5}2$ 无理）。所以 $\Phi-\rho_0$ 这个缺口的本质是\textbf{度的整性}，不是 Cauchy--Schwarz 太粗。
定理~\ref{thm:unattain} 进一步处理例外集：可达性会强制钝角图是\emph{强正则图}（或等团不交并），
而强正则图的重数整性对 $\mathcal E$ 中 $9\le N\le7.1\times10^8$ 的每一个成员都失败，
故 \textbf{$L(N)$ 对一切 $5\le N\le7.1\times10^8$ 不可达，$N=4$ 是唯一取等的阶数}。
（附带一条自我修正：命题~\ref{prop:fracdeg} 的第一版断言「$n^*$ 对一切 $N\ge5$ 无理」，
$n^*(9)=2$ 是反例，见注~\ref{rmk:fracdegfix}。）""")
# honest list fix
rep(r"""命题~\ref{prop:fracdeg}（$n^*$ 的无理性）；命题~\ref{prop:ladder6}（$N=6$ 整数化阶梯的六个闭式）；""",
r"""命题~\ref{prop:fracdeg}（$n^*$ 的整性刻画与 Pell 例外集 $\mathcal E$ 的显式两支递推）；
定理~\ref{thm:unattain}（例外阶数上的强正则图障碍；(b) 支为符号谱计算，(a) 支为\emph{纯整数}判别式检验，
覆盖 $N\le7.1\times10^8$ 的 $13$ 个例外阶数）；命题~\ref{prop:ladder6}（$N=6$ 整数化阶梯的六个闭式）；""")
# re-add the diagnostic remark after rmk:fracdegfix
rep(r"""这是本报告「所有『对一切 $N$』型陈述都必须附带可执行检验」这条纪律的又一次兑现。""",
r"""这是本报告「所有『对一切 $N$』型陈述都必须附带可执行检验」这条纪律的又一次兑现。

\medskip
\noindent\textbf{诊断性读法（这条为什么是本轮最重要的结构性结论）。}
它说明 $\Phi-\rho_0=0.027449$ 这个缺口\emph{不是}因为 Cauchy--Schwarz 太粗：
四条 CS 在冠军处已经同时取等（注~\ref{rmk:quadtight}），把松弛做到「双均匀层内最优」也只能到 $\Phi$。
缺口的真正来源是两条\emph{离散}约束——钝角度必须是整数、钝角图的谱必须是整数——
它们在连续松弛 REL 里完全看不见。这也解释了为什么 \STW\ 的常数偏偏是 $\Phi=L(6)$：
\STW\ 的证明同样没有用到度与谱的整性，它停在的正是行方差机制的连续天花板。
$N=4$ 之所以精确，是因为 $n^*(4)=1$ 与 $A_O=2K_2$ 同时满足了这两条离散约束（注~\ref{rmk:n4exact}）。""")
io.open(main,"w",encoding="utf-8").write(s); print("ok; lines:",len(s.split("\n")))
