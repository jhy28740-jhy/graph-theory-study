import io
tex = r"""
\subsection{冠军的闭式导出：$Z_5$-等变秩三投影恰好是一条代数弧}\label{sec:conefamily}

到目前为止，$\rho_0=\frac{9+\sqrt5+2\sqrt{35}}{18}$ 在报告里的地位是「一个被构造出来、并用高精度数值验证过的下界见证」：
第~\ref{sec:champ}\ 节给出了它的构造，证书脚本核验了它的每一位数字，但\textbf{没有人解释它为什么是这几个根式}。
本小节把这件事变成一个定理。结论是：

\begin{quote}\small
在 $N=6$ 的\emph{全部} $Z_5$-等变秩三投影上（这是一个一维代数弧，不是离散集合），
$\lambda_{\max}(Q_+)$ 是一个显式的一元根式函数 $\rho(s)$，它在弧内部有唯一临界点 $s^*$，
而 $\rho(s^*)$ \textbf{精确等于} $\rho_0$。冠军不是数值搜索的偶然产物，而是一个二次方程的根。
\end{quote}

更意外的是副产品：定理~\ref{thm:zonal} 的见证——正二十面体等角线系（ETF）——
\textbf{就在同一条弧上}，位于 $s=\frac12$；它的值是 $\rho(\frac12)=\frac{5+\sqrt5+\sqrt{30}}{10}=1.271329355255$，
正是第~\ref{sec:kappa} 节数值扫描里那个始终排在冠军之后的「亚军」。
也就是说，本报告里两个最重要的六阶对象（$\Phi$ 的见证与 $\rho_0$ 的见证）是同一条代数弧上的两点，
而 $\rho_0$ 是这条弧的最高点、ETF 是弧的中点。这个统一是我们之前完全没有预料到的。

\subsubsection*{等变分解与分类}

把六个点记作 $\{h\}\cup\{1,2,3,4,5\}$，$Z_5$ 平凡地作用在 $h$ 上、以 $j\mapsto j+1$ 循环作用在锥点上。
一个 $Z_5$-等变的对称矩阵由五个实参数决定：
\[
g_h:=q_{hh},\quad g_c:=q_{jj},\quad A:=q_{hj},\quad d_1:=q_{jk}\ (|j-k|\equiv1),\quad
d_2:=q_{jk}\ (|j-k|\equiv2).
\]
$\R^6$ 分解为三个不可约子空间：平凡分量 $T=\mathrm{span}\{e_h,\ \one_c/\sqrt5\}$（维 $2$），
以及两个二维分量 $\psi_1$（特征标 $\chi_1,\chi_4$）、$\psi_2$（特征标 $\chi_2,\chi_3$）。
等变矩阵在这三个分量上分别作用为
\begin{equation}\label{eq:isoT}
Q\big|_T=\begin{pmatrix} g_h & \sqrt5\,A\\[2pt] \sqrt5\,A & g_c+2d_1+2d_2\end{pmatrix},\qquad
\mu_1=g_c+\tfrac{\sqrt5-1}2d_1-\tfrac{\sqrt5+1}2d_2,\qquad
\mu_2=g_c-\tfrac{\sqrt5+1}2d_1+\tfrac{\sqrt5-1}2d_2 .
\end{equation}

\begin{lemma}[$Z_5$-等变秩三投影的完全分类，已证]\label{lem:z5class}
$Q\in\mathcal P_3(6)$ 是 $Z_5$-等变的，当且仅当存在 $s\in[0,1]$ 与 $\varepsilon\in\{+1,-1\}$，使得
\begin{equation}\label{eq:arc}
g_h=1-s,\qquad g_c=\frac{s+2}5,\qquad A=\varepsilon\sqrt{\frac{s(1-s)}5},\qquad
\{d_1,d_2\}=\Bigl\{\frac{2s-1+\sqrt5}{10},\ \frac{2s-1-\sqrt5}{10}\Bigr\},
\end{equation}
其中 $\{d_1,d_2\}$ 是\emph{无序}对：$d_1>d_2$ 与 $d_1<d_2$ 两种分配都出现，各占一支。
每个这样的 $Q$ 的特征多项式恒等于 $L^3(L-1)^3$（对 $s$ 是恒等式，不是逐点数值事实）。
\end{lemma}

\begin{proof}
$\mathrm{spec}(Q)=\mathrm{spec}(Q|_T)\cup\{\mu_1,\mu_1\}\cup\{\mu_2,\mu_2\}$，而投影的特征值只能是 $0$ 或 $1$，
故 $\mu_1,\mu_2\in\{0,1\}$，且 $Q|_T$ 的特征值也在 $\{0,1\}$ 中。秩条件
$3=\rank(Q|_T)+2\mu_1+2\mu_2$ 只有唯一解：$\rank(Q|_T)=1$，且 $\mu_1,\mu_2$ 恰有一个为 $1$、另一个为 $0$
（$\rank(Q|_T)\le2$ 排除 $0$ 或 $2$ 之外的情形：$0+2+2=4$、$2+0+0=2$、$2+2+0=4$ 都不等于 $3$）。

\textbf{支 I}（$\mu_2=0,\mu_1=1$）。两式相减给 $\mu_1-\mu_2=\sqrt5(d_1-d_2)=1$，即 $d_1-d_2=1/\sqrt5$；
再把 $d_2=d_1-1/\sqrt5$ 代入 $\mu_1=1$：
\[
g_c+d_1\Bigl(\tfrac{\sqrt5-1}2-\tfrac{\sqrt5+1}2\Bigr)+\tfrac{\sqrt5+1}{2\sqrt5}=1
\ \Longrightarrow\ d_1=g_c-\tfrac{5-\sqrt5}{10},\quad d_2=g_c-\tfrac{5+\sqrt5}{10}.
\]
$\rank(Q|_T)=1$ 意味着 $\tr(Q|_T)=1$ 与 $\det(Q|_T)=0$。前者给 $g_h+5g_c=3$（这正是 $\tr Q=3$），
后者给 $5A^2=g_h(g_c+2d_1+2d_2)=(3-5g_c)(5g_c-2)$。令 $s:=5g_c-2$，
则 $g_c=\frac{s+2}5$、$g_h=1-s$、$5A^2=s(1-s)$，代回 $d_1,d_2$ 得 \eqref{eq:arc} 中 $d_1>d_2$ 的那一支；
$A^2\ge0$ 恰好等价于 $s\in[0,1]$。

\textbf{支 II}（$\mu_1=0,\mu_2=1$）。完全平行地得到 $d_2-d_1=1/\sqrt5$，
即 $d_1=g_c-\frac{5+\sqrt5}{10}$、$d_2=g_c-\frac{5-\sqrt5}{10}$，而 $g_h,g_c,A$ 与支 I 逐字相同。
这就是 \eqref{eq:arc} 中 $d_1<d_2$ 的那一支。

最后，特征多项式 $L^3(L-1)^3$ 的恒等式由符号计算直接验证（两支各自展开后所有 $s$ 的幂次相消），
证书脚本 \texttt{证书\_行方差松弛与天花板.py} 的第 (iv) 项给出 \texttt{sympy.factor} 的输出。$\square$
\end{proof}

\begin{remark}[两支在正部意义下是同一条弧]\label{rmk:twobranch}
支 I 满足 $d_1>0>d_2$（对一切 $s\in[0,1]$：$2s-1-\sqrt5\le2-1-\sqrt5<0$），
支 II 满足 $d_2>0>d_1$。于是 $Q_+$ 的做法在两支上都是「保留那个正的锥--锥值、把负的那个抹成 $0$」。
锥点上的置换 $j\mapsto 2j\pmod 5$（固定 $h$）把「距离 $1$」对换成「距离 $2$」对，
因此它是一个置换矩阵 $P$，使 $Q_+^{\mathrm{II}}(s)=P\,Q_+^{\mathrm{I}}(s)\,P^\top$。
\textbf{两支的 $\lambda_{\max}(Q_+)$ 逐点相等}，所以「$\lambda_{\max}(Q_+)$ 的等变弧」实际上只有一条，
参数 $s\in[0,1]$，$\varepsilon$ 与支的选取都不影响它。这条观察把后面所有计算减半。
\end{remark}

\subsubsection*{弧上的 Perron 特征值}

以下固定支 I（由注~\ref{rmk:twobranch} 无损失）。抹去 $d_2$ 之后，$Q_+$ 仍是 $Z_5$-等变的，
它在 $\psi_1,\psi_2$ 上的标量是
\begin{equation}\label{eq:psip}
\mu_1(Q_+)=g_c+\tfrac{\sqrt5-1}2d_1=\frac{7-\sqrt5+(1+\sqrt5)s}{10},\qquad
\mu_2(Q_+)=g_c-\tfrac{\sqrt5+1}2d_1=\frac{2-(\sqrt5-1)s}{10},
\end{equation}
而在平凡分量上是
\[
T_+(s)=\begin{pmatrix}1-s & \sqrt{s(1-s)}\\[2pt] \sqrt{s(1-s)} & \dfrac{3s+1+\sqrt5}{5}\end{pmatrix}
\qquad(\text{右上元}=\sqrt5A=\sqrt{s(1-s)}).
\]
由 \eqref{eq:psip}，$\mu_1(Q_+)\le\frac{7-\sqrt5+1+\sqrt5}{10}=\frac8{10}<1$ 且
\[
T_+(s)_{22}-\mu_1(Q_+)=\frac{(5-\sqrt5)s+3\sqrt5-5}{10}\ \ge\ \frac{3\sqrt5-5}{10}=0.170820>0 ,
\]
故 $\lambda_{\max}(Q_+(s))=\lambda_{\max}(T_+(s))$ 对\emph{一切} $s\in[0,1]$ 成立（不需要任何数值验证）。

\begin{theorem}[弧上的精确 Perron 谱与 $\rho_0$ 的闭式导出，已证]\label{thm:conefam}
定义 $\rho(s):=\lambda_{\max}\bigl(Q_+(s)\bigr)$，$s\in[0,1]$。则
\begin{equation}\label{eq:rhos}
\rho(s)=\frac{6+\sqrt5-2s+\sqrt{P(s)}}{10},\qquad
P(s):=21-8\sqrt5+(36+16\sqrt5)s-36s^2,
\end{equation}
并且：
\begin{enumerate}[label=(\alph*),leftmargin=2.2em]
\item $\rho'(s)=0$ 在 $[0,1]$ 内\textbf{恰有一个}解
$s^*=\dfrac{9+4\sqrt5-\sqrt{35}}{18}=0.668232895938863$；
另一个根 $\frac{9+4\sqrt5+\sqrt{35}}{18}=1.325575094$ 落在弧外。
\item $\rho(s^*)=\dfrac{9+\sqrt5+2\sqrt{35}}{18}=\rho_0$，\textbf{精确相等，无任何数值误差}。
\item 端点值 $\rho(0)=1$、$\rho(1)=\frac{4+\sqrt5}5=1.247213595$，都严格小于 $\rho_0$。
\item 故 $\displaystyle\max_{s\in[0,1]}\lambda_{\max}(Q_+(s))=\rho_0$，且最大化点唯一。
\item Perron 向量在弧的最高点处是 $x=(x_h,x_c,x_c,x_c,x_c,x_c)$，其平方权重\textbf{精确}为
\begin{equation}\label{eq:perronw}
\pi_h=x_h^2=\frac{91-25\sqrt7}{126}=0.197271565264962,\qquad
\pi_c=x_c^2=\frac{7+5\sqrt7}{126}=0.160545686947008 .
\end{equation}
\item $s=\frac12$ 处（支 II）给出 $Q=\frac12I_6+\frac1{2\sqrt5}S$，其中 $S$ 是 $6$ 阶会议矩阵
（$S^2=5I_6$，对角为 $0$、非对角为 $\pm1$），即\textbf{正二十面体等角线系}；
其值为 $\rho(\frac12)=\frac{5+\sqrt5+\sqrt{30}}{10}=1.271329355255145$。
\end{enumerate}
\end{theorem}

\begin{proof}
(a) 由 $2\times2$ 公式，$\rho(s)=\frac12\bigl(1-s+\frac{3s+1+\sqrt5}5\bigr)+\frac12\sqrt{D(s)}$，
$D(s)=\bigl(1-s-\frac{3s+1+\sqrt5}5\bigr)^2+4s(1-s)$。通分化简得
$D(s)=P(s)/25$（$25D(s)$ 展开后与 $P(s)$ 逐项相同），于是 \eqref{eq:rhos} 成立。
求导：$\rho'(s)=-\frac15+\frac{P'(s)}{20\sqrt{P(s)}}$，$P'(s)=36+16\sqrt5-72s$，故 $\rho'=0$ 等价于
$P'(s)=4\sqrt{P(s)}$；两边平方（在 $[0,1]$ 上 $P'(s)\ge36+16\sqrt5-72=35.78>0$，故平方不引入伪根）得
$P'(s)^2=16P(s)$，展开相消后是
\[
18s^2-(18+8\sqrt5)s+(7+4\sqrt5)=0,
\]
判别式 $(18+8\sqrt5)^2-72(7+4\sqrt5)=140$，故 $s=\frac{18+8\sqrt5\pm2\sqrt{35}}{36}=\frac{9+4\sqrt5\pm\sqrt{35}}{18}$，
其中只有取 $-$ 的那个落在 $[0,1]$。

(b) 在临界点上 $\sqrt{P(s^*)}=P'(s^*)/4=9+4\sqrt5-18s^*$，代回 \eqref{eq:rhos}：
\[
\rho(s^*)=\frac{6+\sqrt5-2s^*+9+4\sqrt5-18s^*}{10}=\frac{15+5\sqrt5}{10}-2s^*
=\frac{3+\sqrt5}2-\frac{9+4\sqrt5-\sqrt{35}}9
=\frac{9(3+\sqrt5)-2(9+4\sqrt5-\sqrt{35})}{18}=\frac{9+\sqrt5+2\sqrt{35}}{18}.
\]

(c) $P(0)=21-8\sqrt5=(4-\sqrt5)^2$、$P(1)=21+8\sqrt5=(4+\sqrt5)^2$，
故 $\rho(0)=\frac{6+\sqrt5+4-\sqrt5}{10}=1$、$\rho(1)=\frac{4+\sqrt5}{10}+\frac{4+\sqrt5}{10}=\frac{4+\sqrt5}5$。
两者都小于 $\rho_0=1.281568$。

(d) 由 (a)(c) 与 $\rho$ 在 $[0,1]$ 上连续、在内部只有一个临界点。

(e) $T_+(s^*)$ 的 Perron 向量满足 $(g_h-\rho)x_h+\sqrt5A\,x_c'=0$（这里 $x_c'$ 是 $\one_c/\sqrt5$ 方向的分量），
即 $x_h/x_c=\sqrt5\cdot\frac{\sqrt5A}{\rho-g_h}=\frac{5A}{\rho(s^*)-g_h}$。
用 $\rho(s^*)=\frac{3+\sqrt5}2-2s^*$、$g_h=1-s^*$、$5A=\sqrt{5s^*(1-s^*)}$ 化简得
$x_h/x_c=1.1084925$，再由 $x_h^2+5x_c^2=1$ 得 $x_c^2=\frac1{(x_h/x_c)^2+5}$；
符号计算给出 \eqref{eq:perronw} 的两个根式（证书脚本第 (v) 项逐位核验
$\frac{91-25\sqrt7}{126}$、$\frac{7+5\sqrt7}{126}$ 与数值 Perron 向量之差为 $0$）。

(f) $s=\frac12$ 代入 \eqref{eq:arc}：$g_h=g_c=\frac12$、$|A|=\frac1{2\sqrt5}$、
$\{|d_1|,|d_2|\}=\{\frac1{2\sqrt5}\}$，即所有非对角元的绝对值都是 $\frac1{2\sqrt5}$，
$Q=\frac12I_6+\frac1{2\sqrt5}S$，$S:=2\sqrt5(Q-\frac12I_6)$ 是 $\pm1$ 对称矩阵、零对角；
符号计算给 $S^2=5I_6$，故 $S$ 是 $6$ 阶会议矩阵，$Q$ 是 $6$ 点在 $\R^3$ 中的等角紧标架（正二十面体的 $6$ 条对角轴）。
$\rho(\frac12)=\frac{6+\sqrt5-1+\sqrt{30}}{10}$，其中 $P(\frac12)=21-8\sqrt5+18+8\sqrt5-9=30$。$\square$
\end{proof}

\begin{remark}[这条定理在主线上意味着什么]\label{rmk:coneread}
(i) \textbf{$\rho_0$ 从「见证」升级为「极值」。}
第~\ref{sec:lower}\ 节的定理~\ref{thm:lower} 说 $\rho_0$ 可达；本小节说在整个 $Z_5$-等变层内它\emph{只有}这一个极大值点，
并且这个极大值是二次方程 $18s^2-(18+8\sqrt5)s+(7+4\sqrt5)=0$ 的根的显式函数。
根式 $\sqrt{35}$ 的来源至此清楚了：它是判别式 $140=4\cdot35$ 的平方根，
而 $140$ 来自「弧的参数化」与「$2\times2$ Perron 公式」两次代换的合成，没有任何神秘之处。

(ii) \textbf{它同时解释了 $\Phi$ 的见证与 $\rho_0$ 的见证为什么长得像。}
定理~\ref{thm:zonal} 说 $\Phi$ 是「逐点主函数 $+$ 逐通道线性计费」这一\emph{方案族}的天花板，见证是 ETF；
本小节说 ETF 是等变弧的中点，而弧的最高点比它高 $\rho_0-\frac{5+\sqrt5+\sqrt{30}}{10}=9.9\times10^{-3}$。
换句话说，\STW\ 的 $\Phi$ 与我们的 $\rho_0$ 之间的差距，在等变层内\emph{不是}两个孤立常数的差距，
而是同一条弧上「中点」与「顶点」的差距。

(iii) \textbf{它与 REL 松弛的对照是量化的。}
第~\ref{sec:symlayer} 小节的 $Z_5$-对称数值扫描给出五角星图案类的 REL 天花板 $1.284780995477$，
比 $\rho_0$ 高 $3.2128\times10^{-3}$。本小节证明真正的等变秩三最大值是 $\rho_0$。
所以\textbf{那 $3.2\times10^{-3}$ 完全来自 REL 丢掉的秩三信息}，而不是来自 Cauchy--Schwarz 的粗糙：
REL 保留的是「逐行能量 $a_i+r_i=g_i(1-g_i)$」这一\emph{二次}信息，
丢掉的是引理~\ref{lem:rowid}(iii) 的三角形恒等式 \eqref{eq:tri} 与四点 Gram 行列式。
这正是第~\ref{sec:diag} 小节排在 (B) 位的缺失成分，也是新队首 P1 的靶心。

(iv) \textbf{它给出一条可执行的改进路线。}
弧是一条\emph{一维}对象，秩三约束在它上面是恒等式；要压缩那 $3.2\times10^{-3}$，
必须允许 $Q$ 离开等变层（否则已经做完了）。
最小的下一步是把对称性从 $Z_5$ 降到「hub $+$ 锥点上的一个二面角子群」，
此时等变参数从 $1$ 维变成 $3$--$4$ 维，秩三约束（\eqref{eq:tri} 与 Gram 行列式）第一次成为\emph{非平凡}的方程，
可以期望把 REL 的天花板往下压。这是 P1 的具体形态。
\end{remark}
"""
io.open("_frag_rowvar_4.tex","w",encoding="utf-8").write(tex.lstrip("\n"))
print("written A, lines:",len(tex.split(chr(10))))
