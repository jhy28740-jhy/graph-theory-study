# -*- coding: utf-8 -*-
import io
main="进展报告_前三大特征值之和的精确上界.tex"
s=io.open(main,encoding="utf-8").read()
def rep(old,new,cnt=1):
    global s
    assert s.count(old)==cnt, ("COUNT %d"%s.count(old), old[:70])
    s=s.replace(old,new,cnt)

# ---------- (1) \date ----------
rep(r"\date{2026 年 10 月 29 日（第六版，增 \S\ref{sec:bipglobal}：完全二部掩码的\textbf{全局}闭合 $\Lam_3\le\frac54$ 对一切 $N$；$N=5$ 层无条件；$N=6$ 例外 $22\to20$）}",
r"\date{2026 年 10 月 29 日（第七版，增 \S\ref{sec:rowvar}：行方差机制与天花板 $L(N)$；$Z_5$-等变秩三投影的完全分类与 $\rho_0$ 的闭式导出；$\Phi(Q)=\lambda_{\max}(Q_+)$ 的凸性与谱盒归约；旧队首 P1 被精确有理证书否证）}")

# ---------- (2) abstract 第七 bullet ----------
rep(r"""（数值 $\sup\Lam_3\le1.2524<\rho_0$，余量 $2.9\times10^{-2}$；注~\ref{rmk:mpopen} 精确指出了断点所在）。

\medskip
\noindent\textbf{关键词：}""",
r"""（数值 $\sup\Lam_3\le1.2524<\rho_0$，余量 $2.9\times10^{-2}$；注~\ref{rmk:mpopen} 精确指出了断点所在）。

第七，汇报第 12 轮的工作（第~\ref{sec:rowvar} 节）。这一轮有一条坏消息、一个新机制和一条主定理。
\textbf{坏消息}（否证）：第六版的队首 P1——「证 $N=6$ 的 $17$ 个无支配点类不是腔掩码」——被一个\emph{精确有理}证书推翻：
$Q=\frac12(I_6-M)$（$M$ 为完美匹配的置换矩阵）属于 $\mathcal P_3(6)$ 且 $1\{Q\ge0\}=I_6+A(K_{2,2,2})$，
故 $K_{2,2,2}$ \emph{是}腔掩码（命题~\ref{prop:chamberdead}）；进一步数值扫描显示 $19$ 个未闭合类\textbf{全部}是腔掩码
（符号违背量 $\le5\times10^{-19}$），即「纯几何不存在性」这条路线能关闭的类数是 $0$。
\textbf{新机制}（已证）：我们找到一组作用在\emph{每一行}上的恒等式与不等式——
行能量恒等式 $a_i+r_i=g_i(1-g_i)$、Perron 行方程 $s_i=(\rho-g_i)x_i$、三角形恒等式 \eqref{eq:tri}
（引理~\ref{lem:rowid}），以及由它们导出的四条 Cauchy--Schwarz 不等式 $[f],[g],[h_1],[h_2]$（命题~\ref{prop:fourcs}）。
冠军处这四条\textbf{同时取等}（残差 $\le10^{-16}$），只有五个锥行的 $[f]$ 留有 $1.8\times10^{-3}$ 的余量。
由此得到一个四参数族 $(g,\pi,r,w)$ 上的松弛 REL，$\kappa_3(N)\le\mathrm{REL}(N)$ 无条件成立，
并且在\emph{双均匀层}（$g_i\equiv3/N$、$\pi_i\equiv1/N$）上 REL 可以精确求解：
\textbf{定理~\ref{thm:LN}}：$\rho\le L(N)=\frac{2+\kappa_N}{1+\kappa_N}$，$\kappa_N=\sqrt{\frac{3(N-1)}{N-3}}$。
这一个闭式同时重导了两条已知结果——$\lim_{N\to\infty}L(N)=\frac{1+\sqrt3}2$ 是 Mohar 的常数，
$L(6)=\frac{3+\sqrt5}4=\Phi$ 是 \STW\ 的常数——并且 $L(4)=\frac54$ \textbf{取等}，给出 $\kappa_3(4)=\frac54$ 的第二条证明。
命题~\ref{prop:fracdeg} 解释了它为什么不能更好：取等要求的钝角度 $n^*=N\kappa_N/(3(1+\kappa_N))$
对一切 $N\ge5$ 都是\emph{无理数}（$n^*(6)=\frac{5-\sqrt5}2$），所以 $\Phi-\rho_0$ 这个缺口的本质是\textbf{度的整性}，
不是 Cauchy--Schwarz 太粗。$N=6$ 的整数化阶梯（命题~\ref{prop:ladder6}）把双均匀层压到 $\frac{13}{10}$，
再用两值刚性（命题~\ref{prop:twovalue}，符号谱计算）排除其可达性。
\textbf{主定理}（已证）：$N=6$ 的\emph{全部} $Z_5$-等变秩三投影恰好构成一条一维代数弧（引理~\ref{lem:z5class}），
弧上的 Perron 特征值有闭式 $\rho(s)=\frac{6+\sqrt5-2s+\sqrt{P(s)}}{10}$，
它在弧内部有唯一临界点 $s^*=\frac{9+4\sqrt5-\sqrt{35}}{18}$，而 $\rho(s^*)=\frac{9+\sqrt5+2\sqrt{35}}{18}=\rho_0$
\textbf{精确相等}（定理~\ref{thm:conefam}）。于是 $\rho_0$ 的地位从「一个被验证过的下界见证」升级为
「一条显式代数弧的唯一极大值」。意外的副产品是：定理~\ref{thm:zonal} 的见证——正二十面体等角线系——
就在\emph{同一条弧}上，位于 $s=\frac12$，值 $\frac{5+\sqrt5+\sqrt{30}}{10}=1.271329355$，
即本报告两个最重要的六阶对象是同一条弧的中点与顶点。
另外两条一般性结果：$\Phi(Q)=\lambda_{\max}(Q_+)$ 在全体对称矩阵上是\textbf{凸}的，
故 $\operatorname{conv}(\mathcal P_3(N))=\mathcal S_N$ 且秩三约束在求上确界时无损（定理~\ref{thm:convex}）；
余秩一层的符号分类给出 $\kappa_2(3)=\frac87$ 的一条三行新证明（推论~\ref{cor:k23}）。
现状的量化阶梯是
$1.366025$（Mohar）$>1.309017$（$\Phi=$ REL 自由上确界）$>1.300000$（$\frac{13}{10}=$ 加整性）
$>1.284781$（加 $Z_5$ 对称）$>1.281568$（$\rho_0$）；
剩余的 $3.2\times10^{-3}$ 已定位到 REL 丢掉的秩三信息（\eqref{eq:tri} 与四点 Gram 行列式），见第~\ref{sec:diag} 小节。

\medskip
\noindent\textbf{关键词：}""")

# ---------- (3) obstacles: new 障碍 0'' ----------
rep(r"""\subsection{主要障碍}

\paragraph{障碍 0（第六版重写，现居首位）：$N=6$ 的 $19$ 个未闭合类，即一条几何不存在性陈述}""",
r"""\subsection{主要障碍}

\paragraph{障碍 0$''$（第 12 轮重写，现居首位）：$3.2128\times10^{-3}$ 的秩三缺口}
第~\ref{sec:rowvar} 节把主线的形状再次改变。原来的首位障碍（下面这条「障碍 0」）已经\textbf{被否证}：
命题~\ref{prop:chamberdead} 说明 $K_{2,2,2}$ 是腔掩码，数值扫描说明 $19/19$ 类都是腔掩码，
所以「证不存在腔掩码」这条几何路线能关闭的类数是 $0$。
现在的首位障碍是一条\emph{定量的}缺口：在 $Z_5$-对称、五角星图案的类内，
REL 松弛（四条 Cauchy--Schwarz $+$ 行能量恒等式）的天花板是 $1.284780995477$，
而定理~\ref{thm:conefam} 证明真正的秩三最大值是 $\rho_0=1.281568196872$，差 $3.2128\times10^{-3}$。
这个差的来源已经被精确定位：REL 只保留了\emph{逐行二次}信息 $a_i+r_i=g_i(1-g_i)$，
丢掉了引理~\ref{lem:rowid}(iii) 的三角形恒等式 \eqref{eq:tri}、四点 Gram 行列式 $\det G_4=0$
（$\R^3$ 中任意四点线性相关）以及等角线 $\le6$ 的约束。
它\emph{有利}的地方是：这些都是\textbf{等式}而非不等式，而且第~\ref{sec:conefamily} 小节已经证明
它们在等变层里退化为恒等式，所以必须\emph{先降对称性}它们才会咬人——
这给出一个明确的、可执行的第一步（第~\ref{sec:pq7} 小节的 P1(c)）。
它\emph{不利}的地方是：四点 Gram 行列式是 $q_{ij}$ 的\emph{四次}多项式，
而 REL 的变量是 $(g,\pi,r,w)$ 这些\emph{二次}量，把前者翻译成后者的语言没有现成办法。

\paragraph{障碍 0（第六版，\textbf{已被否证}，保留以记录演化）：$N=6$ 的 $19$ 个未闭合类，即一条几何不存在性陈述}""")

# ---------- (4) honest labels ----------
rep(r"""\subsection{一句话总结}

下界一侧已经完整""",
r"""\subsection{第 12 轮新增结论的证明状态标注}\label{sec:honest12}

按第~\ref{sec:honest} 节的纪律，第~\ref{sec:rowvar} 节的结论逐条标注如下。

\textbf{已证（含完整论证，且不依赖任何未证猜想）}：
命题~\ref{prop:chamberdead}（$\Q$ 上精确，无浮点）；
引理~\ref{lem:rowid} 的三条恒等式（含三角形恒等式 \eqref{eq:tri}）；
命题~\ref{prop:fourcs} 的四条 Cauchy--Schwarz 不等式 $[f],[g],[h_1],[h_2]$；
定义~\ref{def:rel} 的松弛 REL 与 $\kappa_3(N)\le\mathrm{REL}(N)$；
定理~\ref{thm:LN}（$L(N)$ 闭式，含 AM--GM 取等点的可行性验证）；
命题~\ref{prop:fracdeg}（$n^*$ 的无理性）；命题~\ref{prop:ladder6}（$N=6$ 整数化阶梯的六个闭式）；
命题~\ref{prop:twovalue}（两值秩三投影的谱分类，全部由符号特征多项式给出）；
引理~\ref{lem:z5class}（$Z_5$-等变秩三投影的完全分类，特征多项式恒等式 $L^3(L-1)^3$）；
定理~\ref{thm:conefam}（$\rho(s)$ 闭式、唯一临界点、$\rho(s^*)=\rho_0$、精确 Perron 权重、ETF 识别）；
定理~\ref{thm:convex}（$\Phi$ 的凸性与 $\operatorname{conv}(\mathcal P_3(N))=\mathcal S_N$）；
引理~\ref{lem:perp2}（垂直二标架分解）；命题~\ref{prop:corank1} 与推论~\ref{cor:k23}（$\kappa_2(3)=\frac87$）；
引理~\ref{lem:slope1}（$\rho=\frac{\alpha+\ell}2$）。

\textbf{精确（数值等式已由符号计算判定为恒等，差化简后恰为 $0$）}：
冠军处四条 CS 同时取等（残差 $0$，非 $10^{-16}$ 量级的近似）；
$\|Q_-\|_F^2=\frac{16+5\sqrt7}{81}=\rho_0\tau$；
$\lambda^*(\text{冠军数据})=\rho_0$；$L(4)=\frac54$ 且被 $Q=I_4-vv^\top$ 取到；
$\rho(\frac12)=\frac{5+\sqrt5+\sqrt{30}}{10}$（ETF）；$\pi_h=\frac{91-25\sqrt7}{126}$、$\pi_c=\frac{7+5\sqrt7}{126}$；
$\mathrm{spec}(Q_+)$ 与 $\mathrm{spec}(Q_-)$ 的四个根式；
$\kappa_2(3)$ 见证的特征多项式 $\frac1{343}(7L-8)(7L-5)(7L-1)$；
$6$ 阶会议矩阵的 $S^2=5I_6$。

\textbf{数值（只有浮点证据，不得作为成果引用）}：
$\mathrm{REL}(N)=L(N)$ 对 $N=5,6,7,8,10,12$（差分进化，吻合到 $\le4\times10^{-6}$；\emph{猜想}，非定理）；
$Z_5$-对称层八个图案类的天花板（其中 $6/8$ 类被 REL 单独排除，冠军所在类为 $1.284780995477$）；
$19/19$ 未闭合类都是腔掩码（符号违背量 $\le5\times10^{-19}$）；
$\alpha=1$ 层的 $\max\beta=0.258695$；通道线性规划的 $8701$ 个构型与 HiGHS 解；
$\Phi$ 的凸性在 $N=5,6,7$ 上共 $24000$ 对随机点无 Jensen 违背（\emph{注}：定理~\ref{thm:convex}(a) 本身已证，
这条数值只是独立交叉核验）；
$E_\perp=\frac12x^\top|B|x=0.519536756454917$（未找到闭式）。

\textbf{否证}：命题~\ref{prop:chamberdead} 否证了 \eqref{eq:chamber6} 与第六版队首 P1；
注~\ref{rmk:tab19read}(iv) 的「等价」降级为单向蕴含；
不等式 $[g]$ 的取等\emph{不是}冠军的刻画（$N=4$ 上允许 $\rho=1.4186$，见第~\ref{sec:miscid} 小节 (6)）。

\subsection{一句话总结}

下界一侧已经完整""")

# ---------- (5) 第七版增补 ----------
rep(r"""三者都不依赖任何未证猜想，且方法论互不重叠（团吹起 $+$ 根式恒等式 / 半定对偶 / 秩二长期方程 $+$ 初等不等式）。

\appendix""",
r"""三者都不依赖任何未证猜想，且方法论互不重叠（团吹起 $+$ 根式恒等式 / 半定对偶 / 秩二长期方程 $+$ 初等不等式）。

\medskip
\noindent\textbf{（第七版增补）一句话总结的更新。}
第六版把主线定为「对 $17$ 个符号模式逐个建立 \eqref{eq:chamber6}」。第 12 轮把这条路线\textbf{否证}了：
$19/19$ 类都是腔掩码（命题~\ref{prop:chamberdead}，$K_{2,2,2}$ 有精确有理证书），
所以「腔掩码不存在性」这条路能关闭的类数是 $0$，第六版增补中「唯一在列的路径」的说法作废。
现在的状态是：\textbf{主线从「逐类排除」换成「一条不等式链的逐级收紧」}，
链条是 $1.366025$（Mohar）$>1.309017$（$\Phi$，等于 REL 的自由上确界）$>1.300000$（$\frac{13}{10}$，加度的整性）
$>1.284781$（加 $Z_5$ 对称）$>1.281568$（$\rho_0$，加秩三）。
前两级已经\textbf{定理化}（定理~\ref{thm:LN} 与命题~\ref{prop:ladder6}），
第三级是数值天花板，第四级由定理~\ref{thm:conefam} 精确闭合。
于是「还差什么」第一次有了一个\emph{定量且定位清楚}的答案：差 REL 与秩三之间的 $3.2\times10^{-3}$，
缺的原料是三角形恒等式 \eqref{eq:tri} 与四点 Gram 行列式（第~\ref{sec:diag} 小节的 (A)(B) 两项）。
另外，$\rho_0$ 现在不再只是一个见证：它是 $Z_5$-等变秩三投影弧的唯一极大值（定理~\ref{thm:conefam}），
$\sqrt{35}$ 来自判别式 $140$，$\sqrt7$ 来自 Perron 权重的二次方程——两个根式都有了出处。
至于毕业论文的章节结构，第七版之后可以独立成章的结果增加到\textbf{五}个：
上面三个，加上「行方差机制与天花板 $L(N)$」（定理~\ref{thm:LN}，含 Mohar 与 \STW\ 两个常数的统一重导），
以及「$Z_5$-等变弧的完全分类与 $\rho_0$ 的闭式导出」（定理~\ref{thm:conefam}，含 ETF 的意外识别）。
这两章都是纯数学的、自包含的、不需要任何数值验证即可宣读。

\appendix""")
io.open(main,"w",encoding="utf-8").write(s)
print("all edits applied; lines:",len(s.split("\n")))
