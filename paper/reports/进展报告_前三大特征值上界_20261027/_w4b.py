import io
tex = r"""
\subsection{垂直二标架分解：把秩三拆成秩一 $+$ 秩二}\label{sec:perp2}

冠军的 $x$ 不是 $Q$ 的特征向量（$\alpha=0.9766\ne1$），所以 $Qx$ 与 $x$ 之间有一个夹角。
把这件事做成一个分解，就能把 $\mathcal P_3$ 上的问题改写成 $\mathcal P_2$ 上的问题——后者的几何要简单得多
（$\R^3$ 中的 $2$-标架就是一个单位向量加一个正交方向）。

\begin{lemma}[垂直二标架分解，已证]\label{lem:perp2}
设 $Q\in\mathcal P_3(N)$、$x\ge0$、$\|x\|=1$、$\alpha=x^\top Qx>0$。令
\[
p:=\frac{Qx}{\sqrt\alpha},\qquad B:=Q-pp^\top,\qquad q:=p-\sqrt\alpha\,x .
\]
则：
\begin{enumerate}[label=(\roman*),leftmargin=2.2em]
\item $\|p\|=1$，且 $B\in\mathcal P_2(N)$（$B^2=B=B^\top$、$\tr B=2$、$B\succeq0$）；
\item $Bp=0$、$Bx=0$；
\item $\|q\|^2=1-\alpha$，且 $q\perp x$；
\item $Q=pp^\top+B$，从而 $\Lam_3$ 的 Ky Fan 对偶量分裂为「秩一部分 $pp^\top$」与「垂直部分 $B$」之和。
\end{enumerate}
\end{lemma}

\begin{proof}
$Q$ 是投影，故 $Qp=Q^2x/\sqrt\alpha=Qx/\sqrt\alpha=p$，$\|p\|^2=x^\top Q^2x/\alpha=\alpha/\alpha=1$。
于是 $B^2=Q^2-Qpp^\top-pp^\top Q+pp^\top pp^\top=Q-pp^\top-pp^\top+pp^\top=B$，
$\tr B=3-1=2$，$B=Q-pp^\top\succeq0$（因为 $p$ 在 $Q$ 的像空间内且 $\|p\|=1$，$pp^\top\preceq Q$）。
$Bp=Qp-p=0$；$Bx=Qx-p(p^\top x)=\sqrt\alpha p-\sqrt\alpha p=0$（$p^\top x=x^\top Qx/\sqrt\alpha=\sqrt\alpha$）。
$\|q\|^2=1-2\sqrt\alpha\,p^\top x+\alpha=1-2\alpha+\alpha=1-\alpha$，$q^\top x=p^\top x-\sqrt\alpha=0$。$\square$
\end{proof}

\begin{remark}[冠军处的分裂，精确/数值]\label{rmk:perp2read}
冠军处 $\|q\|^2=1-\alpha=0.023404740944563$，
恰好等于命题~\ref{prop:decouple} 里「投影能量损失」那一项——这不是巧合：
$1-\alpha=\|q\|^2$ 就是 $x$ 偏离 $Q$ 的像空间的平方距离。
垂直部分的能量 $E_\perp:=\tfrac12x^\top|B|x$ 满足 $x^\top B_+x=-x^\top B_-x=E_\perp$
（因为 $Bx=0$ 使 $x^\top Bx=0$），冠军处 $E_\perp=0.519536756454917$（数值）。
这个分解的价值在于：$B\in\mathcal P_2(N)$ 且 $Bx=0$，
于是 $\kappa_3$ 可以改写成
\[
\kappa_3=\sup\Bigl\{x^\top\bigl(pp^\top+B\bigr)_+x:\ B\in\mathcal P_2(N),\ Bp=Bx=0,\ \|p\|=1,\ x\ge0,\ \|x\|=1\Bigr\},
\]
其中秩三条件被\emph{完全}吸收进了 $B\in\mathcal P_2$ 与两个正交约束里。
$\mathcal P_2(N)$ 是 $\R^3$ 中的 $2$-标架空间，它的秩二结构（例如四点 Gram 行列式为零）
比秩三情形的对应约束要\textbf{显式}得多。这是 P1 路线的第二种可能形态。
\end{remark}

\subsection{$\Phi(Q)=\lambda_{\max}(Q_+)$ 是凸的：谱盒归约定理}\label{sec:convex}

这是本轮一个短小但\textbf{方向性很强}的定理。它说明「秩三」这个离散约束在求上确界时可以\emph{无损}丢掉，
同时也说明为什么丢掉之后仍然做不出 $\rho_0$。

\begin{theorem}[凸性与谱盒归约，已证]\label{thm:convex}
记 $\mathcal S_N:=\{Q\in\R^{N\times N}_{\rm sym}:0\preceq Q\preceq I,\ \tr Q=3\}$（谱盒），
$\Phi(Q):=\lambda_{\max}(Q_+)$。则
\begin{enumerate}[label=(\alph*),leftmargin=2.2em]
\item $\Phi$ 在\emph{全体}对称矩阵上是凸函数；
\item $\operatorname{conv}\bigl(\mathcal P_3(N)\bigr)=\mathcal S_N$，且 $\mathcal S_N$ 的极点恰是 $\mathcal P_3(N)$；
\item 因此 $\kappa_3(N)=\max_{Q\in\mathcal P_3(N)}\Phi(Q)=\max_{Q\in\mathcal S_N}\Phi(Q)$。
\end{enumerate}
\end{theorem}

\begin{proof}
(a) 因为 $Q_+\ge0$ 逐点非负，对任意单位向量 $y$ 有 $y^\top Q_+y\le|y|^\top Q_+|y|$，故
\[
\Phi(Q)=\max_{\substack{y\ge0\\ \|y\|=1}}y^\top Q_+y=\max_{\substack{y\ge0\\ \|y\|=1}}\sum_{i,j}(q_{ij})_+\,y_iy_j .
\]
固定 $y\ge0$ 时，每一项 $(q_{ij})_+y_iy_j$ 是 $q_{ij}$ 的函数 $t\mapsto c\,(t)_+$ 与 $c=y_iy_j\ge0$ 的复合，
而 $c\,(t)_+$ 在 $c\ge0$ 时是凸函数（它是 $\max(ct,0)$）。有限个凸函数之和仍凸，
故内层函数对 $Q$ 凸；$\Phi$ 是一族凸函数的逐点上确界，因此凸。

(b) $\mathcal P_3(N)$ 是对角阵 $\mathrm{diag}(1,1,1,0,\dots,0)$ 的\emph{酉轨道}。
酉轨道的凸包由 Schur--Horn 定理刻画为 $\{Q:\mathrm{spec}(Q)\prec(1,1,1,0,\dots,0)\}$，
而这正是 $\mathcal S_N$（$0\preceq Q\preceq I$ 且 $\tr Q=3$ 等价于特征值落在 $[0,1]$ 且和为 $3$，
即被 $(1,1,1,0,\dots)$ 优超）。凸紧集的极点：$Q\in\mathcal S_N$ 是极点当且仅当它不能写成两个不同点的中点，
这在谱语言里等价于每个特征值都是 $0$ 或 $1$，即 $Q\in\mathcal P_3(N)$。

(c) 凸函数在凸紧集上的最大值在极点取得，故 $\max_{\mathcal S_N}\Phi=\max_{\mathcal P_3(N)}\Phi$。$\square$
\end{proof}

\begin{remark}[这条定理的两副面孔]\label{rmk:convexread}
(i) \textbf{好的一面}：$\kappa_3(N)$ 可以在谱盒 $\mathcal S_N$ 上计算，秩三约束是\emph{免费}的。
这为「用半定/矩方法给 $\kappa_3$ 配上界」提供了合法性：任何在 $\mathcal S_N$ 上成立的 $\Phi$ 的上界，
自动就是 $\kappa_3(N)$ 的上界，不需要额外论证极值点是投影。
第~\ref{sec:relax} 小节的 REL 松弛与第~\ref{sec:conefamily} 小节的等变弧都落在这个框架里。

(ii) \textbf{坏的一面}：$\Phi$ 是凸的，所以「在 $\mathcal S_N$ 上最大化 $\Phi$」是一个\emph{凸最大化}问题，
不是半定规划；任何把它变成可解问题的做法都必须用一个\emph{凹的或线性的}函数从上方控制 $\Phi$，
而这一步必然损失。两个已经量化过的例子：
\begin{itemize}
\item $\Phi(Q)\le\frac12\bigl(1+\lambda_{\max}(|Q|)\bigr)$（因为 $Q_+=\frac12(Q+|Q|)$）。
右端在 $\mathcal S_N$ 上可以取到 $\Phi=\frac{3+\sqrt5}4$——见证正是 ETF（$\lambda_{\max}(|Q|)=\varphi$）；
再用 $\lambda_{\max}(|Q|)\le\||Q|\|_F=\|Q\|_F=\sqrt3$ 就退化成 Mohar 的 $\frac{1+\sqrt3}2$。
\item 完全丢掉自洽性 $A=1\{Q>0\}$ 之后，命题~\ref{prop:octnogo} 给出八面体测试集上的值 $4/3>\rho_0$。
\end{itemize}
所以凸化路线的\textbf{最好成绩恰好是 \STW\ 的 $\Phi$}——这与定理~\ref{thm:zonal} 的结论完全一致，
而且是同一个现象的两个视角：$\Phi$ 是「线性/凹控制」的天花板，$\rho_0$ 需要非凸信息。
\end{remark}

\subsection{余秩一层：$\kappa_2(3)=\frac87$ 与 $\kappa_3(4)=\frac54$ 的自包含证明}\label{sec:corank1}

$N$ 阶的秩 $N-1$ 投影恰好是 $\{I_N-vv^\top:\|v\|=1\}$，即「余秩一层」。
定理~\ref{thm:convex} 在这一层把 $\kappa_{N-1}(N)$ 化成一个\textbf{球面上的}极值问题，
而球面问题的符号分类是可以穷尽的。这一小节给出两个精确常数。

\begin{proposition}[余秩一层的闭式，已证]\label{prop:corank1}
设 $v\in\R^N$、$\|v\|=1$、$Q=I_N-vv^\top$。
\begin{enumerate}[label=(\alph*),leftmargin=2.2em]
\item 若 $v$ 的坐标有 $m$ 个正、$1$ 个负、其余为零（$N=m+1$），记 $v=(c\one_m,-d)$、$mc^2+d^2=1$，则
\[
\lambda_{\max}(Q_+)=\tfrac12\bigl(1+(m-1)c^2\bigr)+\tfrac12\sqrt{(d^2-c^2)^2+4mc^2d^2},
\]
它在 $c^2=\frac2{3m+1}$、$d^2=\frac{m+1}{3m+1}$ 处取最大值
\begin{equation}\label{eq:fm1}
f_{m,1}=\frac{4m}{3m+1}.
\end{equation}
\item 若 $v$ 有 $m$ 个正、$m$ 个负且两组内部分别相等（$N=2m$，$v=(c\one_m,-c\one_m)$、$c^2=\frac1{2m}$），则
\begin{equation}\label{eq:fmm}
\lambda_{\max}(Q_+)=f_{m,m}=1+\frac{m-1}{2m}=\frac{3m-1}{2m}.
\end{equation}
\end{proposition}

\begin{proof}
(a) $Q_+$ 的非对角正元恰是「正坐标 $\times$ 负坐标」的 $m$ 个交叉位置，值都是 $cd$；对角是 $1-c^2$（$m$ 个）与 $1-d^2$。
由 $m$ 个正坐标的对称性，Perron 向量形如 $(y,y,\dots,y,z)$，于是问题化到二维子空间
$\mathrm{span}\{\one_m/\sqrt m,\ e_{m+1}\}$ 上的矩阵
$\begin{pmatrix}1-c^2&\sqrt m\,cd\\ \sqrt m\,cd&1-d^2\end{pmatrix}$，其特征值即所给公式。
令 $z=c^2$、$d^2=1-mz$，被开方式成为 $1+(2m-2)z+(-3m^2+2m+1)z^2$，
求导并整理得临界点方程 $z\bigl[4m(3m+1)z-8m\bigr]=0$，非零根 $z=\frac2{3m+1}$。
代入后开方式\emph{恰好}等于 $1$（因为 $(m-1)^2+8m(m+1)=(3m+1)^2$），故
$\lambda_{\max}=\frac12\cdot\frac{5m-1}{3m+1}+\frac12=\frac{4m}{3m+1}$。端点 $z=0$ 与 $z=1/m$ 都给 $1<\frac{4m}{3m+1}$。

(b) $Q_+$ 是「对角 $1-\frac1{2m}$ $+$ $\frac1{2m}K_{m,m}$ 的邻接」，
故 $\lambda_{\max}=1-\frac1{2m}+\frac1{2m}\cdot m=1+\frac{m-1}{2m}$。$\square$
\end{proof}

\begin{corollary}[$\kappa_2(3)=\frac87$，已证]\label{cor:k23}
$\kappa_2(3)=\max_{\|v\|=1}\lambda_{\max}\bigl((I_3-vv^\top)_+\bigr)=\frac87$，
极大化点 $v=\bigl(\sqrt{2/7},\sqrt{2/7},-\sqrt{3/7}\bigr)$，
此时 $(I_3-vv^\top)_+$ 的特征多项式\textbf{精确}为 $\frac1{343}(7L-8)(7L-5)(7L-1)$。
\end{corollary}

\begin{proof}
$\mathcal P_2(3)=\{I_3-vv^\top:\|v\|=1\}$（$3$ 维中秩二投影就是补一条直线），
故不需要定理~\ref{thm:convex}，直接对 $v$ 分类。按坐标符号分三种情形：
\emph{全同号}（或含零）时 $Q_+$ 的非对角正元个数为 $0$ 或需要单独讨论，但此时 $Q_+$ 是对角占优的，
$\lambda_{\max}(Q_+)\le\max_i(1-v_i^2)\le1<\frac87$；
\emph{$2$ 正 $1$ 负}（或 $1$ 正 $2$ 负，由 $v\mapsto-v$ 对称）时正是命题~\ref{prop:corank1}(a) 的 $m=2$ 情形，
最大值 $\frac{4\cdot2}{7}=\frac87$，$c^2=\frac27$、$d^2=\frac37$。
特征多项式由符号计算直接展开。$\square$
\end{proof}

\begin{remark}\label{rmk:k23read}
(i) 推论~\ref{cor:k23} 是 $k=2$ 情形常数 $c_2^\Sigma=\frac87$ 的一条\textbf{新的、三行的}证明；
已知证明要走完全二部吹起的谱计算。这条副线值得单独写成一个小节（见第~\ref{sec:pq7} 小节的 P4）。
(ii) $m=3$ 时 \eqref{eq:fm1} 给 $\frac65$，$m=2$ 时 \eqref{eq:fmm} 给 $\frac54$；
两者都在 $N=4$ 的余秩一层内，而 $\frac54>\frac65$，所以 $\kappa_3(4)$ 的余秩一见证是\emph{平衡}的那一个：
$v=(\frac12,\frac12,-\frac12,-\frac12)$，$\lambda_{\max}=\frac54$。这与注~\ref{rmk:n4exact} 的
$Q=I_4-vv^\top$ 是\textbf{同一个}矩阵——于是 $\kappa_3(4)=\frac54$ 现在有了两条独立证明：
行方差机制（定理~\ref{thm:LN} 在 $N=4$ 取等）与余秩一层穷举。
(iii) $f_{m,m}=\frac{3m-1}{2m}\to\frac32$、$f_{m,1}=\frac{4m}{3m+1}\to\frac43$：
余秩一层在 $N$ 大时给出的常数\emph{超过} $\rho_0$，
这提醒我们 $Q=I_N-vv^\top$ 的秩是 $N-1$ 而不是 $3$，它只对 $\kappa_{N-1}(N)$ 有意义，
不能拿来反驳 $\kappa_3\le\rho_0$。
\end{remark}

\subsection{斜率恒等式 $\rho=\frac{\alpha+\ell}2$、$\alpha=1$ 层与通道线性规划}\label{sec:slope1}

\begin{lemma}[斜率恒等式，已证]\label{lem:slope1}
在命题~\ref{prop:perronrig} 的记号下（$\alpha=x^\top Qx$、$\beta=x^\top Q_-x$、$\ell=x^\top|Q|x$、$\rho=\alpha+\beta$），
恒有
\begin{equation}\label{eq:slope1}
\ell=\alpha+2\beta,\qquad\text{从而}\qquad \rho=\frac{\alpha+\ell}2 .
\end{equation}
于是主猜想 $\kappa_3\le\rho_0$ \textbf{逐字等价}于「第四条不等式」$\alpha+\ell\le2\rho_0$，
并且冠军处取等（$\ell=1.586541134688899=2\rho_0-\alpha$，精确到 $10^{-15}$）。
\end{lemma}

\begin{proof}
$|Q|=Q_++Q_-=Q+2Q_-$（逐点），故 $\ell=\alpha+2\beta$；$\rho=\alpha+\beta$ 是命题~\ref{prop:perronrig} 的 (b)。
两式相加除以 $2$ 得 \eqref{eq:slope1} 的第二个式子。$\square$
\end{proof}

\begin{remark}\label{rmk:slope1read}
(i) 引理~\ref{lem:slope1} 把注~\ref{rmk:rigiditylimit} 里那条「缺失的第四不等式」的\emph{形状}钉死了：
它必须是 $\alpha+\ell$ 的线性界，斜率\textbf{只能}是 $1$（任何 $\alpha+c\ell\le C$ 形式在 $c\ne1$ 时
都无法被冠军饱和，因为冠军同时给出 $\rho$ 与 $\alpha,\ell$ 的值）。
(ii) 但它也说明这条路\textbf{没有提供新信息}：$\alpha+\ell\le2\rho_0$ 与 $\rho\le\rho_0$ 是同一句话。
真正的杠杆在于 $\ell=x^\top|Q|x$ 是一个\emph{绝对值}量，
而绝对值量有独立的工具（例如 $|Q|$ 的谱、$\|Q\|_1$、通道分解）。
第~\ref{sec:nogo} 节的定理~\ref{thm:zonal} 正是把「对 $|Q|$ 做逐点主函数」这一整族方案封顶在 $\Phi$；
所以剩下的可能性只有\textbf{耦合}通道（同时用到 $x$ 与 $Q$ 的联合结构）。
\end{remark}

\textbf{$\alpha=1$ 层（数值）。}若 $\alpha=1$，则 $Qx=x$，$x$ 落在 $Q$ 的像空间里，$\rho=1+\beta$。
我们在这一层上做了全局搜索（差分进化 $+$ 多起点，$N=6$），得到 $\max\beta=0.258695<\tau=\rho_0-1=0.281568$，
即 $\alpha=1$ 层给出 $\rho\le1.2587<\rho_0$。这说明\textbf{冠军不在 $\alpha=1$ 层的边界上}，
「$\alpha$ 尽量大」与「$\beta$ 尽量大」之间存在真实的权衡；zonal 见证在这一层不可行。

\textbf{通道线性规划（数值）。}把 $\ell=x^\top|Q|x$ 按「通道」（$q_{ij}$ 的符号与 $x_ix_j$ 的大小组合）拆成
$8701$ 个构型，对每个构型写出线性约束并用 HiGHS 求解，得到：
度 $2$ 通道的最优值 $z=2.7320505=1+\sqrt3$（与 Mohar 常数一致），
度 $4,6,8,10$ 通道相对 $2\rho_0=2.563136$ 的缺口分别是 $+7.6\times10^{-3}$、$+5.2\times10^{-3}$、
$+4.0\times10^{-5}$、$+1.1\times10^{-5}$。
结论：\emph{单独}的 $x$-加权通道族并\textbf{不}被 $\Phi$ 封顶（度 $2$ 通道能到 $1+\sqrt3>2\Phi$），
真正封顶的是\textbf{联合}通道界。这把第~\ref{sec:nogo} 节的否定性结论精确化了一层：
损失不在单个通道，而在通道之间的耦合。

\subsection{若干新的精确恒等式，以及一条被否证的充分性}\label{sec:miscid}

\textbf{(1) $\|Q_-\|_F^2=\rho_0\tau$，精确。}冠军处
\[
\|Q_-\|_F^2=10d_2^2=\frac{16+5\sqrt7}{81}=0.360848846362012=\rho_0(\rho_0-1)=\rho_0\tau ,
\]
符号计算核验差为 $0$。这个值同时等于第~\ref{sec:fourcs} 小节的 $S_O=\sum_ir_i$（不等式 [g] 在冠军处取等），
所以「负部的 Frobenius 能量」$=\rho_0\tau$ 是一个\emph{恒等式级别}的事实，不是近似。

\textbf{(2) $Q_+$ 的完整谱，精确。}由 \eqref{eq:psip} 与定理~\ref{thm:conefam}，冠军的 $Q_+$ 的谱是
\[
\mathrm{spec}(Q_+)=\Bigl\{\rho_0,\ \ \tfrac{7-\sqrt5+(1+\sqrt5)s^*}{10}\ (\text{二重}),\ \
\tfrac{2-(\sqrt5-1)s^*}{10},\ 0,\ 0\Bigr\}
=\{1.281568197,\ 0.692637910^{[2]},\ 0.117401872,\ 0,\ 0\}.
\]
$Q_-$（即 $|Q|$ 的负部支撑）的谱是 $\{0.379920437,\ 0.117401872^{[2]},\ 0,\dots\}$，
其中 $0.117401872$ 与 $Q_+$ 的最小非零特征值\emph{相同}——这是 $Z_5$ 对称性的直接后果
（同一个 $\psi_2$ 分量在 $Q$ 中是 $0$，在 $Q_+$ 与 $Q_-$ 中各留一份）。

\textbf{(3) $|Q|$-路线的一行式与 Mohar 常数的第二次导出。}
$Q_+=\frac12(Q+|Q|)$ 给
\[
\lambda_{\max}(Q_+)\le\tfrac12\bigl(\lambda_{\max}(Q)+\lambda_{\max}(|Q|)\bigr)\le\tfrac12\bigl(1+\|Q\|_F\bigr)=\tfrac{1+\sqrt3}2,
\]
最后一步用 $\lambda_{\max}(|Q|)\le\||Q|\|_F=\|Q\|_F=\sqrt3$。这就是 Mohar 的 $\frac{1+\sqrt3}2$，
与定理~\ref{thm:LN} 的 $\lim_{N\to\infty}L(N)$ \textbf{完全相同}——
「行方差机制的 $N\to\infty$ 极限」与「一次 Cauchy--Schwarz」给同一个常数，
这不是偶然：定理~\ref{thm:LN} 的证明在双均匀层里退化成一次 CS。
冠军处 $\lambda_{\max}(|Q|)=1.602533590389977<\varphi$，故这条路线给 $1.301267$，
而 ETF 处 $\lambda_{\max}(|Q|)=\frac12+\frac{\sqrt5}2=\varphi$ 恰好给 $\frac{1+\varphi}2=\Phi$。
\textbf{\STW\ 的 $\Phi$ 就是 $|Q|$-路线在 ETF 上的取值。}

\textbf{(4) \STW\ 命题 6.1 的一行重导。}
对图 $G$，令 $M=\Sigma\circ Q$（$\Sigma=1\{Q>0\}$ 为掩码，$Q$ 为 $\Lam_3$ 的 Ky Fan 对偶投影），
则 $\Lam_3(G)=\tr(M\cdot A)$ 型的计算给出 $\Lam_3(G)\le n\cdot\lambda_{\max}(Q_+)-3$，
所以 $\kappa_3:=\sup_N\max_{Q\in\mathcal P_3(N)}\lambda_{\max}(Q_+)$ 一旦 $\le\rho_0$ 就得到主定理。
这条重导比 \STW\ 原文短得多，且把「$-3$」的来源讲清楚了（$\tr Q=3$）。

\textbf{(5) KKT 对偶变量（已验证）。}极值问题的 KKT 条件给出一个对偶矩阵 $\Lambda$，
冠军处 $\tr\Lambda=\lambda=\rho_0$ \emph{精确}，$\mathrm{spec}(\Lambda)=\{0.762031,\ 0.259768^{[2]}\}$。
这条恒等式说明冠军的 KKT 乘子没有退化，局部刚性（定理~\ref{thm:champ} 附近的 Hess 认证）是\emph{严格}的局部极大，
不是鞍点。

\textbf{(6) 一条被否证的充分性（否证也是交付物）。}
第~\ref{sec:fourcs} 小节的不等式 [g] $\beta\le\sqrt{S_OW_O}$ 在冠军处取等，
因此很自然地会猜想「$[g]$ 取等 $+$ 秩三 $\Rightarrow$ 冠军」。这是\textbf{假的}：
在 $N=4$ 上，$[g]$ 单独（配上秩三与 Perron 条件）允许 $\rho=1.4186>\rho_0$。
所以 [g] 的取等只是冠军的\emph{必要}特征，不是刻画；四条 CS 同时取等（注~\ref{rmk:quadtight}）才是。
这也再次说明：把不等式一条条加进 REL 时，必须\textbf{同时}检查它们能否被非冠军构型饱和。

\subsection{更新后的方法优先队列（第七版）}\label{sec:pq7}

第~\ref{sec:pq6} 小节的队列（第六版）在本轮之后需要大改：旧队首 P1 被命题~\ref{prop:chamberdead} 否证，
旧 P6 被引理~\ref{lem:slope1} 吸收。新队列如下（完整版见 \texttt{研究方法优先级协议.md} 第七版）：

\begin{center}\small
\begin{tabular}{p{1.1cm}p{3.1cm}p{8.6cm}p{1.6cm}}
\toprule
位次 & 名称 & 内容与判据 & 现状\\
\midrule
\textbf{P1} & 全局：秩三约束 $\times$ REL &
把引理~\ref{lem:rowid}(iii) 的三角形恒等式 \eqref{eq:tri}、四点 Gram 行列式 $\det G_4\ge0$、
以及 $\R^3$ 中等角线 $\le6$ 的约束加进第~\ref{sec:relax} 小节的 REL，
在 $Z_5$-对称五角星图案类内闭合 $3.2128\times10^{-3}$ 的缺口。
判据：REL$^+$ 的天花板降到 $\rho_0$，或找到一个 $>\rho_0$ 的合法构型（后者同样是重大结果）。 & 主攻\\
P2 & 整数化锐化（脱离 $Z_5$） &
把命题~\ref{prop:ladder6} 的整数化阶梯推广到 $N=6$ 的一般图案（不正则的钝角图），
特别是 $n=1$ 子层的精确上确界。判据：双均匀层的上界降到 $\le1.29$。 & 排队\\
P3 & 证明 $\mathrm{REL}(N)=L(N)$ &
定理~\ref{thm:LN} 只对双均匀层成立；REL 的自由上确界在 $N=5,6,7,8,10,12$ 数值上等于 $L(N)$（误差 $\le4\times10^{-6}$）。
若对一切 $N$ 证明相等，则 $\kappa_3(N)\le L(N)$ 成为无条件定理（$N\le6$ 优于或等于 \STW）。 & 排队\\
P4 & 副线：$\kappa_2=\frac87$ 对一切 $N$ &
推论~\ref{cor:k23} 已给 $N=3$。若能把余秩一层 $+$ 定理~\ref{thm:convex} 的组合推到一般 $N$，
就得到 $c_2^\Sigma=\frac87$ 的新证明。判据：$\mathcal P_2(N)$ 上的完整符号分类。 & 副线\\
P5 & 旧队列降级项 &
$k\ge3$ 断点定位、原子化、$\min(a,b)=1$ 族的 $\frac65$、范数改写：全部降为「有余力再做」。 & 降级\\
\midrule
弹出 & 旧 P1（腔掩码不存在性） & 命题~\ref{prop:chamberdead} 否证：$K_{2,2,2}$ 是腔掩码，$19/19$ 未闭合类都是腔掩码。 & 已否证\\
吸收 & 旧 P6（联合 $\alpha+\ell$） & 引理~\ref{lem:slope1} 说明它与主猜想逐字等价，无独立内容。 & 已吸收\\
\bottomrule
\end{tabular}
\end{center}

\textbf{P1 的三种具体形态}（按我们的偏好排序，都是数学方法而非搜索）：
\begin{enumerate}[label=(\alph*),leftmargin=2.2em]
\item \textbf{三角形恒等式路线}。由 \eqref{eq:tri}，$\sum_{i<j<k}q_{ij}q_{jk}q_{ki}$ 被对角线单独决定。
在 $Z_5$-对称层里这个和可以显式算出（只有三种三角形：$(h,j,j{+}1)$、$(h,j,j{+}2)$、$(j,j{+}1,j{+}2)$），
于是 \eqref{eq:tri} 变成 $(g_h,g_c,d_1,d_2,A)$ 的一个\emph{方程}。第~\ref{sec:conefamily} 小节说明这个方程在等变层里是恒等式，
所以必须\textbf{先降低对称性}（见 (c)）它才会咬人。
\item \textbf{四点 Gram 路线}。$Q=VV^\top$、$V\in\mathrm{St}(N,3)$，故任意四行的 Gram 行列式 $\det G_4\ge0$，
且在 $N\ge4$ 时四点必线性相关、$\det G_4=0$ 是\emph{等式}。冠军的六个向量在 $\R^3$ 中，
所有 $\binom64=15$ 个四点行列式都为零，这是 $15$ 个方程——远比 REL 的四条不等式强。
把它们写成 $(g,\pi,r,w)$ 语言的代理形式是 P1 的核心技术工作。
\item \textbf{降对称路线}。把 $Z_5$ 降为 $Z_2$（hub $+$ 锥点上的一个对合），
等变参数从 $1$ 维升到 $4$ 维，秩三约束第一次非平凡。
此时第~\ref{sec:conefamily} 小节的技巧（分类 $\mu$ 特征值）仍然适用，只是分量的重数变了。
\end{enumerate}
"""
io.open("_frag_rowvar_4b.tex","w",encoding="utf-8").write(tex.lstrip("\n"))
print("written B, lines:",len(tex.split(chr(10))))
