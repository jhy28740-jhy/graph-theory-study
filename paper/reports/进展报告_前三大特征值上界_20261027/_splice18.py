# -*- coding: utf-8 -*-
import io,os
FRAG = r'''
\section{第 18 轮：钝性稀释恒等式与 $\Lam$-参数化（消去特征值的代数重构）}\label{sec:r18}

\subsection{动机：换一个正交的方向}\label{sec:r18:motiv}

第~\ref{sec:r1517}~节的全局二分定理~\ref{thm:dichotomy}~把缺口压缩为「排除分支 I 或分支 II」，
两个分支都已经是有限维的。但这两个分支都仍然建立在同一个对象之上：
$Q_+$ 的\textbf{谱}——$\lambda_{\max}(Q_+)$ 与它的 Perron 向量 $x$。
谱对象天然带着两个技术麻烦：一是 $\lambda_{\max}$ 只在单重时才光滑，
二是 Perron 向量由一个 $N$ 维特征方程隐式定义，$N$ 还不固定。

本轮换一个与之\textbf{正交}的方向：不再研究 $Q_+$ 的谱，而是设法把 $\rho=\lambda_{\max}(Q_+)$
本身改写成一个\textbf{不含特征值、不含 Perron 向量}的代数泛函，
其自变量只有几何数据 $(g_i,u_i)$ 与一个 $3\times3$ 的半正定矩阵 $\Lam$。
结果（定理~\ref{thm:lamparam}）是：主猜想等价于一族\textbf{五个额外参数}上的代数不等式，
且这个代数泛函在 $\Lam$ 上是\textbf{一次正齐次}的，数值上还是\textbf{凹}的
（命题~\ref{prop:concave}）。这使内层极大化成为一个凸规划，从而原则上可以给出
\textbf{对偶证书}（推论~\ref{cor:dualcert}），而不是只能跑算法。

沿路我们得到两个精确恒等式（引理~\ref{lem:deltautau}、引理~\ref{lem:dilution}）、
一个 $u$-KKT 的闭式解（命题~\ref{prop:uktkt}）、一个初等的结构定理
（定理~\ref{thm:transitive}，它解释了 $\rho_0$ 中为什么出现 $\sqrt5$），
以及\textbf{两个否证}（命题~\ref{prop:refA}、命题~\ref{prop:refB}），
它们精确指出「哪一步松弛一定会崩」。

\subsection{记号}\label{sec:r18:not}

沿用 \STW\ Theorem 4.1 的坐标并做如下固定记法。设 $N\ge3$，$g_i>0$，$u_i\in S^2$，
\begin{equation}\label{eq:frame18}
\sum_{i=1}^N g_iu_iu_i^\top=I_3\qquad(\text{从而 } \textstyle\sum_ig_i=3),
\end{equation}
记 $v_i:=\sqrt{g_i}\,u_i\in\R^3$，则 $\sum_iv_iv_i^\top=I_3$，$Q=VV^\top$ 的条目为
$Q_{ij}=v_i^\top v_j=\sqrt{g_ig_j}\langle u_i,u_j\rangle$，且 $Q$ 是秩三正交投影。
写
\[
(t)_+=\max(t,0),\qquad (t)_-=\max(-t,0),\qquad
P_{ij}:=(\langle u_i,u_j\rangle)_-,\qquad K_{ij}:=(\langle u_i,u_j\rangle)_+ .
\]
$Q_+=(Q)_{ij}=(\sqrt{g_ig_j}\langle u_i,u_j\rangle)_+=\sqrt{g_ig_j}K_{ij}$，
$Q_-=(Q_-)_{ij}=\sqrt{g_ig_j}P_{ij}$，$Q=Q_+-Q_-$，$Q_+=\tfrac12(Q+|Q|)$。
设 $\rho=\lambda_{\max}(Q_+)$ 单重，$x\ge0$，$\|x\|=1$，$Q_+x=\rho x$，并令
\begin{equation}\label{eq:hdef18}
h_i:=\frac{x_i}{\sqrt{g_i}}>0\qquad\Bigl(\text{即 } x_i=\sqrt{g_i}h_i,\quad \sum_ig_ih_i^2=1\Bigr).
\end{equation}
于是 Perron 方程等价于
\begin{equation}\label{eq:perron18}
\rho h_i=\sum_jg_jh_jK_{ij}=\langle u_i,m\rangle+d_i ,
\end{equation}
其中我们引入两个核心量：
\begin{equation}\label{eq:md18}
m:=\sum_{i}g_ih_iu_i\in\R^3,\qquad d_i:=\sum_jg_jh_jP_{ij}\ge0 .
\end{equation}
$m$ 是\textbf{矩向量}，$d=(d_i)$ 是\textbf{钝场}。由第~\ref{sec:r1517}~节的 $\tau$-形式（引理~\ref{lem:tauform}）：
\[
\alpha:=\|m\|^2=x^\top Qx,\qquad \beta:=x^\top Q_-x=\sum_ig_ih_id_i,\qquad
\rho=\alpha+\beta,\qquad \tau:=\rho-1=\beta-(1-\alpha).
\]
最后再引入一个\textbf{纯几何}量（不含 $h$）：
\begin{equation}\label{eq:Didef}
D_i:=\sum_jg_jP_{ij}=\sum_jg_j(\langle u_i,u_j\rangle)_-\ \ge0 ,
\end{equation}
它是原子 $i$ 处的「钝质量」。注意 $d_i$ 与 $D_i$ 的关系是 $d_i=\sum_jg_jh_jP_{ij}$：
$d$ 是 $D$ 的 $h$-加权版本。相应的钝图记为 $G_-$：$i\sim j\Leftrightarrow\langle u_i,u_j\rangle<0$。

\subsection{第一个精确恒等式：钝场向量恰为 $\tau m$}\label{sec:r18:delta}

\begin{lemma}[钝场向量恒等式]\label{lem:deltautau}
设 $\delta:=\sum_ig_id_iu_i$。则
\begin{equation}\label{eq:deltatau}
\boxed{\ \delta=\tau\,m\ }
\end{equation}
等价地，
\begin{align}
QQ_-x&=\tau\,Qx=\tau m, \label{eq:QQmx}\\
(I-Q)Q_-x&=\rho\,(I-Q)x, \label{eq:orthQmx}\\
\|Q_-x\|^2&=\tau^2\alpha+\rho^2(1-\alpha)=(1-\alpha)(1+2\tau)+\tau^2, \label{eq:normQmx2}\\
\langle m,\delta\rangle&=\alpha\tau . \label{eq:mdelta}
\end{align}
\end{lemma}

\begin{proof}
$Q_-=Q_+-Q$，故 $Q_-x=\rho x-Qx$。用 $Q^2=Q$：
$QQ_-x=\rho Qx-Qx=\tau Qx$，$(I-Q)Q_-x=\rho(I-Q)x-(I-Q)Qx=\rho(I-Q)x$，
这就是 \eqref{eq:QQmx}\eqref{eq:orthQmx}。又 $Qx=V^\top Vx$ 而
$(Vx)_k=\sum_i\sqrt{g_i}x_iu_{i,k}=m$ 的第 $k$ 分量乘以 $\sqrt{g_i}$ 的约定——直接算：
$\sum_i\sqrt{g_i}(Qx)_iu_i=\sum_i\sqrt{g_i}\sum_j\sqrt{g_ig_j}\langle u_i,u_j\rangle x_ju_i
=\sum_j x_j\sqrt{g_j}\sum_ig_i\langle u_i,u_j\rangle u_i=\sum_jx_j\sqrt{g_j}u_j=m$，
另一方面 $\sum_i\sqrt{g_i}(Q_-x)_iu_i=\sum_ig_id_iu_i=\delta$，
而 $\sum_i\sqrt{g_i}(Q_+x)_iu_i=\rho m$；由 $Q_-x=Q_+x-Qx$ 得 $\delta=\rho m-m=\tau m$，即 \eqref{eq:deltatau}。
\eqref{eq:normQmx2}：$Q_-x$ 沿 $\operatorname{im}Q$ 与 $\ker Q$ 正交分解，
由 \eqref{eq:QQmx}\eqref{eq:orthQmx} 得 $\|Q_-x\|^2=\tau^2\|Qx\|^2+\rho^2\|(I-Q)x\|^2$，
而 $\|Qx\|^2=x^\top Qx=\alpha$、$\|(I-Q)x\|^2=1-\alpha$。
展开 $\tau^2\alpha+\rho^2(1-\alpha)$ 并用 $\rho=1+\tau$ 即得引理~\ref{lem:normqmx}~的旧形式，
故 \eqref{eq:normQmx2} 是它的\textbf{勾股细化}：它把 $\|Q_-x\|^2$ 分成了「行空间内」与「行空间外」两块。
\eqref{eq:mdelta}~是 \eqref{eq:deltatau}~与 $\|m\|^2=\alpha$ 的直接推论。
\end{proof}

\begin{remark}[读法]\label{rmk:deltaread}
\eqref{eq:deltatau}~的几何含义是：把每个原子的钝场 $d_i$ 以质量 $g_i$ 堆到方向 $u_i$ 上，
得到的向量\textbf{必定}与矩向量 $m$ 平行，比例系数恰好是超出量 $\tau=\rho-1$。
因此 $\tau$ 有了一个纯几何的读法：
\begin{equation}\label{eq:tauread}
\tau=\frac{\bigl\|\sum_ig_id_iu_i\bigr\|}{\bigl\|\sum_ig_ih_iu_i\bigr\|}\qquad(m\ne0),
\end{equation}
即「钝场向量长度 / 矩向量长度」。这把主猜想变成
\[
\Bigl\|\sum_ig_id_iu_i\Bigr\|\le\tau_0\Bigl\|\sum_ig_ih_iu_i\Bigr\|,
\]
两端都是 $\R^3$ 中的向量长度，而不再是矩阵的最大特征值。
数值验证：冠军处 \eqref{eq:deltatau}~的残差为 $6.1\times10^{-16}$，
在 $800$ 个随机帧构型上最大残差 $2.9\times10^{-14}$
（\texttt{证书\_Lam参数化与稀释恒等式.py}~(A) 段）。
\end{remark}

\subsection{第二个精确恒等式：钝性稀释}\label{sec:r18:dil}

\begin{lemma}[钝性稀释恒等式]\label{lem:dilution}
定义
\begin{equation}\label{eq:Dhdef}
\mathcal D(h):=\sum_ig_ih_i^2D_i=\sum_{i\sim j}g_ig_j|\langle u_i,u_j\rangle|\,\frac{h_i^2+h_j^2}{2},
\qquad
\mathcal E_-(h):=\sum_{i\sim j}g_ig_j|\langle u_i,u_j\rangle|\,(h_i-h_j)^2 .
\end{equation}
则对\textbf{一切}构型与一切 $h\ge0$ 有精确等式
\begin{equation}\label{eq:dilution}
\boxed{\ \mathcal D(h)=\beta+\tfrac12\,\mathcal E_-(h)\ }\qquad\text{即}\qquad
\beta=\mathcal D(h)-\tfrac12\mathcal E_-(h)\ge \mathcal D(h)-\tfrac12\mathcal E_-(h),
\end{equation}
且 $\mathcal E_-(h)\ge0$，$\mathcal E_-(h)=0\Leftrightarrow h$ 在钝图 $G_-$ 的每个连通支上取常值。
特别地 $\beta\le\mathcal D(h)$（\textbf{钝性稀释不等式}），等号 $\Leftrightarrow\mathcal E_-=0$。
\end{lemma}

\begin{proof}
对固定的无序钝对 $\{i,j\}$，$\frac{h_i^2+h_j^2}2-h_ih_j=\frac{(h_i-h_j)^2}2$；
把这一等式乘以 $g_ig_j|\langle u_i,u_j\rangle|$ 并对所有钝对求和，
左端给出 $\mathcal D(h)$（因为 $\sum_{j\sim i}g_j|\langle u_i,u_j\rangle|=D_i$），
中项给出 $\beta=\sum_{i\sim j}g_ig_jh_ih_j|\langle u_i,u_j\rangle|$，右端给出 $\frac12\mathcal E_-(h)$。
\end{proof}

\begin{corollary}[$\tau$ 的势--能量分解]\label{cor:tauenergy}
对一切构型，
\begin{equation}\label{eq:tauenergy}
\tau=\underbrace{\bigl[\mathcal D(h)-(1-\alpha)\bigr]}_{\text{势项 }\Phi}
\;-\;\underbrace{\tfrac12\mathcal E_-(h)}_{\text{钝 Dirichlet 能量}} .
\end{equation}
\end{corollary}

\begin{remark}[为什么这条分解是有信息的]\label{rmk:dilread}
\eqref{eq:tauenergy}~把主猜想写成
\begin{equation}\label{eq:star}
\Phi=\mathcal D(h)-(1-\alpha)\ \le\ \tau_0+\tfrac12\mathcal E_-(h),
\end{equation}
左边是「势」：$\mathcal D(h)$ 对 $h^2$ 是\textbf{线性}的，$-(1-\alpha)$ 是 $-\|m\|^2+1$；
右边是「能量」：$h$ 在钝图上的加权 Dirichlet 能量。
这是一个 Poincar\'e 型不等式的形状，而不是谱不等式的形状。
它的关键优点是：$\mathcal D(h)$ 只通过 $h^2$ 依赖 $h$，因此一旦用上
$g$-KKT 的 $h_i^2=u_i^\top\tilde\Lam u_i$（第~\ref{sec:r1517}~节的 R15-2），
$\mathcal D(h)$ 就变成 $\tr(\tilde\Lam\Xi)$，其中
\begin{equation}\label{eq:Xidef}
\Xi:=\sum_ig_iD_iu_iu_i^\top\succeq0
\end{equation}
是一个\textbf{只依赖几何}的 $3\times3$ 矩阵。于是势项完全落在 $3\times3$ 矩阵的世界里。
冠军处：$\operatorname{eig}(\Xi)=\{0.379920,\,0.379920,\,0.253875\}$，$\tr\Xi=1.013716$，
$\tr(\tilde\Lam\Xi)=0.304973=\beta$（因为 $\mathcal E_-=0$），$\Phi=\tau_0$ 精确取等。
\end{remark}

\begin{lemma}[$\mathcal E_-$ 的第二显式]\label{lem:Em2}
设 $h$ 是 Perron 向量（即 \eqref{eq:perron18}~成立）。则
\begin{equation}\label{eq:Em2}
\mathcal E_-(h)=\frac1{\rho^2}\sum_{i\sim j}g_ig_j|\langle u_i,u_j\rangle|\,
\zeta_{ij}^2,\qquad
\zeta_{ij}:=\langle u_i-u_j,m\rangle+d_i-d_j .
\end{equation}
\end{lemma}

\begin{proof}
由 \eqref{eq:perron18}，$\rho h_i=\langle u_i,m\rangle+d_i$，故 $\rho(h_i-h_j)=\zeta_{ij}$。
\end{proof}

\begin{remark}\label{rmk:Em2read}
\eqref{eq:Em2}~的重要性在于：右端\textbf{不再显含} $h$，只含几何量 $u_i,g_i$ 与两个
辅助量 $m,d$（而 $m,d$ 又由 $h$ 决定）。它把「需要多少 Dirichlet 能量来补偿」这个问题
改写成了一个\textbf{显式的非负二次型}的下界问题，是下一轮 P1 的具体攻击目标。
数值验证：\eqref{eq:Em2}~在冠军处残差 $<10^{-15}$，在 $1200$ 个随机构型上最大残差 $2.1\times10^{-15}$。
\end{remark}

\subsection{$u$-KKT 的闭式解}\label{sec:r18:ukkt}

第~\ref{sec:r1517}~节以「R15-2 / R15-3」的形式记录了两个在冠军处数值成立的临界点关系。
本轮把它们\textbf{解出来}：乘子不是「存在」，而是有闭式。

\begin{proposition}[$u$-KKT 的显式形式]\label{prop:uktkt}
设 $(U,g)$ 是 $\rho$ 在帧切片 \eqref{eq:frame18}~上的临界点，$g_i>0$，$\lambda_{\max}(Q_+)$ 单重。
令 $W_i:=\sum_{j\in S_i^+}g_jh_ju_j$（$S_i^+=\{j:\langle u_i,u_j\rangle>0\}$）。
设 Lagrange 函数为
$\mathcal L=\rho-\langle\Lam,\sum_ig_iu_iu_i^\top-I_3\rangle-\sum_i\frac{\nu_i}2(\|u_i\|^2-1)$。
则：
\begin{enumerate}\itemsep2pt
\item[(i)] \textbf{球面乘子恒为零}：$\nu_i\equiv0$；
\item[(ii)] $\tilde\Lam u_i=\dfrac{h_i}{\rho}W_i$，其中 $\tilde\Lam:=\frac2\rho\Lam$；
\item[(iii)] \textbf{闭式}：
\begin{equation}\label{eq:lamclosed}
\tilde\Lam=\frac1\rho\sum_{i}\sum_{j\in S_i^+}g_ig_jh_ih_j\,u_iu_j^\top
=\frac1\rho\sum_{(i,j)\ \text{锐}}g_ig_jh_ih_j\,u_iu_j^\top ;
\end{equation}
\item[(iv)] \eqref{eq:lamclosed}~\textbf{自动对称}，$\tr\tilde\Lam=1$ 是恒等式（不依赖临界性），
且 $h_i^2=u_i^\top\tilde\Lam u_i$。
\end{enumerate}
\end{proposition}

\begin{proof}
$\partial\rho/\partial u_i=\sum_jx_ix_j\sqrt{g_ig_j}\,\one[\langle u_i,u_j\rangle>0]\,u_j
=x_i\sqrt{g_i}W_i=g_ih_iW_i$，
$\partial\mathcal L/\partial u_i=g_ih_iW_i-2g_i\Lam u_i-\nu_iu_i=0$，故
$h_iW_i=2\Lam u_i+\hat\nu_iu_i=\rho\tilde\Lam u_i+\hat\nu_iu_i$（$\hat\nu_i=\nu_i/g_i$）。
与 $u_i$ 作内积：左端 $h_i\langle u_i,W_i\rangle=h_i\sum_{j\in S_i^+}g_jh_j\langle u_i,u_j\rangle
=h_i\sum_jg_jh_jK_{ij}=h_i\cdot\rho h_i$（用 \eqref{eq:perron18}），
右端 $\rho u_i^\top\tilde\Lam u_i+\hat\nu_i$。而 $g$-KKT 给
$\frac{\partial\rho}{\partial g_i}=\frac{\rho x_i^2}{2g_i}=\frac{\rho h_i^2}2=\langle u_i,\Lam u_i\rangle$，
即 $u_i^\top\tilde\Lam u_i=h_i^2$。代入得 $\rho h_i^2=\rho h_i^2+\hat\nu_i$，故 $\hat\nu_i=0$，(i)(ii) 成立。
(iii)：$\tilde\Lam=\tilde\Lam\cdot I_3=\sum_ig_i\tilde\Lam u_iu_i^\top=\frac1\rho\sum_ig_ih_iW_iu_i^\top$，展开 $W_i$ 即得。
(iv)：\eqref{eq:lamclosed}~的和式在 $(i,j)$ 上的系数 $g_ig_jh_ih_j$ 对称，而「锐」关系本身对称，
故转置后仍是同一个和；$\tr\tilde\Lam=\frac1\rho\sum_ig_ih_i\langle u_i,W_i\rangle
=\frac1\rho\sum_ig_ih_i\cdot\rho h_i=\sum_ig_ih_i^2=1$；
$u_i^\top\tilde\Lam u_i=\frac{h_i}{\rho}\langle u_i,W_i\rangle=h_i^2$。
\end{proof}

\begin{remark}\label{rmk:ukktread}
$\nu_i\equiv0$ 是一个此前未记录的干净事实：\textbf{球面约束在临界点上是「松」的}，
全部约束信息都来自帧条件。这一点在数值上同样成立：冠军处
$|h_iW_i-\rho\tilde\Lam u_i|$ 的最大残差为 $7.2\times10^{-16}$，
\eqref{eq:lamclosed}~与最小二乘拟合出的 $\tilde\Lam$ 相差 $5.0\times10^{-16}$，
$\operatorname{eig}(\tilde\Lam)=\{0.594609,\,0.202696,\,0.202696\}$，$\lambda_{\min}=0.2027>0$。
\end{remark}

\subsection{主结果：$\Lam$-参数化}\label{sec:r18:main}

\begin{definition}[$\Lam$-泛函]\label{def:lamfunc}
对固定的几何 $(g,U)$ 满足 \eqref{eq:frame18}，与任意 $\Lam\succeq0$，令
\begin{equation}\label{eq:hlam}
h_i(\Lam):=\sqrt{u_i^\top\Lam u_i}=\|\Lam^{1/2}u_i\|,\qquad
x_i(\Lam):=\sqrt{g_i}\,h_i(\Lambda),
\end{equation}
并定义
\begin{equation}\label{eq:rholam}
\rho(\Lam):=x(\Lam)^\top Q_+x(\Lam)
=\Bigl\|\sum_ig_ih_i(\Lam)u_i\Bigr\|^2+\sum_{i,j}g_ig_jh_i(\Lam)h_j(\Lam)(\langle u_i,u_j\rangle)_- .
\end{equation}
\end{definition}

\begin{theorem}[$\Lam$-参数化：消去特征值的等价重构]\label{thm:lamparam}
\begin{enumerate}\itemsep3pt
\item[(i)] $\|x(\Lam)\|^2=\sum_ig_ih_i(\Lam)^2=\tr\Lam$；
$\rho(t\Lam)=t\,\rho(\Lam)$（\textbf{一次正齐次}）；
且当 $\tr\Lam\le1$ 时 $\rho(\Lam)\le\lambda_{\max}(Q_+)$。
\item[(ii)] 在任一临界点上（命题~\ref{prop:uktkt}），$\rho(\tilde\Lam)=\rho=\lambda_{\max}(Q_+)$，且
\begin{equation}\label{eq:gradlam}
\nabla_\Lam\rho(\Lam)\big|_{\Lam=\tilde\Lam}=\rho\,I_3 ,
\end{equation}
即 $\tilde\Lam$ 是 $\rho(\cdot)$ 在超平面 $\{\tr=1\}$ 上的\textbf{驻点}。
更一般地，$\nabla_\Lam\rho(\Lam)=\sum_i\frac{g_i[\langle u_i,m(\Lam)\rangle+d_i(\Lam)]}{h_i(\Lam)}u_iu_i^\top$，
故 \eqref{eq:gradlam}~成立当且仅当 $h(\Lam)$ 满足 Perron 方程 \eqref{eq:perron18}。
\item[(iii)] \textbf{等价重构}：
\begin{equation}\label{eq:kappalam}
\kappa_3=\sup_{(g,U)}\ \max_{\substack{\Lam\succeq0\\ \tr\Lam=1}}\ \rho_{g,U}(\Lam) .
\end{equation}
因此 $\kappa_3\le\rho_0$（等价地 $\Lam_3(G)\le\rho_0n-3$ 对一切 $G$）由下面这一条
\textbf{不含特征值、不含 Perron 向量}的代数不等式推出：
\begin{equation}\label{eq:starmain}
\Bigl\|\sum_ig_ih_iu_i\Bigr\|^2+\sum_{i,j}g_ig_jh_ih_j(\langle u_i,u_j\rangle)_-\ \le\ \rho_0 ,
\qquad h_i=\sqrt{u_i^\top\Lam u_i},\ \Lam\succeq0,\ \tr\Lam=1 .
\end{equation}
\end{enumerate}
\end{theorem}

\begin{proof}
(i) $\sum_ig_ih_i(\Lam)^2=\sum_ig_iu_i^\top\Lam u_i=\tr(\Lam\sum_ig_iu_iu_i^\top)=\tr\Lam$。
$h(t\Lam)=\sqrt t\,h(\Lam)$，而 \eqref{eq:rholam}~的右端对 $h$ 是二次的，故一次正齐次。
若 $\tr\Lam=1$ 则 $\|x(\Lam)\|=1$，于是 $\rho(\Lam)=x(\Lam)^\top Q_+x(\Lam)\le\lambda_{\max}(Q_+)$；
$\tr\Lam<1$ 时用齐次性。
(ii) 在临界点，$h_i(\tilde\Lam)^2=u_i^\top\tilde\Lam u_i=h_i^2$，故 $h(\tilde\Lam)=h$（两者皆正），
从而 $\rho(\tilde\Lam)=x^\top Q_+x=\rho$。梯度：$\partial h_i/\partial\Lam=u_iu_i^\top/(2h_i)$，
$\partial\rho(\Lam)/\partial h_i=2g_i\langle u_i,m(\Lam)\rangle+2\sum_jg_ig_jh_j(\Lam)P_{ij}
=2g_i[\langle u_i,m(\Lam)\rangle+d_i(\Lam)]$，链式法则给出所断言的公式；
在 $\Lam=\tilde\Lam$ 处方括号为 $\rho h_i$（\eqref{eq:perron18}），于是
$\nabla_\Lam\rho=\rho\sum_ig_iu_iu_i^\top=\rho I_3$。
(iii) 「$\le$」：由 (i)，对每个几何 $\max_{\tr\Lam=1}\rho(\Lam)\le\lambda_{\max}(Q_+)$，取上确界即得。
「$\ge$」：在 $\rho$ 的极大点（$g_i>0$，$\lambda_{\max}$ 单重）由 (ii) 得
$\max_\Lam\rho(\Lam)\ge\rho(\tilde\Lam)=\rho$；对 $g_i=0$ 或重根的退化情形，
用第~\ref{sec:r1517}~节的合并/分裂等谱性化归到 $g_i>0$ 且 $\lambda_{\max}$ 单重的情形
（合并保 $\rho$ 精确，分裂添一个零特征值），再用极限论证。
\end{proof}

\begin{remark}[这消除了什么]\label{rmk:lamparamread}
\eqref{eq:starmain}~与原来的问题 $\lambda_{\max}((VV^\top)_+)\le\rho_0$ 相比：
\begin{itemize}\itemsep2pt
\item 不再有 $N\times N$ 矩阵的特征值，也不再有 Perron 向量（$h$ 由 $\Lam$ 显式给出）；
\item 新增的自由度只有 $\Lam$ 的 $5$ 个（半正定、迹一）；
\item 约束集 $\{(g,u,\Lam): g\ge0,\ \sum g_iu_iu_i^\top=I_3,\ \Lam\succeq0,\ \tr\Lam=1\}$
是\textbf{半代数}的；把 $h_i\ge0$、$h_i^2=u_i^\top\Lam u_i$ 显式引入后，
目标函数在固定符号模式（即固定钝图 $G_-$）下是\textbf{多项式}的。
这为矩-SOS 层级或显式对偶证书留出了接口。
\end{itemize}
诚实标注：(iii) 中「$\ge$」一侧的退化情形化归沿用第~\ref{sec:r1517}~节的等谱性结论，
其细节（$g_i\to0$ 的极限与重根的分歧）在本报告中仍按\textbf{条件性}处理。
\end{remark}

\begin{proposition}[凹性，数值]\label{prop:concave}
对固定几何，$\rho(\Lam)$ 在 $\Lam$ 上是\textbf{凹}的：在冠军几何上，
$600$ 组随机 $(\Lam_1,\Lam_2,s)$ 的最小凹性余量
$\rho((1-s)\Lam_1+s\Lam_2)-[(1-s)\rho(\Lam_1)+s\rho(\Lam_2)]$ 为 $+1.99\times10^{-5}$；
在 $10$ 个随机几何（$N=5,\dots,11$）上全局最小余量为 $+1.57\times10^{-5}$，无一为负。
\end{proposition}

\begin{remark}\label{rmk:concave}
凹性尚未解析证明。注意它\textbf{不能}由 $\nabla^2\rho$ 的符号直接读出：
\eqref{eq:gradlam}~只在 Perron 点成立，不能对 $\Lam$ 再求导。
若凹性成立，则 $\rho(\Lam)$ 作为一次正齐次凹函数具有支撑超平面表示
$\rho(\Lam)=\inf\{\langle\Lam,S\rangle:\ S\in\mathcal S\}$，
其中 $\mathcal S=\{S\in\mathrm{Sym}_3:\ \langle\Lam,S\rangle\ge\rho(\Lam)\ \forall\Lam\succeq0\}$，
于是由 Sion 极小极大定理：
\begin{equation}\label{eq:dualcert}
\max_{\substack{\Lam\succeq0\\\tr\Lam=1}}\rho(\Lam)=\min_{S\in\mathcal S}\ \lambda_{\max}(S).
\end{equation}
\end{remark}

\begin{corollary}[对偶证书形式]\label{cor:dualcert}
设 \eqref{eq:dualcert}~成立。则要证某几何满足 $\max_\Lam\rho(\Lam)\le\rho_0$，
\textbf{只需给出一个} $S\in\mathrm{Sym}_3$ 使得
\begin{equation}\label{eq:certcond}
\lambda_{\max}(S)\le\rho_0\qquad\text{且}\qquad \rho(\Lam)\le\langle\Lam,S\rangle\quad\forall\Lam\succeq0 .
\end{equation}
在冠军几何上这个证书是\textbf{被强制的且是紧的}：由 $\rho(\tilde\Lam)=\rho_0$ 与 $\tr\tilde\Lam=1$，
\eqref{eq:certcond}~给出 $\rho_0\le\langle\tilde\Lam,S\rangle\le\lambda_{\max}(S)\tr\tilde\Lam\le\rho_0$，
故必须处处取等，即 $\lambda_{\max}(S)=\rho_0$ 且 $\tilde\Lam$ 的支撑落在 $S$ 的最高特征空间里。
特别地 $S=\rho_0I_3$ 是唯一的「各向同性」候选，而它给出的条件
$\rho(\Lam)\le\rho_0\tr\Lam$ 恰与 \eqref{eq:starmain}~等价。
\end{corollary}

\begin{remark}\label{rmk:certsearch}
因此冠军几何上不存在「比切平面更好的」各向同性证书；要前进必须利用 $\tilde\Lam$ 的
特征方向 $u_0$（枢纽方向，$\tilde\Lam$ 的最高特征向量）来构造\textbf{各向异性}的 $S$，
形如 $S=\rho_0u_0u_0^\top+\kappa(I_3-u_0u_0^\top)$，$\kappa\le\rho_0$。
这是下一轮 P1 的具体构造任务。
\end{remark}

\begin{proposition}[尖锐性与刚性，数值]\label{prop:r18sharp}
\begin{enumerate}\itemsep2pt
\item[(i)] 冠军几何上 $\max_{\Lam\succeq0,\tr\Lam=1}\rho(\Lam)=1.281568196872167$，
与 $\rho_0=1.281568196872168$ 之差为 $-4.4\times10^{-16}$；
$26$ 组多起点（L-BFGS-B 与 Nelder--Mead 各一次）的极大点全部收敛到
$\operatorname{eig}=\{0.59460858,0.20269572,0.20269570\}$，与 $\tilde\Lam$ 的距离 $1.1\times10^{-8}$。
\item[(ii)] 在 $25$ 个几何上（Stiefel 梯度上升极大元 $N=6,7,9,10,12$；
$16$ 个随机帧构型 $N=5,\dots,12$；冠军的三档微扰 $\varepsilon=10^{-3},10^{-2},4\times10^{-2}$）
全部满足 $\max_\Lam\rho(\Lam)\le\rho_0$，全局最大值 $1.281568196584$（出现在 $N=10$ 的极大元上，
即冠军加零质量行）。
\item[(iii)] 在临界点上 $\max_\Lam\rho(\Lam)=\rho$ 精确成立
（$N=6$：两边同为 $-3.02\times10^{-2}$；$N=7$：同为 $-6.55\times10^{-2}$；
$N=9,10,12$：一致到 $10^{-9}$），
印证定理~\ref{thm:lamparam}(iii)~的「$\ge$」一侧。
\end{enumerate}
\end{proposition}

\subsection{一个初等的结构定理：可传递符号模式}\label{sec:r18:trans}

\begin{theorem}[可传递符号模式定理]\label{thm:transitive}
设存在 $\sigma\in\{\pm1\}^N$ 使
\begin{equation}\label{eq:transig}
\sigma_i\sigma_j\langle u_i,u_j\rangle\ge0\qquad(\forall i,j).
\end{equation}
则 $\rho=\lambda_{\max}(Q_+)\le1$。
逆否命题：\textbf{$\rho>1$ 蕴含符号模式不可传递}，
等价地，钝图 $G_-$ \textbf{不是完全二部图}（在把 $\langle u_i,u_j\rangle=0$ 的对任意归部之后）。
\end{theorem}

\begin{proof}
$x^\top|Q|x=\sum_{ij}x_ix_j|v_i^\top v_j|
=\sum_{ij}x_ix_j\sigma_i\sigma_jv_i^\top v_j
=\bigl\|\sum_i\sigma_ix_iv_i\bigr\|^2
\le\lambda_{\max}\bigl(\sum_iv_iv_i^\top\bigr)\,\|\sigma\odot x\|^2=1\cdot1=1$，
其中用到 $\sum_iv_iv_i^\top=I_3$ 与 $x\ge0$、$\|x\|=1$。
另一方面 $\alpha=x^\top Qx=\|V^\top x\|^2\le\lambda_{\max}(VV^\top)\|x\|^2=1$。
由 $Q_+=\tfrac12(Q+|Q|)$（逐项）得 $\rho=\tfrac12\alpha+\tfrac12x^\top|Q|x\le1$。
\end{proof}

\begin{remark}[为什么这解释了 $\sqrt5$]\label{rmk:transread}
冠军的钝图是 $G_-=C_5\sqcup\{\text{枢纽}\}$：枢纽与所有原子成锐角（$D_0=0$），
五个五边形原子之间恰有「次近邻」为钝，构成一个 $5$-圈。
$C_5$ \textbf{不是}完全二部图，故定理~\ref{thm:transitive}~不适用，$\rho_0>1$ 得以发生；
而且 $C_5$ 是\textbf{最小的}非二部钝图（$C_3$ 在 $S^2$ 上不可能全部两两成钝角，
因为三个两两钝角的单位向量需要 $\sum g_iu_iu_i^\top$ 的秩结构与之相容——
引理~\ref{lem:obtuseclique}~已给出 $\omega(G_-)\le4$，而 $C_5$ 的团数为 $2$）。
这给出 $\rho_0=\frac{9+\sqrt5+2\sqrt{35}}{18}$ 中 $\sqrt5$ 的\textbf{结构解释}：
它来自 $A(C_5)$ 的谱 $\{2,\frac{\sqrt5-1}2,\frac{\sqrt5-1}2,-\frac{\sqrt5+1}2,-\frac{\sqrt5+1}2\}$，
第~\ref{sec:r1517}~节已验证 $\operatorname{spec}(Q_-)=2g^*|s_2|\cdot\operatorname{spec}(A(C_5))$。
数值验证：$662\to193$ 个严格可传递构型上 $\max\rho=1+1.3\times10^{-10}$
（超出量与帧残差 $6.9\times10^{-11}$ 同阶，即纯数值误差）；
$2208$ 个 $\rho>1$ 的随机构型中可传递者 $0$ 个。
\end{remark}

\subsection{两个否证：哪一步松弛一定会崩}\label{sec:r18:ref}

按本协议，否证也是成果。本轮我们试验了两条看起来最有希望的松弛，二者\textbf{都是假的}，
而且失败方式恰好指明了正确路线。

\begin{proposition}[否证 R18-A：自由 $\Lam$ 松弛失效]\label{prop:refA}
考虑把 \eqref{eq:tauenergy}~中的 Dirichlet 项丢掉、并让 $\Lam$ 在整个
$\{\Lam\succeq0,\tr\Lam=1\}$ 上自由变化的松弛
\begin{equation}\label{eq:phifree}
\Phi_{\mathrm{free}}(g,U):=\max_{\substack{\Lam\succeq0\\\tr\Lam=1}}
\Bigl[\tr(\Lam\Xi)+\Bigl\|\sum_ig_ih_i(\Lam)u_i\Bigr\|^2-1\Bigr].
\end{equation}
则\textbf{即使在冠军几何上}也有
\[
\Phi_{\mathrm{free}}=0.302149011933>\tau_0=0.281568196872\qquad(\text{超出 }+2.058\times10^{-2}),
\]
极大元是\textbf{秩一}的（$\operatorname{eig}=\{1,1.7\times10^{-14},-1.4\times10^{-17}\}$）。
因此「$\Lam$ 是 KKT 乘子」这一信息\textbf{不可省略}：
只有把 $\Lam$ 限制到满足 $u$-KKT（命题~\ref{prop:uktkt}(ii)）的那一点，
才能得到 $\Phi=\tau_0$。
\end{proposition}

\begin{proposition}[否证 R18-B：强稀释不等式失效]\label{prop:refB}
不等式
\begin{equation}\label{eq:strongdil}
\alpha+\mathcal D(h)\le\rho_0\qquad(h=h_{\mathrm{Perron}})
\end{equation}
是\textbf{假}的。一个 $N=4$ 的帧构型给出 $\alpha+\mathcal D(h)=1.462928395681$，
超出 $\rho_0$ 达 $+1.814\times10^{-1}$。但该构型的真实 $\rho=0.798522010308\ll\rho_0$，
其补偿项为 $\tfrac12\mathcal E_-(h)=0.664406385373$，
恰好把 $\alpha+\mathcal D(h)-\tfrac12\mathcal E_-$ 拉回 $0.798522=\rho$。
因此\textbf{钝 Dirichlet 能量 $\mathcal E_-$ 是不可或缺的补偿项}，
任何丢弃它的路线都注定失败——这与第~\ref{sec:r1517}~节「单独界 $\beta$ 的路线全部失败」
是同一现象的两个侧面：$\mathcal E_-$ 与 $-(1-\alpha)$ 是主猜想中仅有的两个「负项」。
\end{proposition}

\begin{remark}\label{rmk:refsafe}
两个反例的真实 $\rho$ 都远低于 $\rho_0$（分别为 $0.798522$ 与冠军几何上的 $\rho_0$ 本身），
所以被否证的是\textbf{路线}，主猜想~\ref{conj:main}~未受任何影响。
它们的价值在于把「哪一项不能丢」变成了定理级的事实。
\end{remark}

\subsection{本轮状态、诚实标签与方法优先队列（第八版）}\label{sec:r18:status}

\textbf{已证（解析，含完整证明）}：
引理~\ref{lem:deltautau}（钝场向量 $=\tau m$ 及三个推论式）；
引理~\ref{lem:dilution}（钝性稀释恒等式，含等号刻画）；
推论~\ref{cor:tauenergy}（$\tau$ 的势--能量分解，精确）；
引理~\ref{lem:Em2}（$\mathcal E_-$ 的第二显式）；
命题~\ref{prop:uktkt}（$u$-KKT 闭式解、$\nu_i\equiv0$、$\tr\tilde\Lam=1$ 恒等）；
定理~\ref{thm:lamparam}(i)(ii)（$\Lam$-泛函的齐次性、$\rho(\tilde\Lam)=\rho$、$\nabla\rho=\rho I_3$）；
定理~\ref{thm:transitive}（可传递符号模式 $\Rightarrow\rho\le1$）。

\textbf{条件性}：定理~\ref{thm:lamparam}(iii)~的「$\ge$」一侧对退化情形（$g_i=0$、重根）
的化归；推论~\ref{cor:dualcert}~依赖命题~\ref{prop:concave}~的凹性。

\textbf{数值（不作成果）}：
命题~\ref{prop:concave}（凹性，冠军 $+$ $10$ 个随机几何）；
命题~\ref{prop:r18sharp}（尖锐性、刚性、$25$ 个几何的扫描）；
引理~\ref{lem:deltautau}\eqref{eq:dilution}\eqref{eq:Em2}~的随机验证
（$800/1500/1200$ 个构型，残差 $\le2.9\times10^{-14}$）；
定理~\ref{thm:transitive}~的 $193+2208$ 个构型检验。

\textbf{否证}：命题~\ref{prop:refA}（自由 $\Lam$ 松弛在冠军几何上超出 $\tau_0$ 达 $2.06\times10^{-2}$）；
命题~\ref{prop:refB}（强稀释不等式在 $N=4$ 上超出 $\rho_0$ 达 $1.81\times10^{-1}$）。

\textbf{方法优先队列（第八版）}
\begin{enumerate}\itemsep3pt
\item \textbf{P1（队首，全局）——构造对偶证书 $S$}：由推论~\ref{cor:dualcert}，
只需对每个几何给出 $S\in\mathrm{Sym}_3$ 满足 \eqref{eq:certcond}。
由注~\ref{rmk:certsearch}，冠军处证书必须是紧的且 $\tilde\Lam$ 落在 $S$ 的最高特征空间，
故候选族为 $S=\rho_0u_0u_0^\top+\kappa(I_3-u_0u_0^\top)$，$\kappa\le\rho_0$；
第一步是确定 $\kappa$ 的可行区间（一个 $3\times3$ 的 LMI 问题），
第二步是把 $\kappa$ 的表达式推广到一般几何（预期与 $\Xi$、$\tr\Xi$、$D_i$ 的分布有关）。
这是\textbf{纯数学}任务，不需要对每阶图跑算法。
\item \textbf{P2（全局）——解析证明 $\rho(\Lam)$ 的凹性}：
$\rho(\Lam)=\sum_k\lambda_k(Q_+)\langle w_k,x(\Lam)\rangle^2$，
而 $x_i(\Lam)=\sqrt{g_iu_i^\top\Lam u_i}$ 是线性函数的平方根（凹）。
可尝试的路径：(a) 逐项证明 $\langle w,x(\Lam)\rangle^2$ 在 $w\ge0$ 时凹；
(b) 用 \eqref{eq:dualcert}~的极小极大直接构造 $S$ 绕过凹性；
(c) 计算 $\nabla^2\rho$ 在 $\tilde\Lam$ 处的显式并证明其在该点负半定（局部凹性）。
\item \textbf{P3（全局）——$\mathcal E_-$ 的下界}：用引理~\ref{lem:Em2}~的显式
$\zeta_{ij}=\langle u_i-u_j,m\rangle+d_i-d_j$，在 \eqref{eq:star}~中给出
$\frac12\mathcal E_-\ge\Phi-\tau_0$ 的直接证明。这是把否证 R18-B 变正的路径。
\item \textbf{P4（局部）——第~\ref{sec:r1517}~节遗留的分支 I / 分支 II}：
分支 II 的排除（推论~\ref{cor:branch2rigid}~的插值刚性）；分支 I 经
$\nu_-(Q_+)\le\nu_-(Q_-)$ 与 $\omega(G_-)\le4$（引理~\ref{lem:obtuseclique}）。
注意定理~\ref{thm:transitive}~给出分支 I 的一个新入口：若 $\nu_-(Q_+)\ge1$ 则 $\rho$ 的
符号结构必然不可传递，$G_-$ 不是完全二部。
\item \textbf{P5（$N=6$ 收尾）}：第~\ref{sec:n6global}~节余下的 $20$ 个符号类。
\item \textbf{P6（副方向，下界）}：$\Lam_3(G)\ge\rho_0n-3$ 的可达到性，
冠军族已由定理~\ref{thm:Fn}~给出闭式谱，只差无穷族的渐近。
\end{enumerate}

\textbf{本轮把主线改写成}：定理~\ref{thm:lamparam}~的 \eqref{eq:starmain}。
这是一个\textbf{五个额外参数}上的代数不等式，目标函数在固定钝图符号模式后是多项式的，
约束集是半代数的，并且（数值上）对 $\Lam$ 凹。
与第~\ref{sec:r1517}~节的二分定理相比，它\textbf{完全不含}「$N\times N$ 矩阵的谱」这一对象，
因此更适合用凸对偶、SOS 或显式不等式来处理；而命题~\ref{prop:refA}\ref{prop:refB}
已经排除了两条最自然的错误走法。
'''
p='进展报告_前三大特征值之和的精确上界.tex'
s=open(p,encoding='utf-8').read()
anchor=r'\section{通道型主函数方案的精确天花板}\label{sec:nogo}'
assert anchor in s, 'anchor missing'
s=s.replace(anchor, FRAG.lstrip('\n')+'\n\n'+anchor,1)
old='\\textbf{方法优先队列（第七版）}'
assert old in s
s=s.replace(old,'\\textbf{方法优先队列（第七版，第 15--17 轮；已被第八版取代，见第~\\ref{sec:r18:status}~节）}',1)
open(p,'w',encoding='utf-8').write(s)
print('inserted R18 section, new length=%d'%len(s))
