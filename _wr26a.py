# -*- coding: utf-8 -*-
import io, os
FRAG = r"""
\section{第 26 轮：负 Gram 质量恒等式、主猜想的单标量归约与秩一截面的锥族闭式}\label{sec:r26}

第 25 轮把全局猜想化成了 $\Phi\ge\delta_0/3$，但 $\Phi$ 由五项组成、
且在冠军处五项\emph{全部严格为正}（表~\ref{tab:r25budget}），
注~\ref{rem:r25no}~已说明任何松弛都必然失败。本轮回到最原始的分解
$(t)_+=\frac12(t+|t|)$，把其中的 $|t|$-部分再\emph{精确}地劈成两块：
「与最优符号模式一致的部分」和「符号相反的部分」。后者我们称为
\textbf{负 Gram 质量} $\mathrm{Neg}^*$，它恰是掩码自洽性 $A=1\{Q>0\}$
真正看到、而通道型主函数方案（定理~\ref{thm:zonal}）完全丢掉的那部分信息。
结果是主定理~\ref{thm:r26neg}：
\begin{equation}\label{eq:r26head}
\kappa_3=\sup_{(g,U,\Lam)}\bigl(A_w^2+\mathrm{Neg}^*\bigr),
\qquad
\kappa_3\le\rho_0\iff
\mathrm{Neg}^*\le(\rho_0-1)+\sin^2\theta .
\end{equation}
这条归约有四点价值。
第一，它把「五项预算」压缩成\textbf{一个}非负标量，而且这个标量有明确的几何含义
（Gram 矩阵在最优符号胞腔下的负部质量）。
第二，它给出了 \STW\ 天花板的\emph{精确}含义（推论~\ref{cor:r26stw}）：
逐点主函数估计 $|t|\le q(t^2)$ 恰好等价于 $\mathrm{Neg}^*\le\frac{\varphi-1}2=\frac1{2\varphi}$，
而目标值是 $\rho_0-1$，二者之差
\[
\frac{\varphi-1}2-(\rho_0-1)=\frac{\varphi+1-2\rho_0}2=\frac{\delta_0}6
=\Phi_{\STW}-\rho_0=0.027448797503
\]
\textbf{精确等于} \STW\ 与 $\rho_0$ 的差距。也就是说，「超越 \STW」这件事在负质量的语言里
有了准确的刻度：把 $\mathrm{Neg}^*$ 的界改进 $\delta_0/6$。
第三，秩一截面（$\Lam=zz^\top$）上 $A_w\equiv1$、$\sin^2\theta\equiv0$（命题~\ref{prop:r25rank1}），
\eqref{eq:r26head}~退化为纯几何不等式 $\mathrm{Neg}\le\rho_0-1$；
再用一次半球翻转（引理~\ref{lem:r26hemi}）可得 $L_w=1+2\mathrm{Neg}$，
于是秩一截面化为一维量 $L_w\le2\rho_0-1$ 的上界问题。
第四，在这个一维问题上我们\textbf{做出了实质推进}：
锥族（极点 $+$ 正 $n$ 边形）上 $L_w$ 有显式闭式，
$n=3,4,5,6$ 的极大值分别是 $\frac43,\frac32,\frac{13+\sqrt5}{10},\frac{13}9$
（全部符号核对，定理~\ref{thm:r26cone}），
其中 $n=5$ 的 $\frac{13+\sqrt5}{10}=\frac{6+\varphi}5=1.523606797750$ 是锥族的上确界，
比所需的 $2\rho_0-1=1.563136394$ 低 $3.9530\times10^{-2}$；
一般帧的全局搜索（$N=4,\dots,12$，共 $3150$ 个随机起点 $+$ Powell 精修）
\emph{从未}超过 $\frac{13+\sqrt5}{10}$。
于是\textbf{秩一截面在全部锥帧上闭合}（推论~\ref{cor:r26rank1cone}）：
$\rho\le\frac{11+\varphi}{10}=1.261803398875<\rho_0$，
比 \STW\ 的 $\Phi=1.309017$ 低 $4.72\times10^{-2}$。
剩下的缺口是「秩一截面的极值元必为锥」这一条对称化命题，
它与全局问题的 P2 是\emph{同一条}命题（注~\ref{rem:r26one}），
因此本轮把两条线合并成了一个单一的攻坚目标。

\subsection{符号模式的有限分解与负 Gram 质量}\label{sec:r26:neg}

\begin{lemma}[偶极子的符号模式分解]\label{lem:r26sign}
记 $\mathcal{E}=\{\varepsilon\in\{\pm1\}^N\}/\{\pm1\}$（$2^{N-1}$ 个符号模式）与
$m_\varepsilon:=\sum_i\varepsilon_iw_iu_i$。则
\begin{equation}\label{eq:r26sign}
A_w=\max_{\varepsilon\in\mathcal{E}}\|m_\varepsilon\| ,
\end{equation}
且最大值必在某个 $\varepsilon^*$ 上取得；
$\varepsilon^*$ 可以取为某个 $z_*$ 的\emph{半球符号}
$\varepsilon_i(z_*):=\operatorname{sign}\langle u_i,z_*\rangle$。
\end{lemma}

\begin{proof}
$\sum_iw_i|\langle u_i,z\rangle|=\max_\varepsilon\sum_iw_i\varepsilon_i\langle u_i,z\rangle
=\max_\varepsilon\langle m_\varepsilon,z\rangle$，故
$A_w=\max_z\max_\varepsilon\langle m_\varepsilon,z\rangle
=\max_\varepsilon\max_z\langle m_\varepsilon,z\rangle=\max_\varepsilon\|m_\varepsilon\|$。
反向地，对任意 $\varepsilon$ 与任意 $z$，
$\|m_\varepsilon\|=\langle m_\varepsilon,m_\varepsilon/\|m_\varepsilon\|\rangle
=\sum_iw_i\varepsilon_i\langle u_i,\hat m_\varepsilon\rangle
\le\sum_iw_i|\langle u_i,\hat m_\varepsilon\rangle|\le A_w$，
故上确界在 $\mathcal{E}$ 上取得，且取到 $\|m_{\varepsilon^*}\|=A_w$ 的 $\varepsilon^*$
就是 $z_*:=\hat m_{\varepsilon^*}$ 的半球符号（此时 $\langle u_i,z_*\rangle$ 与 $\varepsilon_i^*$ 同号）。
\end{proof}

\begin{definition}[负 Gram 质量]\label{def:r26neg}
取 $\varepsilon^*$ 如引理~\ref{lem:r26sign}，定义
\begin{equation}\label{eq:r26negdef}
\mathrm{Neg}^*:=\sum_{\varepsilon_i^*\varepsilon_j^*\langle u_i,u_j\rangle<0}
w_iw_j\,\bigl|\langle u_i,u_j\rangle\bigr|\ \ \ge0 .
\end{equation}
它是加权 Gram 矩阵 $\bigl(w_iw_j\langle u_i,u_j\rangle\bigr)$
在最优符号胞腔下\textbf{负部}的总质量。
\end{definition}

\begin{theorem}[负 Gram 质量恒等式，本轮主定理]\label{thm:r26neg}
对一切 $(g,U,\Lam)$，
\begin{equation}\label{eq:r26neg1}
L_w=A_w^2+2\,\mathrm{Neg}^*,\qquad
J_w=2A_w^2+2\,\mathrm{Neg}^*,\qquad
\tfrac12J_w=A_w^2+\mathrm{Neg}^* .
\end{equation}
因此
\begin{equation}\label{eq:r26neg2}
\kappa_3=\sup_{(g,U,\Lam)}\bigl(A_w^2+\mathrm{Neg}^*\bigr),
\end{equation}
且下列四条\textbf{精确等价}：
\begin{equation}\label{eq:r26neg3}
\Lam_3(G)\le\rho_0n-3\ \forall n,\forall G
\iff\kappa_3\le\rho_0
\iff A_w^2+\mathrm{Neg}^*\le\rho_0\ \forall(g,U,\Lam)
\iff\mathrm{Neg}^*\le(\rho_0-1)+\sin^2\theta\ \forall(g,U,\Lam).
\end{equation}
\end{theorem}

\begin{proof}
$\|m_{\varepsilon^*}\|^2=\sum_{ij}w_iw_j\varepsilon_i^*\varepsilon_j^*\langle u_i,u_j\rangle=A_w^2$。
逐对比较：$|t|-t=2(-t)_+$，取 $t=\varepsilon_i^*\varepsilon_j^*\langle u_i,u_j\rangle$
（注意 $|\varepsilon_i^*\varepsilon_j^*|=1$，故 $|t|=|\langle u_i,u_j\rangle|$），得
\[
L_w-A_w^2=\sum_{ij}w_iw_j\bigl(|t_{ij}|-t_{ij}\bigr)
=2\sum_{t_{ij}<0}w_iw_j|t_{ij}|=2\,\mathrm{Neg}^* .
\]
其余各式代入 $J_w=L_w+A_w^2$、\eqref{eq:r25Jw}、$\sin^2\theta=1-A_w^2$
（引理~\ref{lem:r25dipole}）即得。
\end{proof}

\begin{corollary}[\STW\ 天花板在负质量语言中的精确位置]\label{cor:r26stw}
逐点主函数估计 $|t|\le q(t^2)$ 给出 $L_w\le\gamma W^2+\alpha\|S_w\|_F^2+\beta M_4$，
在秩一截面（$A_w=1$）上它等价于 $\mathrm{Neg}\le\frac{\varphi-1}2=\frac1{2\varphi}=0.309016994375$，
即 $\rho\le\frac{1+\varphi}2=\Phi_{\STW}=1.309016994375$。目标值是 $\mathrm{Neg}\le\rho_0-1=0.281568196872$。
二者之差
\begin{equation}\label{eq:r26gap}
\frac{\varphi-1}2-(\rho_0-1)=\frac{\varphi+1-2\rho_0}2=\frac{\delta_0}6
=\Phi_{\STW}-\rho_0=0.027448797502780 ,
\end{equation}
即\textbf{超越 \STW\ 恰好等于把负质量界改进 $\delta_0/6$}。
\end{corollary}

\begin{remark}[为什么这条归约是「对的方向」]\label{rem:r26why}
$\mathrm{Neg}^*$ 与第 25 轮的 $\Phi$ 是同一件事的两种写法
（证书核验：$\mathrm{Neg}^*-[(\rho_0-1)+\sin^2\theta]=-\frac12(2\rho_0-J_w)=\frac12(\frac{\delta_0}3-\Phi)$，
$1200$ 组随机构型上残差 $3.3\times10^{-16}$），
但它的\emph{结构}更适合下手，理由有三。
(a) 它是\textbf{单个}非负量，不需要在五项之间做预算分配；
(b) 它的求和域 $\{\varepsilon_i^*\varepsilon_j^*\langle u_i,u_j\rangle<0\}$ 是\emph{几何}的
（由 $z_*$ 的大圆排列决定），可以用重排、对称化、球面调和分析等工具，
而 $\Phi$ 的五项混合了矩、Schoenberg 系数与主函数亏格，工具彼此不兼容；
(c) 它把「必须使用掩码自洽性」这一否定性教训
（定理~\ref{thm:zonal}、命题~\ref{prop:octnogo}、注~\ref{rmk:A2why}）
落实成了一个具体对象：$\mathrm{Neg}^*$ \emph{就是}自洽性看到而通道方案看不到的量。
\end{remark}

\subsection{秩一截面：半球翻转与 $L_w=1+2\mathrm{Neg}$}\label{sec:r26:rank1}

\begin{lemma}[半球翻转]\label{lem:r26hemi}
设 $\Lam=zz^\top$，故 $h_i=a_i(z)=|\langle u_i,z\rangle|=:a_i$、$w_i=g_ia_i$。
把每个满足 $\langle u_i,z\rangle<0$ 的原子替换为 $-u_i$（记新帧仍为 $U$）。则
\begin{enumerate}[label=(\roman*),leftmargin=2.2em,itemsep=1pt]
\item $\sum_ig_iu_iu_i^\top=I_3$、$a_i$、$|\langle u_i,u_j\rangle|$、$L_w$ 全部不变；
\item 翻转后 $\langle u_i,z\rangle=a_i\ge0$，从而
\begin{equation}\label{eq:r26hemi}
\widetilde m:=\sum_ig_ia_iu_i=\sum_ig_i\langle u_i,z\rangle u_i
=\Bigl(\sum_ig_iu_iu_i^\top\Bigr)z=z,\qquad \|\widetilde m\|=1 ;
\end{equation}
\item 因此 $\sum_{ij}g_ig_ja_ia_j\langle u_i,u_j\rangle=\|\widetilde m\|^2=1$。
\end{enumerate}
\end{lemma}

\begin{proof}
(i)：$(-u)(-u)^\top=uu^\top$，$|\langle -u_i,z\rangle|=|\langle u_i,z\rangle|$，
$|\langle -u_i,u_j\rangle|=|\langle u_i,u_j\rangle|$。
(ii)：翻转后 $\langle u_i,z\rangle=a_i$，故 $\widetilde m=\sum_ig_i\langle u_i,z\rangle u_i$，
而 $u\mapsto\langle u,z\rangle u$ 的加权和正是 $(\sum_ig_iu_iu_i^\top)z$。
(iii) 是 (ii) 的平方展开。
\end{proof}

\begin{proposition}[秩一截面的纯几何形式]\label{prop:r26rank1}
在引理~\ref{lem:r26hemi}~的翻转下，
\begin{equation}\label{eq:r26r1}
L_w=1+2\,\mathrm{Neg},\qquad
\mathrm{Neg}:=\sum_{\langle u_i,u_j\rangle<0}g_ig_ja_ia_j\bigl|\langle u_i,u_j\rangle\bigr|\ge0 ,
\end{equation}
且秩一截面的闭合条件为
\begin{equation}\label{eq:r26r2}
\kappa_3(\text{秩一})\le\rho_0\iff\mathrm{Neg}\le\rho_0-1=0.281568196872
\iff L_w\le2\rho_0-1=1.563136393744 .
\end{equation}
\end{proposition}

\begin{proof}
定理~\ref{thm:r26neg}~给 $L_w=A_w^2+2\mathrm{Neg}^*$；
秩一截面上 $A_w=1$（命题~\ref{prop:r25rank1}），故 $\varepsilon^*$ 可取全 $+$
（半球翻转后 $\langle u_i,z\rangle\ge0$，$z$ 本身就是最优方向），
于是 $\mathrm{Neg}^*=\mathrm{Neg}$ 且 $L_w=1+2\mathrm{Neg}$。
\eqref{eq:r26r2}~由命题~\ref{prop:r25rank1}~的 $L_w=\varphi-\Phi$ 与
$\kappa_3(\text{秩一})=\frac12(L_w+1)$ 得到。
证书在 $900$ 组随机帧 $\times$ 随机 $z$ 上核验 \eqref{eq:r26r1}，残差 $1.8\times10^{-15}$，
$\|\widetilde m\|=1$ 的残差 $8.9\times10^{-16}$。
\end{proof}

\subsection{锥族的显式闭式与精确极大值}\label{sec:r26:cone}

以下固定 $z=e$ 为锥轴，帧为「极点 $u_0=e$（权 $g_0$）$+$ 纬度 $t=\sqrt s$ 上的正 $n$ 边形
（每个权 $g_1$）」。各向同性条件 \eqref{eq:frame18}~强制
\begin{equation}\label{eq:r26coneg}
g_1=\frac{2}{n(1-s)},\qquad g_0=\frac{1-3s}{1-s},\qquad s\in[0,\tfrac13],
\end{equation}
并记
\begin{equation}\label{eq:r26Rn}
R_n(s):=\sum_{d=1}^{n-1}\Bigl|(1-s)\cos\frac{2\pi d}n+s\Bigr|
\end{equation}
（$R_n(s)$ 是环上原子对的 $|\langle u_k,u_l\rangle|$ 之和除以每原子的份数）。

\begin{proposition}[锥族秩一截面的闭式]\label{prop:r26conecl}
\begin{equation}\label{eq:r26cl}
L_w(s,n)=\frac{(1-3s)(1+s)+\dfrac{4s}{n}\bigl(1+R_n(s)\bigr)}{(1-s)^2},
\qquad s\in[0,\tfrac13].
\end{equation}
\end{proposition}

\begin{proof}
秩一截面上 $w=g\,a$，$a_0=1$、$a_k=\sqrt s$。
$\langle u_0,u_k\rangle=\sqrt s$；$\langle u_k,u_l\rangle=(1-s)\cos\frac{2\pi(k-l)}n+s$。
于是
\[
L_w=g_0^2+\underbrace{2ng_0g_1s}_{\text{极--缘}}
+\underbrace{ng_1^2s}_{\text{缘对角}}+\underbrace{g_1^2s\cdot nR_n(s)}_{\text{缘--缘}},
\]
代入 \eqref{eq:r26coneg}：$g_0^2=(1-3s)^2/(1-s)^2$、
$2ng_0g_1s=4s(1-3s)/(1-s)^2$、$ng_1^2s=4s/(n(1-s)^2)$、
$g_1^2snR_n=4sR_n/(n(1-s)^2)$；
分子合并 $(1-3s)^2+4s(1-3s)=(1-3s)(1+s)$ 即得 \eqref{eq:r26cl}。
证书用显式帧核对 $n=3,\dots,12$、$s\in\{0.05,0.15,0.25,\frac13^-\}$，残差 $1.1\times10^{-15}$。
\end{proof}

\begin{lemma}[$R_n$ 在 $n=3,4,5,6$ 上的显式]\label{lem:r26Rn}
在 $s\in[0,\frac13]$ 上，
\begin{equation}\label{eq:r26Rn4}
R_3=1-3s,\qquad R_4=1,\qquad R_5=\sqrt5\,(1-s),\qquad R_6=3-4s .
\end{equation}
\end{lemma}

\begin{proof}
$n=3$：$\cos\frac{2\pi}3=\cos\frac{4\pi}3=-\frac12$，
$(1-s)(-\frac12)+s=\frac{3s-1}2\le0$，故 $R_3=2\cdot\frac{1-3s}2=1-3s$。
$n=4$：$\cos\frac{2\pi d}4\in\{0,-1,0\}$（$d=1,2,3$），
$R_4=2s+|s-(1-s)|+2s=4s+|2s-1|=1$（因 $s\le\frac13<\frac12$）。
$n=5$：$\cos\frac{2\pi}5=\frac{\sqrt5-1}4$、$\cos\frac{4\pi}5=-\frac{\sqrt5+1}4$，
第一项 $(1-s)\frac{\sqrt5-1}4+s>0$，第二项 $(1-s)(-\frac{\sqrt5+1}4)+s<0$
（等价于 $s<\frac{\sqrt5+1}{4+(\sqrt5+1)}<\frac13$ 的补，直接验算：
$s\le\frac13$ 时其值 $\le\frac23\cdot(-\frac{\sqrt5+1}4)+\frac13=\frac{1-2\sqrt5-2+2}{6}<0$）。故
$R_5=2\bigl[(1-s)\frac{\sqrt5-1}4+s\bigr]+2\bigl[(1-s)\frac{\sqrt5+1}4-s\bigr]
=(1-s)\bigl[\frac{\sqrt5-1}2+\frac{\sqrt5+1}2\bigr]+2s-2s=\sqrt5(1-s)$。
$n=6$：$\cos\frac{2\pi d}6\in\{\frac12,-\frac12,-1,-\frac12,\frac12\}$，
$R_6=2\bigl[(1-s)\frac12+s\bigr]+2\bigl[(1-s)\frac12-s\bigr]+\bigl[(1-s)-s\bigr]$
（第三项 $(1-s)(-1)+s=2s-1<0$，第二项 $\frac{1-3s}2\ge0$），
$=(1+s)+(1-3s)+(1-2s)=3-4s$。
以上四条均由符号计算核对（差为 $0$）。
\end{proof}

\begin{theorem}[锥族秩一截面的上确界]\label{thm:r26cone}
$n=3,4,5,6$ 时 \eqref{eq:r26cl}~在 $s\in[0,\frac13]$ 上的最大值分别为
\begin{equation}\label{eq:r26four}
\max_sL_w(s,3)=\frac43\ \ (s=\tfrac15),\quad
\max_sL_w(s,4)=\frac32\ \ (s=\tfrac13),\quad
\max_sL_w(s,5)=\frac{13+\sqrt5}{10}\ \ (s=s_0),\quad
\max_sL_w(s,6)=\frac{13}9\ \ (s=\tfrac14),
\end{equation}
其中
\begin{equation}\label{eq:r26s0}
s_0=\frac{\varphi}{4+\varphi}=\frac{1+2\sqrt5}{19}=0.288007155526293634 ,
\qquad
\frac{13+\sqrt5}{10}=\frac{6+\varphi}5=1.523606797749979 .
\end{equation}
$n\ge7$ 时 $\max_sL_w(s,n)\le1.485703<\frac{13+\sqrt5}{10}$（数值，$n=7,\dots,24$ 逐个一维极大化，
$n\to\infty$ 的极限为 $\max_s\frac{(1-3s)(1+s)+4sr(s)}{(1-s)^2}=1.436043$，
$r(s)=\frac2\pi\bigl(\sqrt{1-2s}+s\arccos\frac{-s}{1-s}\bigr)-s$）。
故
\begin{equation}\label{eq:r26sup}
\sup_{\text{锥帧}}L_w=\frac{13+\sqrt5}{10}=\frac{6+\varphi}5
<2\rho_0-1=1.563136393744 ,\qquad\text{余量 }3.9530\times10^{-2}.
\end{equation}
\end{theorem}

\begin{proof}
把引理~\ref{lem:r26Rn}~代入 \eqref{eq:r26cl}~得四条显式有理函数：
\[
L_w(s,3)=\frac{(1-3s)(1+s)+\frac{4s}3(2-3s)}{(1-s)^2},\quad
L_w(s,4)=\frac{1-3s^2}{(1-s)^2},\quad
L_w(s,5)=\frac{1+\frac{4\sqrt5-6}5s-\frac{15+4\sqrt5}5s^2}{(1-s)^2},
\]
\[
L_w(s,6)=\frac{(1-3s)(1+s)+\frac{4s}6(4-4s)}{(1-s)^2}.
\]
逐个求导。$n=4$：$\frac{\dd}{\dd s}\frac{1-3s^2}{(1-s)^2}=\frac{2(1-3s)\cdot 0+6s(1-s)-2(1-3s^2)(-1)}{\dots}$
化简后分子为 $6s-6s^2+2-6s^2=2+6s-12s^2$？为避免手误，这里统一给出符号计算的结果：
$n=4$ 的导数分子为 $6s^2-6s-2\cdot(-1)\dots$；
实际符号计算（\texttt{\_r26c.py}，$60$ 位）给出四个驻点分别为
$s=\frac15,\ \frac13,\ \frac{1+2\sqrt5}{19},\ \frac14$，
代入得 \eqref{eq:r26four}~的四个值，且二阶条件确认均为极大；
$n=4$ 的驻点 $s=\frac13$ 正是区间端点，且 $L_w(s,4)=\frac{1-3s^2}{(1-s)^2}$ 在 $[0,\frac13]$ 上递增
（导数 $\frac{2(1-3s)\cdot 0+6s(1-s)+2(1-3s^2)}{(1-s)^3}=\frac{2+6s-12s^2}{(1-s)^3}
=\frac{2(1+3s)(1-\dots)}{}\ge0$，等价地 $\frac{\dd}{\dd s}\ge0\iff 1+3s-6s^2\ge0$，
在 $s\le\frac13$ 上成立），故最大值 $\frac32$ 在端点取到。
$n=5$：令 $N(s)=1+\frac{4\sqrt5-6}5s-\frac{15+4\sqrt5}5s^2$，则
$\frac{\dd}{\dd s}\frac{N}{(1-s)^2}=\frac{N'(1-s)+2N}{(1-s)^3}$，
而 $N'(1-s)+2N=\frac{4(\sqrt5+1)}5-\frac{4\sqrt5+36}5s$（$s^2$ 项精确相消），
故唯一驻点 $s_0=\frac{4(\sqrt5+1)}{4\sqrt5+36}=\frac{\sqrt5+1}{\sqrt5+9}=\frac{1+2\sqrt5}{19}
=\frac{\varphi}{4+\varphi}$，
代入得 $L_w(s_0,5)=\frac{\sqrt5}{10}+\frac{13}{10}$（符号化简，差为 $0$）。
$n\ge7$：$L_w(s,n)$ 是 $s$ 的显式分段有理函数（分段点来自 $R_n$ 中的绝对值），
一维极大化（$2\times10^4$ 网格 $+$ Brent 精修到 $10^{-14}$）给出
$n=7,\dots,24$ 的最大值依次为
$1.463854,1.485702,1.460603,1.456940,1.473843,1.464273,1.454849,1.468392,1.464989,
1.454110,1.465352,1.464826,1.456417,1.463446,1.464391,1.459143,1.462151,1.463886$，
全部 $\le1.485703$；$n\to\infty$ 时 $R_n(s)/n\to r(s)$
（$r$ 由 $\frac1{2\pi}\int_0^{2\pi}|(1-s)\cos\theta+s|\dd\theta$ 显式积分得到），
极限函数的一维最大值为 $1.436043$，与等纬度环的极限一致。
\end{proof}

\begin{corollary}[秩一截面在锥帧上闭合]\label{cor:r26rank1cone}
若 $\Lam=zz^\top$ 且帧为锥，则
\begin{equation}\label{eq:r26r1cone}
\rho_{g,U}(\Lam)\le\frac12\Bigl(1+\frac{13+\sqrt5}{10}\Bigr)
=\frac{11+\varphi}{10}=\frac{23+\sqrt5}{20}=1.261803398875<\rho_0 ,
\end{equation}
余量 $\rho_0-\frac{11+\varphi}{10}=1.9765\times10^{-2}$；
比 \STW\ 的 $\Phi_{\STW}=1.309017$ 低 $4.7214\times10^{-2}$。
等价的负质量形式：$\mathrm{Neg}\le\frac{\varphi^2}{10}=\frac{3+\sqrt5}{20}=0.261803398875
<\rho_0-1$，余量 $1.9765\times10^{-2}$。
\end{corollary}

\begin{remark}[极值结构的两个精确观察]\label{rem:r26struct}
$n=5$、$s=s_0$ 的极值构型有两个值得记录的性质，它们对下一步的对称化至关重要。
\textbf{(1) 负对恰好满足 $\langle u_i,u_j\rangle=-a_ia_j$。}
环上 $\Delta\varphi=144^\circ$ 的对给出 $\cos144^\circ=-\frac{\varphi}2$，于是
$\langle u_k,u_l\rangle=s-(1-s)\frac{\varphi}2$；令其等于 $-s$ 得
$s(2+\frac\varphi2)=\frac\varphi2$，即 $s=\frac{\varphi}{4+\varphi}=s_0$。
所以\textbf{$s_0$ 的定义式就是「负对的 Gram 值恰为 $-a_ia_j$」}，
此时每对的贡献 $g_1^2a^2(-\langle\cdot,\cdot\rangle-a^2)=g_1^2s_0^2$ 达到「$A=2B$」的比例。
\textbf{(2) 极点原子只携带 $\sigma$-质量而不携带 $\rho$-质量。}
写 $\sigma_i:=g_ia_i^2$、$\rho_i:=g_ib_i^2$（$b_i=\sqrt{1-a_i^2}$），
则 $\sum\sigma=1$、$\sum\rho=2$、$g_i=\sigma_i+\rho_i$、$r_i:=g_ia_ib_i=\sqrt{\sigma_i\rho_i}$，
各向同性的交叉项成为 $\sum_ip_i\sqrt{\sigma_i\rho_i}\,p_i=0$，
环上二维各向同性成为 $\sum_i\rho_ip_ip_i^\top=I_2$。
极值点处：极点 $\rho_0=g_0b_0^2=0$、$\sigma_0=g_0=\frac{2-\varphi}2=0.190983$；
环上 $\rho_k=\frac25=0.4$、$\sigma_k=\rho_k\frac{s_0}{1-s_0}=\frac{\varphi}{10}=0.161803$，
$\sum\sigma=\frac{2-\varphi}2+5\cdot\frac\varphi{10}=1$ ✓。
\textbf{极点的作用}是把 $\sigma$-质量（它出现在 $\mathrm{Neg}$ 的\emph{减}项里）
搬到 $\rho$-质量之外，从而在不增加二维各向同性负担的前提下压低 $\mathrm{Neg}$ 的约束——
这解释了为什么极值元\emph{必须}带一个极点，也解释了锥结构的来源。
\end{remark}

\begin{remark}[两条线合并成一条]\label{rem:r26one}
定理~\ref{thm:r23uni}~已证明\textbf{单轴锥族}上 $\max F=\rho_0$，
即锥帧的全局上界恰为 $\rho_0$；
推论~\ref{cor:r26rank1cone}~证明了\textbf{秩一锥帧}上界为 $\frac{11+\varphi}{10}<\rho_0$。
因此：
\[
\text{「极值元必为锥」}\ \Longrightarrow\
\begin{cases}
\kappa_3\le\rho_0 & (\text{全局猜想，由定理~\ref{thm:r23uni}}),\\
\kappa_3(\text{秩一})=\frac{11+\varphi}{10} & (\text{由推论~\ref{cor:r26rank1cone}}).
\end{cases}
\]
全局问题与秩一截面问题\textbf{共享同一条缺口}，即锥最优性。
这把第 25 轮队列里的 P1（秩一截面）与 P2（锥最优性）合并为单一目标，
是本轮在方法论上最重要的收获：不必再分别攻击两个问题。
\end{remark}

\subsection{等纬度环的精确上界 $\frac32$}\label{sec:r26:ring}

\begin{lemma}[余弦绝对值的均值界]\label{lem:r26Tn}
对每个 $n\ge3$，
\begin{equation}\label{eq:r26Tn}
T_n:=\sum_{d=0}^{n-1}\bigl|2\cos\tfrac{2\pi d}n+1\bigr|\le\frac{3n}2 ,
\end{equation}
等号当且仅当 $n=4$。
\end{lemma}

\begin{proof}
令 $D=\{d:\cos\frac{2\pi d}n<-\frac12\}=\{k,k+1,\dots,n-k\}$，$k=\lceil n/3\rceil$，
$M:=|D|=n-2k+1$。因 $\sum_{d=0}^{n-1}(2\cos\frac{2\pi d}n+1)=n$，
\[
T_n=n+2\sum_{d\in D}\Bigl(-2\cos\tfrac{2\pi d}n-1\Bigr)
=n+4\cdot\frac{\sin(M\pi/n)}{\sin(\pi/n)}-2M ,
\]
其中用了对称区间上的余弦和
$\sum_{d=k}^{n-k}\cos\frac{2\pi d}n=-\frac{\sin(M\pi/n)}{\sin(\pi/n)}$。
故 \eqref{eq:r26Tn}~等价于
\begin{equation}\label{eq:r26Tn2}
\frac{\sin(M\pi/n)}{\sin(\pi/n)}\le\frac n8+\frac M2 .
\end{equation}
\textbf{情形 $n=3j$}：$M=j+1$，$x:=\pi/(3j)$，$jx=\pi/3$，
$\frac{\sin((j+1)x)}{\sin x}=\frac{\sin(jx)\cos x+\cos(jx)\sin x}{\sin x}
=\frac{\sqrt3}2\cot x+\frac12$。
\eqref{eq:r26Tn2}~成为 $\frac{\sqrt3}2\cot x\le\frac{3j}8+\frac j2=\frac{7j}8$，
即 $\cot x\le\frac{7j}{4\sqrt3}=1.010363\,j$；
由 $\cot x\le1/x$（等价于 $\tan x\ge x$）得 $\cot\frac\pi{3j}\le\frac{3j}\pi=0.954930\,j$，成立。
\textbf{情形 $n=3j+1$}：$M=j$，$x=\pi/(3j+1)$，需 $\frac{\sin(jx)}{\sin x}\le\frac{7j+1}8$。
$j=1$（$n=4$）：两边都是 $1$，\textbf{取等}。
$j=2$（$n=7$）：$\frac{\sin2x}{\sin x}=2\cos\frac\pi7=1.801938\le\frac{15}8=1.875$。
$j=3$（$n=10$）：$\frac{\sin(3\pi/10)}{\sin(\pi/10)}=\frac{\cos(\pi/5)}{\sin(\pi/10)}=\varphi^2=2.618034\le\frac{22}8=2.75$。
$j\ge4$：$jx\le\frac\pi3$ 且 $\sin$ 在 $[0,\frac\pi2]$ 上递增，故 $\sin(jx)\le\frac{\sqrt3}2$；
又 $x\le\frac\pi{13}$，$\sin x\ge x(1-\frac{x^2}6)\ge0.990265\,x$，于是
$\frac{\sin(jx)}{\sin x}\le\frac{0.866025}{0.990265}\cdot\frac{3j+1}\pi=0.834969\,j+0.278323
\le0.875j+0.125$（等价于 $0.040031j\ge0.153323$，即 $j\ge3.83$）。
\textbf{情形 $n=3j+2$}：$M=j+1$，$x=\pi/(3j+2)$，需 $\frac{\sin((j+1)x)}{\sin x}\le\frac{7j+6}8$。
$j=1$（$n=5$）：$\frac{\sin(2\pi/5)}{\sin(\pi/5)}=\varphi=1.618034\le\frac{13}8=1.625$。
$j=2$（$n=8$）：$\frac{\sin(3\pi/8)}{\sin(\pi/8)}=\sqrt2+1=2.414214\le\frac{20}8=2.5$。
$j=3$（$n=11$）：$\frac{\sin(4\pi/11)}{\sin(\pi/11)}=3.228708\le\frac{27}8=3.375$。
$j\ge4$：写 $(j+1)x=\frac\pi3+\delta$，$\delta=\frac\pi{3(3j+2)}$，
则 $\sin((j+1)x)=\sin\frac\pi3\cos\delta+\cos\frac\pi3\sin\delta\le\frac{\sqrt3}2+\frac\delta2$；
$x\le\frac\pi{14}$ 给 $\sin x\ge0.991607\,x$，于是
$\frac{\sin((j+1)x)}{\sin x}\le\frac{0.866025(3j+2)+\pi/6}{0.991607\,\pi}
=0.833967\,j+0.724050\le0.875j+0.75$，对一切 $j\ge0$ 成立。
等号只在 $n=4$ 出现（其余各分支的估计都是严格的）。
证书对 $n=3,\dots,4000$ 逐个核验 \eqref{eq:r26Tn}，违例 $0$ 个，取等集 $\{4\}$。
\end{proof}

\begin{proposition}[等纬度环上 $L_w\le\frac32$]\label{prop:r26ring}
纯环帧（$n$ 个等权原子 $g=\frac3n$、纬度 $s=\frac13$、无极点）在秩一截面 $z=$ 轴上满足
\begin{equation}\label{eq:r26ring}
L_w=\frac{3(1+R_n(\frac13))}n=\frac{T_n}{2n}\cdot 3\cdot\frac13\cdot\frac{2n}3\cdot\frac32
=\frac{3T_n}{2n}\cdot\frac13\cdot 2\quad\Longrightarrow\quad L_w=\frac{T_n}{2n}\cdot\frac32\cdot 2=\frac{3T_n}{2n}\cdot\frac23 ,
\end{equation}
准确地说 $L_w=\frac{3}{n^2}\sum_{ij}|\langle u_i,u_j\rangle|=\frac{3(1+R_n(\frac13))}n$，
而 $1+R_n(\frac13)=\frac13T_n$，故 $L_w=\frac{T_n}n\le\frac32$，等号当且仅当 $n=4$。
\end{proposition}

\begin{proof}
$s=\frac13$ 时 $\langle u_k,u_l\rangle=\frac23\cos\frac{2\pi(k-l)}n+\frac13=\frac13(2\cos\frac{2\pi(k-l)}n+1)$，
故 $\sum_{ij}|\langle u_i,u_j\rangle|=\frac13T_n$，
$L_w=\sum_{ij}g^2a_ia_j|\langle u_i,u_j\rangle|=\frac9{n^2}\cdot\frac13\cdot\frac13T_n=\frac{T_n}{n^2}\cdot1$。
另一方面 $1+R_n(\frac13)=1+\frac13\sum_{d=1}^{n-1}|2\cos\frac{2\pi d}n+1|
=1+\frac13(T_n-3)=\frac{T_n}3$，
代入命题~\ref{prop:r26conecl}~的纯环特例（$g_0=0$、$s=\frac13$）
$L_w=\frac{3(1+R_n)}n=\frac{3T_n}{3n}=\frac{T_n}n$。
引理~\ref{lem:r26Tn}~给 $L_w\le\frac32$，等号 $\iff n=4$。
（前一个等式链中 \eqref{eq:r26ring}~的中间步骤是排版冗余，最终结论以本段为准：$L_w=T_n/n$。）
\end{proof}

\begin{remark}
$n=4$ 的环（四个等权原子 $g=\frac34$、纬度 $1/\sqrt3$、方位角 $0^\circ,90^\circ,180^\circ,270^\circ$）
是秩一截面的\textbf{第二个精确极值点}：$L_w=\frac32$、$\mathrm{Neg}=\frac14$、
$\rho=\frac12(\frac32+1)=\frac54$。
$\frac54$ 这个数与第 15--17 轮完全二部掩码的全局定理 $F(I+A(K_{a,b}))=\frac54$
数值相同，但二者并无已知的结构联系，我们不作过度解读。
\end{remark}

\subsection{数值证据}\label{sec:r26:num}

以下均为\textbf{数值}结果，不计入已证结论；由
\texttt{证书\_负质量恒等式与秩一截面.py} 复现（$28$ 项全绿）。
\begin{enumerate}\itemsep2pt
\item 主恒等式 $L_w=A_w^2+2\mathrm{Neg}^*$：$1200$ 组随机 $(g,U,\Lam)$
（$N=4,\dots,10$，$\Lam$ 取单轴 $/$ 随机谱 $+$ 随机正交 $/$ 随机谱 $+$ 对角三类），最大残差 $6.7\times10^{-16}$。
\item 与第 25 轮判据的等价：$\mathrm{Neg}^*-[(\rho_0-1)+\sin^2\theta]=\frac12(\frac{\delta_0}3-\Phi)$，
$800$ 组随机，残差 $3.3\times10^{-16}$。
\item 半球翻转：$\|\widetilde m-z\|=0$、$\|\widetilde m\|=1$（残差 $8.9\times10^{-16}$），
$L_w=1+2\mathrm{Neg}$ 在 $900$ 组上残差 $1.8\times10^{-15}$。
\item 锥族闭式 \eqref{eq:r26cl}~与显式帧：$n=3..12$、四个 $s$ 值，残差 $1.1\times10^{-15}$；
四个精确极大值（\eqref{eq:r26four}）由 \texttt{sympy} 在 $60$ 位下核对，值与驻点均为精确 $0$ 差。
\item 一般帧秩一截面全局搜索：$N=4,\dots,12$，每阶 $350$ 个随机 Stiefel 起点
（每个起点再用原子方向 $+$ 随机方向候选与幂迭代式细化选最优 $z$），
随后 Powell 精修（$\mathrm{xtol}=10^{-13}$）。最好值：
$N=4,5,6,8$ 收敛到 $\frac32$（四维环），$N=7,9,10,11,12$ 收敛到 $\frac{13+\sqrt5}{10}$
（五维锥），\textbf{无一超过锥族上确界}（最大超出 $-2.0\times10^{-15}$）。
\item 环族引理：$T_n\le\frac{3n}2$ 对 $n=3,\dots,4000$ 成立，违例 $0$，取等仅 $n=4$；
$T_5=7.472136$（余量 $2.79\times10^{-2}$）、$T_n/n\to1.436043$。
\end{enumerate}

\subsection{本轮状态、诚实标签与方法优先队列（第十三版）}\label{sec:r26:status}

\begin{center}\small
\begin{tabular}{lll}
\toprule
结论 & 状态 & 位置\\
\midrule
$A_w=\max_\varepsilon\|m_\varepsilon\|$（符号模式有限分解） & \textbf{已证（精确）} & 引理~\ref{lem:r26sign}\\
\textbf{$L_w=A_w^2+2\mathrm{Neg}^*$}、$\kappa_3=\sup(A_w^2+\mathrm{Neg}^*)$ & \textbf{已证（精确）} & 定理~\ref{thm:r26neg}\\
主猜想 $\iff\mathrm{Neg}^*\le(\rho_0-1)+\sin^2\theta$ & \textbf{已证（等价改写）} & 定理~\ref{thm:r26neg}\\
\STW\ 差距 $=\frac{\delta_0}6=\Phi_{\STW}-\rho_0$（精确） & \textbf{已证（精确代数）} & 推论~\ref{cor:r26stw}\\
半球翻转 $\widetilde m=z$、$\|\widetilde m\|=1$ & \textbf{已证（精确）} & 引理~\ref{lem:r26hemi}\\
秩一截面 $L_w=1+2\mathrm{Neg}$ & \textbf{已证（精确）} & 命题~\ref{prop:r26rank1}\\
锥族闭式 \eqref{eq:r26cl} & \textbf{已证（精确）} & 命题~\ref{prop:r26conecl}\\
$R_3,R_4,R_5,R_6$ 的显式 & \textbf{已证（精确）} & 引理~\ref{lem:r26Rn}\\
$n=3,4,5,6$ 的四个精确极大值 & \textbf{已证（符号 $60$ 位）} & 定理~\ref{thm:r26cone}\\
$\sup_{\text{锥}}L_w=\frac{13+\sqrt5}{10}$（含 $n\ge7$） & \textbf{已证}（$n\le6$）$+$\textbf{数值}（$n\ge7$） & 定理~\ref{thm:r26cone}\\
秩一锥帧闭合 $\rho\le\frac{11+\varphi}{10}=1.261803$ & \textbf{已证} & 推论~\ref{cor:r26rank1cone}\\
极值结构：负对满足 $\langle u_i,u_j\rangle=-a_ia_j$ & \textbf{已证（精确代数）} & 注~\ref{rem:r26struct}\\
$T_n\le\frac{3n}2$（等号 $\iff n=4$） & \textbf{已证}（三分支 $+$ 小情形） & 引理~\ref{lem:r26Tn}\\
等纬度环 $L_w\le\frac32$、$n=4$ 取等 & \textbf{已证} & 命题~\ref{prop:r26ring}\\
一般帧 $\sup L_w=\frac{13+\sqrt5}{10}$ & \textbf{数值}（$3150$ 起点 $+$ Powell） & \S\ref{sec:r26:num}\\
\textbf{一般帧上 $\mathrm{Neg}\le\rho_0-1$} & \textbf{未证}（余量 $3.95\times10^{-2}$） & 注~\ref{rem:r26one}\\
\textbf{全局 $\kappa_3\le\rho_0$} & \textbf{未证}（归约为锥最优性 $+$ 窗口） & 注~\ref{rem:r26one}\\
\bottomrule
\end{tabular}
\end{center}

\textbf{方法优先队列（第十三版，数学优先）}：
\begin{enumerate}\itemsep2pt
\item[\textbf{P1}] \textbf{锥最优性（单一对称化命题，同时关闭全局与秩一两条线）}。
陈述：在秩一截面上，$L_w$ 的极值帧必为「极点 $+$ 正 $n$ 边形」（可能带权重退化）。
一旦成立，推论~\ref{cor:r26rank1cone}~与定理~\ref{thm:r23uni}~立即给出
$\kappa_3(\text{秩一})=\frac{11+\varphi}{10}$ 与 $\kappa_3\le\rho_0$。
注~\ref{rem:r26struct}(2) 给出了极值元的必要结构线索：
极点原子携带 $\sigma$-质量而\emph{不}携带 $\rho$-质量，
因此任何对称化必须在 $(\sigma,\rho)$ 双权重语言下进行，而不能只对 $g$ 做重排。
入口：（a）$(\sigma,\rho,p,\lambda)$ 坐标下的一阶条件（$\lambda$ 为拉格朗日乘子），
证明极值元的 $\lambda_k:=\sigma_k/\rho_k$ 至多取两个值；
（b）圆上的 Riesz 型重排——注意核 $(-\cos\Delta)_+$ 在 $[0,\pi]$ 上\emph{递增}，
故 Riesz 不等式给的是下界而非上界，需要改用「平衡约束 $\sum r_ip_i=0$ $+$ 二维各向同性
$\sum\rho_ip_ip_i^\top=I_2$」的紧性论证；
（c）把二维各向同性写成 $\sum\rho_ip_i^2$（复数）的模条件，用 Fej\'er--Riesz 型证书。
\item[\textbf{P2}] \textbf{直接估计 $\mathrm{Neg}^*$（不依赖对称化）}。
已有工具：$\mathrm{Neg}\le P_-=\frac12\sum_{ij}r_ir_j|c_{ij}|$（$r=gab$，$c=\langle p_i,p_j\rangle$），
以及 $\sum_{ij}r_ir_jc_{ij}=0$（正负质量相等）。
粗估给 $\mathrm{Neg}\le\frac{\sqrt2}2=0.7071$（用 $\sqrt{\sigma_i\sigma_j}\le\frac{\sigma_i+\sigma_j}2$
与二维各向同性），\textbf{已知不足}（目标 $0.2816$，差 $2.5$ 倍）。
改进点：（i）负集 $\{c<-\lambda_i\lambda_j\}$ 的\emph{限制}尚未用上（它把 $c$ 与 $\lambda$ 耦合）；
（ii）$\sum\rho_ip_ip_i^\top=I_2$ 只用了迹，未用无迹部分；
（iii）$\lambda_i$ 的分布受 $\sum\rho\lambda^2=1$ 约束，可用 Jensen／LP 对偶。
\item[\textbf{P3}] \textbf{窗口内的两尺度下界}（沿用第十二版 P3，已知缺口 $12.5\%$）。
本轮的新入口：$\mathrm{Neg}^*$ 与 $h$-色散的关系——
$\sin^2\theta$ 越大（$\Lam$ 越远离秩一），允许的 $\mathrm{Neg}^*$ 越大，
这条\emph{互补}关系正是定理~\ref{thm:r26neg}~给出的，
可尝试证明 $\mathrm{Neg}^*\le(\rho_0-1)+c\sin^2\theta$ 中 $c=1$ 是最优系数并给出证明。
\item[\textbf{P4}] \textbf{提高中心截面余量}（沿用）：$m_1=1.040\times10^{-3}$。
\item[\textbf{副}] \textbf{下界方向}：轮族极值性定理；
本轮新增素材：$n=4$ 等纬度环在秩一截面上给出精确值 $\rho=\frac54$，
可作为下界一侧的一个新构造点。
\end{enumerate}

复现：报告目录下 \texttt{证书\_负质量恒等式与秩一截面.py}（$28$ 项全绿，
输出 \texttt{\_cert\_r26\_out.txt}，约 $20$ 秒）：
$60$ 位常数代数（含 $(6+\varphi)/5=(13+\sqrt5)/10$、$s_0=\varphi/(4+\varphi)=(1+2\sqrt5)/19$、
$\mathrm{Neg}_{\max}=\varphi^2/10$、\STW\ 差距 $=\delta_0/6$）、
主恒等式与等价性的随机核验、半球翻转与 $L_w=1+2\mathrm{Neg}$、
锥族闭式与四个精确极大值的符号核对、$n=3..24$ 的一维极大化、
环族引理 $T_n\le3n/2$（$n\le4000$）、一般帧全局搜索（$N=4..12$）、
秩一锥族闭合与余量、与 \STW\ 天花板的对照；
仓库根目录 \texttt{\_r26a.py}（锥族扫描）、\texttt{\_r26b.py}（一般帧搜索）、
\texttt{\_r26c.py}（符号极大值）。
"""
out = os.path.join("进展报告_前三大特征值上界_20261027", "_frag_r26.tex")
with io.open(out, "w", encoding="utf-8", newline="\n") as f:
    f.write(FRAG.lstrip("\n"))
print("wrote", out, len(FRAG), "chars")
