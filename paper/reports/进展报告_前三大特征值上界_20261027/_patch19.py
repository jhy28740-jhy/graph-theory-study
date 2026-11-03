# -*- coding: utf-8 -*-
import io
P = '进展报告_前三大特征值之和的精确上界.tex'
s = io.open(P, encoding='utf-8').read()
orig_len = len(s)
frag = io.open('_frag_r19.tex', encoding='utf-8').read()

# ---- 0. self-contained alpha/beta inside the FRAGMENT
a = r"""记 $Q=VV^\top$、$Q_\pm=(\pm Q)_+$、$\alpha=\|V^\top x\|^2$、$\beta=x^\top Q_-x$，
则 $\rho=\alpha+\beta$，并且只用下面四条数据"""
b = r"""记 $Q=VV^\top$、$Q_\pm=(\pm Q)_+$（故 $Q_+=Q+Q_-$），
$\alpha:=x^\top Qx=\|V^\top x\|^2$、$\beta:=x^\top Q_-x\ge0$，
则 $\rho=x^\top Q_+x=\alpha+\beta$。只用下面四条数据"""
assert a in frag, 'frag alpha/beta anchor missing'
frag = frag.replace(a, b)
a = r"$\rho=\alpha+\beta$ 是第~\ref{sec:r18}~节的分解。由 \eqref{eq:bootdata}，"
b = r"$\rho=\alpha+\beta$ 来自 $Q_+=Q+Q_-$ 的逐项分解。由 \eqref{eq:bootdata}，"
assert a in frag
frag = frag.replace(a, b)

# ---- 1. prop:concave -> superseded pointer
a = r"""\begin{proposition}[凹性，数值]\label{prop:concave}
对固定几何，$\rho(\Lam)$ 在 $\Lam$ 上是\textbf{凹}的：在冠军几何上，
$600$ 组随机 $(\Lam_1,\Lam_2,s)$ 的最小凹性余量
$\rho((1-s)\Lam_1+s\Lam_2)-[(1-s)\rho(\Lam_1)+s\rho(\Lam_2)]$ 为 $+1.99\times10^{-5}$；
在 $10$ 个随机几何（$N=5,\dots,11$）上全局最小余量为 $+1.57\times10^{-5}$，无一为负。
\end{proposition}"""
b = r"""\begin{proposition}[凹性，数值；\textbf{已在第~\ref{sec:r19}~节升级为定理~\ref{thm:geo}}]\label{prop:concave}
对固定几何，$\rho(\Lam)$ 在 $\Lam$ 上是\textbf{凹}的：在冠军几何上，
$600$ 组随机 $(\Lam_1,\Lam_2,s)$ 的最小凹性余量
$\rho((1-s)\Lam_1+s\Lam_2)-[(1-s)\rho(\Lam_1)+s\rho(\Lam_2)]$ 为 $+1.99\times10^{-5}$；
在 $10$ 个随机几何（$N=5,\dots,11$）上全局最小余量为 $+1.57\times10^{-5}$，无一为负。
\end{proposition}

\begin{remark}[第 19 轮回填]\label{rmk:concavefill}
本命题现在是\textbf{定理}。写 $L_i=g_iu_i^\top\Lam u_i=x_i(\Lam)^2$，则
$\rho(\Lam)=\sum_{ij}(Q_+)_{ij}\sqrt{L_iL_j}$，而 $\sqrt{ab}=\min_{t>0}\frac{ta+b/t}2$
是仿射族的下确界，故凹性由「线性 $\to$ 几何平均 $\to$ 非负组合」三步解析得到，
见定理~\ref{thm:geo}。同一表示还给出闭式梯度（推论~\ref{cor:gradpsi}）、
精确对偶（定理~\ref{thm:exactdual}）与内层规划的全局极大点刻画（推论~\ref{cor:argmax}）。
\end{remark}"""
assert a in s, 'prop:concave anchor missing'
s = s.replace(a, b)

# ---- 2. rmk:concave
a = r"""\begin{remark}\label{rmk:concave}
凹性尚未解析证明。注意它\textbf{不能}由 $\nabla^2\rho$ 的符号直接读出：
\eqref{eq:gradlam}~只在 Perron 点成立，不能对 $\Lam$ 再求导。
若凹性成立，则 $\rho(\Lam)$ 作为一次正齐次凹函数具有支撑超平面表示"""
b = r"""\begin{remark}\label{rmk:concave}
（第 19 轮更新：凹性\textbf{已}解析证明，见定理~\ref{thm:geo}。下面保留当时的判断；
其中「不能由 $\nabla^2\rho$ 读出」仍然正确，正确做法是把 $\rho$ 写成
线性函数的几何平均的非负组合。）
注意凹性\textbf{不能}由 $\nabla^2\rho$ 的符号直接读出：
\eqref{eq:gradlam}~只在 Perron 点成立，不能对 $\Lam$ 再求导。
$\rho(\Lam)$ 作为一次正齐次凹函数具有支撑超平面表示"""
assert a in s, 'rmk:concave anchor missing'
s = s.replace(a, b)

# ---- 3. cor:dualcert unconditional
a = r"""\begin{corollary}[对偶证书形式]\label{cor:dualcert}
设 \eqref{eq:dualcert}~成立。则要证某几何满足 $\max_\Lam\rho(\Lam)\le\rho_0$，"""
b = r"""\begin{corollary}[对偶证书形式]\label{cor:dualcert}
（第 19 轮更新：\eqref{eq:dualcert}~已由定理~\ref{thm:exactdual}~无条件证明，
且其中的 $\min$ 可限制在显式对偶族 $\{S(T):T>0\}$ 上，最小值在
$t^*_{ij}=x_j/x_i$ 处取到，此时 $S(T^*)=\rho I_3$。本推论的结论因此无条件成立；
但它同时被命题~\ref{prop:P1dead}~判为\textbf{无用}——冠军处可行证书只有 $\rho_0I_3$ 一个。）
则要证某几何满足 $\max_\Lam\rho(\Lam)\le\rho_0$，"""
assert a in s, 'cor:dualcert anchor missing'
s = s.replace(a, b)

# ---- 4. sec:r18:status
a = r"\textbf{方法优先队列（第八版）}"
b = (r"\textbf{方法优先队列（第八版；P2 已由定理~\ref{thm:geo}~完成，"
     r"P1 已由命题~\ref{prop:P1dead}~否证并弹出，现行队列为第~\ref{sec:r19:status}~节的第九版）}")
assert a in s
s = s.replace(a, b)

a = r"""\textbf{条件性}：定理~\ref{thm:lamparam}(iii)~的「$\ge$」一侧对退化情形（$g_i=0$、重根）
的化归；推论~\ref{cor:dualcert}~依赖命题~\ref{prop:concave}~的凹性。"""
b = r"""\textbf{条件性}：定理~\ref{thm:lamparam}(iii)~的「$\ge$」一侧对退化情形（$g_i=0$、重根）
的化归；推论~\ref{cor:dualcert}~依赖命题~\ref{prop:concave}~的凹性。
\emph{（第 19 轮更新：凹性已解析证明，故推论~\ref{cor:dualcert}~无条件；
「$\ge$」一侧在 $x>0$ 时亦无条件，见推论~\ref{cor:argmax}。）}"""
assert a in s
s = s.replace(a, b)

# ---- 5. insert Round 19
anchor = r"\section{通道型主函数方案的精确天花板}\label{sec:nogo}"
assert anchor in s, 'sec:nogo anchor missing'
s = s.replace(anchor, frag + anchor, 1)

io.open(P, 'w', encoding='utf-8', newline='\n').write(s)
labs = ['thm:geo','lem:geo','cor:gradpsi','thm:exactdual','cor:argmax','lem:lampos',
        'prop:P1dead','prop:ceilweyl','prop:ceilharm','thm:section','cor:sectspecial',
        'prop:sectmap','cor:sectreduce','sec:r19','sec:r19:status','rmk:concavefill']
print('tex patched: %d -> %d chars (+%d)' % (orig_len, len(s), len(s) - orig_len))
print('missing new labels:', [x for x in labs if ('\\label{'+x+'}') not in s])
print('amsthm addendum not used:', 'addendum' not in s)
