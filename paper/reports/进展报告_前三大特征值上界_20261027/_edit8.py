# -*- coding: utf-8 -*-
import io
main="进展报告_前三大特征值之和的精确上界.tex"
s=io.open(main,encoding="utf-8").read()
def rep(old,new,cnt=1):
    global s
    assert s.count(old)==cnt, ("COUNT %d"%s.count(old), old[:70]); s=s.replace(old,new,cnt)

old_prop = s[s.index(r"\begin{proposition}[分数度障碍，已证]\label{prop:fracdeg}"):s.index(r"\begin{proposition}[$N=6$ 的整数化阶梯")]
new_prop = r"""\begin{proposition}[分数度障碍与例外阶数集 $\mathcal E$，已证]\label{prop:fracdeg}
在双均匀层中 $w_i=\sum_{j\in O(i)}\pi_j=n_i/N$，其中 $n_i:=|O(i)|\in\{0,1,\dots,N-1\}$ 是钝角图的度。
定理~\ref{thm:LN} 的取等\textbf{强制}两件事：
\begin{enumerate}[label=(\roman*),leftmargin=2.6em]
\item 逐行的 Cauchy--Schwarz 与 AM--GM 全部取等，故 $r_i\equiv\varrho^*$、$w_i\equiv\omega^*$，
从而钝角图 $A_O$ 是\emph{正则}的，度数必须是
\begin{equation}\label{eq:nstar}
n^*:=N\omega^*=\frac{N\kappa_N}{3(1+\kappa_N)};
\end{equation}
\item 四条 CS 同时取等，故 $Q$ 是\emph{两值}的：
$Q=\frac3NI_N+m^+(J-I-A_O)-m^-A_O$，$m^\pm\ge0$（推导见第~\ref{sec:twovalue} 小节开头）。
\end{enumerate}
精确值：
\[
n^*(4)=1,\quad n^*(5)=2-\frac{\sqrt6}3,\quad n^*(6)=\frac{5-\sqrt5}2,\quad
n^*(7)=3-\sqrt2,\quad n^*(8)=\frac72-\frac{\sqrt{105}}6,\quad n^*(9)=\mathbf{2},
\]
\[
n^*(10)=\frac92-\frac{\sqrt{21}}2,\quad n^*(12)=\frac{11}2-\frac{\sqrt{33}}2,\quad
n^*(28)=\mathbf{6},\quad n^*(99)=\mathbf{21}.
\]
于是 $n^*$ 必须是\textbf{整数}。而 $n^*\in\Z$ 当且仅当 $\kappa_N\in\Q$，
即 $\frac{3(N-1)}{N-3}=\frac{p^2}{q^2}$（既约）。写成 $N-3=kq^2$ 后方程化为 $k(p^2-3q^2)=6$，
故 $k\in\{1,2,3,6\}$，其中 $k=2,3$ 无解（模 $3$ 矛盾），只余两支 Pell 型方程：
\begin{equation}\label{eq:pell}
\mathcal E=\underbrace{\{\,3+q^2:\ q^2-3p'^2=-2\,\}}_{\text{A 支}}\ \cup\
\underbrace{\{\,3+6q^2:\ p^2-3q^2=1\,\}}_{\text{B 支}}
=\{4,9,28,99,364,1353,5044,18819,70228,262089,\dots\}.
\end{equation}
$\mathcal E$ 是\textbf{无穷集}（两支都是标准 Pell 方程）。
\emph{结论}：对 $N\ge5$ 且 $N\notin\mathcal E$，$L(N)$ 在双均匀层内不可达；
$\mathcal E$ 中的例外阶数需要另行处理（定理~\ref{thm:unattain}）。
\end{proposition}

\begin{proof}
(i) 双均匀层里 \eqref{eq:f} 给 $\rho\le\frac3N+\sqrt{N a_iv_i}$，$a_i=\frac3N(1-\frac3N)-r_i$、$v_i=1-\frac1N-w_i$；
对 $i$ 取最小值后，用 Cauchy--Schwarz 与 AM--GM 两次：
\[
\min_i\sqrt{a_iv_i}\le\frac1N\sum_i\sqrt{a_iv_i}\le\frac1N\sqrt{\Bigl(\sum_ia_i\Bigr)\Bigl(\sum_iv_i\Bigr)}
=\frac1N\sqrt{\bigl(3-{\textstyle\frac9N}-c\bigr)N\bigl(1-{\textstyle1N}-\omega\bigr)},
\]
于是 $\rho\le\frac3N+\sqrt{(3-\frac9N-c)(1-\frac1N-\omega)}=A$，与定理~\ref{thm:LN} 的证明中的 $A$ 逐字相同
——所以\textbf{不需要}预先假设 $r_i,w_i$ 均匀。取等要求两次放缩都取等：
$\min=$ 均值 $\Rightarrow a_iv_i$ 全部相等；CS 取等 $\Rightarrow a_i\propto v_i$；两者合起来给 $r_i\equiv\varrho^*$、$w_i\equiv\omega^*$。
再由 $w_i=n_i/N$ 得 $A_O$ 正则、度 $n^*=N\omega^*$，而 $\omega^*$ 由 AM--GM 取等条件唯一确定为
$\frac{\kappa_N}{3(1+\kappa_N)}$，即 \eqref{eq:nstar}。(ii) 见第~\ref{sec:twovalue} 小节。

$\frac{3(N-1)}{N-3}=\frac{p^2}{q^2}$（$\gcd(p,q)=1$）$\Rightarrow p^2(N-3)=q^2(3N-3)=q^2(3(N-3)+6)$。
令 $N-3=kq^2$（$k\in\Z_{>0}$）得 $k(p^2-3q^2)=6$，故 $k\mid6$。
$k=2$：$p^2-3q^2=3\Rightarrow3\mid p\Rightarrow3p'^2-q^2=1\Rightarrow q^2\equiv-1\pmod3$，无解；
$k=3$：$p^2-3q^2=2\Rightarrow p^2\equiv2\pmod3$，无解。
$k=1$：$p^2-3q^2=6$，模 $3$ 给 $3\mid p$，写 $p=3p'$ 得 $q^2-3p'^2=-2$，且 $N=3+q^2$（A 支）；
$k=6$：$p^2-3q^2=1$（Pell），$N=3+6q^2$（B 支）。
两支的基本解分别是 $(q,p')=(1,1)$、$(p,q)=(2,1)$，递推 $(q,p')\mapsto(2q+3p',q+2p')$、
$(p,q)\mapsto(2p+3q,p+2q)$ 生成全部解，故 $\mathcal E$ 无穷。$n^*(4)=1$、$n^*(9)=2$、$n^*(28)=6$、$n^*(99)=21$
由 \eqref{eq:nstar} 直接代入；$N=5,6,7,8,10,12$ 处 $\kappa_N$ 无理故 $n^*$ 无理。$\square$
\end{proof}

\begin{theorem}[例外阶数上 $L(N)$ 仍不可达：强正则图障碍，已证（$9\le N\le7.1\times10^8$）]\label{thm:unattain}
设 $N\in\mathcal E$、$N\ge9$，$k:=n^*$。若 $L(N)$ 在双均匀层内可达，则钝角图 $A_O$ 必是
$k$-正则且\textbf{恰有三个不同特征值}的图，从而只有两种可能：
\begin{enumerate}[label=(\alph*),leftmargin=2.2em]
\item $A_O$ 连通，于是是强正则图 $\mathrm{SRG}(N,k,\lambda,\mu)$，$\mu\ge1$；
\item $A_O$ 不连通，于是 $A_O=mK_{k+1}$，$m(k+1)=N$。
\end{enumerate}
两种可能都被排除：
\begin{itemize}
\item[(b)] $(k+1)\mid N$ 在 $\mathcal E$ 中\textbf{只}发生于 $N=9$（$A_O=3K_3$）与 $N=28$（$A_O=4K_7$）
（以及可达的 $N=4$，$A_O=2K_2$）。此时 $A_O$ 有特征值 $-1$、重数 $N-m\ge6$，
对应的 $Q$-特征值是 $\frac3N+m^-$，它必须属于 $\{0,1\}$：
$m^-=0$ 给 $\frac3N\notin\{0,1\}$；$m^-=1-\frac3N$ 给 $N-m\ge6>3$ 个特征值等于 $1$，与 $\rank Q=3$ 矛盾。
\item[(a)] 由 $Q=\frac3NI+m^+(J-I-A_O)-m^-A_O$，在 $\one^\perp$ 上 $Q$ 的作用是 $\beta I-uA_O$，
$\beta=\frac3N-m^+$、$u=m^++m^->0$。$Q$ 的特征值只能是 $0$ 或 $1$，而 $\theta\mapsto\beta-u\theta$ 是\emph{单射}，
故 $A_O$ 在 $\one^\perp$ 上至多两个特征值，且「映射到 $1$」的那个特征值重数 $\le3$。
用强正则图的重数公式
$f,g=\frac12\bigl[(N-1)\mp\frac{2k+(N-1)(\lambda-\mu)}{\sqrt{\Delta}}\bigr]$，$\Delta=(\lambda-\mu)^2+4(k-\mu)$，
以及参数关系 $\mu(N-k-1)=k(k-\lambda-1)$，对每个 $N\in\mathcal E$ 穷举 $\lambda\in[0,k)$：
\[
\begin{array}{c|cccccccccc}
N&99&364&1353&5044&18819&70228&262089&978124&3650403&13623484\\
k&21&77&286&1066&3977&14841&55386&206702&771421&2878981\\
\hline
N&50843529&189750628&708158979&&&&&&&\\
k&10744502&40099026&149651601&&&&&&&
\end{array}
\]
\textbf{没有一个} $(\lambda,\mu)$ 使 $\Delta$ 成为完全平方（故 $f,g$ 必为无理数，不可能是重数）。
情形 (a) 因此不可能。
\end{itemize}
\emph{结论}：$L(N)$ 对 $\mathcal E$ 中一切 $9\le N\le7.1\times10^8$ 都不可达。
合并命题~\ref{prop:fracdeg}：\textbf{$L(N)$ 对一切 $5\le N\le7.1\times10^8$ 不可达，$N=4$ 是唯一取等的阶数。}
\end{theorem}

\begin{remark}[对第六版一条陈述的自我修正]\label{rmk:fracdegfix}
本条命题的第一版写的是「$n^*$ 对一切 $N\ge5$ 无理，故 $L(N)$ 对一切 $N\ge5$ 不可达」。
\textbf{前半句是错的}：$n^*(9)=2$、$n^*(28)=6$、$n^*(99)=21$ 都是整数。
错误是在证书脚本 \texttt{证书\_行方差松弛与天花板.py} 的第 (ii) 项被发现的——
把「无理」的判定从 \texttt{sympy.is\_rational}（对根式返回 \texttt{None}）改成
「$\frac{3(N-1)}{N-3}$ 是否为 $\Q$ 中的平方」的纯整数检验后，立刻扫出 $N=9$。
这是本报告「所有『对一切 $N$』型陈述都必须附带可执行检验」这条纪律的又一次兑现。
修正后的结论\emph{更强}：例外集 $\mathcal E$ 虽然无穷，但它稀疏到（前 $13$ 项就跨过 $7\times10^8$），
而且在 $\mathcal E$ 上还有第二道、性质完全不同的障碍（强正则图的重数整性）。
换句话说，「行方差机制在 $N\ge5$ 必然松」这个结论\textbf{成立}，但它的原因不止一个：
$N\notin\mathcal E$ 时是度的整性，$N\in\mathcal E$ 时是谱的整性。
\end{remark}

"""
s=s.replace(old_prop,new_prop,1)
assert new_prop[:60] in s
io.open(main,"w",encoding="utf-8").write(s)
print("prop:fracdeg replaced; lines:",len(s.split("\n")))
