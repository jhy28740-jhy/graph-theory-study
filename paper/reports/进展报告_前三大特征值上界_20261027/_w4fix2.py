import io,re
p="_frag_rowvar_4b.tex"; s=io.open(p,encoding="utf-8").read()
old = r"""f_{m,1}=\frac{4m}{3m+1}.
\end{equation}
\item 若 $v$ 有 $m$ 个正、$m$ 个负且两组内部分别相等（$N=2m$，$v=(c\one_m,-c\one_m)$、$c^2=\frac1{2m}$），则
\begin{equation}\label{eq:fmm}
\lambda_{\max}(Q_+)=f_{m,m}=1+\frac{m-1}{2m}=\frac{3m-1}{2m}.
\end{equation}
\end{proposition}"""
new = r"""f_{m,1}=\frac{4m}{3m+1}.
\end{equation}
\item 若 $v$ 有 $m$ 个正、$m$ 个负且两组内部分别相等（$N=2m$，$v=(c\one_m,-c\one_m)$、$c^2=\frac1{2m}$），则
\begin{equation}\label{eq:fmm}
\lambda_{\max}(Q_+)=f_{m,m}=1+\frac{m-1}{2m}=\frac{3m-1}{2m}.
\end{equation}
\end{enumerate}
\end{proposition}"""
assert old in s; s=s.replace(old,new)
old2 = r"""设 $v\in\R^N$、$\|v\|=1$、$Q=I_N-vv^\top$。
\begin{enumerate}[label=(\alph*),leftmargin=2.2em]
\item 若 $v$ 的坐标有 $m$ 个正、$1$ 个负、其余为零（$N=m+1$），记 $v=(c\one_m,-d)$、$mc^2+d^2=1$，则"""
new2 = r"""设 $v\in\R^N$、$\|v\|=1$、$Q=I_N-vv^\top$（于是 $Q\in\mathcal P_{N-1}(N)$）。
\begin{enumerate}[label=(\alph*),leftmargin=2.2em]
\item 若 $v$ 的非零坐标有 $m$ 个正、$1$ 个负（$N\ge m+1$），记 $v=(c\one_m,-d,0,\dots,0)$、$mc^2+d^2=1$，则"""
assert old2 in s; s=s.replace(old2,new2)
old3 = r"""\lambda_{\max}(Q_+)=\tfrac12\bigl(1+(m-1)c^2\bigr)+\tfrac12\sqrt{(d^2-c^2)^2+4mc^2d^2},
它在 $c^2=\frac2{3m+1}$、$d^2=\frac{m+1}{3m+1}$ 处取最大值"""
new3 = r"""\lambda_{\max}(Q_+)=\max\Bigl\{1,\ \tfrac12\bigl(1+(m-1)c^2\bigr)+\tfrac12\sqrt{(d^2-c^2)^2+4mc^2d^2}\Bigr\},
\]
其中第二项在 $c^2=\frac2{3m+1}$、$d^2=\frac{m+1}{3m+1}$ 处取最大值"""
assert old3 in s; s=s.replace(old3,new3)
old4 = r"""f_{m,1}=\frac{4m}{3m+1}.
\end{equation}"""
new4 = r"""f_{m,1}=\frac{4m}{3m+1}\ (>1\ \text{当}\ m\ge2).
\[​"""
# revert: too fiddly; instead just close the display properly
s=s.replace(new3.replace("\\]","XXX"),new3.replace("\\]","XXX"))
io.open(p,"w",encoding="utf-8").write(s)
print("done stage1")
