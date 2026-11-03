import io
p="_frag_rowvar_4b.tex"
L=io.open(p,encoding="utf-8").read().split("\n")
L[103]=r"设 $v\in\R^N$、$\|v\|=1$、$Q=I_N-vv^\top$（于是 $Q\in\mathcal P_{N-1}(N)$，即 $Q$ 在余秩一层）。"
L[105]=r"\item 若 $v$ 的非零坐标有 $m$ 个正、$1$ 个负（$N\ge m+1$），记 $v=(c\one_m,-d,0,\dots,0)$、$mc^2+d^2=1$。"
L[106]=r"此时 $Q_+$ 在剩下的 $N-m-1$ 个坐标方向上是孤立的对角元 $1$，而在前 $m+1$ 个方向上给出"
L[107]=r"\lambda_{\max}(Q_+)=\max\Bigl\{1,\ \tfrac12\bigl(1+(m-1)c^2\bigr)+\tfrac12\sqrt{(d^2-c^2)^2+4mc^2d^2}\Bigr\},"
L[109]=r"其中第二项在 $c^2=\frac2{3m+1}$、$d^2=\frac{m+1}{3m+1}$ 处取最大值"
L[111]=r"f_{m,1}=\frac{4m}{3m+1}\qquad(m\ge2\ \text{时}\ f_{m,1}>1).".replace("\\qquad(","\\qquad(")
L[116]=L[116]+"\n"+r"\end{enumerate}"
# fix cor:k23 proof classification
for i,l in enumerate(L):
    if "\\emph{全同号}" in l:
        L[i]=r"\emph{三个坐标同号，或某个坐标为零且其余两个同号}时 $Q_+$ 的非对角元全部被抹成 $0$，于是它是对角阵，"
        L[i+1]=r"$\lambda_{\max}(Q_+)=\max_i(1-v_i^2)\le1<\frac87$；"
        print("patched proof lines",i,i+1)
io.open(p,"w",encoding="utf-8").write("\n".join(L))
print("ok")
