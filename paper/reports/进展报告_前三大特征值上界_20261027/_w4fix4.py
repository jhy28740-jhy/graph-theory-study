import io
p="_frag_rowvar_4b.tex"
L=io.open(p,encoding="utf-8").read().split("\n")
L[106]=r"此时 $Q_+$ 在剩下的 $N-m-1$ 个坐标方向上是孤立的对角元 $1$，而在前 $m+1$ 个方向上给出"+"\n"+r"\["
L[142]=r"故不需要定理~\ref{thm:convex}，直接对 $v$ 的符号模式分类（三类）："
L[144]=r"$\lambda_{\max}(Q_+)=\max_i(1-v_i^2)\le1<\frac87$；\emph{恰有一个坐标为零、另两个异号}时是 (a) 的 $m=1$ 情形，$f_{1,1}=1$；"
io.open(p,"w",encoding="utf-8").write("\n".join(L))
print("ok")
import re
s="\n".join(L)
for env in ["proposition","enumerate","equation","proof","remark","lemma","theorem","corollary","tabular","center","itemize"]:
    b=len(re.findall(r"\\begin\{%s\}"%env,s)); e=len(re.findall(r"\\end\{%s\}"%env,s))
    print(env,b,e,"OK" if b==e else "*** MISMATCH")
print("braces:",s.count("{")-s.count("}"),"  [:",s.count("[")-s.count("]"))
