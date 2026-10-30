# -*- coding: utf-8 -*-
import io, shutil
tex = "进展报告_前三大特征值之和的精确上界.tex"
frag = "_frag_r28.tex"
shutil.copyfile(tex, tex + ".bak_r28")
with io.open(tex, encoding="utf-8") as f:
    src = f.read()
with io.open(frag, encoding="utf-8") as f:
    fr = f.read()
if not fr.endswith("\n"):
    fr += "\n"

anchor = "\\section{通道型主函数方案的精确天花板}\\label{sec:nogo}"
assert src.count(anchor) == 1, ("anchor", src.count(anchor))
src = src.replace(anchor, fr + "\n" + anchor, 1)

d_old_head = "（第十一版，在第十版基础上增一节（\\S\\ref{sec:r27}，第 27 轮）。"
d_new_head = "（第十二版，在第十一版基础上再增一节（\\S\\ref{sec:r28}，第 28 轮）。"
assert src.count(d_old_head) == 1, ("head", src.count(d_old_head))
src = src.replace(d_old_head, d_new_head, 1)

tail_marker = "第十版原有的 \\S\\ref{sec:r25}、\\S\\ref{sec:r26} 与第九版的 \\S\\ref{sec:r24} 内容全部保留。）}"
r28sum = ("\\S\\ref{sec:r28}（第 28 轮）：确认第 27 轮链条中唯一的松弛在 "
          "$P_{\\text{champ}}=\\frac{3+\\sqrt5}5<\\frac{10}{3\\pi}$；用 SLSQP 数值证实 "
          "$\\sup P=\\frac{3+\\sqrt5}5$（$N=3..16$，$5$-环最优），并\\textbf{证明}环族极值 "
          "$\\max_nP(n)=P(5)=\\frac{3+\\sqrt5}5$（命题~\\ref{prop:r28ring}）与刚性 "
          "$|\\widehat\\mu(5)|=2\\Rightarrow|\\widehat\\mu(5m)|=2$（引理~\\ref{lem:r28rigid5}，解释了冠军值的机制来源）；"
          "提出精确值猜想~\\ref{conj:r28P}（$\\sup P=\\frac{3+\\sqrt5}5$），其后果是 "
          "$\\kappa_3(\\text{秩一})=\\frac{11+\\varphi}{10}=1.261803$~精确无条件（推论~\\ref{cor:r28exact}），"
          "但严格证明是符号不定的矩问题、本轮未解；同时\\textbf{精确定位}全局障碍——增稠因子 "
          "$\\kappa_i=h_i/\\sqrt{s_i}>1$~使刚性恒等式 $\\sum\\rho\\rho c^2=2$~失效，只能覆盖单轴薄层 "
          "$a\\in[0.876,1]$、够不到冠军 $a^*=0.5946$，故 $\\kappa_3\\le\\rho_0$~的唯一剩余缺口是锥最优性 P1"
          "（注~\\ref{rem:r28path}）。")
new_tail = r28sum + tail_marker
assert src.count(tail_marker) == 1, ("tail", src.count(tail_marker))
src = src.replace(tail_marker, new_tail, 1)

with io.open(tex, "w", encoding="utf-8") as f:
    f.write(src)
print("R28 spliced OK; lines =", src.count("\n") + 1)
