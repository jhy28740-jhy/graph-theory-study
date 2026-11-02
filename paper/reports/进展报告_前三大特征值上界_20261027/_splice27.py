# -*- coding: utf-8 -*-
import io, shutil, sys
tex = "进展报告_前三大特征值之和的精确上界.tex"
frag = "_frag_r27.tex"
shutil.copyfile(tex, tex + ".bak_r27")
with io.open(tex, encoding="utf-8") as f:
    src = f.read()
with io.open(frag, encoding="utf-8") as f:
    fr = f.read()
if not fr.endswith("\n"):
    fr += "\n"

anchor = "\\section{通道型主函数方案的精确天花板}\\label{sec:nogo}"
assert src.count(anchor) == 1, ("anchor count", src.count(anchor))
src = src.replace(anchor, fr + "\n" + anchor, 1)

# fix rem:r26struct typo: remove leftover -a^2
old = "g_1^2a^2(-\\langle\\cdot,\\cdot\\rangle-a^2)=g_1^2s_0^2"
new = "g_1^2a^2(-\\langle\\cdot,\\cdot\\rangle)=g_1^2s_0^2"
assert src.count(old) == 1, ("typo count", src.count(old))
src = src.replace(old, new, 1)

# bump \date: 第十版 -> 第十一版, and append R27 summary before the closing
d_old_head = "（第十版，在第九版基础上增两节。"
d_new_head = "（第十一版，在第十版基础上增一节（\\S\\ref{sec:r27}，第 27 轮）。"
assert src.count(d_old_head) == 1, ("date head", src.count(d_old_head))
src = src.replace(d_old_head, d_new_head, 1)

d_old_tail = "第九版原有的 \\S\\ref{sec:r24} 内容全部保留。）}"
r27sum = ("\\S\\ref{sec:r27}（第 27 轮）：在秩一截面上用「Young 不等式 $+$ 刚性恒等式 "
          "$\\sum_{ij}\\rho_i\\rho_jc_{ij}^2=\\|I_2\\|_F^2=2$ $+$ 圆上 Fourier」的链条，\\textbf{无条件}证明 "
          "$\\mathrm{Neg}\\le\\frac5{6\\pi}=0.265258$，从而 $\\kappa_3(\\text{秩一})\\le1+\\frac5{6\\pi}=1.265258<\\rho_0$"
          "（定理~\\ref{thm:r27main}），把第 26 轮\\emph{条件性}的 $\\frac{11+\\varphi}{10}=1.261803$~升级为\\emph{无条件}结果、"
          "并较 \\STW\\ 的 $\\Phi=1.309017$~改进 $4.38\\times10^{-2}$；关键机制是二维各向同性 $\\sum\\rho_ie^{2i\\varphi_i}=0$~"
          "自动消灭钝角余弦平方核 $K(u)=(\\cos u)_-^2$~的第二 Fourier 系数 $\\widehat K_2=\\frac18$"
          "（除常数项外最大的正系数），使 $P\\le\\frac{10}{3\\pi}$（定理~\\ref{thm:r27P}；闭式 $S=\\frac5{24}-\\frac\\pi{16}$，"
          "引理~\\ref{lem:r27S}）；另附秩一刚性 $A_w\\equiv1$~的一行新证明（$\\|Z_S\\|^2\\le\\sigma_S$，引理~\\ref{lem:r27zs}、"
          "命题~\\ref{prop:r27aw}），并指出本轮 Fourier 机器为秩一\\emph{专属}、一般 $\\Lam$~的 Young 余项权重 "
          "$g_ih_ib_i^2/a_i\\ne\\rho_i$~使刚性恒等式失效，故全局闭合仍须锥最优性（注~\\ref{rem:r27global}）。")
d_new_tail = (r27sum + "第十版原有的 \\S\\ref{sec:r25}、\\S\\ref{sec:r26} 与第九版的 \\S\\ref{sec:r24} 内容全部保留。）}")
assert src.count(d_old_tail) == 1, ("date tail", src.count(d_old_tail))
src = src.replace(d_old_tail, d_new_tail, 1)

with io.open(tex, "w", encoding="utf-8") as f:
    f.write(src)
print("spliced OK; new line count =", src.count("\n") + 1)
