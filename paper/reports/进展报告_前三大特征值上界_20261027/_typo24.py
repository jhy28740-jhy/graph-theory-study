# -*- coding: utf-8 -*-
import io
p = '进展报告_前三大特征值之和的精确上界.tex'
s = io.open(p, encoding='utf-8-sig').read()
reps = [
# --- (1) prop:r24sharp: split the wide display in two
(r"""\[
m_{\rm pole}=\delta,\quad m_j=\frac{3-\delta}n\ (j=1,\dots,n),\qquad
\ell_{\rm pole}=e_3,\quad
\ell_j=\Bigl(\sqrt{1-u_2}\cos\tfrac{2\pi j}n,\ \sqrt{1-u_2}\sin\tfrac{2\pi j}n,\ \sqrt{u_2}\Bigr),
\]""",
 r"""\[
m_{\rm pole}=\delta,\qquad m_j=\frac{3-\delta}n\quad(j=1,\dots,n),\qquad \ell_{\rm pole}=e_3,
\]
\[
\ell_j=\Bigl(\sqrt{1-u_2}\cos\tfrac{2\pi j}n,\ \sqrt{1-u_2}\sin\tfrac{2\pi j}n,\ \sqrt{u_2}\Bigr),
\]"""),
# --- (2) thm:r24centre step 1: split the wide display
(r"""\[
6\rho_0-J=6\rho_0-L-A^2=B\Bigl(m_4-\frac15\Bigr)+\frac V2+D-\delta_0,
\qquad \delta_0=3\varphi+3-6\rho_0=\frac{9+7\sqrt5-4\sqrt{35}}6=0.1646927850166773 ,
\]""",
 r"""\[
\begin{gathered}
6\rho_0-J=6\rho_0-L-A^2=B\Bigl(m_4-\frac15\Bigr)+\frac V2+D-\delta_0,\\[2pt]
\delta_0=3\varphi+3-6\rho_0=\frac{9+7\sqrt5-4\sqrt{35}}6=0.1646927850166773 .
\end{gathered}
\]"""),
# --- (3) branch 1: pull Phi>=Q(V) out of the display, drop the cases block
(r"""\[
\Phi(V)\ge c_2V^2+c_1V+c_0=:Q(V),\qquad
\begin{cases}
c_2=\frac{35B}{72}R_{\rm w}^2=0.4130017665424713,\\[1pt]
c_1=\frac12-\frac{7B}{27}R_{\rm w}=-0.0237916619419962,\\[1pt]
c_0=B\bigl(\frac{19}{81}-\frac15\bigr)-\delta_0=\frac{14B}{405}-\delta_0=0.0013825973096503.
\end{cases}
\]""",
 r"""$\Phi(V)\ge c_2V^2+c_1V+c_0=:Q(V)$，其中
\[
c_2=\frac{35B}{72}R_{\rm w}^2=0.4130017665424713,\qquad
c_1=\frac12-\frac{7B}{27}R_{\rm w}=-0.0237916619419962,
\]
\[
c_0=B\Bigl(\frac{19}{81}-\frac15\Bigr)-\delta_0=\frac{14B}{405}-\delta_0=0.0013825973096503 .
\]"""),
# --- (4) lem:r24mom proof: the long inline chain -> display
(r"""第二个等号：$(A_4)_{ijkl}=\sum_am_a\xi_{ai}\xi_{aj}\xi_{ak}\xi_{al}$，故
$\|A_4\|_F^2=\sum_{a,b}m_am_b\bigl(\sum_i\xi_{ai}\xi_{bi}\bigr)^4=\sum_{a,b}m_am_b\langle\ell_a,\ell_b\rangle^4=9m_4$
（与 $m_4$ 的定义一致，也与命题~\ref{prop:r23harm}~的记号一致）。""",
 r"""第二个等号：$(A_4)_{ijkl}=\sum_am_a\xi_{ai}\xi_{aj}\xi_{ak}\xi_{al}$，故
\[
\|A_4\|_F^2=\sum_{a,b}m_am_b\Bigl(\sum_i\xi_{ai}\xi_{bi}\Bigr)^4
=\sum_{a,b}m_am_b\langle\ell_a,\ell_b\rangle^4=9m_4
\]
（与 $m_4$ 的定义一致，也与命题~\ref{prop:r23harm}~的记号一致）。"""),
# --- (5) rem:r24drop: the three moment relations -> display
(r"""原因是这六块恰好携带了两个被各向同性完全决定的量：
$\mathcal{M}_{400}=\sum_im_iu_i^2=\mathcal{X}$ 与
$\mathcal{M}_{220}+\mathcal{M}_{202}=\sum_im_iu_i(1-u_i)=1-\mathcal{X}$，
以及 $\mathcal{M}_{040}+\mathcal{M}_{004}+2\mathcal{M}_{022}=\sum_im_i(1-u_i)^2=1+\mathcal{X}$。""",
 r"""原因是这六块恰好携带了三个被各向同性完全决定的量：
\[
\mathcal{M}_{400}=\sum_im_iu_i^2=\mathcal{X},\qquad
\mathcal{M}_{220}+\mathcal{M}_{202}=\sum_im_iu_i(1-u_i)=1-\mathcal{X},
\]
\[
\mathcal{M}_{040}+\mathcal{M}_{004}+2\mathcal{M}_{022}=\sum_im_i(1-u_i)^2=1+\mathcal{X}.
\]"""),
# --- (6) cor:r24Rw: h(V) inline -> display
(r"""由 $A^2=3-\frac V2$，$A$ 随 $V$ 递减，故 $h(V):=R(A(V))\,V
=\frac{(1+1/\sqrt3)^2V}{\sqrt3(\sqrt3+\sqrt{3-V/2})}$ 关于 $V$ 严格递增，且
$\sigma_u^2\le h(V)$。""",
 r"""由 $A^2=3-\frac V2$，$A$ 随 $V$ 递减，故
\[
h(V):=R(A(V))\,V=\frac{(1+1/\sqrt3)^2V}{\sqrt3\,\bigl(\sqrt3+\sqrt{3-V/2}\bigr)}
\]
关于 $V$ 严格递增，且 $\sigma_u^2\le h(V)$。"""),
]
for old, new in reps:
    assert old in s, 'NOT FOUND:\n' + old[:120]
    s = s.replace(old, new, 1)
io.open(p, 'w', encoding='utf-8').write(s)
print('all 6 typography fixes applied; lines =', s.count('\n')+1)
