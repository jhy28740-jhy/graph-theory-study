# 论文文件夹：前三大特征值之和的精确上界

本文件夹是这篇论文的**唯一完整入口**，论文本体、参考文献、汇报材料、审核报告全部在这里。

- 题目：**A joint moment bound for the three largest graph eigenvalues**
- 作者：Haoyang Ji 日期：November 2, 2026
- 篇幅：20 页 A4，编译器 **pdflatex**（不是 xelatex）
- 在线 PDF：https://github.com/jhy28740-jhy/graph-theory-study/blob/main/paper/top3_joint_moment_bound.pdf

## 主结果

对任意 n ≥ 3 阶简单图 G：

    λ₁(G) + λ₂(G) + λ₃(G) ≤ c·n − 3

    c = 106106902139519/107310156800000 + 41√3/240 − 53/9000000
      = 1.284673255723918505…

改进 STW（arXiv:2608.02429）Proposition 6.1 的系数 (3+√5)/4 = 1.309016994374947…
系数改进量 0.024343738651028919（相对改进 1.86%），常数项 −3 不变。

结构定理：每个极大 graphon 的典范本征函数律至多有 **24** 个支撑点。

诚实标注：**最优系数与极图均未确定**。
已知数值天花板 (29+10√5)/40 = 1.284016994374948（比 c 小 0.000656261）；
数值锚点 ρ₀ = (9+√5+2√35)/18 = 1.281568196872168（K₁∨C₅ 吹起族的极限）。

## 文件夹内容

| 路径 | 说明 |
|---|---|
| `top3_joint_moment_bound.tex` | 论文源码，1442 行 |
| `top3_joint_moment_bound.pdf` | 论文成品，20 页 |
| `参考文献链接.md` | **全部链接清单**：12 条参考文献 + 已下载未引用文献 + 汇报材料 + 待办 R1 |
| `literature/` | arXiv 原文 PDF 与检索工作区（STW、KLMPT、Tang λ₃、Huang–Wei 谱和等） |
| `reports/汇报_前三大特征值上界_20261102/` | 24 页自编译汇报 PPT（pptx + pdf）、16 页讲稿与排版建议（tex + pdf）、PPT 生成脚本与图表资源 |
| `reports/审核报告_top3_joint_moment_bound_20261102/` | 33 页详细审核报告（tex + pdf）：逐节复核、符号验证、数值验证、问题清单 G1/G2、表述建议 S1–S7、松弛量八项分解、天花板分析、优先级排序 |
| `reports/进展报告_前三大特征值上界_20261027/` | 早期进展报告与其工作脚本 |

## 重新编译

    set PATH=%LOCALAPPDATA%\Programs\MiKTeX\miktex\bin\x64;%PATH%
    pdflatex top3_joint_moment_bound.tex
    pdflatex top3_joint_moment_bound.tex

两遍以解析交叉引用。当前日志干净：0 error、0 overfull/underfull、
12 条 bibitem 与全部 \cite 闭合、65 个 label 无重复无悬空引用。

中文文档（讲稿、审核报告）用 xelatex + ctexart，同样跑两遍。

## 版本状态

本次修订（2026-11-02 16:35）新增三处：

1. `\author{Haoyang Ji}`、`\date{November 2, 2026}`（并写入 PDF 元数据 pdfauthor）
2. 第 113 行紧性依据修正：原来写「T_W 在无穷维空间上自伴故紧」是不成立的理由，
   改为「W ∈ L²([0,1]²) ⇒ T_W 为 Hilbert–Schmidt 故紧；又自伴，谱为实且仅在 0 处聚集」
3. 上一轮审核的 G1、G2、S1–S7 共 9 项已全部落实（落实位置见审核报告与本文件末节）

SHA256(tex) = `516C37E0A826AF959B1BE61DAD16D18FB8AB917193C0A81636E903115EF62360`
SHA256(pdf) = `34CD29471CE520E6BF148F80D2A8274515508F34822D56EC4C7C3760DEAF2280`

## 仍未结的一项

**R1**：Ebrahimi Boroojeni 2008 SFU 硕士论文定理 3.4.1 尚未获取，
新颖性声明的最终确认悬置。检索入口 https://summit.sfu.ca/ （需校园网）。
