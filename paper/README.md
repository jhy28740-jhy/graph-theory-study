# 论文文件夹：前三大特征值之和的精确上界

本文件夹是这篇论文的**唯一完整入口**，论文本体、参考文献、汇报材料、审核报告全部在这里。

- 题目：**A joint moment bound for the three largest graph eigenvalues**
- 作者：Haoyang Ji 日期：November 3, 2026
- 篇幅：22 页 A4，编译器 **pdflatex**（不是 xelatex）
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
**无条件缺口 c − ρ₀ = 0.003105058851750628**（论文第 10 节 Further questions 已明确写出）。

### 查新状态（2026-11-03 第二次复核，详见 `参考文献链接.md` §七）

- **主定理新颖性安全**：主流文献中无人改进过 STW Proposition 6.1 的 (3+√5)/4，c 仍是公开可查的最好系数。
- **必须引用的先行记录**：Zenodo 上 T. Newton 的两条自存档
  （doi:10.5281/zenodo.21785906，2026-08-04；doi:10.5281/zenodo.22082842，2026-08-24）
  以「pole–pentagon law」提出了**完全相同**的常数 α₃ = 1/2+√5/18+√35/9 = ρ₀
  与**完全相同**的权重 (7+5√7)/126，并声称 Λ₃ ≤ α₃n − 3；但该声称**明确是条件性的**
  （定理标题为 *Candidate*，依赖机器证书链，作者自请独立审计），且未经同行评审。
  论文已新增 `Remark 2.3` 如实说明，本文**不主张该常数的优先权**，
  Theorem 1.1 的无条件性与自足性不受影响。
- 本次修订新增 4 条引用：`HW`（Huang–Wei，S₂ 精确解）、`Paulin`（二次不等式）、`Newton1`、`Newton2`。

## 文件夹内容

| 路径 | 说明 |
|---|---|
| `top3_joint_moment_bound.tex` | 论文源码，1534 行 |
| `top3_joint_moment_bound.pdf` | 论文成品，22 页 |
| `参考文献链接.md` | **全部链接清单**：16 条参考文献 + 已下载未引用文献 + 汇报材料 + R1 结论 + §七 第二次查新复核 |
| `literature/` | arXiv 原文 PDF 与纯文本抽取（STW、KLMPT、Tang λ₃、Huang–Wei 谱和、Paulin、Ebrahimi 2008 等），共 18 个文件 |
| `literature/weighted_spectral_literature_audit/` | 加权谱极值的查新审计报告（tex + pdf + 4 张复核截图 + 3 篇笔记 + 检索脚本） |
| `reports/汇报_前三大特征值上界_20261102/` | 24 页自编译汇报 PPT（pptx + pdf）、16 页讲稿与排版建议（tex + pdf）、PPT 生成脚本与图表资源 |
| `reports/审核报告_top3_joint_moment_bound_20261102/` | 33 页详细审核报告（tex + pdf）：逐节复核、符号验证、数值验证、问题清单 G1/G2、表述建议 S1–S7、松弛量八项分解、天花板分析、优先级排序 |
| `reports/进展报告_前三大特征值上界_20261027/` | 早期进展报告（tex + pdf）与 24 个 `证书_*.py` 可复现证书脚本 |
| `reports/查新复核报告_前三大特征值上界_20261103/` | **7 页第二次查新复核报告（tex + pdf）**：文献比对表、Newton Zenodo 记录逐项核验与原文引文、时间线与优先权判断、风险评级、λ₁+2λ₂ 副线查新、给导师的一句话说明 |
| `literature/novelty_20261103/` | 本次查新的原始抓取数据（5 个 json/log）与 4 个检索脚本 |
| `literature/Newton2026a.pdf` / `Newton2026b.pdf` | 两条 Zenodo 记录的 PDF 原件与纯文本抽取 |

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

SHA256(tex) = `DF8E1AC03C8EE9024B2182B977B1EDC76AB03513335A83BFC3D7796806621E08`
SHA256(pdf) = `A8579BE83ED85D76B66DED5F3A559B590E5D6D8D3152C77F3AF17ABD4DD0D1F3`
（PDF 构建于 2026-11-03 09:39，20 页，首页已含 `Haoyang Ji` 与 `November 2, 2026`。）

## R1 已解决（2026-11-03）

Ebrahimi 2008 SFU 硕士论文已直接下载（Summit 为开放获取，无需校园网）：
`literature/Ebrahimi2008_SFU_etd3462.pdf`。其 Theorem 3.4.1 给出
τ₃ ≤ 3/√5 = 1.341640786… < 1.3417，**弱于** (3+√5)/4，
故本论文新颖性不受影响。显式系数链：3/√5 > (3+√5)/4 > c。
详见 `参考文献链接.md` 第五节。

## 本文件夹的整理状态（2026-11-03 第二次收尾）

`paper/` 现在是 **85 个文件 / 13.3 MB**，只含论文本体、参考文献、汇报材料、审核报告四类内容。
第一次整理（2026-11-02）只归了根目录，本文件夹内部仍混着三类东西，本次已全部移出到
`archive/_cleanup_20261103/`（**只移动，未删除**，逐条明细见该目录下的 `_move_log.txt`）：

1. **已弃用的 dashi PPT 产物** `_dashi_弃用/`（209 个文件、11.24 MB，含整套模板运行时、
   woff2 字体、vendor JS）。汇报 PPT 已改用自编译的 python-pptx 版式引擎，dashi 不再使用。
2. **查新检索的原始抓取数据** `weighted_spectral_literature_audit/raw/` 与
   `search_log.json`、`forward_log.json`（72 个文件、7.3 MB 的 openalex / crossref /
   arXiv / bing / 万方 / CNKI 原始响应）。审计报告本体（tex + pdf）保留在 `literature/`。
3. **LaTeX 构建垃圾与一次性脚本**：51 个 `.aux` / `.log` / `.out` / `.toc`，
   进展报告目录下 128 个 `_*.py`、`_frag_*.tex`、`*.bak_*`、`_cert_*.txt`，
   以及 `literature/.tmp_lit/`、`_paulin.html`、与 `.txt` 重复的 `*.lay.txt` / `*.pdftxt`。

另删除了一份与 `literature/_2608.02429.pdf` 字节完全相同的重复 STW PDF（SHA256 一致，
已移至 archive 留存）。