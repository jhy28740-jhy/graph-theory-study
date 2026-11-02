# 图论仓库：谱极值研究 + 代数图论课程资料

本仓库有两条线：一条是**前三大特征值之和的精确上界**研究（核心成果，见 `paper/`），
一条是《代数图论》课程整理资料（见 `notes/`）。

2026-11-02 做过一次全量分类整理：根目录原有 884 个散落文件 + 26 个目录，
现按类别归入 8 个顶层目录。**整理只做移动，未删除任何文件**；
移动明细见 `archive/_reorg_log_20261102.txt`。

## 顶层结构

| 目录 | 内容 |
|---|---|
| `paper/` | **主论文及其全部相关材料**：论文 tex/pdf、参考文献链接清单、arXiv 原文、24 页汇报 PPT、讲稿与排版建议、33 页审核报告、早期进展报告。要发给人看，取这一个文件夹即可 |
| `thesis_spectral_extremes/` | 学位论文「谱极值」主线工作区 |
| `notes/` | 课程讲义与知识点梳理、研究进展与恢复入口（含 24/43/53/61/115 页存档）、研究方法优先级协议、文献综述笔记、迁移记录 |
| `scripts/` | 数值验证与证书检验脚本（约 500 个 Python，另有 PowerShell / Node） |
| `data/` | 脚本输出：`.npy` / `.npz` / `.json` / `.bin` / `.txt` 中间结果 |
| `tex_fragments/` | LaTeX 片段（论文拼装用） |
| `assets/` | 图片（论文页面裁剪图等） |
| `archive/` | 历史阶段归档：`s13_stage2`–`s13_stage12`、`s13_checkpoint`、三个阶段性 zip、`migration_20261020`、`_pdfbuild`、`_chapters`、`_clean`、`output`、`tmp`、`__pycache__` |

## 一、主论文（核心成果）

**A joint moment bound for the three largest graph eigenvalues**
Haoyang Ji，November 2, 2026，20 页，pdflatex。

在线 PDF：https://github.com/jhy28740-jhy/graph-theory-study/blob/main/paper/top3_joint_moment_bound.pdf

主结果：对任意 n ≥ 3 阶简单图 G，

    λ₁(G) + λ₂(G) + λ₃(G) ≤ c·n − 3,
    c = 1.284673255723918505…

改进 Sivashankar–Tang–Wakhare（arXiv:2608.02429）Proposition 6.1 的系数 (3+√5)/4，
相对改进 1.86%。另证每个极大 graphon 的典范本征函数律至多 24 个支撑点。
最优系数与极图均未确定。

详见 `paper/README.md` 与 `paper/参考文献链接.md`。

## 二、课程资料（原仓库主体）

| 文件 | 说明 | 规模 |
| --- | --- | --- |
| `notes/代数图论课程知识点梳理.md` / `.pdf` | 全书 13 章知识点梳理：每章按「定义 → 命题/定理 → 证明思路 → 与前后的连贯性」展开，并给出每章小结 | 约 43 页 PDF |
| `notes/第8章课程讲义_强正则图与设计.md` / `.pdf` | 第 8 章授课讲义（12 学时）：6 讲 + 习题课，含时间分配、板书设计、作业、Q&A 预演 | 约 58 页 PDF |
| `scripts/_ch8_verify.py` | 讲义附录的数值验证脚本，可复现讲义中全部「实测数据」 | — |

**Part I 谱图论**：1 图论基本概念 / 2 矩阵理论工具 / 3 图的关联矩阵 / 4 谱界 /
5 谱极值问题 / 6 谱图论近期进展选讲。
**Part II 代数图论**：7 对称性与 Cayley 图 / 8 强正则图与设计 / 9 特征标与有限 Fourier 分析 /
10 代数图的谱计算 / 11 有限域与有限几何 / 12 二次特征与 Paley 图 / 13 割、流与临界群。

第 8 章六讲：强正则图（§8.1）/ 邻接代数（§8.2）/ 设计图与部分设计图（§8.3）/
距离正则图（§8.4）/ 区组设计与关联矩阵（§8.5）/ 全章总结与习题课。

课程 PDF 由 `markdown-it` + `markdown-it-texmath`（KaTeX 服务端渲染）+ `mermaid` +
Playwright/Chromium 生成，构建目录已归档至 `archive/_pdfbuild/`。

## 三、旧路径映射（整理后的对照）

| 整理前 | 整理后 |
|---|---|
| `_26*.pdf`、`_ref*.pdf`、`_paulin*`、`.tmp_lit/`、`weighted_spectral_literature_audit/` | `paper/literature/` |
| `汇报_前三大特征值上界_20261102/`、`审核报告_top3_joint_moment_bound_20261102/`、`进展报告_前三大特征值上界_20261027/` | `paper/reports/` 下同名目录 |
| `代数图论课程知识点梳理.*`、`第8章课程讲义_*`、`研究进展与恢复入口*.md`、`研究方法优先级协议.md`、`review_spectral_moore.*`、`迁移完成记录_*.md` | `notes/` |
| `*.py`、`*.ps1`、`*.mjs` | `scripts/` |
| `*.npy`、`*.npz`、`*.json`、`*.bin`、其余 `*.txt` | `data/` |
| `*.tex`（片段，`review_spectral_moore.tex` 除外） | `tex_fragments/` |
| `*.png` | `assets/` |
| `s13_stage*/`、`s13_checkpoint/`、`s13_*.zip`、`migration_20261020/`、`_pdfbuild/`、`_chapters/`、`_clean/`、`output/`、`tmp/`、`__pycache__/` | `archive/` |

注意：`notes/研究进展与恢复入口*.md` 等文档中记载的相对路径可能仍指向整理前的位置，
按上表换算即可。脚本内部若使用相对路径读取 `data/` 下的中间结果，需要相应调整。
