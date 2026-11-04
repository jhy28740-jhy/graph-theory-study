# 图论仓库：谱极值研究 + 代数图论课程资料

本仓库有两条线：一条是**前三大特征值之和的精确上界**研究（核心成果，见 `paper/`），
一条是《代数图论》课程整理资料（见 `notes/`）。

毕业论文成稿见 `thesis_main/`（全英文 60 页，2026-11-04 入库，commit `8f79fb6`）。

2026-11-02 做过一次全量分类整理：根目录原有 884 个散落文件 + 26 个目录，
现按类别归入 8 个顶层目录。**整理只做移动，未删除任何文件**；
移动明细见 `archive/_reorg_log_20261102.txt`。

2026-11-03 做了第二次收尾：第一次只归了根目录，`paper/` 内部仍混着
已弃用的 dashi PPT 产物（11.2 MB）、查新检索的原始抓取数据（7.3 MB）、
51 个 LaTeX 构建垃圾与 128 个一次性脚本，本次全部移入 `archive/_cleanup_20261103/`
（明细见该目录下 `_move_log.txt`）。`paper/` 因此从 465 个被追踪文件降到 **85 个 / 13.3 MB**。
同时补入了第一次整理时漏提交的 9 个文件（`notes/` 8 个进展记录 + `scripts/_FINAL_RERUN.txt`）。

## 顶层结构

| 目录 | 内容 |
|---|---|
| `paper/` | **主论文及其全部相关材料**：论文 tex/pdf、参考文献链接清单、arXiv 原文、24 页汇报 PPT、讲稿与排版建议、33 页审核报告、早期进展报告。要发给人看，取这一个文件夹即可 |
| `thesis_main/` | **毕业论文（全英文，60 页）**：在主论文 22 页核心之外扩写——新增文献综述、预备知识与桥引理、加权 [0,1] 适用范围、第 7 章数值验证与可复现脚本；含编译好的 `main.pdf` 与 preprint↔thesis 编号对照表 |
| `thesis_spectral_extremes/` | 学位论文「谱极值」主线工作区（约 1.88 GB，整体不入库；其中 `chapter_F2/ch_F2_classification.tex` 是待入库的新章节） |
| `notes/` | 课程讲义与知识点梳理、研究进展与恢复入口（含 24/43/53/61/115 页存档）、研究方法优先级协议、文献综述笔记、迁移记录，共 15 个文件 |
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

**毕业论文（全英文，60 页）**：`thesis_main/main.tex`，pdflatex 三遍编译，零未定义引用、零 >1pt 溢出。

在线 PDF：https://github.com/jhy28740-jhy/graph-theory-study/blob/main/thesis_main/main.pdf

相对主论文新增的内容：Ch1 文献综述（Hong、Nikiforov、Mohar、EMNA、STW、Paulin、HW、KLMPT 等）；
Ch2 预备知识（Ky Fan 最大原理、graphon 的 Σ_k 与割度量、blow-up 谱、投影/各向同性律桥引理 Lemma 2.3）；
§6.4 加权 [0,1] 适用范围（Prop 6.8 加权形式 + Prop 6.10「加权不改变渐近常数」的否定性结果，并如实披露其非新颖性）；
Ch7 数值验证协议与表 7.1–7.4（n≤7 全枚举 2 131 016 个标号图零反例、等权/最优权 blow-up、加权随机对照）；
Ch8 总结、缺口与开放问题。preprint↔thesis 编号对照见 `thesis_main/numbering-cross-table.md`，
可复现脚本与 JSON 输出见 `thesis_main/scripts/` 与 `thesis_main/data/`。
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

## 四、`notes/` 明细

| 文件 | 说明 |
|---|---|
| `代数图论课程知识点梳理.md` / `.pdf` | 全书 13 章知识点梳理（约 43 页 PDF） |
| `第8章课程讲义_强正则图与设计.md` / `.pdf` | 第 8 章授课讲义，12 学时（约 58 页 PDF） |
| `README_课程资料_原始.md` | 课程资料部分的原始说明 |
| `研究进展与恢复入口.md` | **当前研究状态的唯一入口**：结论清单、证明状态、证书位置、下一步 |
| `研究进展与恢复入口_24页存档.md` / `_43页` / `_53页` / `_61页` / `_115页存档.md` | 上述入口文件的历史版本存档 |
| `研究方法优先级协议.md` | 方法优先级协议：数学方法优先于平凡算法搜索，附已尝试方法的排队记录 |
| `迁移完成记录_20261020.md` | 工作区迁移记录 |
| `review_spectral_moore.md` / `.tex` | 谱 Moore 问题的文献复核笔记 |

## 五、两条研究线的当前状态

| 线 | 目标 | 状态 | 位置 |
|---|---|---|---|
| **主线（毕业论文）** | Λ₃ = λ₁+λ₂+λ₃ ≤ c·n − 3 的精确上界 | 已证 c = 1.284673255723918505，改进 STW 的 (3+√5)/4 约 1.86%；最优系数与极图未确定。**2026-11-03 第二次查新：主定理未被抢先；但下界常数 c\*=(9+√5+2√35)/18 与 K₁∨C₅ 模板已被 Zenodo 自存档 doi:10.5281/zenodo.21785906（2026-08-04）以条件性方式提出，论文已补注 2.3 与 4 条引用如实披露** | `paper/`（预印本）与 `thesis_main/`（毕业论文成稿）；查新报告见 `paper/reports/查新复核报告_前三大特征值上界_20261103/` |
| **副线（增页数的章节）** | λ₁+2λ₂ 的极图分类 | 正在撰写独立章节，已完成 §1–§10（问题与常数、谱区域、三候选归约、全阶夹逼、图核唯一性、稳定性、算术与 Pell 平局、n≤12 精确分类、阶梯图类全阶封闭、单点扩张族无条件封闭）；尚缺 §11 剩余缺口与 §12 小结。**该目录被 .gitignore 排除，成稿后需 `git add -f` 入库** | `thesis_spectral_extremes/chapter_F2/ch_F2_classification.tex` |

论文在线链接汇总见 `paper/参考文献链接.md`。