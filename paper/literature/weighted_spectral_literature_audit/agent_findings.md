# λ₁+kλ₂ 文献尽调：独立核验记录

检索与核验日期：2026 年 10 月 14 日（Asia/Shanghai）。研究量均按邻接特征值的代数降序理解。本文没有把关键词未命中解释为不存在文献；对最新预印本只报告已核对的陈述，没有声称逐行审计全部证明。

## 1. 最重要的选题判断

一般简单图上，任意固定正权重的渐近常数已由 Roland Paulin 的完整二维谱区域决定。因而“求 λ₁+2λ₂ 或一般 λ₁+kλ₂ 的渐近最优系数”不能作为尚未解决的新方向。

λ₁+λ₂ 的有限阶精确极值已有 Huang–Wei 的 2026 年预印本。树类的一般正权组合也已有 Kumar–Mohar–Pragada–Zhan 的直接论文，且 k=2 的充分大阶精确唯一极图已被覆盖。

在本次实际取得并检查的来源中，尚未找到“所有简单图上固定 k>1 的每个有限 n 精确极值和全部极图”或其结构稳定性定理的直接论文。但是 Google Scholar 不可达、Semantic Scholar 限流，且普通搜索引擎的公式索引很弱，所以只能写“本次未找到”，不能宣称已经证明是文献空白。

## 2. 直接决定任意 k 渐近常数的来源

Roland Paulin, *Quadratic inequalities between the largest eigenvalues of a graph*, arXiv:2609.37980v1, 2026-09-29。

- [摘要](https://arxiv.org/abs/2609.37980v1)
- [全文](https://arxiv.org/html/2609.37980v1)
- [PDF](https://arxiv.org/pdf/2609.37980v1)
- 本地已核对全文：`s13_checkpoint/sources/2609.37980.full.txt`。

原文 Corollary 2.1 定义 x=(λ₁+1)/n、y=(λ₂+1)/n，并给出闭包及极限区域的精确等式

    R₁,₂=S₁,₂={(x,y): 0≤y≤x, x²+xy+2y²≤x+y}。

因此，对 k≥0，以下公式是已知区域的直接支持函数计算，而非新的图论定理：

    h(k)=[k+3+2√(2k²−2k+4)]/7，0≤k≤3；
    h(k)=(k+1)/2，k≥3。

每个 n≥2 的简单图满足

    λ₁+kλ₂ ≤ n h(k)−(k+1)，

且归一化极大值的极限恰为 h(k)。特别地，k=2 的最优渐近系数为

    h(2)=(5+4√2)/7，
    λ₁+2λ₂ ≤ [(5+4√2)/7]n−3。

推导复核：椭圆中心为 (3/7,1/7)，二次型矩阵为 [[1,1/2],[1/2,2]]。全椭圆在方向 (1,k) 的支持值为第一式；其支撑点满足 x≥y 当且仅当 0≤k≤3。在 k≥3 分支，约束 y≤x 截断后的支撑点为 (1/2,1/2)。闭包精确刻画保证渐近可达，但不自动给出任意固定 n 的整数参数最优图、全部等号图或结构稳定性。

## 3. 不加权的有限 n 最大值已存在直接结果

Jing Huang, Wei Wei, *The Exact Maximum of the Spectral Sum of Graphs*, arXiv:2607.23081v1, 2026-07-25。

- [摘要](https://arxiv.org/abs/2607.23081v1)
- [全文](https://arxiv.org/html/2607.23081v1)
- [PDF](https://arxiv.org/pdf/2607.23081v1)
- 本地已核对全文：`s13_checkpoint/sources/2607.23081.full.txt`。

Theorem 1.1 明确声称：每个 n≥5，λ₁+λ₂ 的唯一最大图为

    K*n = K(n,ceil(s_n/2),floor(s_n/2))，
    K(n,p,q)=K_(n−p−q) ∨ (K_p ⊔ K_q)。

若 n=7a+r，r=0,…,6，s_n 依次为 4a、4a、4a+1、4a+2、4a+2、4a+3、4a+4。原文还给出 S₂(K*n)≤8n/7−2，等号恰在 7|n。它处理的是 k=1；全文中的 weighted Ferrers quotients 指证明所用商矩阵，不是宣称覆盖任意 λ₁+kλ₂。

较早来源：Hitesh Kumar, Lele Liu, Hermie Monterde, Shivaramakrishna Pragada, Michael Tait, *Maximum spectral sum of graphs*, [arXiv:2604.00512v2](https://arxiv.org/abs/2604.00512v2)，v2 2026-05-03。其 Theorem 1.2 是 λ₁+λ₂≤8n/7；不能把 −2 修正或有限阶唯一性归给该文。其 Conjecture 5.1 是上述有限阶问题。

## 4. 标题相近但实际为极小值的期刊文献

Shaowei Sun, Yaping Min, Kinkar Chandra Das, *Extremal graphs for the sum of two largest eigenvalues*, **AIMS Mathematics** 11(5) (2026), 15028–15036，DOI 10.3934/math.2026618；正式发表日期 2026-05-28。

- [DOI](https://doi.org/10.3934/math.2026618)
- [期刊页面](https://www.aimspress.com/article/doi/10.3934/math.2026618)
- [正式 PDF](https://www.aimspress.com/aimspress-data/math/2026/5/PDF/math-11-05-618.pdf)
- [Crossref 元数据](https://api.crossref.org/works/10.3934/math.2026618)

已通过 Crossref 和期刊页核对标题、三位作者、卷期页码及摘要，并由独立子核验查正式 PDF：

- Theorem 3.1：期刊页 15033（PDF 第 6 页），连通 n≥467 图的 S₂<4 恰对应文中四族 Y_n、W_n、P_n、C_n。
- Theorem 3.2：期刊页 15034（PDF 第 7 页），连通 n≥467 时 S₂(G)≥S₂(P_n)，等号当且仅当 G≅P_n。

期刊摘要的原句为 “the path graph **minimizes** this sum among all connected graphs of order n≥467”。所以该文不解决所有图的最大 λ₁+λ₂，更不解决 k>1 的最大加权谱和。不能仅见标题 “Extremal” 就认定它是最大值论文。

证据保存于 `s13_stage12/literature_search_raw/aims_2026618.txt`、`aims_2026618.result.json`、`crossref_2026618.txt`、`crossref_2026618.result.json`。

## 5. 树类已经有直接一般权重论文

Hitesh Kumar, Bojan Mohar, Shivaramakrishna Pragada, Hanmeng Zhan, *Convex combination of first and second eigenvalues of trees*, [arXiv:2601.10036v2](https://arxiv.org/abs/2601.10036v2)，初稿 2026-01-15，v2 2026-08-07。

- [PDF](https://arxiv.org/pdf/2601.10036v2)
- 元数据检索证据：`weighted_spectral_literature_audit/raw/arxiv_query2.json`。

目标函数为 Ψ(T,α)=αλ₁(T)+(1−α)λ₂(T)。任意固定 k>1 可写为 λ₁+kλ₂=(k+1)Ψ(T,1/(k+1))，因此是直接相同的加权问题，限制是图类为树。

独立核读原文 Theorem 1.4（PDF 第 2–3 页）：固定 α∈(0,1/2)，对充分大 n（可依赖 α），给出唯一最大树。特别地 k=2 即 α=1/3 落在已解决分支，极图为

    DC(floor((n−3)/2),ceil((n−3)/2),3)。

这里沿用原文 DC 记号：长为 3 个顶点的路径两端分别附加给定数量的叶。由 Proposition 2.3（PDF 第 4 页）及 Theorem 1.4，k=2 的充分大 n 极大值为

    odd n: √((n+1)/2)+2√((n−3)/2)；
    even n: √((n−1+√5)/2)+2√((n−1−√5)/2)。

Theorem 1.7（PDF 第 3 页）还给出固定 k≥1 的树类渐近常数：以 √(n−1) 归一化后极大值趋于 (1+k)/√2。Theorem 1.5 的全 n≥5 精确结果是 α=1/2，即不加权谱和；不能把这一“全 n”直接转移到 k=2。

因此不宜再把“树上的 λ₁+2λ₂ 渐近或充分大 n 极图”提议为未有人做的新方向。有限数量的小阶情形不等于一个未经研究的完整方向。

## 6. 其他真实相关来源及适用边界

Vladimir Nikiforov, *Linear combinations of graph eigenvalues*, **Electronic Journal of Linear Algebra** 15 (2006), 329–336，[DOI](https://doi.org/10.13001/1081-3810.1242)，[arXiv:math/0608198v2](https://arxiv.org/abs/math/0608198v2)。Theorem 1 给固定有限极端特征值线性组合的归一化最大值极限存在性，直接覆盖 λ₁+kλ₂，但不确定任意 k 的极值常数。原文 Proposition 5、6 讨论 λ₁+λ₂ 的早期线性界和构造。不能把“极限存在”写成“极值问题已经精确求解”。

Javad Ebrahimi B., Bojan Mohar, Vladimir Nikiforov, Azhvan Sheikh Ahmady, *On the sum of two largest eigenvalues of a symmetric matrix*, **Linear Algebra and its Applications** 429 (2008), 2781–2787，[DOI](https://doi.org/10.1016/j.laa.2008.06.016)。本次核对 Crossref 书目信息，并由上述 2026 年全文交叉核对其 8/7 猜想及 7a 顶点的 K_(3a)∨(K_(2a)⊔K_(2a)) 构造。ScienceDirect 全文访问受挑战页限制，本次没有冒称独立核读全文全部定理。

Bojan Mohar, *On the sum of k largest eigenvalues of graphs and symmetric matrices*, **Journal of Combinatorial Theory, Series B** 99 (2009), 306–313，[DOI](https://doi.org/10.1016/j.jctb.2008.07.001)。这里 k 是求和项数，不能与“第二特征值的权重 k”混同。本次仅确认书目和它被新文献用于一般部分谱和背景，不引用未亲自核实的具体定理。

Evans M. Harrell II, Joachim Stubbe, *On sums of graph eigenvalues*, [arXiv:1308.5340](https://arxiv.org/abs/1308.5340)，期刊 DOI [10.1016/j.laa.2014.05.001](https://doi.org/10.1016/j.laa.2014.05.001)。摘要包含多种矩阵的谱和变分估计和格嵌入界，不是所有简单图 k>1 加权前两特征值的有限阶极图定理。

Lele Liu, *Two conjectures in spectral graph theory involving the linear combinations of graph eigenvalues*, [arXiv:2206.03723v1](https://arxiv.org/abs/2206.03723v1)。摘要中的两个问题分别是 λ₁(G)+λ₁(补图) 与 signless-Laplacian Q-spread；均不是同一个邻接矩阵的 λ₁+kλ₂。

Shaowei Sun, Mengyao Guo, Hongyan Ge, Kinkar Chandra Das, *Hoffman-type Results for the Sum of k Largest Eigenvalues of a Graph*, [arXiv:2609.23707](https://arxiv.org/abs/2609.23707)。独立摘要/全文筛查确认它与小谱和、极小问题有关，本次未将其当作加权前两特征值最大值解。

## 7. 实际检索记录与访问限制

本次实际执行了下列查询，不仅按题名推断：

- arXiv：`all:"weighted sum" AND all:graph AND all:eigenvalues`；`all:"linear combinations" AND all:"graph eigenvalues"`；`all:"lambda_1" AND all:"lambda_2" AND all:"spectral"`；父任务另执行 `all:"spectral sum" AND all:"graph"`。前三条记录在 `s13_stage12/literature_search_raw/arxiv_*.result.json`。公式查询主要命中无关的天体物理文章，表明公式词法检索不足。
- [arXiv 加权查询入口](https://export.arxiv.org/api/query?search_query=all%3A%22weighted%20sum%22%20AND%20all%3Agraph%20AND%20all%3Aeigenvalues&start=0&max_results=30)。
- [arXiv 线性组合查询入口](https://export.arxiv.org/api/query?search_query=all%3A%22linear%20combinations%22%20AND%20all%3A%22graph%20eigenvalues%22&start=0&max_results=30)。
- Crossref：`sum of two largest adjacency eigenvalues graph`、`weighted sum largest eigenvalues adjacency graph`、`weighted linear combination largest second eigenvalue graph stability`，并按 DOI 精确核验。广义查询混入大量 Laplacian、signless Laplacian、distance matrix 结果，均未当作直接命中。
- [Crossref 加权与稳定性查询入口](https://api.crossref.org/works?query.bibliographic=weighted%20linear%20combination%20largest%20second%20eigenvalue%20graph%20stability&rows=30)。
- OpenAlex：`sum of two largest adjacency eigenvalues`、`maximum sum of largest two eigenvalues adjacency graph`、`weighted sum adjacency eigenvalues graph`、`sum of largest eigenvalues graph adjacency`、`spectral sum adjacency graph`、`lambda1 k lambda2 graph`、`lambda1 2lambda2 graph`。前 5 项各取前 50 条并筛选；后两项几乎全无关。原始响应存于 `s13_stage12/literature_search_raw/openalex0.txt` 至 `openalex6.txt`。
- [OpenAlex 查询入口](https://api.openalex.org/works?search=sum%20of%20two%20largest%20adjacency%20eigenvalues&per-page=50)。
- Semantic Scholar：请求 `largest two adjacency eigenvalues graph`，接口返回 HTTP 429 Too Many Requests；未宣称完成该库检索。
- [Semantic Scholar 查询入口](https://api.semanticscholar.org/graph/v1/paper/search?query=largest%20two%20adjacency%20eigenvalues%20graph&limit=10&fields=title,authors,year,url,externalIds,abstract)。
- Google Scholar：两次请求分别为 `"weighted sum" "eigenvalues" graph largest` 和 `"lambda_1" "lambda_2" "graph" "weighted"`，均连接超时（WinError 10060）。失败日志已保存 `google_scholar_weighted.result.json` 与 `google_scholar_formula.result.json`。
- [Google Scholar 加权查询入口](https://scholar.google.com/scholar?q=%22weighted%20sum%22%20%22eigenvalues%22%20graph%20largest)。
- Bing：`"lambda_1" "2lambda_2" graph` 和 `"weighted sum" "largest" "eigenvalues" "graph"` 页面均明确“没有与此相关的结果”。这只是检索结果，不是不存在相关论文的证据。

本轮可支持的谨慎结论是：一般加权渐近问题已有直接决定性来源；不加权有限阶和树类加权问题已有直接结果；一般图 k>1 的有限阶精确分类或结构稳定性在已检查来源内未找到，仍需进一步针对性核验后才能提出明确的新颖性主张。
