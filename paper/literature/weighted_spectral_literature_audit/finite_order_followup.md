# 固定 k>1 的有限阶极图与查新记录

检索日期：2026-10-15（Asia/Shanghai）。目标为

    M_k(n)=max_{|V(G)|=n}{lambda_1(G)+k lambda_2(G)}, k>1.

特征值按邻接矩阵的代数降序排列。下文严格区分一般简单图、树、渐近结果、最大值和最小值；“未找到”只表示本次实际访问并核对的来源内未找到。

## 直接涉及所有简单图

Roland Paulin, *Quadratic inequalities between the largest eigenvalues of a graph*, [arXiv:2609.37980v1](https://arxiv.org/abs/2609.37980v1), Corollary 2.1。令 x=(lambda_1+1)/n、y=(lambda_2+1)/n，则所有有限简单图满足

    0 <= y <= x,  x^2 + xy + 2y^2 <= x+y,

且该区域正是图序列的渐近闭包。它决定固定 k>=1 的支持函数，但没有给出每个固定 n 的整数优化、全部等号图或稳定性。

由该区域的一元支持函数计算

    C_k = [k+3+2 sqrt(2(k^2-k+2))]/7,  1<=k<=3;
    C_k = (k+1)/2,                              k>=3.

因此每个有限简单图满足 lambda_1+k lambda_2 <= C_k n-(k+1)，并且 lim M_k(n)/n=C_k。特别地 C_2=(5+4 sqrt(2))/7。这是 k=2 的最佳渐近系数和统一上界，不是 M_2(n) 对每个 n 的精确公式。

## k=1 的一般图有限阶精确结果

Jing Huang, Wei Wei, *The Exact Maximum of the Spectral Sum of Graphs*, [arXiv:2607.23081v1](https://arxiv.org/abs/2607.23081v1)，Theorem 1.1 声称对每个 n>=5，lambda_1+lambda_2 的唯一极图为

    K_n^*=K_(n-s_n) join (K_ceil(s_n/2) disjoint-union K_floor(s_n/2)),

其中 s_n 按 n mod 7 分段确定；并给出 M_1(n)<=8n/7-2，且 7|n 时取等。它只处理等权 k=1，不能扩展成 k>1 的结果；这是预印本结论。

较早的 Hitesh Kumar、Lele Liu、Hermie Monterde、Shivaramakrishna Pragada、Michael Tait，*Maximum spectral sum of graphs*, [arXiv:2604.00512v2](https://arxiv.org/abs/2604.00512v2)，Theorem 1.2，只证明 lambda_1+lambda_2<=8n/7；有限阶同余分类和唯一性是 Huang--Wei 的后续预印本结论。

## 题名相近但研究最小值

Shaowei Sun, Yaping Min, Kinkar Chandra Das, *Extremal graphs for the sum of two largest eigenvalues*, AIMS Mathematics 11(5) (2026), 15028-15036，DOI [10.3934/math.2026618](https://doi.org/10.3934/math.2026618)。[期刊页](https://www.aimspress.com/article/doi/10.3934/math.2026618)，[正式 PDF](https://www.aimspress.com/aimspress-data/math/2026/5/PDF/math-11-05-618.pdf)。摘要以及 Theorem 3.1、3.2（PDF 第6-7页；期刊页15033-15034）处理连通大阶图的最小 lambda_1+lambda_2：n>=467 时路径 P_n 唯一最小化该和。它不是一般图最大值，也没有 lambda_1+2lambda_2 定理。

## 已有的树类固定权重结果

Hitesh Kumar, Bojan Mohar, Shivaramakrishna Pragada, Hanmeng Zhan, *Convex combination of first and second eigenvalues of trees*, [arXiv:2601.10036v2](https://arxiv.org/abs/2601.10036v2)，Theorem 1.4。目标为 Psi(T,alpha)=alpha lambda_1(T)+(1-alpha)lambda_2(T)，故 lambda_1+k lambda_2=(k+1)Psi(T,1/(k+1))。该定理只对树、固定 0<alpha<1/2、充分大 n 给出唯一极值树：奇数 n 为平衡长度3双彗星；偶数 n 在 alpha_0=(sqrt(5)-1)/(2sqrt(5)) 两侧分别为长度4平衡双彗星和长度3双彗星。取 alpha=1/3（k=2）时，树类已有充分大阶唯一极图；这不能推广为一般简单图结果。Theorem 1.7 也只是树类归一化渐近结论。

## 相近但不直接解决的来源

Vladimir Nikiforov, *Linear combinations of graph eigenvalues*, EJL A 15 (2006), 329-336，[DOI](https://doi.org/10.13001/1081-3810.1242)：保证固定极端特征值线性组合的归一化最大值极限存在，并非 M_k(n) 的有限阶精确分类。

Bojan Mohar, *On the sum of k largest eigenvalues of graphs and symmetric matrices*, JCTB 99 (2009), 306-313，[DOI](https://doi.org/10.1016/j.jctb.2008.07.001)：这里的 k 是前 k 个特征值求和的项数，目标是 lambda_1+...+lambda_k，不是 lambda_1+k lambda_2。

带有 A_alpha、signless Laplacian、distance eigenvalue 或 weighted matrix 的命中项均不是同一邻接矩阵上的 lambda_1+k lambda_2 一般图有限阶最大值，未列为直接结果。

## 定向检索与限制

实际查询了 arXiv、OpenAlex、Crossref，关键词包括 finite order、exact maximum、extremal graph、uniqueness、stability、weighted largest two eigenvalues、lambda_1+2lambda_2、convex combination of first and second eigenvalues；原始响应保存在 weighted_spectral_literature_audit/raw/、s13_stage12/literature_search_raw/ 及 agent_findings.md。

Google Scholar 请求超时，Semantic Scholar 返回 HTTP 429，zbMATH 返回 403，MathSciNet 公共入口需要机构登录；因此没有把检索失败写成不存在相关工作。

**结论。** 一般图 k=1 的有限阶精确最大值已有 Huang--Wei；任意固定 k>=1 的渐近支持函数已有 Paulin；树类 k=2 的充分大阶极图已有 Kumar--Mohar--Pragada--Zhan。除此之外，本次没有找到一般简单图固定 k=2 或 k>1 的有限阶精确最大值、极图唯一性或近极值稳定性直接论文。这个方向可作为候选研究问题，但开题时应保留“本次检索未找到”的限定。
