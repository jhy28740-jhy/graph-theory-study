# 发表状态核查报告（2026-11-09）

核查对象：`paper/top3_joint_moment_bound.pdf` ——《A joint moment bound for the three largest
graph eigenvalues》(Haoyang Ji, 手稿日期 2026-11-03)。

## 结论

- **截至 2026-11-09（Asia/Shanghai），该结果尚未公开发表，也未上任何预印本平台。**
  arXiv、Crossref、OpenAlex、Zenodo 四个结构化索引均无匹配记录；手稿 PDF 本身不带
  arXiv 戳或投稿回执。
- 唯一的先行公开记录是 Tony Newton（Newton Astro Labs, London）的两条 Zenodo 自存档：
  doi:10.5281/zenodo.21785906（2026-08-04）与 doi:10.5281/zenodo.22082842（2026-08-24）。
  它们宣布同一常数 c_w = (9+√5+2√35)/18 = 1.281568196872168… 与上界 Λ₃(G) ≤ c_w·n − 3，
  但**按其自述为条件性**（最终正性层为机器证书，接受须以独立审计为前提），
  不构成发表，也不覆盖本论文的无条件系数界与最大化子结构定理。
- 本结论与 2026-11-04 的本地《查新复核报告_top3与lam12lam2.pdf》一致：该报告独立核查过
  arXiv math.CO、OpenAlex、Zenodo、OSF，结论为"没有任何工作改进过 STW (arXiv:2608.02429)
  Proposition 6.1 的系数 (3+√5)/4"。

## 检索明细（2026-11-09 实测）

| 渠道 | 检索式 | 结果 |
|---|---|---|
| arXiv API | `ti:"joint moment bound"` | 0 条 |
| arXiv API | `all:"three largest graph eigenvalues"` | 0 条 |
| arXiv API | `au:"Haoyang Ji"` | 6 条，均为动力系统/粒子物理方向（同名不同人），无谱图论文 |
| arXiv API（对照） | `ti:"maximum spectral sum"` | 命中 Kumar–Liu–Monterde–Pragada–Tait（arXiv:2604.00512），证明索引覆盖 2026 年文献、检索式有效 |
| Crossref | `query.bibliographic=joint moment bound three largest graph eigenvalues` | 前 10 条无匹配 |
| OpenAlex | `search="joint moment bound"` | count=0（另一次宽检索 count=2，命中均无关） |
| Zenodo | `q="joint moment bound"` | total=0 |
| Semantic Scholar | paper/search | HTTP 429 限流，未获结果（缺失路） |
| DuckDuckGo / Bing HTML | 精确标题 | 反爬拦截，无响应（缺失路） |

## 含义与建议

- 新颖性窗口仍然开放：已发表/预印本文献中最强上界仍是 STW 的 (3+√5)/4 ≈ 1.3090170；
  本论文的无条件系数 c = 1.2846732557… 没有任何公开竞争记录（Newton 的条件性宣称除外）。
- 若准备投稿或上 arXiv，当前无撞车风险；优先权覆盖按论文 Remark 2.3 的写法引用 Newton
  两条 Zenodo 记录即可。
- 两路缺失（S2 限流、搜索引擎反爬）不影响结论：四个结构化数据库已覆盖数学预印本与
  期刊/学位论文的主要索引面。

## 复现命令

```
curl -s "https://export.arxiv.org/api/query?search_query=ti:%22joint+moment+bound%22&max_results=15"
curl -s "https://export.arxiv.org/api/query?search_query=all:%22three+largest+graph+eigenvalues%22&max_results=15"
curl -s "https://export.arxiv.org/api/query?search_query=au:%22Haoyang+Ji%22&max_results=50"
curl -s "https://api.crossref.org/works?query.bibliographic=joint+moment+bound+three+largest+graph+eigenvalues&rows=10"
curl -s "https://api.openalex.org/works?search=joint%20moment%20bound&per-page=5"
curl -s "https://zenodo.org/api/records?q=%22joint%20moment%20bound%22&size=5"
```