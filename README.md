# 图论学习 · 代数图论课程资料

本仓库收录《代数图论》(Algebraic Graph Theory) 课程的整理资料，源文档为 Markdown，同时提供排版好的 PDF 便于打印与投影。

## 目录

| 文件 | 说明 | 规模 |
| --- | --- | --- |
| `代数图论课程知识点梳理.md` / `.pdf` | 全书 13 章知识点梳理：每章按「定义 → 命题/定理 → 证明思路 → 与前后的连贯性」展开，并给出每章小结 | 约 43 页 PDF |
| `第8章课程讲义_强正则图与设计.md` / `.pdf` | 第 8 章授课讲义（12 学时）：6 讲 + 习题课，含时间分配、板书设计、作业、Q&A 预演 | 约 58 页 PDF |
| `_ch8_verify.py` | 讲义附录的数值验证脚本，可复现讲义中全部「实测数据」 | — |

## 一、知识点梳理（13 章）

**Part I 谱图论**
1. 图论基本概念
2. 矩阵理论工具
3. 图的关联矩阵
4. 谱界
5. 谱极值问题
6. 谱图论近期进展选讲

**Part II 代数图论**
7. 对称性与 Cayley 图
8. 强正则图与设计
9. 特征标与有限 Fourier 分析
10. 代数图的谱计算
11. 有限域与有限几何
12. 二次特征与 Paley 图
13. 割、流与临界群

文末附「全书脉络与教学建议」。

## 二、第 8 章讲义结构

| 讲次 | 主题 | 对应教材 |
| --- | --- | --- |
| 第一讲 | 强正则图：把「正则」升级到「一对点一致」 | §8.1 |
| 第二讲 | 邻接代数：一个矩阵恒等式统治全章 | §8.2 |
| 第三讲 | 设计图与部分设计图：二部版本 | §8.3 |
| 第四讲 | 距离正则图：把一致性推到所有距离层 | §8.4 |
| 第五讲 | 区组设计与关联矩阵 | §8.5 |
| 第六讲 | 全章总结与习题课（题组 A–E） | — |

另含：常见疑问与课堂 Q&A 预演、板书总模板与 12 学时时间分配表、延伸阅读、数值验证脚本附录。

## 三、复现 PDF

PDF 由 `markdown-it` + `markdown-it-texmath`（KaTeX 服务端渲染公式）+ `mermaid` + Playwright/Chromium 生成：

```bash
mkdir _pdfbuild && cd _pdfbuild
npm i markdown-it markdown-it-texmath katex mermaid playwright
# 将 build.mjs 与 vendor/ 放入本目录后
node build.mjs
```

数学公式为服务端 KaTeX 渲染，图表为 mermaid v11（注意：`dist/mermaid.min.js` 的 API 位于
`window.__esbuild_esm_mermaid_nm.mermaid.default`，而非 `window.mermaid`）。

## 四、运行验证脚本

```bash
python _ch8_verify.py
```

输出强正则图参数检验、拉丁方图 / Shrikhande 图、Petersen 与 Tutte–Coxeter 图、超立方体距离正则性、
区组设计关联矩阵、Fano 平面等全部实测数据，与讲义正文一一对应。
