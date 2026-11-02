import json, os
goal = {
  "title": "图前三大特征值之和的上界：常数改进与余量分析",
  "goal": "向导师汇报毕业论文主线进展：在 STW 命题 6.1 基础上把前三大特征值之和的上界常数从 (3+√5)/4 降到 c=1.284673，说明两条证明主线、独立审核结论、松驰归属、框架天花板与 P0-P3 优先级",
  "audience": "指导教师 / 论文答辩委员会",
  "owner": "谱极值图论课题组",
  "randomSeed": "top3-bound-20261102-k7q",
  "pageCount": 15,
  "themePack": "theme07",
  "slides": [
    {"layout": "theme07_page005", "props": {
      "eyebrow": "谱极值图论 · 毕业论文",
      "reportTag": "中期汇报 · 2026.11",
      "titleL1": "图前三大特征值之和",
      "titleL2Em": "上界常数 c",
      "titleL2Rest": "的改进与余量分析",
      "brandLabel": "毕业论文中期汇报",
      "brandSub": "导师审阅版",
      "sub": "用联合矩恒等式把 (3+√5)/4 降到 c=1.284673",
      "thesis": "主定理已证并通过独立审核，两处承重缺口已修补",
      "closing": "较 STW 常数改进 0.024344",
      "featureNumber": "1.284673",
      "featureLabel": "本文上界常数 c",
      "specs": [["主定理", "已证 · 独立审核"], ["改进幅度", "0.024344"], ["框架天花板", "1.284017"], ["穷举反例", "0 / 33864 图"]]
    }},
    {"layout": "theme07_page006", "props": {
      "eyebrow": "本页速览",
      "title": "摘要",
      "sub": "问题 · 方法 · 结果 · 优先级",
      "lead": ["对 n 阶简单图，前三大特征值之和满足 ", "Λ₃ ≤ 1.284673 n − 3", "，把 STW 命题 6.1 的系数 ", "(3+√5)/4", " 降低了 0.024344，且已通过独立审核。"],
      "keywords": ["联合矩证书", "有限支撑变分", "分离观测", "符号证书复验", "余量量化"],
      "closing": "下一步：证小分支为空，把 c 推到 1.284017",
      "chartTitle": "冠军松驰来源占比",
      "shareLabel": "占比",
      "tracks": [
        {"label": "径向迹项 r_tr", "pct": 61.5},
        {"label": "四次型 Q̂", "pct": 26.6},
        {"label": "正则项 ℓ·r_reg", "pct": 7.4},
        {"label": "R₀ 项", "pct": 3.7},
        {"label": "离对角 r_off/9", "pct": 0.8}
      ],
      "metrics": [
        {"label": "本文常数 c", "lead": "1.284673", "tail": "", "unit": "已证"},
        {"label": "较 STW 改进", "lead": "0.024344", "tail": "", "unit": "绝对值"},
        {"label": "框架天花板", "lead": "1.284017", "tail": "", "unit": "(29+10√5)/40"},
        {"label": "数值上确界", "lead": "1.281568", "tail": "", "unit": "ρ₀"}
      ],
      "cardCount": 4
    }},
    {"layout": "theme07_page007", "props": {
      "eyebrow": "汇报结构",
      "title": "汇报路径",
      "sub": "从问题设定到后续优先级",
      "lead": "先交代问题与 STW 基准，再给出主定理与两条证明主线，随后汇报独立审核结论、余量量化与框架天花板，最后落到 P0-P3 优先级。",
      "closing": "先建立框架，再看数字与判断",
      "chapters": [
        {"no": "01", "zh": "问题与 STW 基准", "en": "Setup & Baseline"},
        {"no": "02", "zh": "主定理与常数", "en": "Main Theorem"},
        {"no": "03", "zh": "两条证明主线", "en": "Two Proof Lines"},
        {"no": "04", "zh": "独立审核结论", "en": "Independent Audit"},
        {"no": "05", "zh": "余量与天花板", "en": "Slack & Ceiling"},
        {"no": "06", "zh": "优先级与后续", "en": "Priorities"}
      ],
      "cardCount": 6
    }},
    {"layout": "theme07_page008", "props": {
      "eyebrow": "方法总览",
      "title": "两条主线",
      "sub": "结构线收敛，证书线闭合",
      "lead": "结构线把无穷维变分压到有限支撑的原子测度；证书线用一个联合矩恒等式给出全局不等式，两条线在冠军测度处对接取值。",
      "closing": "全程用推理，不对每阶图跑算法",
      "layers": [
        {"no": "L1", "zh": "结构线", "en": "Structure", "desc": "变分 → 有限支撑 → 分离观测", "q": "支撑有多大？"},
        {"no": "L2", "zh": "证书线", "en": "Certificate", "desc": "联合矩恒等式 + 带符号四次式", "q": "八项都非负吗？"},
        {"no": "L3", "zh": "对接", "en": "Matching", "desc": "冠军测度处证书取值恰为 c", "q": "唯一吗？"}
      ],
      "steps": [
        {"label": "归一化", "detail": "把 Λ₃ 写成三投影的联合矩"},
        {"label": "变分", "detail": "对径向分布求上确界，导出极值方程"},
        {"label": "离散化", "detail": "证极值测度支撑有限，原子数 ≤24"},
        {"label": "证书", "detail": "构造八项非负分解，逐项给闭式常数"}
      ],
      "cardCount": 3
    }},
    {"layout": "theme07_page025", "props": {
      "eyebrow": "主定理",
      "title": "上界常数",
      "sub": "对一切 n 阶简单图成立",
      "numLead": "1.284673",
      "numTail": "",
      "numUnit": "Λ₃ ≤ c n − 3",
      "numCaption": "本文定理 1.2 的常数",
      "note": "STW 为 (3+√5)/4=1.309017，降低 0.024344",
      "closing": "−3 的位移由小阶穷举坐实",
      "aux": [
        {"label": "较 STW 改进", "value": "0.024344", "unit": "绝对值"},
        {"label": "框架天花板", "value": "1.284017", "unit": "(29+10√5)/40"},
        {"label": "数值上确界 ρ₀", "value": "1.281568", "unit": "(9+√5+2√35)/18"}
      ],
      "auxCount": 3
    }},
    {"layout": "theme07_page048", "props": {
      "eyebrow": "冠军松驰构成",
      "title": "八项分解",
      "titleTail": "谁在吃掉余量",
      "sub": "证书在冠军测度处的松驰按项归一化",
      "lead": "松驰集中在径向迹项与四次型两块，合计约 88%；其余六项不到 12%。要继续收紧必须动主项，微调末项几乎无收益。",
      "closing": "先攻 r_tr，再重拟合 Q̂",
      "panelTitle": "松驰来源占比 · 归一化",
      "stats": [
        {"value": "1.284673", "unit": "c", "label": "本文已证常数"},
        {"value": "0.024344", "unit": "Δ", "label": "较 STW 改进"}
      ],
      "components": [
        {"name": "径向迹项", "en": "r_tr", "pct": 61.5, "note": "主项 −E[RR′z_rad·g]，源于 9J−V>0"},
        {"name": "四次型", "en": "Q̂", "pct": 26.6, "note": "其中 Q3 系数松 14.98 倍"},
        {"name": "正则项", "en": "ℓ·r_reg", "pct": 7.4, "note": "与径向方差耦合，可部分吸收"},
        {"name": "其余五项", "en": "R₀ · r_off", "pct": 4.5, "note": "R₀ 3.68 · r_off/9 0.78"}
      ],
      "componentCount": 4
    }},
    {"layout": "theme07_page023", "props": {
      "eyebrow": "常数阶梯",
      "title": "STW → ρ₀",
      "sub": "总差 0.027449 的三段归属",
      "lead": "每一段都对应一个明确的数学任务，不是数值调参",
      "closing": "已证的只有第一段，其余两段是 P1 与 P2",
      "chartLabel": "常数下降量 · 绝对值",
      "totalLabel": "STW − ρ₀ 总差",
      "totalUnit": "0.027449",
      "steps": [
        {"name": "已证改进", "value": 0.024344, "note": "本文主定理，(3+√5)/4 → c"},
        {"name": "小分支回收", "value": 0.000656, "note": "P1：c → 框架天花板 1.284017"},
        {"name": "框架外余量", "value": 0.002449, "note": "P2：天花板 → ρ₀=1.281568"}
      ],
      "stepCount": 3
    }},
    {"layout": "theme07_page059", "props": {
      "eyebrow": "第二部分",
      "index": "02",
      "title": "审核与余量",
      "sub": "从独立复核到可回收空间",
      "lead": "这一部分交代独立审核的结论与缺口处置，量化冠军处松驰的归属，并给出本框架的天花板。",
      "closing": "先看清余量，再谈优化",
      "keywords": ["独立审核", "G1 / G2 修补", "松驰量化", "数值证据", "天花板"],
      "keywordCount": 5
    }},
    {"layout": "theme07_page060", "props": {
      "eyebrow": "独立审核",
      "marker": "33 页审核报告",
      "segment": "结论 · 缺口处置",
      "title": "主定理成立",
      "titleTail": "两处缺口已修补",
      "lead": "审核覆盖全部定理、引理与 16 张符号证书，并交叉复核穷举与随机实验；主定理证明有效，发现两处承重缺口 G1、G2，均已在论文正文修补并重新编译复验。",
      "statLine": "G1 配对论证 · G2 常值性 已补",
      "closing": "R1 文献未读为残余风险",
      "funnelTitle": "复核 → 定位 → 修补",
      "stages": [
        {"name": "符号证书复核", "value": 16, "note": "残差全部 ≤1.7e-16"},
        {"name": "承重缺口", "value": 2, "note": "引理 2.6 / 引理 4.2", "conv": "已定位"},
        {"name": "修补完成", "value": 2, "note": "写入正文并复验编译", "conv": "100%"}
      ],
      "metrics": [
        {"lbl": "穷举反例", "val": "0", "unit": "n≤6 共 33864 图"},
        {"lbl": "证书残差", "val": "1.7e-16", "unit": "机器精度上界"},
        {"lbl": "随机律反例", "val": "0", "unit": "4000 组分布"},
        {"lbl": "残余风险", "val": "1", "unit": "Ebrahimi 2008 未读"}
      ],
      "stageCount": 3,
      "metricCount": 4
    }},
    {"layout": "theme07_page070", "props": {
      "eyebrow": "数值与机器验证",
      "marker": "证据口径",
      "segment": "验证 · 数据来源",
      "title": "数值证据清单",
      "titleTail": "验证 · 完备度评级",
      "lead": "所有实验都只做交叉验证，不作为证明依据；完备度 100 表示穷举或完备性论证，90 表示抽样或离散网格。",
      "statLine": "全部实验零反例 · 最小余量 +0.0224818",
      "anchorValue": "0",
      "anchorLabel": "全部实验的反例总数",
      "closing": "口径透明，结论才可复核",
      "colHeads": {"dim": "验证项", "level": "完备度", "note": "口径与结论"},
      "rows": [
        {"dim": "小阶穷举", "level": 100, "note": "n=3..6 全图 33864 个，零反例"},
        {"dim": "符号证书", "level": 100, "note": "16 张证书残差 ≤1.7e-16"},
        {"dim": "闭式常数", "level": 100, "note": "c、ρ₀、天花板三值独立复算一致"},
        {"dim": "k 网格全局最大", "level": 90, "note": "k=6..20 上确界恰为 ρ₀，从不超过"},
        {"dim": "吹起族渐近", "level": 90, "note": "n=126..2016，Λ₃/n → ρ₀"},
        {"dim": "随机律抽样", "level": 90, "note": "4000 组分布，零反例"}
      ],
      "rowCount": 6,
      "flowTitle": "验证流程",
      "flow": ["构造符号证书", "高精度复验", "小阶穷举", "结构化吹起", "随机律抽样"]
    }},
    {"layout": "theme07_page062", "props": {
      "eyebrow": "框架天花板",
      "title": "还能降多少",
      "sub": "本方法内可达到的最优常数",
      "kicker": "小分支为空",
      "numLead": "1.284017",
      "numTail": "",
      "numUnit": "(29+10√5)/40",
      "numCaption": "较当前 c 再降 0.000656",
      "ceilingLabel": "框架最优 c",
      "gaugeCap": "STW → ρ₀ 总差归属",
      "gaugeAnnotation": "占 0.027449 比例",
      "gauge": [
        {"label": "已证改进", "value": 88.7},
        {"label": "P1 小分支", "value": 2.4},
        {"label": "P2 框架外", "value": 8.9}
      ],
      "note": "达到天花板等价于证小分支为空",
      "closing": "P2 需八项同时取零",
      "aux": [
        {"label": "当前 c", "value": "1.284673", "unit": "本文已证"},
        {"label": "框架天花板", "value": "1.284017", "unit": "(29+10√5)/40"},
        {"label": "数值上确界", "value": "1.281568", "unit": "ρ₀ 猜想值"}
      ],
      "auxCount": 3
    }},
    {"layout": "theme07_page012", "props": {
      "eyebrow": "优先级矩阵",
      "title": "难度 × 收益",
      "sub": "四条后续路线的定位",
      "lead": "横轴是证明难度，纵轴是预期收益；左上是最该先做的一步，右下两条路线已经论证过收益极低，明确排除。",
      "closing": "先做左上，再评估右上",
      "axisX": "证明难度",
      "axisY": "预期收益",
      "axisYNote": "可回收松驰",
      "axisXNote": "所需新工具",
      "quadrants": [
        {"name": "P1 小分支为空", "pos": "低难度 · 高收益", "dirs": ["保留 q² 符号", "径向方差互斥"], "col": 1, "row": 1},
        {"name": "P2 精确值 ρ₀", "pos": "高难度 · 高收益", "dirs": ["八项全零证书", "24 原子降到 6"], "col": 2, "row": 1},
        {"name": "P3 吹起渐近", "pos": "低难度 · 低收益", "dirs": ["Weyl 误差 O(1/n)", "坐实 −3 位移"], "col": 1, "row": 2},
        {"name": "已排除路线", "pos": "高难度 · 低收益", "dirs": ["逐阶算法搜索", "阈值微调 0.16%"], "col": 2, "row": 2}
      ]
    }},
    {"layout": "theme07_page064", "props": {
      "eyebrow": "后续优先级",
      "segment": "路线 · P0-P3",
      "title": "先全局后局部",
      "titleTail": "路线 · 排序结果",
      "lead": "P0 已完成论文加固；P1 攻小分支为空以触及框架天花板；P2 追求精确值，需先降原子数；P3 补严格渐近，风险最低但收益也最小。",
      "closing": "以数学推理为主，算法只作交叉验证",
      "criteriaTitle": "排序依据",
      "criteria": [
        {"name": "数学优先", "note": "能推理就不靠逐阶算法"},
        {"name": "收益难度比", "note": "优先可回收松驰大的台阶"},
        {"name": "承重性", "note": "先补影响主定理的缺口"},
        {"name": "可验证性", "note": "每一步都要有机器证书"}
      ],
      "cardsTitle": "四个优先级 · P0-P3",
      "cards": [
        {"name": "P0 论文加固", "en": "Done", "note": "G1/G2 修补 + 附录 A + 记号统一", "certainty": 100},
        {"name": "P1 小分支为空", "en": "Ceiling", "note": "回收 0.000656，c 降至 1.284017", "certainty": 45},
        {"name": "P2 精确值 ρ₀", "en": "Exact Value", "note": "需八项全零证书，先把 24 原子降到 6", "certainty": 15},
        {"name": "P3 吹起渐近", "en": "Asymptotics", "note": "Weyl 给 O(1/n)，坐实 −3 位移", "certainty": 70}
      ],
      "itemCount": 4,
      "cardCount": 4
    }},
    {"layout": "theme07_page015", "props": {
      "eyebrow": "结论",
      "kicker": "三句话",
      "quoteLead": "上界已从",
      "quoteEm1": "(3+√5)/4",
      "quoteMid": "降到",
      "quoteEm2": "1.284673",
      "quoteTail": "，并完成独立审核。",
      "conclusions": [
        {"idx": "01", "title": "方法", "note": "联合矩证书取代逐项放缩"},
        {"idx": "02", "title": "结果", "note": "改进 0.024344，主定理已证"},
        {"idx": "03", "title": "边界", "note": "框架天花板 1.284017，非终值"}
      ],
      "itemCount": 3,
      "source": "详见 top3_joint_moment_bound.pdf 与独立审核报告"
    }},
    {"layout": "theme07_page067", "props": {
      "eyebrow": "结语",
      "kicker": "下一步",
      "quoteLead": "主攻",
      "quoteEm": "小分支为空",
      "quoteTail": "，副线补严格渐近。",
      "tagline": "以数学推理为主，算法只作交叉验证。",
      "brandLabel": "毕业论文中期汇报",
      "brandSub": "欢迎导师指正"
    }}
  ]
}
out = os.path.join("汇报_前三大特征值上界_20261102", "deck", "goal.json")
with open(out, "w", encoding="utf-8") as fh:
    json.dump(goal, fh, ensure_ascii=False, indent=2)
print("written", out, os.path.getsize(out), "slides", len(goal["slides"]))
