# ============ 1 封面 ============
s = prs.slides.add_slide(BLANK)
rect(s, 0, 0, SW, SH, fill=WHITE)
rect(s, 0, 0, SW, 2.05, fill=BLUE)
rect(s, 0, 2.05, SW, 0.06, fill=AMBER)
_, tf = tb(s, 0.85, 0.52, 11.6, 0.4)
line(tf, "谱极值图论 · 硕士学位论文汇报 · 导师审阅版", 14, RGBColor.from_string("C9D3E4"), True, first=True, space_after=0)
_, tf = tb(s, 0.85, 0.98, 11.6, 0.9)
line(tf, "图前三大特征值之和的上界", 40, WHITE, True, first=True, space_after=0)
_, tf = tb(s, 0.85, 2.42, 11.6, 0.5)
line(tf, "常数改进、余量诊断与后续路线：从 STW 的 (3+√5)/4 到 c = 1.284673", 19, BLUE, True, first=True, space_after=0)
_, tf = tb(s, 0.85, 3.02, 11.6, 0.4)
line(tf, "Λ₃(G) = λ₁+λ₂+λ₃ 的显式上界 · 联合矩证书 · 独立审核 · 框架天花板与 P0–P3 路线", 14, GREY, first=True, space_after=0)
x = 0.85
x = chip(s, x, 3.62, "[已证] 主定理 Λ₃ ≤ cn−3", GREEN, GREENL, w=2.62) + 0.18
x = chip(s, x, 3.62, "[已核] 33 页独立审核", BLUE2, BLUEL, w=2.30) + 0.18
x = chip(s, x, 3.62, "[已完成] P0 论文加固", GREEN, GREENL, w=2.30) + 0.18
chip(s, x, 3.62, "[未完成] 精确上界 ρ₀", RED, REDL, w=2.20)
cards = [("1.284673", "本文上界常数 c", "显式代数数，精确算术验证"),
         ("0.024344", "较 STW 的改进", "(3+√5)/4 − c，常数项 −3 不变"),
         ("1.284017", "框架天花板", "(29+10√5)/40，需证小分支为空"),
         ("1.281568", "数值上确界 ρ₀", "(9+√5+2√35)/18，加权轮图下界")]
cw = 2.86
for i, (v, t, sub) in enumerate(cards):
    x0 = 0.85 + i * (cw + 0.20)
    rect(s, x0, 4.35, cw, 1.62, fill=GREYL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.06)
    rect(s, x0, 4.35, 0.07, 1.62, fill=[BLUE, GREEN, AMBER, BLUE2][i])
    _, tf = tb(s, x0 + 0.24, 4.52, cw - 0.4, 0.55)
    line(tf, v, 26, [BLUE, GREEN, AMBER, BLUE2][i], True, first=True, space_after=0)
    _, tf = tb(s, x0 + 0.24, 5.12, cw - 0.4, 0.32)
    line(tf, t, 13.5, INK, True, first=True, space_after=0)
    _, tf = tb(s, x0 + 0.24, 5.46, cw - 0.4, 0.42)
    line(tf, sub, 10.5, GREY, first=True, space_after=0)
_, tf = tb(s, 0.85, 6.32, 11.6, 0.6)
line(tf, "2026 年 11 月 9 日　|　期刊版：top3_joint_moment_bound.tex（22 页纠正版）　|　毕业论文：thesis_main/main.tex（60 页，Master’s）", 11.5, GREY, first=True, space_after=0)

# ============ 2 一页结论 ============
s, y = base("一页结论：做了什么、还差什么", "EXECUTIVE SUMMARY")
_, tf = tb(s, 0.58, y, 7.55, 4.6)
line(tf, "已经拿到手的", 16, BLUE, True, first=True, space_after=6)
line(tf, "• 主定理：Λ₃(G) ≤ cn − 3 对一切 n ≥ 3 成立，c = 1.284673255723918505…，较 STW 命题 6.1 的 (3+√5)/4 降低 0.0243437；等价地，一切二阶矩有限的各向同性律满足 E|X·X′| + |EX|² ≤ 2c。", 13.5, space_after=6, indent=0.22)
line(tf, "• 结构定理：每个极大 graphon 的典范律至多 24 个原子，原核几乎处处是二值阶梯核，正点团至多 6 点。", 13.5, space_after=6, indent=0.22)
line(tf, "• 证明经 33 页独立审核：两处承重跳步（引理 2.6、引理 4.2）已补进正文，七条表述建议 S1–S7 全部落实，新增附录 A 交代常数来源与验证协议。", 13.5, space_after=10, indent=0.22)
line(tf, "还没有拿到手的", 16, RED, True, space_after=6)
line(tf, "• 精确上界 M = ρ₀：需要一个在冠军律上逐项取零的证书，而目前八项无一项为零。", 13.5, space_after=6, indent=0.22)
line(tf, "• 框架天花板 (29+10√5)/40：等价于证明小分支为空，只能回收当前缺口的 21.1%。", 13.5, space_after=6, indent=0.22)
line(tf, "• 下界的严格渐近：ρ₀ 目前只有 graphon 构造与数值吹起族，Weyl 严格化未写。", 13.5, space_after=0, indent=0.22)
rect(s, 8.42, y, 4.36, 4.6, fill=BLUEL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.04)
_, tf = tb(s, 8.68, y + 0.24, 3.86, 4.1)
line(tf, "缺口在哪里", 15, BLUE, True, first=True, space_after=8)
line(tf, "冠军律处证书松驰 6.210×10⁻³，集中在两处：", 12.5, space_after=8)
line(tf, "61.5%", 24, RED, True, space_after=2)
line(tf, "径向迹项 r_tr（其中 87% 来自 9J−V = 0.055 > 0）", 11.5, GREY, space_after=10)
line(tf, "26.5%", 24, AMBER, True, space_after=2)
line(tf, "四次型 Q̂（内含一处 14.98 倍的最坏值放缩 Q₃）", 11.5, GREY, space_after=10)
line(tf, "其余六项合计 < 12%，微调末项无收益。", 12, INK, True, space_after=0)

# ============ 3 目录 ============
s, y = base("汇报路径", "CONTENTS")
items = [("01", "问题设定与已知基准", "Λ₃ 的定义、STW 命题 6.1、为什么分别放缩没有余量"),
         ("02", "主定理与证明架构", "常数 c 的闭式、两条证明主线与唯一接力点"),
         ("03", "结构线：从变分到 24 个原子", "二阶变分、符号类秩控制、分离观测、极小依赖"),
         ("04", "证书线：联合矩恒等式", "八项非负分解、接力链、分支阈值"),
         ("05", "验证、余量与天花板", "机器验证六组实验、松驰七组归属、(29+10√5)/40"),
         ("06", "审核、阻碍与路线", "G1/G2 修补、五条阻碍、优先级矩阵与 P0–P3")]
for i, (no, t, sub) in enumerate(items):
    col, row = i % 2, i // 2
    x0 = 0.58 + col * 6.22
    y0 = y + 0.10 + row * 1.42
    rect(s, x0, y0, 5.98, 1.24, fill=GREYL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.06)
    rect(s, x0, y0, 0.07, 1.24, fill=BLUE2)
    _, tf = tb(s, x0 + 0.26, y0 + 0.14, 0.8, 0.5)
    line(tf, no, 24, BLUE2, True, first=True, space_after=0)
    _, tf = tb(s, x0 + 1.05, y0 + 0.16, 4.75, 0.4)
    line(tf, t, 15.5, INK, True, first=True, space_after=0)
    _, tf = tb(s, x0 + 1.05, y0 + 0.58, 4.75, 0.55)
    line(tf, sub, 11, GREY, first=True, space_after=0)

# ============ 4 背景与桥 ============
s, y = base("背景与桥：从 Nikiforov 图极限到各向同性律", "00 · 背景")
_, tf = tb(s, 0.58, y, 6.0, 4.6)
line(tf, "已具备的底座：graphon", 16, BLUE, True, first=True, space_after=6)
line(tf, "• 对称可测 W:[0,1]²→[0,1]；有限图 G ↔ 阶梯 graphon W_G；blow-up 不改变 W。", 13, space_after=5, indent=0.22)
line(tf, "• 割度量给出紧性：极值序列有极限 graphon，极大值可取到。", 13, space_after=5, indent=0.22)
line(tf, "• Nikiforov 解析方法：谱极值量写成 graphon 泛函，取极限后变分刻画。", 13, space_after=5, indent=0.22)
line(tf, "• 但 λ₂,λ₃ 不是同态密度，不能直接套密度型紧性——本文第一个难点。", 13, AMBER, space_after=0, indent=0.22)
rect(s, 6.95, y, 5.83, 4.6, fill=BLUEL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.04)
_, tf = tb(s, 7.2, y + 0.24, 5.3, 4.1)
line(tf, "本文的第一块砖：Ky Fan 桥", 15.5, BLUE, True, first=True, space_after=7)
line(tf, "• Ky Fan：λ₁+λ₂+λ₃ = max tr(P_U A)（dim U=3），迹型泛函，极限下可转移。", 12.5, space_after=6, indent=0.22)
line(tf, "• 极值 graphon 的三维最大化子空间给特征函数 f₁,f₂,f₃；点 x ↦ v(x)∈R³。", 12.5, space_after=6, indent=0.22)
line(tf, "• 结构被律 μ=law(v(x)) 捕获；归一化后 μ 各向同性，E[vvᵀ]∝I。", 12.5, space_after=6, indent=0.22)
line(tf, "• 谱问题 → 概率/矩问题：对所有各向同性律控制一个 μ-决定的量。", 12.5, INK, True, space_after=0, indent=0.22)
_, tf = tb(s, 0.58, y + 4.78, 12.2, 0.6)
line(tf, "这座桥即毕业论文 Lemma 2.3：把 Nikiforov 的图极限方法从 λ₁ 推广到 Λ₃ 的关键一步。", 12.5, GREY, first=True)

# ============ 4 问题设定 ============
s, y = base("问题设定：只给定 n，Λ₃ 最大能到多少", "01 · 问题与基准")
_, tf = tb(s, 0.58, y, 7.7, 4.5)
line(tf, "• Λ₃(G) = λ₁(G)+λ₂(G)+λ₃(G)，邻接特征值按代数序排列、保留符号；极值问题只以阶数 n 为约束。", 14, first=True, space_after=8, indent=0.22)
line(tf, "• 出处：Ebrahimi 2008 硕士论文最早系统提出；前 k 大特征值之和的一般界见 Mohar 2009。", 14, space_after=8, indent=0.22)
line(tf, "• 语言：graphon。Σ_k(W) = sup_{EFFᵀ=I_k} ∫∫ W(x,y)F(x)·F(y)，M = sup_W Σ₃(W) 可达；有限图的归一邻接矩阵是 step graphon。", 14, space_after=8, indent=0.22)
line(tf, "• 约定：Σ_k 对 graphon 取前 k 个正特征值之和（不足补零）；Λ₃(G) 不做此约定，符号保留。", 14, space_after=8, indent=0.22)
line(tf, "• 为什么 k = 3 是最有信息量的情形：秩三绝对投影常数 γ(3) = φ 是精确值，分别放缩已无余量，任何改进都必须来自耦合（下一页）。", 14, space_after=0, indent=0.22)
rect(s, 8.62, y, 4.16, 4.5, fill=BLUEL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.04)
_, tf = tb(s, 8.86, y + 0.22, 3.7, 4.0)
line(tf, "归一化语言", 15, BLUE, True, first=True, space_after=8)
line(tf, "X ∈ ℝ³，E[XXᵀ] = I₃，X′ 为独立副本：", 12.5, space_after=10)
pic(s, "eq_AB.png", 8.86, y + 1.15, w=3.62)
_, tf = tb(s, 8.86, y + 2.35, 3.7, 1.9)
line(tf, "STW 的界 = 两次独立放缩之和；", 12.5, space_after=6)
line(tf, "本文的界 = 对 (A, B) 同时取大值的定量限制。", 12.5, BLUE, True, space_after=0)

# ============ 5 已知基准 STW ============
s, y = base("已知基准：STW 命题 6.1 与它的余量结构", "01 · 问题与基准")
pic(s, "eq_stw.png", 0.72, y + 0.05, w=3.5)
_, tf = tb(s, 4.6, y, 8.2, 1.0)
line(tf, "Sivashankar–Tang–Wakhare（2026）命题 6.1，证明用秩三绝对投影常数 γ(3) = (1+√5)/2 = φ；", 13.5, first=True, space_after=4)
line(tf, "所需投影常数界由 Deręgowska–Lewandowska 给出。", 13.5, GREY, space_after=0)
rect(s, 0.58, y + 1.25, 12.2, 0.02, fill=LINE)
_, tf = tb(s, 0.58, y + 1.5, 6.0, 3.4)
line(tf, "余量结构：两个量天然耦合", 15.5, BLUE, True, first=True, space_after=8)
line(tf, "• 分别估计 A ≤ φ 与 B ≤ 1，相加即得 STW 界。", 13.5, space_after=7, indent=0.22)
line(tf, "• γ(3) = φ 是精确的：在 A 取等的方向上没有任何余量可捞。", 13.5, space_after=7, indent=0.22)
line(tf, "• 因此改进的唯一来源是对 A 与 B 同时取大值的定量限制。", 13.5, space_after=7, indent=0.22)
line(tf, "• 这个“混合投影常数”问题 STW 在其 §6.2 已提出但未解决；本文直接采用它的各向同性律表述。", 13.5, space_after=0, indent=0.22)
rect(s, 6.95, y + 1.5, 5.83, 3.4, fill=AMBERL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.04)
_, tf = tb(s, 7.2, y + 1.72, 5.35, 3.0)
line(tf, "这意味着什么", 15, AMBER, True, first=True, space_after=8)
line(tf, "改进必须是结构性的：不是把某个常数调紧一点，而是换一种同时控制 A 与 B 的恒等式。", 13, space_after=8)
line(tf, "本文的答案：一条联合的三阶 + 四阶矩恒等式，让“均值大”与“绝对能量大”两种情形的误差互相支付。", 13, space_after=0)

# ============ 6 主定理 ============
s, y = base("主定理：显式常数 c 与两条等价表述", "02 · 主定理")
pic(s, "eq_c.png", 0.66, y + 0.02, w=8.6)
_, tf = tb(s, 0.58, y + 1.05, 12.2, 0.4)
line(tf, "定理 1.1（main）[已证]　对一切 n ≥ 3 阶有限简单图：", 14.5, INK, True, first=True, space_after=0)
pic(s, "eq_main.png", 1.1, y + 1.5, w=2.6)
pic(s, "eq_law.png", 5.2, y + 1.5, w=4.3)
_, tf = tb(s, 9.7, y + 1.55, 3.1, 0.9)
line(tf, "右式对一切满足 E[XXᵀ]=I₃ 的律成立；只假设二阶矩有限，允许零向量与无界半径。", 11.5, GREY, first=True, space_after=0)
rect(s, 0.58, y + 2.75, 12.2, 0.02, fill=LINE)
_, tf = tb(s, 0.58, y + 2.98, 6.0, 2.3)
line(tf, "• 系数较 STW 降低 0.024343738651028919…，常数项 −3 不变。", 13.5, first=True, space_after=7, indent=0.22)
line(tf, "• 图论版是律版的转移：归一邻接矩阵是 step graphon，转移产生精确的 −3（来源 ΣᵢQᵢ = 3）。", 13.5, space_after=7, indent=0.22)
line(tf, "• 不声称 c 最优；不与涉及额外图参数的估计比较。", 13.5, space_after=0, indent=0.22)
rect(s, 6.95, y + 2.98, 5.83, 2.3, fill=GREENL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.04)
_, tf = tb(s, 7.2, y + 3.18, 5.35, 1.95)
line(tf, "常数是“算出来的”，不是“拟合出来的”", 14, GREEN, True, first=True, space_after=8)
line(tf, "全部系数在精确有理算术与二次域 ℚ(√3)、ℚ(√5) 中验证，从不用浮点；16 张符号证书残差 ≤ 1.7×10⁻¹⁶。", 12.5, space_after=0)

# ============ 7 证明架构总览 ============
s, y = base("证明架构：两条主线，一个接力点", "02 · 主定理")
rect(s, 0.58, y, 5.85, 2.5, fill=BLUEL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.05)
_, tf = tb(s, 0.85, y + 0.18, 5.35, 2.2)
line(tf, "结构线（负责收敛）", 16, BLUE, True, first=True, space_after=6)
line(tf, "二阶变分半负定 → 符号类秩控制 → 矩保持变分 → 分离二元观测 → 极小依赖", 12.5, space_after=6)
line(tf, "产出：m ≤ 24 个原子；正点团 |J| ≤ 6；原核 = 二值阶梯核 a.e.", 12.5, BLUE, True, space_after=0)
rect(s, 6.92, y, 5.86, 2.5, fill=AMBERL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.05)
_, tf = tb(s, 7.19, y + 0.18, 5.35, 2.2)
line(tf, "证书线（负责闭合）", 16, AMBER, True, first=True, space_after=6)
line(tf, "保留 e、V、回归误差、J → 带符号四次控制多项式 → 方向 Hermite + 径向比较 → 八项非负分解", 12.5, space_after=6)
line(tf, "产出：D ≥ d₀（大分支）、D ≥ d₁（小分支），合成 M ≤ c", 12.5, AMBER, True, space_after=0)
ar = rect(s, 6.05, y + 2.62, 1.25, 0.5, fill=GREY, shape=MSO_SHAPE.DOWN_ARROW)
_, tf = tb(s, 0.58, y + 3.25, 12.2, 0.5)
line(tf, "唯一接力点：结构线的 m ≤ 24 是证书线第八项 χ − C²/24 ≥ 0 的非负依据；除此之外两线独立。", 14, RED, True, first=True, space_after=0, align=PP_ALIGN.CENTER)
rect(s, 3.4, y + 3.85, 6.55, 0.72, fill=BLUE, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.12)
_, tf = tb(s, 3.4, y + 3.98, 6.55, 0.5, anchor=MSO_ANCHOR.MIDDLE)
line(tf, "主定理 1.1：Λ₃(G) ≤ cn − 3", 17, WHITE, True, first=True, space_after=0, align=PP_ALIGN.CENTER)
_, tf = tb(s, 0.58, y + 4.85, 12.2, 0.6)
line(tf, "先例：变分/克隆/卦限归约见 KMMP 与 KLMPT；线性形式幂次的极小依赖一般界属 Li–Liu。本文不声称新的依赖一般界，只在自己的具体假设下证分离观测引理。", 11.5, GREY, first=True, space_after=0)

# ============ 8 结构线 I ============
s, y = base("结构线 I：从二阶变分到“有限”", "03 · 结构线")
steps = [("定理 3.1", "在实际概率空间上做二阶变分，极大元处的核半负定"),
         ("引理 3.2 / 3.3", "四个符号类上控制正部 Gram 核的秩；|T| ≤ 2^{d−1} 维数下降"),
         ("定理 3.4", "典范支撑有限：|S| ≤ 1204（“有限”的定量版本）")]
for i, (t, d) in enumerate(steps):
    x0 = 0.58 + i * 4.14
    rect(s, x0, y, 3.94, 1.28, fill=GREYL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.07)
    rect(s, x0, y, 3.94, 0.36, fill=BLUE2)
    _, tf = tb(s, x0 + 0.16, y + 0.05, 3.6, 0.3)
    line(tf, t, 13, WHITE, True, first=True, space_after=0)
    _, tf = tb(s, x0 + 0.16, y + 0.46, 3.62, 0.75)
    line(tf, d, 11.5, INK, first=True, space_after=0)
    if i < 2:
        rect(s, x0 + 3.94, y + 0.5, 0.2, 0.28, fill=GREY, shape=MSO_SHAPE.RIGHT_ARROW)
_, tf = tb(s, 0.58, y + 1.62, 12.2, 0.4)
line(tf, "支撑的“形状”：三个行级引理", 15.5, BLUE, True, first=True, space_after=0)
rows = [["引理", "内容", "审核状态"],
        ["2.4 ellipse-structure", "椭球恒等式：FᵀLF = M a.e.", "[已核] 独立重推通过"],
        ["2.5 row-insertion", "行插入不等式：rᵀL⁻¹r ≤ M", "[已核] 独立重推通过"],
        ["2.6 contact-structure", "阈值公式 Σⱼ pⱼBⱼv = Lv", "[缺口→已证] G1 修补（约十行）"]]
table(s, 0.58, y + 2.05, 12.2, rows, [3.1, 5.6, 3.5], size=12, rowh=0.42)
_, tf = tb(s, 0.58, y + 3.9, 12.2, 1.0)
line(tf, "注意：1204 只是“有限”的定量版本；真正的 24 必须靠下一页的分离观测。二阶变分在“实际”概率空间上做，而不是抽象参数空间，这是后续所有刚性结论的源头。", 12.5, GREY, first=True, space_after=0)

# ============ 9 结构线 II ============
s, y = base("结构线 II：分离观测与 m ≤ 24", "03 · 结构线")
_, tf = tb(s, 0.58, y, 6.1, 4.4)
line(tf, "定理 1.2（structure-main）[已证]", 15.5, BLUE, True, first=True, space_after=7)
line(tf, "(a) 典范律支撑至多 24 点；", 13.5, space_after=5, indent=0.22)
line(tf, "(b) W 在正测度纤维上几乎处处等于二值阶梯核 1{vᵢ·v>0}；", 13.5, space_after=5, indent=0.22)
line(tf, "(c) Σpᵢvᵢvᵢᵀ = I₃，vᵢᵀLv = M，vᵢ₁ > 0，vᵢ·v ≠ 0，ΣⱼpBᵢvⱼ = Lvᵢ；", 13.5, space_after=5, indent=0.22)
line(tf, "(d) 每个正点团 J 上 {vᵢvᵀ} 线性无关，故 |J| ≤ 6。", 13.5, space_after=8, indent=0.22)
line(tf, "断言对每一个原始极大元成立，不是“存在某个小支撑极大元”。", 12, GREY, True, space_after=0)
rect(s, 6.95, y, 5.83, 4.4, fill=GREYL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.04)
_, tf = tb(s, 7.2, y + 0.2, 5.35, 4.0)
line(tf, "从 1204 到 24 的三步", 15, BLUE, True, first=True, space_after=8)
line(tf, "① 矩保持变分 ⇒ 典范向量对极大核的每一个二值行都满足向量关系；", 12.5, space_after=7)
line(tf, "② 引理 4.2：这些二值行分离不同的支撑点（G2 所在，已补一句）；", 12.5, space_after=7)
line(tf, "③ 极小依赖论证：正点团上不可能残存二阶矩关系 ⇒ {vᵢvᵢᵀ} 线性无关 ⇒ |J| ≤ 6；四个符号类各 ≤ 6 ⇒ m ≤ 24。", 12.5, space_after=8)
line(tf, "承重性：m ≤ 24 ⇒ 证书第八项 χ − C²/24 ≥ 0。", 12.5, RED, True, space_after=0)

# ============ 10 证书线 I ============
s, y = base("证书线 I：保留四类误差，一个联合恒等式", "04 · 证书线")
_, tf = tb(s, 0.58, y, 5.5, 2.0)
line(tf, "保留量（单个各向同性律）", 15, BLUE, True, first=True, space_after=6)
line(tf, "e 均值亏损　·　V 径向方差", 13, space_after=4)
line(tf, "回归误差　·　J 三/四阶加权调和矩", 13, space_after=6)
line(tf, "带符号四次控制多项式 f：在 w = −3/8 与 1/2 双接触、在 w = 1 接触；五个条件定五个系数，故 41503 = 7³·11² 由插值唯一确定。", 12, GREY, space_after=0)
pic(s, "eq_D.png", 6.3, y + 0.1, w=6.4)
rows = [["编号", "非负项", "含义"],
        ["T1", "r_tr = (8/9)P − E[RR′z_rad g]", "径向迹项"],
        ["T2", "Q̂（含 Q₃ 与四次残量）", "带符号四次型净贡献"],
        ["T3", "ℓ·r_reg，ℓ = 937/1400", "回归误差项"],
        ["T4", "R₀", "四次型松驰中的储备"],
        ["T5", "r_off / 9", "离对角项"],
        ["T6", "K²/1125", "二次储备"],
        ["T7", "(κ/9)(χ − C²/24)，κ = 3/125000", "正质量自对修正，用 m ≤ 24"]]
table(s, 0.58, y + 2.15, 12.2, rows, [1.1, 6.3, 4.8], size=11.5, rowh=0.36)
_, tf = tb(s, 0.58, y + 5.0, 12.2, 0.5)
line(tf, "八项非负量归并展示为七组（T2 内含两项）。恒等式是精确的：没有使用任何未量化的稳定性缺口。", 12, GREY, True, first=True, space_after=0)

# ============ 11 证书线 II ============
s, y = base("证书线 II：接力链与分支阈值", "04 · 证书线")
chain = [("引理 6.1", "σ 的两条点态不等式"), ("引理 6.3", "径向支付"), ("命题 7.2", "大分支 D ≥ d₀"),
         ("引理 8.2", "z*² ≤ (3/2)σ"), ("命题 8.5", "小分支 D ≥ d₁"), ("引理 9.1", "谱储备"), ("命题 9.3", "M ≤ c")]
bw = 1.62
for i, (t, d) in enumerate(chain):
    x0 = 0.58 + i * (bw + 0.13)
    fill = GREEN if i == len(chain) - 1 else BLUEL
    rect(s, x0, y, bw, 1.05, fill=fill, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.08)
    _, tf = tb(s, x0 + 0.08, y + 0.12, bw - 0.16, 0.32)
    line(tf, t, 12.5, BLUE if i < len(chain) - 1 else WHITE, True, first=True, space_after=0, align=PP_ALIGN.CENTER)
    _, tf = tb(s, x0 + 0.08, y + 0.5, bw - 0.16, 0.5)
    line(tf, d, 10, INK if i < len(chain) - 1 else WHITE, first=True, space_after=0, align=PP_ALIGN.CENTER)
    if i < len(chain) - 1:
        rect(s, x0 + bw, y + 0.42, 0.13, 0.2, fill=GREY, shape=MSO_SHAPE.RIGHT_ARROW)
_, tf = tb(s, 0.58, y + 1.45, 12.2, 0.4)
line(tf, "分支分界（式 smallbranch）：e < 1/20、V < 21/104、C > 5/3，由 α ≥ 26/105、β ≥ 6/7、ρ ≥ 1/2 三个粗下界推出。", 13, INK, True, first=True, space_after=0)
rect(s, 0.58, y + 2.05, 5.95, 2.5, fill=GREYL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.05)
_, tf = tb(s, 0.83, y + 2.25, 5.45, 2.1)
line(tf, "两处已知的非紧常数（如实标注）", 14, BLUE, True, first=True, space_after=7)
line(tf, "• 引理 8.2 的 3/2 不紧：真实上确界 z*²/σ 在两个分支上是 1.200224 与 0.615384（注 rem-theta）。", 12.5, space_after=6, indent=0.22)
line(tf, "• 半径阈值 4/3 → 32/27 可再得 12.5% 收益，但当前证明不需要，故不改（注 rem-radii）。", 12.5, space_after=0, indent=0.22)
rect(s, 6.83, y + 2.05, 5.95, 2.5, fill=AMBERL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.05)
_, tf = tb(s, 7.08, y + 2.25, 5.45, 2.1)
line(tf, "这条链的设计意图", 14, AMBER, True, first=True, space_after=7)
line(tf, "径向与方向误差不被分别丢掉，而是留在同一个恒等式里相互抵消；", 12.5, space_after=6)
line(tf, "因此可以在不改进任何单个投影常数的前提下改进它们的和——这就是超过 STW 的全部来源。", 12.5, space_after=0)

# ============ 12 为什么超过 STW ============
s, y = base("为什么这条路线能超过 STW", "04 · 证书线")
reasons = [("1", "STW 的界来自两次独立放缩", "A ≤ φ 与 B ≤ 1 相加；而 γ(3) = φ 精确，分别放缩这条路已经没有余量。"),
           ("2", "本文把 A 与 B 放进同一个恒等式", "让“均值大”与“绝对能量大”两种情形的误差互相支付，于是不改进任何单个投影常数也能改进它们的和。"),
           ("3", "证书是定量的，不是定性的", "一条精确恒等式 + 八项显式非负量，没有未量化的稳定性缺口；因此每一项的松紧都可以被计算、被指认、被改进。")]
for i, (no, t, d) in enumerate(reasons):
    y0 = y + 0.1 + i * 1.42
    rect(s, 0.58, y0, 12.2, 1.26, fill=GREYL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.06)
    rect(s, 0.58, y0, 0.9, 1.26, fill=BLUE)
    _, tf = tb(s, 0.58, y0 + 0.36, 0.9, 0.6)
    line(tf, no, 30, WHITE, True, first=True, space_after=0, align=PP_ALIGN.CENTER)
    _, tf = tb(s, 1.75, y0 + 0.16, 10.8, 0.4)
    line(tf, t, 16, INK, True, first=True, space_after=0)
    _, tf = tb(s, 1.75, y0 + 0.6, 10.8, 0.6)
    line(tf, d, 12.5, GREY, first=True, space_after=0)
_, tf = tb(s, 0.58, y + 4.5, 12.2, 0.5)
line(tf, "推论：第 3 条使“缺口诊断”成为可能——下一页开始量化缺口到底在哪里、有多大、能不能回收。", 13, BLUE, True, first=True, space_after=0)

# ============ 13 机器验证 ============
s, y = base("机器验证与数值证据：只交叉验证，不作证明", "05 · 验证与余量")
rows = [["实验", "口径", "结果", "完备度"],
        ["有限图转移", "3 ≤ n ≤ 6 全部标号图 2³+2⁶+2¹⁰+2¹⁵ = 33864", "零反例；最小余量 +0.0224818（K₁∨C₅，n=6）", "100"],
        ["符号证书", "16 张证书的恒等式与配方", "残差全部 ≤ 1.7×10⁻¹⁶", "100"],
        ["闭式常数", "c、ρ₀、天花板三值", "独立复算一致", "100"],
        ["k 网格全局最大", "k = 6,…,20 离散网格上 (A+B)/2", "全局最大值恰为 ρ₀，从不超过", "90"],
        ["吹起族渐近", "K₁∨C₅ 吹起，n = 126…2016", "Λ₃/n 从 1.25770 单调升到 1.28008 → ρ₀", "90"],
        ["随机律抽样", "4000 组各向同性分布", "零反例", "90"]]
table(s, 0.58, y, 12.2, rows, [2.2, 4.3, 4.5, 1.2], size=11.5, rowh=0.5)
_, tf = tb(s, 0.58, y + 3.35, 12.2, 1.2)
line(tf, "完备度 100 = 穷举或有完备性论证；90 = 抽样或离散网格。所有“>0”与恒等式断言另在精确有理算术与 ℚ(√3)、ℚ(√5) 中复核，点态不等式额外在 2×10⁶+1 个网格点上扫描。", 12.5, GREY, first=True, space_after=6)
line(tf, "口径提示：33864 是标号图数目；未标号图只有 205 个。论文附录建议写明 labelled（已列入排版建议 T3）。", 12.5, AMBER, True, space_after=0)

# ============ 14 松驰归属 ============
s, y = base("冠军律处的松驰归属：缺口集中在两处", "05 · 验证与余量")
pic(s, "chart_slack.png", 0.62, y - 0.05, w=8.5)
rect(s, 9.35, y, 3.43, 4.55, fill=GREYL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.05)
_, tf = tb(s, 9.6, y + 0.2, 2.95, 4.2)
line(tf, "两条直接结论", 14.5, BLUE, True, first=True, space_after=8)
line(tf, "① r_tr 与 Q̂ 合计 88%：要继续收紧必须动主项，微调末四项毫无意义。", 12.5, space_after=9)
line(tf, "② Q̂ 中 15.27% 来自一处最坏值放缩（Q₃），系数松 14.98 倍——全文最容易识别的靶子。", 12.5, space_after=9)
line(tf, "松驰总量 6.210×10⁻³ = 2(c−ρ₀)。", 11.5, GREY, True, space_after=0)

# ============ 15 Q3 的 15 倍 ============
s, y = base("最便宜的靶子：Q₃ 的一处最坏值放缩", "05 · 验证与余量")
cards = [("0.0027055", "被代入的下界", "系数 (183/200 − BIG·q²) 在 q² ≤ 20/57 上的最小值", AMBER),
         ("0.040527", "冠军律上的真值", "同一系数在冠军律处的取值", GREEN),
         ("14.98×", "放缩倍数", "两者之比，即此处白白丢掉的精度", RED),
         ("15.27%", "占全部松驰", "在 6.210×10⁻³ 中的份额", BLUE)]
cw = 2.95
for i, (v, t, sub, col) in enumerate(cards):
    x0 = 0.58 + i * (cw + 0.13)
    rect(s, x0, y, cw, 1.75, fill=GREYL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.06)
    rect(s, x0, y, cw, 0.09, fill=col)
    _, tf = tb(s, x0 + 0.2, y + 0.24, cw - 0.4, 0.5)
    line(tf, v, 24, col, True, first=True, space_after=0)
    _, tf = tb(s, x0 + 0.2, y + 0.78, cw - 0.4, 0.32)
    line(tf, t, 13, INK, True, first=True, space_after=0)
    _, tf = tb(s, x0 + 0.2, y + 1.12, cw - 0.4, 0.55)
    line(tf, sub, 10.5, GREY, first=True, space_after=0)
_, tf = tb(s, 0.58, y + 2.1, 12.2, 2.6)
line(tf, "为什么不能“直接代冠军值”", 15.5, RED, True, first=True, space_after=8)
line(tf, "• Q 对 q² 的净依赖是符号不定的：保留 q² 为变量 ≠ 代入冠军值就变好。", 13.5, space_after=7, indent=0.22)
line(tf, "• 正确做法：把 Q 看成 (√e, √V, √J, q) 的、系数含 q² 的二次型，在小分支约束集 {e<1/20, V<21/104, q²≤20/57, J≤V/3} 上做一维（q²）分段有理配方；", 13.5, space_after=7, indent=0.22)
line(tf, "• 或者用 SOS 重新拟合。两种做法都可以用精确有理算术完成并逐条验证，不需要浮点 SDP 的结论进入证明。", 13.5, space_after=7, indent=0.22)
line(tf, "• 这是 P1 的第一条子路线：属于“重新配方”，不依赖新想法，信心较高。", 13.5, BLUE, True, space_after=0, indent=0.22)

# ============ 16 常数阶梯 ============
s, y = base("常数阶梯：总差 0.027449 的三段归属", "05 · 验证与余量")
pic(s, "chart_ladder.png", 0.62, y - 0.1, w=8.6)
rect(s, 9.42, y, 3.36, 4.5, fill=BLUEL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.05)
_, tf = tb(s, 9.66, y + 0.2, 2.9, 4.1)
line(tf, "三段占比", 14.5, BLUE, True, first=True, space_after=8)
line(tf, "88.7%", 22, GREEN, True, space_after=2)
line(tf, "已证改进 0.024344", 11.5, GREY, space_after=9)
line(tf, "2.4%", 22, BLUE2, True, space_after=2)
line(tf, "P1 可回收 0.000656", 11.5, GREY, space_after=9)
line(tf, "8.9%", 22, AMBER, True, space_after=2)
line(tf, "框架外余量 0.002449（P2）", 11.5, GREY, space_after=9)
line(tf, "形状本身就是结论：最难的部分已经过去，但最后一公里需要换工具。", 11.5, INK, True, space_after=0)

# ============ 17 框架天花板 ============
s, y = base("本框架的天花板：(29+10√5)/40", "05 · 验证与余量")
pic(s, "eq_ceiling.png", 0.72, y + 0.02, w=5.6)
_, tf = tb(s, 6.7, y, 6.1, 1.2)
line(tf, "达到它的充要条件：小分支为空，即任何达到 M 的律都不落在 e<1/20、V<21/104、C>5/3 区域。", 13.5, first=True, space_after=0)
rect(s, 0.58, y + 1.35, 12.2, 0.02, fill=LINE)
_, tf = tb(s, 0.58, y + 1.6, 6.0, 3.2)
line(tf, "成功后的三个好处", 15.5, BLUE, True, first=True, space_after=8)
line(tf, "• 改进量从 0.0243437 变成干净的 1/40 = 0.025；", 13.5, space_after=7, indent=0.22)
line(tf, "• c 从“明显凑出来的”表达式变成有闭式意义的数，可读性与可信度大幅提升；", 13.5, space_after=7, indent=0.22)
line(tf, "• 两个分支给出同一个常数，是自洽性的强检验。", 13.5, space_after=0, indent=0.22)
rect(s, 6.95, y + 1.6, 5.83, 3.2, fill=AMBERL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.04)
_, tf = tb(s, 7.2, y + 1.8, 5.35, 2.8)
line(tf, "必须如实说的比例", 15, AMBER, True, first=True, space_after=8)
line(tf, "可回收 c − 天花板 = 0.000656261，", 13, space_after=5)
line(tf, "只占当前缺口 c − ρ₀ = 0.003105059 的 21.1%。", 13, space_after=8)
line(tf, "即：P1 回收五分之一，剩下 78.9% 在框架之外，属于 P2。", 13, RED, True, space_after=0)
pic(s, "eq_rho0.png", 0.72, y + 4.95, w=4.6)
_, tf = tb(s, 5.6, y + 5.0, 7.2, 0.6)
line(tf, "下界端点来自加权轮图（K₁∨C₅ 六块阶梯核）；其全局最优性尚未证明。", 12, GREY, first=True, space_after=0)

# ============ 18 审核 G1/G2 ============
s, y = base("独立审核：两处承重跳步与修补", "06 · 审核与路线")
rect(s, 0.58, y, 6.0, 3.1, fill=REDL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.05)
_, tf = tb(s, 0.83, y + 0.18, 5.5, 2.8)
line(tf, "G1 · 引理 2.6（contact-structure）", 14.5, RED, True, first=True, space_after=7)
line(tf, "缺：极大行在零集意义下唯一的论证。原文只写 “Uniqueness forces ∫_A F = 0”。", 12.5, space_after=7)
line(tf, "补（约十行）：把 A₀ = {v·F=0} 上的候选行 h_v + 1_A 与已知行配对，二者同为极大元；由椭球恒等式得 ∫_A F = 0 对一切 A ⊂ A₀ 成立 ⇒ F = 0 a.e. ⇒ A₀ 零测。", 12.5, space_after=7)
line(tf, "承重：阈值公式 ⇒ A + B = 2M 的入口。", 12, RED, True, space_after=0)
rect(s, 6.78, y, 6.0, 3.1, fill=REDL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.05)
_, tf = tb(s, 7.03, y + 0.18, 5.5, 2.8)
line(tf, "G2 · 引理 4.2（binary，分离观测）", 14.5, RED, True, first=True, space_after=7)
line(tf, "缺：dim H = 2、0 ∈ H 一支中“所有观测在 I 上常值”的情形未处理。", 12.5, space_after=7)
line(tf, "补（一句）：若每个分离观测在 I 上常值，则 I 上原子在该观测族下不可区分，与 I 取自正团且 |I| ≥ 2 矛盾，故该情形不出现。", 12.5, space_after=7)
line(tf, "承重：m ≤ 24 ⇒ 证书第七项 χ − C²/24 ≥ 0。", 12, RED, True, space_after=0)
_, tf = tb(s, 0.58, y + 3.35, 12.2, 1.4)
line(tf, "审核方式与结论", 15, BLUE, True, first=True, space_after=6)
line(tf, "33 页报告逐条复核定理、引理与 16 张符号证书，并独立重推引理 2.5 等关键步骤；结论：主定理成立，两处跳步均可在五行以内补齐且无连带影响；另给七条表述建议 S1–S7 与残余风险 R1（Ebrahimi 2008 全文未获取）。", 12.5, space_after=0)

# ============ 19 S1-S7 ============
s, y = base("S1–S7 与附录 A：已落实的表述加固", "06 · 审核与路线")
rows = [["编号", "位置", "改动"],
        ["S1", "引理 2.3", "补一句：块均值零子空间使 Σ₃ 等于商矩阵的谱和"],
        ["S2", "引理 2.5", "显式写出 AB 试探标架的单位化归一"],
        ["S3", "§5.1", "补半句 “R_t = Y 故 E[R_t] = x”，不留推导空档"],
        ["S4", "引理 8.2", "新增注 rem-theta：3/2 不紧，真实上确界 1.200224 / 0.615384"],
        ["S5", "新增附录 A", "常数分“唯一确定 / 选定后验证”两类 + 验证协议（精确算术、网格扫描、33864 图）"],
        ["S6", "定理 1.2", "重排为 (a)–(d) 四条，支撑数与二值核不再混在一句"],
        ["S7", "§9", "新增注 rem-radii：√8<3 相容，4/3→32/27 的 12.5% 收益不需要"]]
table(s, 0.58, y, 12.2, rows, [1.0, 2.2, 9.0], size=11.5, rowh=0.44)
_, tf = tb(s, 0.58, y + 3.7, 12.2, 1.0)
line(tf, "效果：论文 19 页 → 20 页，pdflatex 两遍编译无警告、无未解析引用；附录 A 直接回应“系数是调出来的”这一类质疑。", 12.5, GREEN, True, first=True, space_after=0)

# ============ 20 阻碍清单 ============
s, y = base("当前阻碍：五条必须明说的问题", "06 · 审核与路线")
rows = [["#", "阻碍", "性质", "对路线的影响"],
        ["1", "r_tr 在冠军律上 = 3.8×10⁻³ ≠ 0（9J−V = 0.055 > 0）", "根本性", "紧证书不该含 r_tr：径向支付路线很可能整体放弃（P2 起点）"],
        ["2", "Q 对 q² 净依赖符号不定", "技术性", "P1 第 1 条需参数化配方或 SOS 重拟合，不能直接代值"],
        ["3", "收紧小分支阈值仅得 0.16%", "收益极低", "零风险但只腾空间，不能单独解决问题"],
        ["4", "八项无一为零 ⇒ 需要新等号情形", "本质性", "P2 不宜承诺；近期唯一动作是“24 降到 6”"],
        ["5", "Ebrahimi 2008 全文未获取（R1）", "文献风险", "不影响正确性，影响问题归属表述；建议馆际互借"]]
table(s, 0.58, y, 12.2, rows, [0.5, 4.6, 1.5, 5.6], size=11.5, rowh=0.62)
_, tf = tb(s, 0.58, y + 3.6, 12.2, 1.2)
line(tf, "方法上的自我约束：极图无法由算法推出，逐阶跑算法对本问题没有意义。所有数值与优化只用于发现常数与交叉验证；最终证书必须用精确有理算术与二次域算术逐条验证。", 13, BLUE, True, first=True, space_after=0)

# ============ 21 优先级矩阵 ============
s, y = base("优先级矩阵：难度 × 收益", "06 · 审核与路线")
qx, qy, qw, qh = 1.5, y + 0.35, 5.4, 2.0
rect(s, qx, qy, qw, qh, fill=GREENL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.04)
rect(s, qx + qw + 0.15, qy, qw, qh, fill=REDL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.04)
rect(s, qx, qy + qh + 0.15, qw, qh, fill=BLUEL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.04)
rect(s, qx + qw + 0.15, qy + qh + 0.15, qw, qh, fill=GREYL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.04)
def qtext(x0, y0, title, sub, col):
    _, tf = tb(s, x0 + 0.25, y0 + 0.18, qw - 0.5, 0.4)
    line(tf, title, 16, col, True, first=True, space_after=0)
    _, tf = tb(s, x0 + 0.25, y0 + 0.62, qw - 0.5, qh - 0.8)
    line(tf, sub, 12, INK, first=True, space_after=0)
qtext(qx, qy, "P1 小分支为空（先做）", "低难度 · 高收益：保留 q² 参数化配方；证“高度集中与小径向方差不能同时”", GREEN)
qtext(qx + qw + 0.15, qy, "P2 精确值 ρ₀", "高难度 · 高收益：八项全零证书；第一步把 24 降到 6", RED)
qtext(qx, qy + qh + 0.15, "P3 吹起渐近（并行）", "低难度 · 低收益：Weyl 给 O(1/n)，坐实 −3 位移与线性缺口", BLUE2)
qtext(qx + qw + 0.15, qy + qh + 0.15, "已明确排除", "逐阶算法搜索（极图推不出）；阈值微调（仅 0.16%）", GREY)
_, tf = tb(s, 0.62, qy + 1.6, 0.8, 1.0)
line(tf, "收\n益\n↑", 13, GREY, True, first=True, space_after=0, align=PP_ALIGN.CENTER)
_, tf = tb(s, qx, qy + 2 * qh + 0.28, 2 * qw + 0.15, 0.35)
line(tf, "证明难度 →", 13, GREY, True, first=True, space_after=0, align=PP_ALIGN.CENTER)

# ============ 22 P0-P3 路线图 ============
s, y = base("P0–P3 路线图与时间预算", "06 · 审核与路线")
names = ["P0 论文加固", "P1 小分支为空", "P2 精确值 ρ₀", "P3 吹起渐近"]
cols = [GREEN, BLUE2, RED, AMBER]
for i, (nm, col) in enumerate(zip(names, cols)):
    x0 = 0.58 + i * 3.13
    rect(s, x0, y, 2.95, 0.62, fill=col, shape=MSO_SHAPE.CHEVRON)
    _, tf = tb(s, x0 + 0.35, y + 0.14, 2.5, 0.4)
    line(tf, nm, 14, WHITE, True, first=True, space_after=0)
rows = [["优先级", "目标", "把握度（主观）", "时间预算", "状态"],
        ["P0", "G1/G2 修补 + S1–S7 + 附录 A", "100%", "1–2 天", "已完成"],
        ["P1", "c = (29+10√5)/40，回收 21.1% 缺口", "45%", "2 天评估 + 2–4 周", "主攻"],
        ["P2", "M = ρ₀：八项全零的新证书", "15%", "长期，不承诺", "仅做“24→6”"],
        ["P3", "Λ₃(Gₙ)/n → ρ₀ 严格化，误差 O(1/n)", "70%", "3–5 天", "与 P1 并行"]]
table(s, 0.58, y + 0.95, 12.2, rows, [1.2, 4.6, 2.0, 2.4, 2.0], size=12, rowh=0.5)
_, tf = tb(s, 0.58, y + 3.6, 12.2, 1.3)
line(tf, "排序依据：数学优先（能推理就不靠逐阶算法）· 收益/难度比 · 承重性（先补影响主定理的缺口）· 可验证性（每步都要有机器证书）。", 13, BLUE, True, first=True, space_after=6)
line(tf, "决策点：先花两天做 P1 第 3 条（零风险）并画出小分支下 (e, V, J, q²) 可行域，用图形判断第 2 条是否有希望，再决定投入。", 12.5, GREY, space_after=0)

# ============ 23 结论 + 排版建议 ============
s, y = base("结论三句话与论文排版建议要点", "06 · 审核与路线")
_, tf = tb(s, 0.58, y, 5.9, 2.2)
line(tf, "结论", 16, BLUE, True, first=True, space_after=8)
line(tf, "① 方法：联合矩证书取代逐项放缩；", 13.5, space_after=6)
line(tf, "② 结果：改进 0.024344，主定理已证并通过独立审核；", 13.5, space_after=6)
line(tf, "③ 边界：框架天花板 1.284017，c 不是终值，精确上界仍未解决。", 13.5, space_after=0)
rows = [["级", "排版建议要点（详见讲稿 §5）"],
        ["T1", "§1 末路线图编号化；命题 9.3 八项编号表；附录 A 常数一览表；§10 拆三个小标题"],
        ["T2", "加 Notation 段；c_w 改名避免与 c 混淆；§5.1 长段拆三；矩预算逐条加用途；8 处引理加“用途”注"],
        ["T3", "33864 写明标号口径；摘要补 |J|≤6；定理 1.1 假设前移；BIG 等记号给前向引用"],
        ["T4", "学位论文版改模板与行距；参考文献格式统一；引入 cleveref；补两张图（支撑点正团、吹起族收敛）"],
        ["T5", "期刊短文 → 学位论文章节映射（预备/结构线/证书线/转移/验证/展望六章）"]]
table(s, 6.7, y, 6.08, rows, [0.7, 5.38], size=10.5, rowh=0.62)
_, tf = tb(s, 0.58, y + 2.6, 5.9, 2.2)
line(tf, "为什么这些建议重要", 14.5, BLUE, True, first=True, space_after=7)
line(tf, "T1/T2 决定审稿人能否在十分钟内看懂依赖结构；T3 的四条都是“动手核对就会卡住”的口径问题；T4/T5 决定这份材料能否直接长成一章学位论文。", 12.5, GREY, space_after=0)

# ============ 25 适用范围与披露 ============
s, y = base("适用范围与新颖性披露", "06 · 审核与路线")
_, tf = tb(s, 0.58, y, 6.0, 4.6)
line(tf, "加权 [0,1] 适用范围（毕业论文 §6.4）", 16, BLUE, True, first=True, space_after=6)
line(tf, "• Prop 6.8：主证书对边权 a_ij∈[0,1]、a_ii=0 的加权图同样成立。", 13, space_after=5, indent=0.22)
line(tf, "• Prop 6.10：加权不改变渐近常数（graphon 密度论证），故无渐近增益。", 13, space_after=5, indent=0.22)
line(tf, "• 定位：适用范围结果，非新定理；与主线仅共享 Ky Fan 桥。", 13, space_after=0, indent=0.22)
rect(s, 6.95, y, 5.83, 4.6, fill=BLUEL, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.04)
_, tf = tb(s, 7.2, y + 0.24, 5.3, 4.1)
line(tf, "新颖性与优先权披露", 15.5, AMBER, True, first=True, space_after=7)
line(tf, "• 加权形式不新颖：STW Prop 6.1 证明本身只用 a_ii=0、0≤a_ij≤1；Mohar 2009 已覆盖对称矩阵。", 12.5, space_after=6, indent=0.22)
line(tf, "• 下界常数 c_w 与 K₁∨C₅ 模板：Zenodo 自存档（2026-08-04）条件性先提出，论文注 3.2 已让渡。", 12.5, space_after=6, indent=0.22)
line(tf, "• 主定理 Λ₃≤cn−3 与 24-原子结构：二次查新未被抢先。", 12.5, space_after=0, indent=0.22)
_, tf = tb(s, 0.58, y + 4.78, 12.2, 0.6)
line(tf, "以上均已在毕业论文 Declaration 与 §6.4 如实披露，汇报口径与论文一致。", 12.5, GREY, first=True)

# ============ 24 结语 ============
s, y = base("结语：下一步与交付物", "CLOSING")
rect(s, 0.58, y, 12.2, 1.05, fill=BLUE, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.08)
_, tf = tb(s, 0.9, y + 0.22, 11.6, 0.7)
line(tf, "主攻：证明小分支为空，把 c 推到 (29+10√5)/40。　副线：吹起族的严格渐近。", 17, WHITE, True, first=True, space_after=4)
line(tf, "方法上继续坚持：数学推理优先，算法只作交叉验证。", 13, RGBColor.from_string("C9D3E4"), space_after=0)
_, tf = tb(s, 0.58, y + 1.4, 6.0, 3.2)
line(tf, "交付物", 15.5, BLUE, True, first=True, space_after=7)
line(tf, "• 论文：top3_joint_moment_bound.tex / .pdf（20 页，P0 加固后）", 12.5, space_after=5, indent=0.22)
line(tf, "• 审核：审核报告_top3_joint_moment_bound.pdf（33 页）", 12.5, space_after=5, indent=0.22)
line(tf, "• 讲稿：讲稿与排版建议.pdf（16 页，含逐页口播与 T1–T5）", 12.5, space_after=5, indent=0.22)
line(tf, "• 本幻灯片：自编译 .pptx（python-pptx 版式引擎 + matplotlib 公式/图表）", 12.5, space_after=5, indent=0.22)
line(tf, "• git：github.com:jhy28740-jhy/graph-theory-study.git（main）", 12.5, space_after=0, indent=0.22)
_, tf = tb(s, 6.95, y + 1.4, 5.83, 3.2)
line(tf, "请导师重点指正三处", 15.5, BLUE, True, first=True, space_after=7)
line(tf, "① P1 的投入比例：是否值得为 21.1% 的回收投入 2–4 周；", 12.5, space_after=6, indent=0.22)
line(tf, "② P2 的定位：写成“目标与部分进展”是否合适；", 12.5, space_after=6, indent=0.22)
line(tf, "③ 学位论文的章节映射（讲稿 §5 T5）是否按此展开。", 12.5, space_after=0, indent=0.22)

OUT = os.path.join(HERE, "汇报PPT_自编译_前三大特征值上界_20261102.pptx")
prs.save(OUT)
print("saved", OUT, "slides:", len(prs.slides.__iter__.__self__._sldIdLst))
