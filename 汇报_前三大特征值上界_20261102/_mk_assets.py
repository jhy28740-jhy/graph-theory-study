import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, os

OUT = "assets"
plt.rcParams.update({
    "font.family": ["Microsoft YaHei"],
    "axes.unicode_minus": False,
    "mathtext.fontset": "cm",
    "figure.dpi": 200,
    "savefig.transparent": True,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.06,
})
INK = "#1A1A1A"
BLUE = "#1F3864"
BLUE2 = "#2F5597"
GREEN = "#2E7D32"
AMBER = "#B26A00"
RED = "#B3261E"
GREY = "#8A9099"

def eq(name, tex, size=24, color=INK):
    fig = plt.figure(figsize=(0.1, 0.1))
    t = fig.text(0, 0, tex, fontsize=size, color=color)
    fig.savefig(os.path.join(OUT, name))
    plt.close(fig)

eq("eq_stw.png", r"$\Lambda_3(G)\;\leq\;\dfrac{3+\sqrt{5}}{4}\,n-3$", 26)
eq("eq_main.png", r"$\Lambda_3(G)\;\leq\;c\,n-3$", 28)
eq("eq_law.png", r"$\mathbb{E}\,|X\cdot X'|\;+\;|\mathbb{E}\,X|^{2}\;\leq\;2c$", 26)
eq("eq_c.png", r"$c=\dfrac{106106902139519}{107310156800000}+\dfrac{41\sqrt{3}}{240}-\dfrac{53}{9000000}=1.284673255723918505\ldots$", 21)
eq("eq_ceiling.png", r"$\dfrac{29+10\sqrt{5}}{40}=\dfrac{3+\sqrt{5}}{4}-\dfrac{1}{40}=1.2840169943749\ldots$", 24)
eq("eq_rho0.png", r"$\rho_0=\dfrac{9+\sqrt{5}+2\sqrt{35}}{18}=1.2815681968\ldots$", 24)
eq("eq_D.png", r"$D=\dfrac{3+\sqrt{5}}{2}-A-B,\qquad d=\dfrac{3+\sqrt{5}}{2}-2c=0.0486874773\ldots$", 22)
eq("eq_AB.png", r"$A=\mathbb{E}\,|X\cdot X'|,\qquad B=|\mathbb{E}\,X|^{2},\qquad A\leq\varphi=\dfrac{1+\sqrt{5}}{2},\quad B\leq 1$", 22)

# ---- chart 1: slack composition (7 groups) ----
names = ["$r_{\\rm tr}$ 径向迹项", "$\\widehat{Q}$ 四次型", "$\\ell\\,r_{\\rm reg}$ 回归",
         "$\\mathcal{R}_0$ 储备", "$r_{\\rm off}/9$ 离对角", "$K^2/1125$", "$\\frac{\\kappa}{9}(\\chi-\\frac{C^2}{24})$"]
vals = np.array([3.8e-3, 1.6e-3, 4.6e-4, 2.3e-4, 4.9e-5, 5.7e-6, 9.9e-7])
pct = np.array([61.46, 26.54, 7.44, 3.68, 0.78, 0.09, 0.02])
fig, ax = plt.subplots(figsize=(9.6, 4.3))
y = np.arange(len(names))[::-1]
cols = [RED, AMBER, BLUE2, BLUE2, GREY, GREY, GREY]
ax.barh(y, pct, color=cols, height=0.62)
for yy, p, v in zip(y, pct, vals):
    ax.text(p + 1.0, yy, f"{p:.2f}%   ({v:.1e})", va="center", fontsize=10.5, color=INK)
ax.set_yticks(y)
ax.set_yticklabels(names, fontsize=11.5)
ax.set_xlim(0, 78)
ax.set_xlabel("占冠军律处证书松驰的比例（松驰总量 $6.210\\times10^{-3}=2(c-\\rho_0)$）", fontsize=11)
ax.spines[["top", "right"]].set_visible(False)
ax.grid(axis="x", color="#E1E4E8", lw=0.8)
ax.set_axisbelow(True)
ax.set_title("证书八项（归并七组）在冠军律上的松驰归属", fontsize=13.5, color=BLUE, pad=10)
fig.savefig(os.path.join(OUT, "chart_slack.png"), transparent=False, facecolor="white")
plt.close(fig)

# ---- chart 2: constant ladder waterfall ----
labels = ["STW\n$(3+\\sqrt{5})/4$", "已证改进\n$-0.024344$", "本文 $c$", "P1 回收\n$-0.000656$", "天花板", "P2 余量\n$-0.002449$", "$\\rho_0$"]
start = 1.309016994374948
steps = [None, -0.024343738651029, None, -0.000656261348971, None, -0.002448797502780, None]
levels = [start, None, 1.284673255723919, None, 1.284016994374948, None, 1.281568196872168]
fig, ax = plt.subplots(figsize=(10.2, 4.4))
bottoms, heights, colors = [], [], []
cur = start
for i, s in enumerate(steps):
    if s is None:
        bottoms.append(levels[i] - 0.004)
        heights.append(0.004)
        colors.append(BLUE if i in (0, 2) else (GREEN if i == 4 else AMBER))
    else:
        bottoms.append(cur + s)
        heights.append(-s)
        colors.append("#C9D3E4")
        cur = cur + s
x = np.arange(len(labels))
ax.bar(x, heights, bottom=bottoms, color=colors, width=0.62)
for i, (b, h) in enumerate(zip(bottoms, heights)):
    top = b + h
    if steps[i] is None:
        ax.text(i, top + 0.0006, f"{levels[i]:.6f}", ha="center", fontsize=10.5, color=INK)
    else:
        ax.text(i, top + h / 2, f"{steps[i]:+.6f}", ha="center", va="center", fontsize=9.5, color="#4A5568")
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=10.5)
ax.set_ylim(1.278, 1.313)
ax.set_ylabel("上界常数", fontsize=11)
ax.spines[["top", "right"]].set_visible(False)
ax.grid(axis="y", color="#E1E4E8", lw=0.8)
ax.set_axisbelow(True)
ax.set_title("常数阶梯：$\\mathrm{STW}\\to\\rho_0$ 总差 $0.027449$ 的三段归属（88.7% / 2.4% / 8.9%）", fontsize=13, color=BLUE, pad=10)
fig.savefig(os.path.join(OUT, "chart_ladder.png"), transparent=False, facecolor="white")
plt.close(fig)

# ---- chart 3: 8-term values at champion, log scale ----
fig, ax = plt.subplots(figsize=(9.6, 3.9))
short = ["$r_{\\rm tr}$", "$\\widehat{Q}$", "$\\ell r_{\\rm reg}$", "$\\mathcal{R}_0$", "$r_{\\rm off}/9$", "$K^2/1125$", "$\\frac{\\kappa}{9}(\\chi-\\frac{C^2}{24})$"]
ax.bar(np.arange(len(short)), vals, color=[RED, AMBER, BLUE2, BLUE2, GREY, GREY, GREY], width=0.6)
ax.set_yscale("log")
ax.set_xticks(np.arange(len(short)))
ax.set_xticklabels(short, fontsize=12)
for i, v in enumerate(vals):
    ax.text(i, v * 1.25, f"{v:.1e}", ha="center", fontsize=9.5, color=INK)
ax.set_ylabel("冠军律上的取值（对数轴）", fontsize=11)
ax.spines[["top", "right"]].set_visible(False)
ax.grid(axis="y", color="#E1E4E8", lw=0.8)
ax.set_axisbelow(True)
ax.set_title("紧证书的障碍：八项在冠军律上无一项为零", fontsize=13, color=RED, pad=8)
fig.savefig(os.path.join(OUT, "chart_terms.png"), transparent=False, facecolor="white")
plt.close(fig)

print("assets:", sorted(os.listdir(OUT)))
