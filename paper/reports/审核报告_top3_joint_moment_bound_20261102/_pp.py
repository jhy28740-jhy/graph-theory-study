import io
lines = io.open("审核报告_top3_joint_moment_bound.tex", encoding="utf-8").read().split("\n")
for i in range(1439, 1560):
    print(lines[i][:175])
