import io, re
s = io.open("审核报告_top3_joint_moment_bound.tex", encoding="utf-8").read()
lines = s.split("\n")
for i, ln in enumerate(lines):
    if re.search(r"\\subsection\*\{S\d|\\subsubsection\*\{S\d|\\subsection\{S\d|S[1-7]（", ln):
        print(str(i+1).rjust(5), ln.strip()[:160])
print("==== P blocks ====")
for i, ln in enumerate(lines):
    if re.search(r"P[0-3]（", ln):
        print(str(i+1).rjust(5), ln.strip()[:170])
