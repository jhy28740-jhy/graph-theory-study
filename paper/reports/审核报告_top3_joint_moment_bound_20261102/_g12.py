import io, re
s = io.open("审核报告_top3_joint_moment_bound.tex", encoding="utf-8").read()
lines = s.split("\n")
hits = [i for i, ln in enumerate(lines) if re.search(r"G1|G2", ln)]
print("hits:", len(hits))
seen = set()
for i in hits[:14]:
    a = max(0, i-1); b = min(len(lines), i+10)
    if a in seen: continue
    seen.add(a)
    print("----- line", i+1)
    for j in range(a, b):
        print("   ", lines[j][:170])
