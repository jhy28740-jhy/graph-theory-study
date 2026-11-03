import io, re
s = io.open("讲稿与排版建议.tex", encoding="utf-8").read()
for m in re.finditer(r"^\\(section|subsection)\*?\{(.{0,40})", s, re.M):
    print(m.group(0)[:60])
print("chars", len(s))
