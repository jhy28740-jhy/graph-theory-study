import io
p = "证书_N6全局精确认证.py"
s = io.open(p, encoding="utf-8").read()
frag = io.open("_fragDG.py", encoding="utf-8").read()
marker = 'print()\nprint("D--G 段见下（另一次运行）")\nprint("=" * 88)\nprint("A--C 段总计 %d 项检查，全部通过: %s" % (len(ok), all(ok)))\nsys.exit(0 if all(ok) else 1)\n'
assert marker in s, "marker not found"
s = s.replace(marker, frag)
io.open(p, "w", encoding="utf-8").write(s)
print("spliced, len =", len(s))
