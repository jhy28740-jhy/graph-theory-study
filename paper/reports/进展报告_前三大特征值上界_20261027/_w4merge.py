import io,re
a=io.open("_frag_rowvar_4.tex",encoding="utf-8").read()
b=io.open("_frag_rowvar_4b.tex",encoding="utf-8").read()
s=a.rstrip("\n")+"\n\n"+b.lstrip("\n")
io.open("_frag_rowvar_4.tex","w",encoding="utf-8").write(s)
for env in ["proposition","enumerate","equation","proof","remark","lemma","theorem","corollary","tabular","center","itemize","quote","subsubsection*"]:
    bb=len(re.findall(r"\\begin\{%s\}"%re.escape(env),s)); ee=len(re.findall(r"\\end\{%s\}"%re.escape(env),s))
    if bb!=ee: print("*** MISMATCH",env,bb,ee)
print("merged lines:",len(s.split("\n")),"braces bal:",s.count("{")-s.count("}"))
labels=re.findall(r"\\label\{([^}]*)\}",s)
print("labels:",len(labels),"dupes:",[x for x in set(labels) if labels.count(x)>1])
refs=re.findall(r"\\(?:ref|eqref)\{([^}]*)\}",s)
print("refs used:",sorted(set(refs)))
