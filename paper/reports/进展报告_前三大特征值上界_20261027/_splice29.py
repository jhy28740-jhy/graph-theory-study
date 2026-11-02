import io,sys,shutil
tex="进展报告_前三大特征值之和的精确上界.tex"
shutil.copyfile(tex, tex+".bak_r29")
s=open(tex,encoding="utf-8").read()
frag=open("_frag_r29.tex",encoding="utf-8").read()
anchor="\\section{通道型主函数方案的精确天花板}\\label{sec:nogo}"
assert s.count(anchor)==1, ("anchor count", s.count(anchor))
assert "\\label{sec:r29}" not in s, "r29 already present"
s2=s.replace(anchor, frag+"\n\n"+anchor, 1)
open(tex,"w",encoding="utf-8").write(s2)
print("SPLICE OK. old len",len(s),"new len",len(s2),"delta",len(s2)-len(s))
print("anchor now at char", s2.index(anchor))