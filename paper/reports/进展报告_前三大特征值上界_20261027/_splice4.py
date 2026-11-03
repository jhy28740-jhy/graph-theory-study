import io,re
main="进展报告_前三大特征值之和的精确上界.tex"
s=io.open(main,encoding="utf-8").read()
orig_len=len(s.split("\n"))
anchor="\\section{本轮新结果：通道型主函数方案的精确天花板}\\label{sec:nogo}"
assert s.count(anchor)==1, s.count(anchor)
frags="".join(io.open("_frag_rowvar_%d.tex"%i,encoding="utf-8").read().rstrip("\n")+"\n\n" for i in (1,2,3,4))
newsec="\\section{本轮新结果（第 12 轮）：正部的行方差机制、天花板 $L(N)$，与 $Z_5$-等变弧的完全分类}\\label{sec:rowvar}"
frags=frags.replace("\\section{本轮新结果：正部的行方差机制、天花板 $L(N)$，与旧队首 P1 的否证}\\label{sec:rowvar}",newsec,1)
assert newsec in frags
replacement=frags+"\\section{通道型主函数方案的精确天花板}\\label{sec:nogo}"
s=s.replace(anchor,replacement,1)
io.open(main,"w",encoding="utf-8").write(s)
print("spliced: %d -> %d lines"%(orig_len,len(s.split("\n"))))
