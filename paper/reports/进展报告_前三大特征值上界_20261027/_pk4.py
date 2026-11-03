import io
for f in ["_frag_r25.tex","_frag_r26.tex"]:
    s=io.open(f,encoding="utf-8").read().split("\n")
    print("="*25,f)
    for i,l in enumerate(s,1):
        if "证书" in l or "_cert_r2" in l or "项全绿" in l:
            print("  ",i,":",l[:170])
