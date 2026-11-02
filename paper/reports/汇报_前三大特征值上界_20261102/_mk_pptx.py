import io, os
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
g = {"__file__": os.path.join(os.getcwd(), "_mk_pptx.py")}
for f in ("_ppt_lib.py", "_ppt_slides.py"):
    src = io.open(f, encoding="utf-8").read()
    exec(compile(src, f, "exec"), g)
