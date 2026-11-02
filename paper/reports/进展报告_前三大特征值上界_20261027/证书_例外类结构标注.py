# -*- coding: utf-8 -*-
r"""证书（数据标注）：Omega > rho0 的例外类的结构标注 —— 完全多部 / 含诱导 C5 / graph6 解码
      (2026-10-29, 第六版新增; 报告 命题 prop:tab19 与 注 rmk:tab19read 的唯一数据来源)

做什么
    读入 _danger23.json（23 个 Omega > rho0 的同构类：1 个 N=5 + 22 个 N=6），
    用【标准 graph6 解码】(ord(ch)-63, 不依赖 networkx 的读入) 还原邻接矩阵，然后判定：
      (a) 是否完全多部  <=> 补图的每个连通分支都是团（并输出各部大小，降序）；
      (b) 是否含诱导 C5 <=> 某 5 元子集上恰有 5 条边且每个顶点度 2；
      (c) 度序列（降序）。
    再扣除本轮已闭合的三类（N=6 冠军、K_{2,4}、K_{3,3}），输出剩下 19 类的 LaTeX 表格行。

产物（写在课题根目录）
    _exc23_struct.json : 23 类全表（含 N=5 那一类）
    _exc19_open.json   : 19 个未闭合类
    _tab19.tex         : 命题 prop:tab19 的表格行（graph6 特殊字符已转义）

运行（在课题根目录）：
  C:\Python314\python.exe -X utf8 -u 进展报告_前三大特征值上界_20261027\证书_例外类结构标注.py
"""
import os, io, json, itertools

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NFAIL = []
NOK = 0


def chk(ok, msg):
    global NOK
    print(("  [OK]   " if ok else "  [FAIL] ") + msg)
    if ok:
        NOK += 1
    else:
        NFAIL.append(msg)


def decode_g6(g6):
    """标准 graph6 解码：首字符 ord-63 = 顶点数；其后每字符 6 bit（高位在前），
       按 (j,i) i<j 的上三角列优先顺序填充。返回 (n, A)。"""
    s = g6.split(">>graph6<<")[-1]
    n = ord(s[0]) - 63
    bits = []
    for ch in s[1:]:
        v = ord(ch) - 63
        bits += [(v >> k) & 1 for k in (5, 4, 3, 2, 1, 0)]
    A = [[0] * n for _ in range(n)]
    k = 0
    for j in range(1, n):
        for i in range(j):
            A[i][j] = A[j][i] = bits[k]
            k += 1
    m = n * (n - 1) // 2
    assert k == m, ("graph6 位数不符", k, m)
    # graph6 把位串补齐到 6 的倍数, 故 len(bits) 可以比 m 多出至多 5 个填充位
    assert 0 <= len(bits) - m < 6, ("graph6 填充位异常", len(bits), m)
    return n, A


def components(n, A):
    seen = [False] * n
    comps = []
    for s0 in range(n):
        if seen[s0]:
            continue
        st = [s0]
        seen[s0] = True
        c = []
        while st:
            u = st.pop()
            c.append(u)
            for v in range(n):
                if A[u][v] and not seen[v]:
                    seen[v] = True
                    st.append(v)
        comps.append(sorted(c))
    return comps


def is_clique(A, c):
    return all(A[i][j] for i, j in itertools.combinations(c, 2))


def complement(n, A):
    return [[0 if i == j else 1 - A[i][j] for j in range(n)] for i in range(n)]


def has_induced_C5(n, A):
    for S in itertools.combinations(range(n), 5):
        sub = [[A[i][j] for j in S] for i in S]
        if sum(sub[i][j] for i, j in itertools.combinations(range(5), 2)) != 5:
            continue
        if all(sum(sub[i]) == 2 for i in range(5)):
            return True
    return False


def esc(g):
    """graph6 码里的 LaTeX 特殊字符转义（^ ~ { } _ # & % $）。"""
    out = ""
    for ch in g:
        if ch == "^":
            out += "\\^{}"
        elif ch == "~":
            out += "\\textasciitilde{}"
        elif ch == "{":
            out += "\\{"
        elif ch == "}":
            out += "\\}"
        elif ch == "_":
            out += "\\_"
        elif ch in "#&%$":
            out += "\\" + ch
        else:
            out += ch
    return out


def main():
    rec = json.load(io.open(os.path.join(ROOT, "_danger23.json"), encoding="utf-8"))
    chk(len(rec) == 23, "_danger23.json 含 %d 条记录（应为 23）" % len(rec))
    rows = []
    for idx, r in enumerate(rec):
        n, A = decode_g6(r["g6"])
        deg = sorted((sum(A[i]) for i in range(n)), reverse=True)
        if "deg[" in r["desc"]:
            rec_deg = sorted((int(t) for t in
                              r["desc"].split("deg[")[1].split("]")[0].replace(" ", "").split(",")),
                             reverse=True)
            chk(deg == rec_deg,
                "idx %2d: graph6 解码的度序列 %s 与原始记录 %s 一致" % (idx, deg, rec_deg))
        chk(sum(deg) // 2 == r["edges"], "idx %2d: 边数 %d 与原始记录一致" % (idx, sum(deg) // 2))
        C = complement(n, A)
        comps = components(n, C)
        mp = all(is_clique(C, c) for c in comps)
        sizes = sorted((len(c) for c in comps), reverse=True) if mp else []
        c5 = has_induced_C5(n, A)
        rows.append(dict(idx=idx, N=n, g6=r["g6"].split(">>graph6<<")[-1], deg=deg,
                         Omega=r["Omega"], F=r["F"], sc=r["sc"], edges=r["edges"],
                         multipartite=mp, sizes=sizes, C5=c5, desc=r["desc"]))
    io.open(os.path.join(ROOT, "_exc23_struct.json"), "w", encoding="utf-8").write(
        json.dumps(rows, ensure_ascii=False, indent=1))

    n6 = [r for r in rows if r["N"] == 6]
    n5 = [r for r in rows if r["N"] == 5]
    chk(len(n5) == 1 and len(n6) == 22, "N=5 例外 %d 个、N=6 例外 %d 个" % (len(n5), len(n6)))
    chk(n5[0]["multipartite"] and n5[0]["sizes"] == [3, 2],
        "N=5 的唯一例外是完全二部 K_{%s}（= K_{2,3}）=> 定理 thm:bipglobal 无条件关闭"
        % ",".join(map(str, n5[0]["sizes"])))
    mp6 = [(r["idx"], r["sizes"]) for r in n6 if r["multipartite"]]
    chk(len(mp6) == 4, "N=6 的 22 个例外中完全多部共 %d 个: %s" % (len(mp6), mp6))
    bip6 = [(r["idx"], r["sizes"]) for r in n6 if r["multipartite"] and len(r["sizes"]) == 2]
    chk(sorted(bip6) == [(16, [4, 2]), (17, [3, 3])],
        "其中完全【二部】恰 2 个: %s（= K_{2,4}, K_{3,3}）=> 本轮认证 F=5/4 精确, 例外 22->20" % bip6)
    c56 = [r["idx"] for r in n6 if r["C5"]]
    chk(len(c56) == 3, "N=6 例外中含诱导 C5 的共 %d 个: %s（含冠军 idx 0）" % (len(c56), c56))

    # 已闭合的三类：冠军（thm:champglobal）、K_{2,4}、K_{3,3}
    HANDLED = {0, 16, 17}
    open19 = sorted([r for r in n6 if r["idx"] not in HANDLED], key=lambda r: -r["Omega"])
    chk(len(open19) == 19, "扣除冠军与两个完全二部类后, 未闭合 %d 类" % len(open19))
    chk(all(r["sc"] == 0 for r in open19), "19 类在数值 F-极大化点处全部【不自洽】(sc=0)")
    chk(all(r["F"] <= (9 + 6 ** 0.5 * 4) / 15 + 1e-9 for r in open19),
        "19 类数值 max F 全部 <= (9+4sqrt6)/15 = kappa_3(5) < rho0")
    dom = [r for r in open19 if max(r["deg"]) == 5]
    chk(len(dom) == 2, "19 类中有支配点（度 5）的 %d 个: %s；其余 %d 个无支配点"
        % (len(dom), [r["g6"] for r in dom], 19 - len(dom)))
    mp19 = [r for r in open19 if r["multipartite"]]
    c519 = [r for r in open19 if r["C5"]]
    chk(len(mp19) == 2 and len(c519) == 2,
        "19 类中完全多部 %d 个 %s、含诱导 C5 %d 个 %s、两者都不是 %d 个"
        % (len(mp19), [r["g6"] for r in mp19], len(c519), [r["g6"] for r in c519],
           19 - len(set([r["idx"] for r in mp19]) | set([r["idx"] for r in c519]))))
    from collections import Counter
    cnt = Counter("%.9f" % r["F"] for r in open19)
    print("      数值 max F 的分布: " + ", ".join("%s x%d" % (k, v) for k, v in sorted(cnt.items(), reverse=True)))

    io.open(os.path.join(ROOT, "_exc19_open.json"), "w", encoding="utf-8").write(
        json.dumps(open19, ensure_ascii=False, indent=1))
    lines = []
    for i, r in enumerate(open19, 1):
        tag = ("$K_{%s}$" % ",".join(map(str, r["sizes"]))) if r["multipartite"] else "--"
        lines.append("%d & \\texttt{%s} & $%s$ & $%.6f$ & $%.9f$ & %s & %s \\\\"
                     % (i, esc(r["g6"]), ",".join(map(str, r["deg"])), r["Omega"], r["F"],
                        tag, "是" if r["C5"] else "--"))
    io.open(os.path.join(ROOT, "_tab19.tex"), "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    chk(len(lines) == 19, "生成 _tab19.tex 共 %d 行（命题 prop:tab19 的表格）" % len(lines))

    print("\n  总计 %d 项检查, 全部通过: %s" % (NOK, not NFAIL))
    print("  -> 写入 _exc23_struct.json / _exc19_open.json / _tab19.tex")
    if NFAIL:
        print("!!! 有 FAIL 项:")
        for m in NFAIL:
            print("   -", m)
        raise SystemExit(1)
    print("=== 全部 [OK] ===")


if __name__ == "__main__":
    main()
