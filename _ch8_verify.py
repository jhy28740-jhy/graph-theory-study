"""第 8 章 数值验证脚本（配合《第8章课程讲义_强正则图与设计.md》）"""
import itertools as it
from collections import Counter, defaultdict
from fractions import Fraction
from math import isqrt

import networkx as nx
import numpy as np


def adjacency(G):
    return nx.adjacency_matrix(G).toarray().astype(int)


def spectrum(G):
    uq, cnt = np.unique(np.round(np.linalg.eigvalsh(adjacency(G)), 6), return_counts=True)
    return dict(zip(uq.tolist(), cnt.tolist()))


def gates(n, d, a, c, name=""):
    """四道闸门 + 谱."""
    g1 = d * (d - a - 1) == (n - d - 1) * c
    disc = (a - c) ** 2 + 4 * (d - c)
    sq = isqrt(disc)
    t = (n - 1) * (a - c) + 2 * d
    a_bar, c_bar = n - 2 * d + c - 2, n - 2 * d + a
    print(f"  {name}({n},{d},{a},{c})")
    print(f"   gate1 Prop8.1.1: {d*(d-a-1)} vs {(n-d-1)*c} -> {'PASS' if g1 else 'FAIL'}")
    if sq * sq == disc and sq > 0:
        r, s = Fraction(a - c + sq, 2), Fraction(a - c - sq, 2)
        f = Fraction(-d - Fraction((n - 1) * (a - c - sq), 2), sq)
        g = Fraction(d + Fraction((n - 1) * (a - c + sq), 2), sq)
        ok = f.denominator == 1 and g.denominator == 1
        print(f"   gate2 disc={disc} 完全平方, r={r}, s={s}")
        print(f"   gate3 f={f}, g={g}, f+g={f+g} (n-1={n-1}) -> {'PASS' if ok else 'FAIL'}")
    else:
        print(f"   gate2 disc={disc} 非完全平方; (n-1)(a-c)+2d={t} -> "
              f"{'quadratic type (f=g=(n-1)/2)' if t == 0 else 'FAIL'}")
    print(f"   gate4 补图=({n},{n-d-1},{a_bar},{c_bar}) -> "
          f"{'PASS' if min(a_bar, c_bar, n - d - 1) >= 0 else 'FAIL'}")


def srg_report(G, name):
    A = adjacency(G); A2 = A @ A
    n = len(A)
    deg = A.sum(1)
    if len(set(deg.tolist())) != 1:
        print(f"{name}: 不正则"); return
    d = int(deg[0])
    av = {int(A2[u][v]) for u in range(n) for v in range(n) if u != v and A[u][v]}
    cv = {int(A2[u][v]) for u in range(n) for v in range(n) if u != v and not A[u][v]}
    is_srg = len(av) == 1 and len(cv) == 1
    print(f"{name}: n={n} d={d} a值={sorted(av)} c值={sorted(cv)} 强正则={is_srg}")
    if is_srg:
        gates(n, d, av.pop(), cv.pop(), "   参数")
        print(f"   谱={spectrum(G)} 直径={nx.diameter(G)} 三角形数="
              f"{sum(1 for t in it.combinations(G.nodes, 3) if all(G.has_edge(*e) for e in it.combinations(t, 2)))}")


def link_graph(G, v):
    return G.subgraph(list(G.neighbors(v)))


def intersection_array(G):
    src = list(G.nodes)[0]
    dist = nx.single_source_shortest_path_length(G, src)
    D = max(dist.values())
    k = [0] * (D + 1)
    for dv in dist.values():
        k[dv] += 1
    a, b, c, ok = [], [], [], True
    for i in range(D + 1):
        cnt = {(Counter(dist.get(y, -1) for y in G.neighbors(x))[i - 1],
                Counter(dist.get(y, -1) for y in G.neighbors(x))[i],
                Counter(dist.get(y, -1) for y in G.neighbors(x))[i + 1])
               for x in (v for v, dv in dist.items() if dv == i)}
        if len(cnt) != 1:
            ok = False; a.append(None); b.append(None); c.append(None); continue
        ci, ai, bi = cnt.pop()
        a.append(ai); b.append(bi); c.append(ci)
    print(f"   距离正则={ok} D={D} k_i={k} sum={sum(k)} |V|={G.number_of_nodes()}")
    print(f"   a_i={a} b_i={b[:D]} c_i={c[1:]} 交数组={{{','.join(map(str,b[:D]))};{','.join(map(str,c[1:]))}}}")
    if ok:
        print(f"   (8.5) 成立={all(k[i]*b[i] == k[i+1]*c[i+1] for i in range(D))}")
    print(f"   不同特征值个数={len(spectrum(G))} (D+1={D+1})")


def design_report(B, name):
    A = adjacency(B); A2 = A @ A
    nodes = list(B)
    idx = {v: i for i, v in enumerate(nodes)}
    col = nx.bipartite.color(B)
    vals = Counter()
    for p in (0, 1):
        for u, v in it.combinations([w for w in nodes if col[w] == p], 2):
            vals[int(A2[idx[u]][idx[v]])] += 1
    m, d = len(nodes) // 2, int(A.sum(1)[0])
    cs = sorted(vals)
    kind = f"设计图 (m={m},d={d},c={cs[0]})" if len(cs) == 1 else f"部分设计图 (m={m},d={d},c1={cs[0]},c2={cs[-1]})"
    print(f"{name}: |V|={len(nodes)} {kind} 公共邻居值分布={dict(sorted(vals.items()))}")
    if len(cs) == 1:
        c = cs[0]
        print(f"   c(m-1)={c*(m-1)} vs d(d-1)={d*(d-1)} -> {c*(m-1) == d*(d-1)}; "
              f"极值(c=1)={c == 1}; Prop8.3.1 预测围长={'4' if c > 1 else '6'}")
    print(f"   直径={nx.diameter(B)} 围长={nx.girth(B)} 谱={spectrum(B)}")


def bip_double(G):
    B = nx.Graph()
    B.add_nodes_from([(v, 0) for v in G] + [(v, 1) for v in G])
    for u, v in G.edges:
        B.add_edge((u, 0), (v, 1)); B.add_edge((u, 1), (v, 0))
    return B


def rook(n):
    G = nx.Graph(); G.add_nodes_from([(i, j) for i in range(n) for j in range(n)])
    for x, y in it.combinations(G.nodes, 2):
        if x[0] == y[0] or x[1] == y[1]:
            G.add_edge(x, y)
    return G


def shrikhande():
    G = nx.Graph(); G.add_nodes_from([(i, j) for i in range(4) for j in range(4)])
    for (x, y) in G.nodes:
        for dx, dy in [(1, 0), (0, 1), (1, 1)]:
            G.add_edge((x, y), ((x + dx) % 4, (y + dy) % 4))
            G.add_edge((x, y), ((x - dx) % 4, (y - dy) % 4))
    return G


def X_G(mods):
    """Thm 8.1.1: Cay(G×G, {(s,1),(1,s),(s,s): s≠1})."""
    G = list(it.product(*[range(m) for m in mods]))
    zero = tuple(0 for _ in mods)
    add = lambda x, y: tuple((p + q) % m for p, q, m in zip(x, y, mods))
    S = [s for s in G if s != zero]
    X = nx.Graph(); X.add_nodes_from([(x, y) for x in G for y in G])
    for (x, y) in X.nodes:
        for s in S:
            X.add_edge((x, y), (add(x, s), y))
            X.add_edge((x, y), (x, add(y, s)))
            X.add_edge((x, y), (add(x, s), add(y, s)))
    return X, G, zero, S, add


def tutte_coxeter():
    V = list(range(6))
    def pm(vs):
        if not vs:
            yield frozenset(); return
        a = vs[0]
        for b in vs[1:]:
            for rest in pm([x for x in vs if x not in (a, b)]):
                yield rest | {frozenset((a, b))}
    M = sorted(pm(V))
    T = nx.Graph()
    T.add_nodes_from([("e", tuple(sorted(e))) for e in it.combinations(V, 2)] +
                     [("m", i) for i in range(len(M))])
    for i, m in enumerate(M):
        for e in m:
            T.add_edge(("m", i), ("e", tuple(sorted(e))))
    return T


def incidence_graph(blocks, v):
    B = nx.Graph()
    B.add_nodes_from([("p", i) for i in range(v)] + [("B", j) for j in range(len(blocks))])
    for j, blk in enumerate(blocks):
        for x in blk:
            B.add_edge(("p", x), ("B", j))
    return B


def design_check(blocks, v, k, lam, name):
    b = len(blocks)
    N = np.zeros((v, b), int)
    for j, blk in enumerate(blocks):
        for x in blk:
            N[x, j] = 1
    r = int(N.sum(1)[0])
    pairs = Counter(tuple(sorted(p)) for blk in blocks for p in it.combinations(blk, 2))
    print(f"{name}: v={v} k={k} λ={lam} b={b} r={r} 点对覆盖值={sorted(set(pairs.values()))}")
    print(f"   (8.6): vr==bk -> {v*r == b*k}; r(k-1)==λ(v-1) -> {r*(k-1) == lam*(v-1)}")
    Iv, Jv = np.eye(v, dtype=int), np.ones((v, v), int)
    print(f"   (8.7): NN^T==(r-λ)I+λJ -> {np.array_equal(N @ N.T, (r-lam)*Iv + lam*Jv)}")
    if b == v:
        Ib, Jb = np.eye(b, dtype=int), np.ones((b, b), int)
        print(f"   对称 (b=v): r==k -> {r == k}; (S3) N^T N==(k-λ)I+λJ -> "
              f"{np.array_equal(N.T @ N, (k-lam)*Ib + lam*Jb)}; "
              f"(S2) 两块交 λ 点 -> "
              f"{set(len(set(x) & set(y)) for x, y in it.combinations(blocks, 2)) == {lam}}")
    G = incidence_graph(blocks, v)
    print(f"   关联图: |V|={G.number_of_nodes()} 直径={nx.diameter(G)} 围长={nx.girth(G)}")
    print(f"   Cor8.5.1 预测: ±{k}, ±sqrt({k-lam})={np.sqrt(k-lam):.6f}, 重数 v-1={v-1}")
    print(f"   实测谱: {spectrum(G)}")
    design_report(G, "   作为设计图")


def projective_plane(q):
    """PG(2,q) 的点 = F_q^3 的 1-维子空间, 线 = 2-维子空间 (q 为素数时可用)."""
    vecs = [p for p in it.product(range(q), repeat=3) if any(p)]
    def canon(p):
        nz = next(i for i in range(3) if p[i])
        inv = pow(p[nz], -1, q)
        return tuple((x * inv) % q for x in p)
    pts = sorted({canon(p) for p in vecs})
    idxp = {p: i for i, p in enumerate(pts)}
    lines = set()
    for a, b in it.combinations(pts, 2):
        line = {canon(tuple((s * a[i] + t * b[i]) % q for i in range(3)))
                for s, t in it.product(range(q), repeat=2) if (s, t) != (0, 0)}
        lines.add(frozenset(idxp[p] for p in line))
    lines = [sorted(l) for l in lines]
    return pts, lines


print("=" * 72); print("1. Ex 8.1.1-8.1.5 例子群"); print("=" * 72)
for G, nm in [(nx.petersen_graph(), "Petersen"), (nx.cycle_graph(5), "C5"),
              (nx.cycle_graph(4), "C4"), (nx.complete_bipartite_graph(4, 4), "K4,4"),
              (rook(4), "rook K4xK4"), (shrikhande(), "Shrikhande")]:
    srg_report(G, nm)

print("\n" + "=" * 72); print("2. 双胞胎: link graph 作为二阶不变量"); print("=" * 72)
R, S = rook(4), shrikhande()
print(f"rook ~ Shrikhande ? {nx.is_isomorphic(R, S)}")
for G, nm in [(R, "rook K4xK4"), (S, "Shrikhande")]:
    L = link_graph(G, list(G.nodes)[0])
    tri = sum(1 for t in it.combinations(L.nodes, 3) if all(L.has_edge(*e) for e in it.combinations(t, 2)))
    print(f"  {nm}: link |V|={L.number_of_nodes()} |E|={L.number_of_edges()} 度={sorted({d for _,d in L.degree})} "
          f"分支={nx.number_connected_components(L)} 三角形={tri} 是C6={nx.is_isomorphic(L, nx.cycle_graph(6))}")

print("\n" + "=" * 72); print("3. Thm 8.1.1 / Cor 8.1.1: 三岛与桥三角形"); print("=" * 72)
graphs = {}
for mods, nm in [((4,), "Z4"), ((2, 2), "Z2xZ2")]:
    X, G, zero, S, add = X_G(mods); graphs[nm] = X
    n = len(G)
    invol = [s for s in G if s != zero and add(s, s) == zero]
    print(f"{nm}: |G|={n} 2阶元={len(invol)} 期望参数=({n*n},{3*n-3},{n},6)")
    srg_report(X, f"  X_{nm}")
    L = link_graph(X, (zero, zero)); nb = list(X.neighbors((zero, zero)))
    S1, S2, S3 = {(s, zero) for s in S}, {(zero, s) for s in S}, {(s, s) for s in S}
    isl = lambda v: 1 if v in S1 else (2 if v in S2 else 3)
    tris = [t for t in it.combinations(nb, 3) if all(L.has_edge(*e) for e in it.combinations(t, 2))]
    within = sum(1 for t in tris if len({isl(x) for x in t}) == 1)
    print(f"  link: |V|={L.number_of_nodes()} |E|={L.number_of_edges()} 度={sorted({d for _,d in L.degree})}")
    print(f"  link 三角形: 总={len(tris)} 岛内={within} (预测 3C(n-1,3)={3*len(list(it.combinations(range(n-1),3)))}) "
          f"桥型={len(tris)-within} (预测 #2阶元={len(invol)}) -> {len(tris)-within == len(invol)}")
print(f"X_Z4 ~ X_Z2xZ2 ? {nx.is_isomorphic(graphs['Z4'], graphs['Z2xZ2'])}")
CR = nx.complement(rook(4))
print(f"complement(rook4) ~ X_Z2xZ2 ? {nx.is_isomorphic(CR, graphs['Z2xZ2'])} | ~ X_Z4 ? {nx.is_isomorphic(CR, graphs['Z4'])}")

print("\n" + "=" * 72); print("4. 四道闸门"); print("=" * 72)
for t in [(7, 4, 2, 2, "BAD"), (50, 7, 0, 1, "Hoffman-Singleton"), (3250, 57, 0, 1, "open"),
          (15, 6, 1, 3, "compl T(6)"), (21, 10, 1, 6, "gate1-fail"), (28, 12, 6, 4, "ok?"),
          (36, 15, 6, 6, "ok?"), (16, 9, 4, 6, "X_G |G|=4"), (16, 6, 2, 2, "twins")]:
    gates(*t[:4], t[4] + " ")

print("\n" + "=" * 72); print("5. §8.3 设计图与部分设计图"); print("=" * 72)
fano = [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]]
design_report(incidence_graph(fano, 7), "Heawood = Fano 关联图")
design_report(nx.cycle_graph(6), "C6")
design_report(bip_double(rook(4)), "B(rook4)")
design_report(bip_double(shrikhande()), "B(Shrikhande)")
print(f"B(rook4) ~ B(Shrikhande) ? {nx.is_isomorphic(bip_double(rook(4)), bip_double(shrikhande()))}")
design_report(bip_double(nx.petersen_graph()), "Desargues = B(Petersen)")
for n in (3, 4, 5):
    design_report(nx.hypercube_graph(n), f"Q{n}")
design_report(tutte_coxeter(), "Tutte-Coxeter")

print("\n" + "=" * 72); print("6. 设计图参数枚举 (m<=20, d<=6) 与 Moore 界"); print("=" * 72)
for m in range(2, 21):
    for d in range(2, min(6, m - 1) + 1):
        if d * (d - 1) % (m - 1) == 0 and d * (d - 1) // (m - 1) >= 1:
            c = d * (d - 1) // (m - 1)
            print(f"  (m={m},d={d},c={c}) d-c={d-c} 极值(c=1)={c == 1} m=d^2-d+1? {m == d*d-d+1}")
for d in range(2, 8):
    m = d * d - d + 1
    print(f"  d={d}: c=1 -> m={m}, |V|=2m={2*m}; Moore(围长6)={2*(1+(d-1)+(d-1)**2)} 相等={2*m == 2*(1+(d-1)+(d-1)**2)}")

print("\n" + "=" * 72); print("7. §8.4 距离正则图与交数组"); print("=" * 72)
for G, nm in [(nx.cycle_graph(6), "C6"), (incidence_graph(fano, 7), "Heawood"),
              (nx.petersen_graph(), "Petersen"), (nx.complete_bipartite_graph(4, 4), "K4,4"),
              (rook(4), "rook4"), (shrikhande(), "Shrikhande"), (graphs["Z2xZ2"], "X_Z2xZ2"),
              (nx.hypercube_graph(3), "Q3"), (nx.hypercube_graph(4), "Q4"),
              (nx.hypercube_graph(5), "Q5"), (nx.desargues_graph(), "Desargues"),
              (tutte_coxeter(), "Tutte-Coxeter")]:
    print(f"-- {nm}"); intersection_array(G)

print("\n" + "=" * 72); print("8. §8.5 区组设计与关联矩阵"); print("=" * 72)
design_check(fano, 7, 3, 1, "Fano 平面 2-(7,3,1)")
QR = sorted({(x * x) % 11 for x in range(1, 11)})
paley = [sorted({(y + t) % 11 for y in QR}) for t in range(11)]
print(f"QR mod 11 = {QR}; 差集检验 = {dict(sorted(Counter((a-b) % 11 for a, b in it.permutations(QR, 2)).items()))}")
design_check(paley, 11, 5, 2, "Paley 双平面 2-(11,5,2)")

print("\n" + "=" * 72); print("9. Ex 8.5.2: PG(2,q) 射影平面"); print("=" * 72)
for q in (2, 3, 4, 5):
    try:
        pts, lines = projective_plane(q)
    except Exception as e:
        print(f"  q={q}: 跳过 ({e})"); continue
    v, k, lam = q * q + q + 1, q + 1, 1
    print(f"  q={q}: 点数={len(pts)} 线数={len(lines)} 每线点数={sorted({len(l) for l in lines})} "
          f"期望 v=q²+q+1={v}, k=q+1={k}")
    if len(pts) == v and len(lines) == v and all(len(l) == k for l in lines):
        pairs = Counter(tuple(sorted(p)) for l in lines for p in it.combinations(l, 2))
        print(f"     点对覆盖值={sorted(set(pairs.values()))} -> 2-({v},{k},{lam}) 设计: {set(pairs.values()) == {lam}}")
        design_check(lines, v, k, lam, f"     PG(2,{q})")

