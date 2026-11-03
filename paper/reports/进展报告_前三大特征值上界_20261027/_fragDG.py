# ---------- 共用工具 ----------
rng = np.random.default_rng(20261029)
def flip_target(Sig, w):
    M = (w[:, None] * Sig) * w[None, :]
    ev, U = np.linalg.eigh(M)
    P = U[:, -3:] @ U[:, -3:].T
    S2 = (P >= 0).astype(float); np.fill_diagonal(S2, 1.0)
    return S2, P, float(ev[-1] + ev[-2] + ev[-3]), M
def g6of(Sig):
    return nx.to_graph6_bytes(nx.from_numpy_array(Sig - np.eye(Sig.shape[0]))).decode()
def is_sc(Sig, ns=200, seed=13):
    v, p, ev, Q, z = _core.Fmax(Sig, ns=ns, seed=seed)
    w = np.sqrt(np.clip(p, 0, None)); w = w / np.linalg.norm(w)
    S2, P, lam3, M = flip_target(Sig, w)
    return int(np.array_equal(S2, Sig)), v, float(w.min()), P, M

print()
print("=" * 88)
print("D. lem:flip（翻转优占）：Lam_3(D_w Sigma D_w) <= F(Sigma')；w>0 且 Sigma'!=Sigma 时严格")
print("=" * 88)
worst = 1e9; nstrict = 0; nstrictneed = 0; ntest = 0
for _ in range(120):
    N = int(rng.integers(5, 8))
    H = nx.gnp_random_graph(N, float(rng.uniform(0.25, 0.75)), seed=int(rng.integers(1, 10 ** 6)))
    Sig = np.eye(N) + nx.to_numpy_array(H)
    w = rng.random(N); w = np.sqrt(w / w.sum())
    S2, P, lam3, M = flip_target(Sig, w)
    F2 = _core.Fmax(S2, ns=60, seed=5)[0]
    ntest += 1; worst = min(worst, F2 - lam3)
    if (not np.array_equal(S2, Sig)) and w.min() > 0:
        nstrictneed += 1; nstrict += int(F2 - lam3 > 1e-12)
chk("%d 个随机掩码 x 随机 w：Lam_3(D_w Sigma D_w) <= F(1{Pi_3>=0}) 全部成立，最小余量 %.3e >= 0"
    % (ntest, worst), worst > -1e-9)
chk("其中 Sigma'!=Sigma 且 w>0 的 %d 例全部严格（>1e-12）：%d/%d，与 lem:flip 注记的严格性判据一致"
    % (nstrictneed, nstrict, nstrictneed), nstrict == nstrictneed and nstrictneed >= 25)
print("  翻转链（从 21 个 Omega>rho0 的六原子例外出发；每步在自己的极大化点处翻转）：")
chain_ok = True; chain_rows = []
for g6 in BAD21:
    g = nx.from_graph6_bytes(g6.encode()); Sig = np.eye(6) + nx.to_numpy_array(g)
    cur = Sig; Fprev = -1.0; path = [g6]; mono = True
    for _ in range(10):
        vcur, p, ev, Q, z = _core.Fmax(cur, ns=200, seed=13)
        w = np.sqrt(np.clip(p, 0, None)); w = w / np.linalg.norm(w)
        S2, P, lam3, M = flip_target(cur, w)
        if vcur < Fprev - 1e-12: mono = False
        Fprev = vcur
        if np.array_equal(S2, cur): break
        if vcur >= Fprev - 1e-9 and (S2 != cur).any(): pass
        cur = S2; path.append(g6of(cur))
    scflag, vfin, wmin, Pf, Mf = is_sc(cur)
    incert = g6of(cur) in CERT6 or any(nx.is_isomorphic(nx.from_graph6_bytes(g6of(cur).encode()),
                    nx.from_graph6_bytes(c.encode())) for c in CERT6)
    chain_ok &= (scflag == 1)
    chain_rows.append(dict(start=g6, path=path, final=g6of(cur), F=vfin, sc=scflag, incert=incert, mono=mono))
    print("    %s -> %s   链长 %d   F(末)=%.9f  末端自洽=%d  末端已认证=%s"
          % (g6.replace(">>graph6<<", ""), " -> ".join(x.replace(">>graph6<<", "") for x in path[1:]) or "(不动)",
             len(path) - 1, vfin, scflag, incert))
chk("21 条翻转链全部终止于\\textbf{自洽}掩码（cor:flipsc 的有限降链在数值上完全复现），"
    "且链上 F 单调不减：%d/21" % sum(1 for r in chain_rows if r["sc"] == 1 and r["mono"]),
    all(r["sc"] == 1 and r["mono"] for r in chain_rows))
chk("21 条链的末端全部落在已认证的 134 类之中（=> 每条链给出 F(起点) <= F(末端) <= Omega(末端) <= rho0，"
    "唯一未被证明的环节是「数值求得的极大化点确为真极大化点」）：%d/21"
    % sum(1 for r in chain_rows if r["incert"]), all(r["incert"] for r in chain_rows))

print()
print("=" * 88)
print("E. lem:isolated（边界刚性）：w_i=0 => 翻转像的第 i 行除对角外全为 0")
print("=" * 88)
nbad = 0; ntest2 = 0
for _ in range(120):
    N = int(rng.integers(5, 8))
    H = nx.gnp_random_graph(N, float(rng.uniform(0.3, 0.8)), seed=int(rng.integers(1, 10 ** 6)))
    Sig = np.eye(N) + nx.to_numpy_array(H)
    i = int(rng.integers(0, N))
    w = rng.random(N); w[i] = 0.0; w = np.sqrt(w / w.sum())
    S2, P, lam3, M = flip_target(Sig, w); ntest2 += 1
    row = S2[i].copy(); row[i] = 0.0
    if np.abs(row).max() > 0: nbad += 1
chk("%d 次试验：把某个 w_i 置零后，1{Pi_3(D_w Sigma D_w)>=0} 的第 i 行除对角外\\textbf{全为 0}"
    "（违例 %d 次）——即 lem:isolated 的「零权原子必是孤立点」" % (ntest2, nbad), nbad == 0)
wbd = np.ones(6) / 6 ** 0.5; wbd[2] = 0.0; wbd = wbd / np.linalg.norm(wbd)
S2b, Pb, l3b, Mb = flip_target(Ac, wbd)
chk("冠军掩码上把第 3 个原子权置零：翻转像第 3 行除对角外全 0（%s），"
    "而 K1 v C5 无孤立点，故冠军的任何自洽见证必为正权" % ("是" if np.abs(S2b[2] - np.eye(6)[2]).max() == 0 else "否"),
    np.abs(S2b[2] - np.eye(6)[2]).max() == 0)

print()
print("=" * 88)
print("F. lem:chamber：Sigma 是 Parseval 框架的非锐角图 <=> Sigma=1{Q>=0}, Q in P_3(N)")
print("=" * 88)
vch, pch, evch, Qch, zch = _core.Fmax(Ac, ns=400, seed=20261029)
G = np.eye(6) - Qch                       # Gram(w_1..w_6)
evG = np.linalg.eigvalsh(G)
W = np.linalg.eigh(G)[1][:, 3:] @ np.diag(np.sqrt(np.clip(np.linalg.eigvalsh(G)[3:], 0, None)))  # 6x3
frame = W.T @ W
offd = G - np.diag(np.diag(G))
chk("冠军处 Q=Pi_3(M*) 满足 Q^2=Q（残差 %.2e）、tr Q=%.12f、Gram=I-Q 的特征值 %s（三个 0、三个 1）"
    % (np.abs(Qch @ Qch - Qch).max(), np.trace(Qch), np.array2string(evG, precision=2)),
    np.abs(Qch @ Qch - Qch).max() < 1e-12 and abs(np.trace(Qch) - 3) < 1e-12)
chk("由此得到的 6x3 矩阵 W 满足 Parseval 框架条件 sum_i w_i w_i^T = I_3（残差 %.2e），"
    "且 sign(-Gram) 的非对角部分与冠军掩码 Sigma0 逐条目相符（%d 处不符）"
    % (np.abs(frame - np.eye(3)).max(), int(((((offd <= 0).astype(float)) != (Ac - np.eye(6))).sum() / 2)),
    np.abs(frame - np.eye(3)).max() < 1e-11
    and int((((offd <= 0).astype(float)) != (Ac - np.eye(6))).sum() / 2) == 0)
nrm = np.diag(G)
chk("框架条件等价于 I3/3 in conv{u_i u_i^T}：权重 ||w_i||^2/3 = %s（和 = %.12f），"
    "且 sum_i (||w_i||^2/3) u_i u_i^T - I3/3 的范数 = %.2e"
    % (np.array2string(nrm / 3, precision=6), (nrm / 3).sum(),
       np.abs((W / np.sqrt(np.clip(nrm, 1e-300, None))[:, None]).T
              @ np.diag(nrm / 3) @ (W / np.sqrt(np.clip(nrm, 1e-300, None))[:, None]) - np.eye(3) / 3).max()),
    abs((nrm / 3).sum() - 1) < 1e-12)

print()
print("=" * 88)
print("G. prop:kappa6red：21 行例外表 + 球面采样下的翻转像符号模式扫描")
print("=" * 88)
print("  %-7s %8s %8s %3s %3s %3s %9s %9s %-8s %-22s" % (
    "g6", "F", "Omega", "al", "om", "dl", "min p", "signmgn", "flip Om", "flip in cert?"))
def alpha_omega(g):
    n = g.number_of_nodes(); a = w = 0
    for k in range(1, n + 1):
        for S in itertools.combinations(g.nodes(), k):
            if not any(g.has_edge(x, y) for x, y in itertools.combinations(S, 2)): a = max(a, k)
            if all(g.has_edge(x, y) for x, y in itertools.combinations(S, 2)): w = max(w, k)
    return a, w
tab = []
for g6 in BAD21:
    g = nx.from_graph6_bytes(g6.encode()); Sig = np.eye(6) + nx.to_numpy_array(g)
    a, w = alpha_omega(g); ds = sorted(d for _, d in g.degree())
    v, p, ev, Q, z = _core.Fmax(Sig, ns=250, seed=17)
    om = Omega_sdp(Sig)[0]
    pr = [(i, j) for i in range(6) for j in range(i + 1, 6)]
    sm = float(min(abs(Q[i, j]) for i, j in pr))
    S2, P, l3, M = flip_target(Sig, np.sqrt(np.clip(p, 0, None)) / np.linalg.norm(np.sqrt(np.clip(p, 0, None))))
    om2 = Omega_sdp(np.eye(6) + nx.to_numpy_array(nx.from_numpy_array(S2 - np.eye(6))))[0]
    g2 = nx.from_numpy_array(S2 - np.eye(6))
    incert = any(nx.is_isomorphic(g2, nx.from_graph6_bytes(c.encode())) for c in CERT6)
    # 球面采样：翻转像的符号模式集合
    pats = {}
    rs = np.random.default_rng(99)
    samples = [np.ones(6) / 6 ** 0.5, np.sqrt(np.clip(p, 0, None)) / np.linalg.norm(np.sqrt(np.clip(p, 0, None)))]
    for i in range(6):
        e = np.ones(6); e[i] = 0; samples.append(e / np.linalg.norm(e))
    samples += list(rs.dirichlet(np.ones(6) * 0.7, size=400) ** 0.5)
    samples += list((rs.dirichlet(np.ones(6) * 0.05, size=200)) ** 0.5)
    for ww in samples:
        S3, P3, l33, M3 = flip_target(Sig, ww)
        key = None
        for k in pats:
            if nx.is_isomorphic(nx.from_numpy_array(S3 - np.eye(6)), nx.from_numpy_array(k - np.eye(6))):
                key = k; break
        if key is None:
            pats[S3.copy()] = Omega_sdp(np.eye(6) + nx.to_numpy_array(nx.from_numpy_array(S3 - np.eye(6))))[0]
    pOm = list(pats.values())
    tab.append(dict(g6=g6, F=v, Om=om, alpha=a, omega=w, delta=min(ds), degseq=ds,
                    minp=float(p.min()), signmargin=sm, flipOm=om2, flipcert=incert,
                    npat=len(pats), maxpatOm=max(pOm)))
    print("  %-7s %8.5f %8.5f %3d %3d %3d %9.1e %9.1e %8.5f %-22s" % (
        g6.replace(">>graph6<<", ""), v, om, a, w, min(ds), p.min(), sm, om2,
        ("yes" if incert else "NO")))
    print("        球面采样 %d 个点 -> 翻转像的同构类 %d 个，其 Omega 的最大值 = %.6f (%s rho0)"
          % (len(samples), len(pats), max(pOm), "<=" if max(pOm) <= RHO0 else ">"))
chk("21 个例外全部满足：min p <= %.2e（极大化点在数值上落在单纯形的边界面上），"
    "且 delta(H)>=1 的有 %d 个（唯一 delta=0 的是 %s）"
    % (max(t["minp"] for t in tab), sum(1 for t in tab if t["delta"] >= 1),
       [t["g6"].replace(">>graph6<<", "") for t in tab if t["delta"] == 0][0]),
    max(t["minp"] for t in tab) < 1e-14 and sum(1 for t in tab if t["delta"] >= 1) == 20)
chk("21 个例外在各自数值极大化点处的翻转像\\textbf{全部}落在已认证的 134 类中（Omega<=rho0），"
    "最大 Omega(翻转像) = %.6f，比 rho0 低 %.3e" % (max(t["flipOm"] for t in tab),
    RHO0 - max(t["flipOm"] for t in tab)), all(t["flipcert"] for t in tab))
chk("球面采样（每类 608 个点，含顶点/面心/均匀/Dirichlet 两档）下，21 个例外的翻转像同构类共 %d 种，"
    "其 Omega 的最大值 = %.6f <= rho0（余量 %.3e）——即「翻转像恒为已认证类」这条足以推出 F<=rho0 的"
    "前提在数值上成立，剩下的只是对 S^5 的严格覆盖"
    % (sum(t["npat"] for t in tab), max(t["maxpatOm"] for t in tab), RHO0 - max(t["maxpatOm"] for t in tab)),
    max(t["maxpatOm"] for t in tab) <= RHO0)
chk("21 个例外的 max F = %.12f = kappa_3(5) = (9+4sqrt6)/15（差 %.2e），比 rho0 低 %.3e"
    % (max(t["F"] for t in tab), abs(max(t["F"] for t in tab) - float((9 + 4 * sp.sqrt(6)) / 15)),
       RHO0 - max(t["F"] for t in tab)),
    abs(max(t["F"] for t in tab) - float((9 + 4 * sp.sqrt(6)) / 15)) < 1e-11)
json.dump(tab, open(os.path.join(ROOT, "_bad21table.json"), "w"), indent=1)

print()
print("=" * 88)
print("总计 %d 项检查，全部通过: %s" % (len(ok), all(ok)))
print("=" * 88)
sys.exit(0 if all(ok) else 1)
