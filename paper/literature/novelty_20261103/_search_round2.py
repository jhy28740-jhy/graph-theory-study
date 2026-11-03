import re, time, json, sys, urllib.parse, urllib.request
UA = {"User-Agent": "Mozilla/5.0 (novelty-check)"}
def get(url, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            sys.stderr.write("retry %s\n" % e); time.sleep(5+4*i)
    return ""

def absx(aid):
    t = get("https://export.arxiv.org/abs/" + aid)
    m = re.search(r'<blockquote class="abstract[^"]*">(.*?)</blockquote>', t, re.S)
    ti = re.search(r'<h1 class="title[^"]*">(.*?)</h1>', t, re.S)
    dt = re.search(r'\[Submitted on(.*?)\]', t, re.S)
    hist = re.findall(r'<b>\[v\d+\]</b>\s*([^<]*)', t)
    clean = lambda s: re.sub(r"<[^>]+>", "", s).strip() if s else ""
    return {"id": aid, "title": clean(ti.group(1) if ti else ""),
            "abstract": re.sub(r"\s+", " ", clean(m.group(1) if m else "")),
            "submitted": clean(dt.group(1) if dt else ""), "versions": hist}

for aid in ["2601.10036", "2610.22838", "2607.23081", "2609.23707", "1602.02069"]:
    d = absx(aid)
    print("\n=== %s (%s) ===" % (d["id"], d["submitted"][:40]))
    print("TITLE:", d["title"])
    print("VERS:", d["versions"])
    print("ABS:", d["abstract"][:1600])
    time.sleep(4)

print("\n\n########## extra arXiv queries ##########")
qs = [
 'abs:"lambda_1" AND abs:"lambda_2" AND abs:"extremal" AND cat:math.CO',
 'ti:"spectral sum" AND cat:math.CO',
 'ti:"eigenvalues" AND abs:"three-block" AND cat:math.CO',
 'abs:"Pell equation" AND abs:"eigenvalue" AND cat:math.CO',
 'abs:"stair" AND abs:"threshold graph" AND abs:"eigenvalues"',
 'abs:"largest eigenvalue" AND abs:"twice the second"',
 'abs:"weighted spectral" AND cat:math.CO',
 'ti:"sum of the largest" AND cat:math.CO',
]
for q in qs:
    u = ("https://export.arxiv.org/api/query?search_query=" + urllib.parse.quote(q) +
         "&start=0&max_results=15&sortBy=submittedDate&sortOrder=descending")
    t = get(u)
    print("\n## " + q)
    for m in re.finditer(r"<entry>(.*?)</entry>", t, re.S):
        e = m.group(1)
        g = lambda n: (re.search(r"<%s>(.*?)</%s>" % (n, n), e, re.S) or [None, ""])[1]
        print("  - %s %s %s" % (g("published")[:10], g("id").split("/abs/")[-1], re.sub(r"\s+"," ",g("title"))[:105]))
    time.sleep(4)
