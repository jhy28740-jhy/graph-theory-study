import json, time, urllib.parse, urllib.request, sys

UA = {"User-Agent": "Mozilla/5.0 (research-novelty-check; mailto:research@example.org)"}

def get(url, tries=3):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            sys.stderr.write("retry %d %s : %s\n" % (i, url[:90], e))
            time.sleep(4 + 4 * i)
    return ""

def openalex(q, extra=""):
    u = ("https://api.openalex.org/works?search=" + urllib.parse.quote(q) +
         "&per-page=15&select=id,doi,title,publication_year,primary_location,cited_by_count,authorships" + extra)
    t = get(u)
    if not t:
        return []
    try:
        d = json.loads(t)
    except Exception:
        return []
    out = []
    for w in d.get("results", []):
        src = ((w.get("primary_location") or {}).get("source") or {})
        out.append({
            "title": w.get("title"),
            "year": w.get("publication_year"),
            "doi": w.get("doi"),
            "venue": src.get("display_name"),
            "cites": w.get("cited_by_count"),
            "auth": [a.get("author", {}).get("display_name") for a in (w.get("authorships") or [])][:5],
        })
    return out

def arxiv(q, mx=20):
    u = ("https://export.arxiv.org/api/query?search_query=" + urllib.parse.quote(q) +
         "&start=0&max_results=%d&sortBy=submittedDate&sortOrder=descending" % mx)
    t = get(u)
    out = []
    if not t:
        return out
    import re
    for m in re.finditer(r"<entry>(.*?)</entry>", t, re.S):
        e = m.group(1)
        def tag(name):
            mm = re.search(r"<%s>(.*?)</%s>" % (name, name), e, re.S)
            return re.sub(r"\s+", " ", mm.group(1)).strip() if mm else ""
        out.append({"id": tag("id"), "title": tag("title"), "pub": tag("published")[:10],
                    "auth": re.findall(r"<name>(.*?)</name>", e)[:4]})
    return out

queries_oa = [
 'lambda_1 + 2 lambda_2 graph eigenvalues extremal',
 'weighted sum of the largest and second largest adjacency eigenvalue',
 'maximum of lambda_1 + 2 lambda_2 over graphs',
 'spectral sum lambda1 plus two lambda2 extremal graph',
 'convex combination first and second adjacency eigenvalues',
 'extremal graphs for weighted spectral sums threshold graphs',
 'three-block graph spectral radius second eigenvalue exact maximum',
]
queries_ax = [
 'all:"lambda_1+2\\lambda_2"',
 'abs:"second largest eigenvalue" AND abs:"weighted" AND abs:"extremal graph"',
 'abs:"threshold graph" AND abs:"eigenvalues" AND abs:"extremal"',
 'abs:"spectral sum" AND cat:math.CO',
 'abs:"lambda_1" AND abs:"lambda_2" AND abs:"exact maximum" AND cat:math.CO',
]

res = {"openalex": {}, "arxiv": {}}
for q in queries_oa:
    res["openalex"][q] = openalex(q)
    time.sleep(1.5)
for q in queries_ax:
    res["arxiv"][q] = arxiv(q)
    time.sleep(4)

import os
outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(outdir, exist_ok=True)
with open(os.path.join(outdir, "lam12lam2_hits.json"), "w", encoding="utf-8") as f:
    json.dump(res, f, ensure_ascii=False, indent=1)

print("=== OpenAlex ===")
for q, hits in res["openalex"].items():
    print("\n## " + q)
    for h in hits[:8]:
        print("  - [%s] %s | %s | doi=%s | c=%s" % (h["year"], (h["title"] or "")[:110], h["venue"], h["doi"], h["cites"]))
print("\n\n=== arXiv ===")
for q, hits in res["arxiv"].items():
    print("\n## " + q)
    for h in hits[:12]:
        print("  - %s %s %s" % (h["pub"], h["id"].split("/abs/")[-1], (h["title"] or "")[:110]))
