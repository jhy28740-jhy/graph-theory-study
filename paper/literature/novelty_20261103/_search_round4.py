import json, time, sys, urllib.parse, urllib.request
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/126 Safari/537.36"}
def get(u, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            sys.stderr.write("retry %s\n" % e); time.sleep(3+3*i)
    return ""

print("########## Zenodo record metadata ##########")
for rid in ["21785905", "21785906", "22082841", "22082842"]:
    t = get("https://zenodo.org/api/records/" + rid)
    if not t:
        print(rid, "-> no response"); continue
    try:
        d = json.loads(t)
    except Exception:
        print(rid, "-> not json:", t[:120]); continue
    md = d.get("metadata", {})
    print("\n--- %s ---" % rid)
    print(" doi:", md.get("doi"), "| related doi:", d.get("doi"))
    print(" title:", (md.get("title") or "")[:130])
    print(" pubdate:", md.get("publication_date"), "| created:", d.get("created"), "| updated:", d.get("updated"))
    print(" creators:", [(c.get("name"), c.get("affiliations")) for c in md.get("creators", [])])
    print(" stats:", d.get("stats"))
    print(" versions:", json.dumps(d.get("links", {}).get("versions", ""))[:80])
    print(" conceptrecid:", (d.get("links") or {}).get("latest"), "|", md.get("relations"))
    time.sleep(1.5)

print("\n\n########## OpenAlex title.search (lam1+2lam2 variants) ##########")
qs = ["first and second largest eigenvalues graphs",
      "eigenvalue sum extremal graph exact",
      "weighted eigenvalue sum graph maximum",
      "spectral extremal problem weighted eigenvalues",
      "three block graph eigenvalues",
      "stair graph eigenvalues extremal",
      "second largest eigenvalue maximum graphs exact",
      "Pell equation spectral graph theory",
      "pole pentagon spectral sum",
      "third largest eigenvalue sum graphs"]
for q in qs:
    u = ("https://api.openalex.org/works?filter=title.search:" + urllib.parse.quote(q) +
         "&per-page=10&select=id,doi,title,publication_year,primary_location")
    t = get(u)
    print("\n## " + q)
    if not t:
        print("   (no response)"); continue
    d = json.loads(t)
    print("   count:", d.get("meta", {}).get("count"))
    for w in d.get("results", []):
        src = ((w.get("primary_location") or {}).get("source") or {})
        print("   - [%s] %s | %s" % (w.get("publication_year"), (w.get("title") or "")[:95], src.get("display_name")))
    time.sleep(1.2)
