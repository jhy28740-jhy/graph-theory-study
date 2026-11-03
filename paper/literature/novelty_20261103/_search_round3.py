import json, time, sys, urllib.parse, urllib.request
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/126 Safari/537.36"}
def get(u, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            sys.stderr.write("retry %s\n" % e); time.sleep(4+3*i)
    return ""

print("########## ZENODO ##########")
for q in ['"spectral sum" graph eigenvalue', '"lambda_1" "lambda_2" graph extremal',
          'third spectral sum graphs', 'threshold graph eigenvalues extremal weighted',
          '"Newton Astro Labs"']:
    u = "https://zenodo.org/api/records?size=20&sort=newest&q=" + urllib.parse.quote(q)
    t = get(u)
    print("\n## " + q)
    if not t:
        print("   (no response)"); continue
    try:
        d = json.loads(t)
    except Exception as e:
        print("   parse fail", e); continue
    print("   total:", d.get("hits", {}).get("total"))
    for h in d.get("hits", {}).get("hits", []):
        md = h.get("metadata", {})
        cr = md.get("creators", [])
        print("   - %s | %s | %s | %s" % (h.get("id"), md.get("publication_date"),
              (md.get("title") or "")[:100], [c.get("name") for c in cr][:3]))
    time.sleep(2)

print("\n\n########## OPENALEX title-only ##########")
for q in ['weighted spectral sum graphs', 'lambda1 plus two lambda2',
          'exact maximum weighted combination adjacency eigenvalues',
          'third spectral sum graphs exact bound']:
    u = ("https://api.openalex.org/works?filter=title.search:" + urllib.parse.quote(q) +
         "&per-page=12&select=id,doi,title,publication_year,primary_location")
    t = get(u)
    print("\n## " + q)
    if not t:
        print("   (no response)"); continue
    d = json.loads(t)
    print("   count:", d.get("meta", {}).get("count"))
    for w in d.get("results", []):
        src = ((w.get("primary_location") or {}).get("source") or {})
        print("   - [%s] %s | %s | %s" % (w.get("publication_year"), (w.get("title") or "")[:100],
              src.get("display_name"), w.get("doi")))
    time.sleep(1.5)
