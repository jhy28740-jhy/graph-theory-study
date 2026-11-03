"""Public scholarly literature audit. Raw responses and failures are retained."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import urllib.request, urllib.parse, urllib.error
import json, hashlib, datetime, re, xml.etree.ElementTree as ET

BASE = Path(__file__).resolve().parent
RAW = BASE / 'raw'
RAW.mkdir(exist_ok=True)
q = urllib.parse.quote
jobs = []
def add(key, source, query, url, kind):
    jobs.append(dict(key=key, source=source, query=query, url=url, kind=kind))
add('arxiv_seed', 'arXiv', 'id_list=2609.37980,2607.23081,2604.00512', 'https://export.arxiv.org/api/query?id_list=2609.37980,2607.23081,2604.00512', 'atom')
for i, query in enumerate([
    'all:"weighted sum" AND all:"eigenvalues" AND all:"graph"',
    'all:"linear combination" AND all:"eigenvalues" AND all:"graph"',
    'all:"spectral sum" AND all:"graph"',
]):
    add(f'arxiv_query{i}', 'arXiv', query, 'https://export.arxiv.org/api/query?search_query='+q(query)+'&start=0&max_results=100&sortBy=submittedDate&sortOrder=descending', 'atom')
for i, query in enumerate([
    'weighted sum largest two eigenvalues graph',
    'linear combination first second adjacency eigenvalues',
    'Quadratic inequalities between the largest eigenvalues of a graph',
    'The Exact Maximum of the Spectral Sum of Graphs',
    'spectral sum stability graph',
]):
    add(f'openalex{i}', 'OpenAlex', query, 'https://api.openalex.org/works?search='+q(query)+'&per-page=50', 'openalex')
for key, doi in [('crossref_2008','10.1016/j.laa.2008.06.016'), ('crossref_2015','10.1016/j.laa.2015.05.016'), ('crossref_2026618','10.3934/math.2026618')]:
    add(key, 'Crossref DOI metadata', doi, 'https://api.crossref.org/works/'+q(doi,safe=''), 'crossref')
for i, query in enumerate(['Quadratic inequalities between the largest eigenvalues of a graph','The Exact Maximum of the Spectral Sum of Graphs']):
    add(f'crossref_title{i}', 'Crossref title search', query, 'https://api.crossref.org/works?query.title='+q(query)+'&rows=10', 'crossref')
query = '"weighted sum" "eigenvalues" graph'
add('scholar', 'Google Scholar', query, 'https://scholar.google.com/scholar?q='+q(query), 'html')
query = 'ti:"Quadratic inequalities between the largest eigenvalues"'
add('zbmath', 'zbMATH Open', query, 'https://zbmath.org/?q='+q(query), 'html')
query = 'Quadratic inequalities between the largest eigenvalues'
add('mathscinet', 'MathSciNet public search', query, 'https://mathscinet.ams.org/mathscinet/search/publications.html?pg1=TI&s1='+q(query), 'html')
query = '图 邻接特征值 加权和'
add('cnki', 'CNKI public search', query, 'https://kns.cnki.net/kns8s/defaultresult/index?kw='+q(query), 'html')
add('wanfang', 'Wanfang public search', query, 'https://s.wanfangdata.com.cn/paper?q='+q(query), 'html')
add('semantic_scholar', 'Semantic Scholar', 'weighted sum largest eigenvalues graph', 'https://api.semanticscholar.org/graph/v1/paper/search?query='+q('weighted sum largest eigenvalues graph')+'&limit=20&fields=title,authors,year,url,externalIds', 'json')

def parse(body, kind):
    if kind == 'atom':
        ns={'a':'http://www.w3.org/2005/Atom','x':'http://arxiv.org/schemas/atom','o':'http://a9.com/-/spec/opensearch/1.1/'}
        root=ET.fromstring(body)
        entries=[]
        for e in root.findall('a:entry',ns):
            entries.append({k: ' '.join((e.findtext(path,'',ns)).split()) for k,path in [('id','a:id'),('title','a:title'),('published','a:published'),('updated','a:updated'),('summary','a:summary'),('journal_ref','x:journal_ref'),('doi','x:doi')] } | {'authors':[a.findtext('a:name','',ns) for a in e.findall('a:author',ns)]})
        return {'total':root.findtext('o:totalResults','',ns),'items':entries}
    if kind in ['openalex','crossref','json']:
        obj=json.loads(body)
        if kind=='json': return obj
        records=obj.get('results',[]) if kind=='openalex' else obj['message'].get('items',[obj['message']])
        out=[]
        for r in records:
            title=r.get('title')
            if isinstance(title,list): title='; '.join(title)
            out.append(dict(id=r.get('id'), title=title, doi=r.get('doi',r.get('DOI')), year=r.get('publication_year',r.get('published')), authors=[a.get('author',{}).get('display_name') for a in r.get('authorships',[])] or [(a.get('given','')+' '+a.get('family','')).strip() for a in r.get('author',[])], cited_by_count=r.get('cited_by_count'), source=r.get('container-title'), type=r.get('type'), locations=r.get('locations'), abstract=r.get('abstract')))
        return {'total':obj.get('meta',{}).get('count',obj.get('message',{}).get('total-results')), 'items':out}
    txt=body.decode('utf-8','replace')
    title=re.search(r'<title[^>]*>(.*?)</title>',txt,re.S|re.I)
    # HTML retrieval is only an access probe until search results are inspected.
    clean=re.sub('<[^>]+>',' ',re.sub(r'<(script|style)\b[^>]*>.*?</\1>','',txt,flags=re.S|re.I))
    return {'title':title.group(1) if title else None, 'text_excerpt':' '.join(clean.split())[:3000], 'interpretation':'Access probe; do not treat HTTP success as completed database search.'}

def fetch(job):
    record=dict(job, retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    try:
        req=urllib.request.Request(job['url'], headers={'User-Agent':'AcademicLiteratureAudit/1.0'})
        with urllib.request.urlopen(req, timeout=18) as r:
            body=r.read(); record.update(http_status=r.status, final_url=r.url, content_type=r.headers.get('Content-Type'))
        raw=RAW/(job['key']+'.raw')
        raw.write_bytes(body)
        record.update(raw_file=str(raw.relative_to(BASE)), bytes=len(body), sha256=hashlib.sha256(body).hexdigest(), result=parse(body,job['kind']))
    except Exception as e:
        record.update(error=repr(e))
        if isinstance(e,urllib.error.HTTPError): record['http_status']=e.code
    (RAW/(job['key']+'.json')).write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
    return record

if __name__=='__main__':
    results=[]
    with ThreadPoolExecutor(max_workers=4) as ex:
        for f in as_completed([ex.submit(fetch,j) for j in jobs]):
            r=f.result();results.append(r)
            items=r.get('result',{}).get('items',[])
            print(json.dumps({'key':r['key'],'status':r.get('http_status'), 'error':r.get('error'), 'total':r.get('result',{}).get('total'), 'top_titles':[x.get('title') for x in items[:8]], 'html':r.get('result',{}).get('text_excerpt','')[:400]},ensure_ascii=False),flush=True)
    (BASE/'search_log.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
