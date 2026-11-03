from search_audit import fetch, q, BASE
from concurrent.futures import ThreadPoolExecutor, as_completed
import json

jobs=[]
for key, wid in [('paulin','W7214931116'),('huang_wei','W7171495516'),('kumar','W7148251236')]:
    jobs.append(dict(key='citing_'+key, source='OpenAlex forward citations', query='cites:'+wid, url='https://api.openalex.org/works?filter=cites:'+wid+'&per-page=100', kind='openalex'))
query='all:"convex combination" AND all:"eigenvalues" AND all:"graph"'
jobs.append(dict(key='arxiv_convex',source='arXiv',query=query,url='https://export.arxiv.org/api/query?search_query='+q(query)+'&max_results=100&sortBy=submittedDate&sortOrder=descending',kind='atom'))
query='all:"stability" AND all:"spectral sum"'
jobs.append(dict(key='arxiv_stability',source='arXiv',query=query,url='https://export.arxiv.org/api/query?search_query='+q(query)+'&max_results=100',kind='atom'))
query='"weighted sum" "adjacency eigenvalues"'
jobs.append(dict(key='bing_semantic',source='Bing supplementary web search',query=query,url='https://www.bing.com/search?q='+q(query),kind='html'))
records=[]
with ThreadPoolExecutor(max_workers=3) as ex:
    for f in as_completed([ex.submit(fetch,j) for j in jobs]):
        r=f.result();records.append(r)
        print(json.dumps({'key':r['key'],'status':r.get('http_status'),'error':r.get('error'),'total':r.get('result',{}).get('total'),'titles':[x.get('title') for x in r.get('result',{}).get('items',[])], 'excerpt':r.get('result',{}).get('text_excerpt','')[:1800]},ensure_ascii=False),flush=True)
(BASE/'forward_log.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
