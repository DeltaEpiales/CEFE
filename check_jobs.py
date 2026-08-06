import urllib.request
import json
req1 = urllib.request.Request('https://api.github.com/repos/DeltaEpiales/CEFE/actions/runs?per_page=1', headers={'User-Agent': 'Mozilla/5.0'})
data1 = json.loads(urllib.request.urlopen(req1).read())
latest_run_id = data1['workflow_runs'][0]['id']
req = urllib.request.Request(f'https://api.github.com/repos/DeltaEpiales/CEFE/actions/runs/{latest_run_id}/jobs', headers={'User-Agent': 'Mozilla/5.0'})
try:
    data = json.loads(urllib.request.urlopen(req).read())
    for j in data['jobs']:
        print(f"{j['name']}: status={j['status']}, conclusion={j['conclusion']}")
        if j['conclusion'] == 'failure':
            print(f"FAILED JOB: {j['name']}")
            print(f"Log URL: {j['url']}")
except Exception as e:
    print(e)
