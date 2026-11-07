"""Recompute every displayed/missing regional window from the included CpG counts.

Run from any directory. This reads repository data and prints a verification
receipt; it never modifies inputs or inserts drawing zeros into the estimator.
"""
from pathlib import Path
from collections import defaultdict,Counter
from bisect import bisect_left
import csv,gzip,json,math,sys

R=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(R/'code/figures'))
from panel_helpers import summarize
D=R/'data/regional_methylation'
with (D/'windows.tsv').open() as f:
    expected={(x['array_id'],x['run'],int(x['bin_index'])):x for x in csv.DictReader(f,delimiter='\t')}
arrays=json.loads((D/'arrays.json').read_text())
observed=Counter();checked=0;sites=0
for path in sorted((D/'per_run').glob('ERR*/CpGs.tsv.gz')):
    run=path.parent.name;byref=defaultdict(list)
    with gzip.open(path,'rt') as f:
        for x in csv.DictReader(f,delimiter='\t'):
            assert int(x['valid_calls'])==sum(int(x[k]) for k in ['C','m','h'])
            byref[x['reference']].append(x);sites+=1
    positions={ref:[int(x['reference_C_position0']) for x in rr] for ref,rr in byref.items()}
    for r in arrays:
        a,b=int(r['plot_start0']),int(r['plot_end0']);ref=r['seqid']
        pp=positions.get(ref,[]);rows=byref.get(ref,[])[bisect_left(pp,a):bisect_left(pp,b)]
        for x in summarize(rows,a,b,10000):
            e=expected[r['array_id'],run,x['bin_index']]
            for k in ['start0','end0','observed_CpGs','valid_observations','five_mC_observations']:
                assert int(e[k])==x[k],(r['array_id'],run,k,e[k],x[k])
            assert x['status']==e['status']
            for k in ['mean_CpG_5mC_pct','call_weighted_5mC_pct']:
                assert (x[k] is None and e[k]=='') or (x[k] is not None and math.isclose(x[k],float(e[k]),abs_tol=1e-9))
            assert float(e['plot_5mC_pct'])==(x['mean_CpG_5mC_pct'] if x['status']=='display' else 0.)
            assert (e['zero_placeholder']=='True')==(x['status']!='display')
            observed[x['status']]+=1;checked+=1
assert checked==len(expected)==52074 and sites==2774105
print(json.dumps({'passed':True,'sample_reference_CpG_rows':sites,'array_sample_windows':checked,'statuses':dict(observed),'coordinates_and_both_weighting_schemes_verified':True},indent=2))
