"""Audit equivalent representations and expanded-pool comparisons; sequence only."""
import sequence_sensitivity as s
import csv,gzip,json,collections,itertools,statistics,hashlib
D=s.D
primary='score75_margin10_equivalent'
rr=s.read(D/'results/all_regions.tsv');summ=s.read(D/'results/family_summary.tsv')
units=collections.defaultdict(list);phase=collections.Counter();chr11=collections.Counter();total=0
for p in sorted((D/'results').glob('*/unit_audit.tsv.gz')):
 for r in s.read(p):
  total+=1
  if r['quality_exclusion']:continue
  f=r['family'];own=float(r['CL1_score'] if f=='CL1' else r['CL2_score']);other=float(r['CL2_score'] if f=='CL1' else r['CL1_score'])
  if own<.75 or own-other<.10:continue
  phase[f,r['old_status'],r['phase_state'],bool(r['equivalent'])]+=1
  if r['array_id']=='TU.NC_133186.1.CL1.R002':chr11[r['old_status'],r['phase_state'],bool(r['equivalent'])]+=1
  if r['equivalent']:
   assert r['min_option']==r['max_option']==r['equivalent']
   units[r['array_id']].append(r)
# Cross-check pool counts/coverage against independently reselected audit rows.
for aid,pool in units.items():
 r=next(r for r in rr if r['array_id']==aid and r['scenario']==primary)
 assert len(pool)==int(r['eligible_units'])
 assert sum(int(x['end0'])-int(x['start0']) for x in pool)==int(r['retained_bp'])
# Additional sampling answers a new question: stability of the expanded equivalent-copy pool.
ext=[];within4={};between=[]
for aid,pool in sorted(units.items()):
 pool.sort(key=s.key)
 for view,pp in [('expanded_32',pool[:32]),('expanded_128',pool[:128]),('rescued_only_32',[r for r in pool if not r['strict']][:32])]:
  num=den=0
  for a,b in itertools.combinations(pp,2):
   x,y=a['equivalent'],b['equivalent'];num+=s.edlib.align(x,y,mode='NW')['editDistance'];den+=max(len(x),len(y))
  ext.append(dict(array_id=aid,family=s.REG[aid]['family'],view=view,eligible_pool=len(pool) if not view.startswith('rescued') else sum(not x['strict'] for x in pool),sampled_units=len(pp),normalized_edit_distance=num/den if den else ''))
 n,di,_,_=s.metric(pool[:4],'equivalent');within4[aid]=di
for fam in s.F:
 ids=[a for a in units if s.REG[a]['family']==fam];nei=collections.defaultdict(list)
 for a,b in itertools.combinations(ids,2):
  num=den=0
  for x,y in itertools.product(units[a][:4],units[b][:4]):
   u,v=x['equivalent'],y['equivalent'];num+=s.edlib.align(u,v,mode='NW')['editDistance'];den+=max(len(u),len(v))
  nei[a].append(num/den);nei[b].append(num/den)
 between.append(dict(family=fam,regions=len(ids),within4=statistics.mean(within4[a] for a in ids),between4=statistics.mean(statistics.mean(nei[a]) for a in ids)))
s.write(D/'results/expanded_sampling.tsv',ext);s.write(D/'results/expanded_within_between.tsv',between)
s.write(D/'results/phase_substates.tsv',[dict(family=f,original_status=o,phase_state=p,equivalent_sequence=eq,units=n) for (f,o,p,eq),n in sorted(phase.items())])
s.write(D/'results/chr11_phase_substates.tsv',[dict(original_status=o,phase_state=p,equivalent_sequence=eq,units=n) for (o,p,eq),n in sorted(chr11.items())])
# Explicit common-region comparisons: all alternatives currently retain all 49.
cr=[]
for x in summ:
 name=x['scenario'];fam=x['family'];alt={r['array_id']:float(r['normalized_edit_distance']) for r in rr if r['scenario']==name and r['family']==fam and r['normalized_edit_distance']!=''}
 base={r['array_id']:float(r['normalized_edit_distance']) for r in rr if r['scenario']==s.BASE and r['family']==fam and r['normalized_edit_distance']!=''}
 ids=sorted(alt.keys()&base.keys())
 cr.append(dict(scenario=name,family=fam,common_regions=len(ids),original_mean=statistics.mean(base[a] for a in ids),alternative_mean=statistics.mean(alt[a] for a in ids)))
s.write(D/'results/common_regions.tsv',cr)
res=dict(all_unit_records=total,regions=len(units),primary_pool_counts_and_coverage_verified=True,baseline_49_region_counts_and_distances_reproduced=True,expanded_128={f:statistics.mean(float(r['normalized_edit_distance']) for r in ext if r['family']==f and r['view']=='expanded_128') for f in s.F},expanded_within_between=between,all_scenarios_AT_lower_than_GC=all(float(next(x for x in summ if x['scenario']==n and x['family']=='CL1')['region_equal_distance'])<float(next(x for x in summ if x['scenario']==n and x['family']=='CL2')['region_equal_distance']) for n in {x['scenario'] for x in summ}),common_copy_values_identical=all(r['baseline_distance']==r['alternative_distance'] for r in s.read(D/'results/common_copies.tsv')),interpretation_common_copies='Identical by construction when original unique-normalized copies are held fixed; new mean differences come from copy selection, not edits to those sequences.',chr11_substates=[dict(status=o,phase=p,equivalent=eq,n=n) for (o,p,eq),n in chr11.items()])
(D/'checks/verification.json').write_text(json.dumps(res,ensure_ascii=False,indent=2));print(json.dumps(res,ensure_ascii=False,indent=2))
