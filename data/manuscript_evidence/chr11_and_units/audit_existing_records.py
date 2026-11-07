from pathlib import Path
import csv,gzip,json,collections,hashlib,statistics,shutil
D=Path(__file__).resolve().parents[1];ROOT=D.parents[1];W=ROOT/'output/danio_array_revision_2026-11-06';E=ROOT/'output/danio_cl1_cl2_execution_2026-11-03';O=D/'evidence/round2';O.mkdir(exist_ok=True)
def rows(p):
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rr):
 with (gzip.open(p,'wt') if p.suffix=='.gz' else p.open('w')) as f:
  w=csv.DictWriter(f,rr[0].keys(),delimiter='\t');w.writeheader();w.writerows(rr)
aid='TU.NC_133186.1.CL1.R002';members=[x for x in rows(W/'array_catalog/membership.tsv') if x['array_id']==aid];fids={x['source_id'] for x in members}
p=E/'genomes/TU/catalog_units/NC_133186.1/units.tsv.gz';group=collections.defaultdict(list)
for n,r in enumerate(rows(p),2):
 if r['fragment'] in fids:group[r['fragment']].append(dict(source_line=n,**r))
records=[]
for fid,rr in group.items():
 rr.sort(key=lambda r:(int(r['start0']),int(r['end0'])));end=-1
 for i,r in enumerate(rr):
  a,b=int(r['start0']),int(r['end0']);overlap=a<end or (i+1<len(rr) and int(rr[i+1]['start0'])<b);end=max(end,b)
  records.append(dict(array_id=aid,**r,overlap_excluded=overlap,final_status='overlap' if overlap else r['pass_status']))
assert len(records)==8077
counts=collections.Counter((r['pass_status'],r['normalization'] or '(not_attempted)',r['final_status']) for r in records)
summ=[dict(array_id=aid,original_status=a,normalization_substatus=b,final_status=c,units=n) for (a,b,c),n in sorted(counts.items())]
write(O/'chr11_units_with_substates.tsv.gz',records);write(O/'chr11_substate_counts.tsv',summ)
kept=[r for r in records if r['final_status']=='source_family_only'];assert len(kept)==7698 and {r['normalization'] for r in kept}=={'cycle_offset_or_span_tie'}
alltie=[r for r in records if r['normalization']=='cycle_offset_or_span_tie'];assert len(alltie)==7725
rep=rows(D/'evidence/TU_copy_representation.tsv');coverage=[]
for fam in ['CL1','CL2']:
 vals=[float(r['eligible_fraction_of_satellite'])*100 for r in rep if r['family']==fam]
 coverage.append(dict(family=fam,regions=len(vals),min_percent=min(vals),median_percent=statistics.median(vals),max_percent=max(vals),regions_below_1pct=sum(v<1 for v in vals),regions_below_10pct=sum(v<10 for v in vals)))
write(O/'eligible_coverage_distribution.tsv',coverage)
# Read the existing plotted-window table; no new methylation calculation.
wp=W/'publication_merged_2026-11-06/checks/methylation_windows.tsv';wins=rows(wp)
assert all((r['status']=='display')==(int(r['observed_CpGs'])>=10) for r in wins)
assert all(float(r['plot_5mC_pct'])==0 for r in wins if r['status']!='display')
from Bio import SeqIO
seq={r.id:str(r.seq) for r in SeqIO.parse(E/'inputs/three.fa','fasta')};A=seq['CL2_half1'];B=seq['CL2_half2'];diag=[dict(position1=i+1,A=a,B=b) for i,(a,b) in enumerate(zip(A,B)) if a!=b]
assert diag==[dict(position1=6,A='G',B='A'),dict(position1=49,A='T',B='C')]
source_paths=[p,E/'code/catalog_units.py',E/'code/phase.py',E/'inputs/three.fa',D/'code/panel_helpers.py',wp]
shutil.copy2(E/'code/catalog_units.py',O/'catalog_units.py');shutil.copy2(E/'code/phase.py',O/'phase.py');shutil.copy2(E/'inputs/three.fa',O/'three.fa')
res={'chr11_array':aid,'source_rows':8077,'raw_tie_rows':7725,'tie_rows_reclassified_as_overlap':27,'prior_chr11_wide_tie_count':7726,'prior_count_scope_difference':'One tie row belongs to a shorter chr11 AT region outside the large array; not part of these 8077 rows.','final_source_family_only':7698,'final_eligible':15,'normalization_code':'catalog_units.py: normalize compares both directions with Edlib HW against three tandem templates; stores (start modulo template length, aligned target span); >1 best pair yields cycle_offset_or_span_tie.','confirmed':'All 7698 non-overlap family-only rows carry this narrower stored substatus.','unresolved':'Combined status does not separate start-only, span-only or both; it does not establish a biological cause. No alignment or threshold/phase rerun performed.','coverage':coverage,'AB_diagnostic_positions':diag,'AB_scope':'Two diagnostic-base classes, not two unique full-length sequences; phase.rawcall uses 82..102 bp, ACGT, both global alignments concordant.','track_unit':'Distinct reference C positions within each individual and 10 kb bin, each with at least one callable observation; count is not read observations. panel_helpers.summarize asserts unique positions.','existing_window_records_checked':len(wins),'window_states':dict(collections.Counter(r['status'] for r in wins)),'sources':[dict(path=str(q),sha256=hashlib.sha256(q.read_bytes()).hexdigest()) for q in source_paths]}
(O/'unit_audit.json').write_text(json.dumps(res,ensure_ascii=False,indent=2));print(json.dumps({k:v for k,v in res.items() if k!='sources'},ensure_ascii=False,indent=2))
