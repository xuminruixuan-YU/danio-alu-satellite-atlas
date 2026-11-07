"""Frozen-region sequence-only sensitivity. No methylation inputs or job entrypoints."""
import csv,gzip,json,hashlib,time,os,sys,itertools,collections,statistics,resource
from pathlib import Path
from functools import lru_cache
from concurrent.futures import ProcessPoolExecutor,as_completed
import edlib
from Bio import SeqIO
D=Path(__file__).resolve().parents[1]
for sub in ['results','checks']: (D/sub).mkdir(exist_ok=True)
W=D/'inputs/frozen_catalog'
UNIT_ROOT=D/'inputs/frozen_units'
T={r.id:str(r.seq) for name in ['two','three'] for r in SeqIO.parse(D/'inputs'/f'{name}.fa','fasta')}
F={'CL1':T['CL1_186nt'],'CL2':T['CL2_half1']}
CONDS=[(.70,.10),(.75,.10),(.80,.10),(.90,.10),(.75,.05),(.75,.15)]
PHASES=['strict','equivalent','min_option','max_option']
BASE='score75_margin10_strict'
def sid(c,p):return f'score{round(c[0]*100)}_margin{round(c[1]*100)}_{p}'
def read(p):
 with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rr):
 if not rr:return
 with (gzip.open(p,'wt') if p.suffix=='.gz' else p.open('w')) as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def reciprocal(x,y):
 return max(1-max(edlib.align(x,t*(2+len(x)//len(t)),mode='HW')['editDistance']/len(x),edlib.align(t,x*(2+len(t)//len(x)),mode='HW')['editDistance']/len(t)) for t in {y,rc(y)})
@lru_cache(maxsize=60000)
def scores(s):
 sims={}
 for k,t in T.items():
  ratio=max(len(s),len(t))/min(len(s),len(t))
  sims[k]=reciprocal(s,t) if ratio<=5.2 and abs(ratio-round(ratio))/ratio<=.15 else 0.
 return sims['CL1_186nt'],max(v for k,v in sims.items() if k.startswith('CL2'))
def rotate(q,t,start,span):
 if start==0:return q
 ref=(t*3)[start:start+span];a=edlib.align(q,ref,mode='NW',task='path');z=edlib.getNiceAlignment(a,q,ref)
 ti=start;qi=0
 for x,y in zip(z['query_aligned'],z['target_aligned']):
  if ti==len(t) and y!='-':return q[qi:]+q[:qi]
  ti+=y!='-';qi+=x!='-'
 return ''
@lru_cache(maxsize=60000)
def norm(s,f):
 t=F[f];opts=[]
 for strand,q in [('+',s),('-',rc(s))]:
  a=edlib.align(q,t*3,mode='HW',task='locations')
  loc=sorted({(l%len(t),h-l+1) for l,h in a['locations']})
  opts.append((a['editDistance'],strand,q,loc))
 best=[o for o in opts if o[0]==min(v[0] for v in opts)]
 if len(best)>1:return ('direction_tie',0,0,0,'','','','')
 _,strand,q,loc=best[0];seq=[rotate(q,t,a,b) for a,b in loc]
 starts=len({a for a,b in loc});spans=len({b for a,b in loc})
 status='unique' if len(loc)==1 else 'offset_and_span_tie' if starts>1 and spans>1 else 'offset_tie' if starts>1 else 'span_tie'
 if not all(seq):return ('origin_not_covered',len(loc),starts,spans,'','','',strand)
 unique=len(set(seq))
 return (status,len(loc),starts,spans,seq[0] if unique==1 else '',seq[0],seq[-1],strand)
def key(x):return hashlib.sha256(('copy20261103|'+x['copy_coordinates']).encode()).digest()[:8]
def metric(rr,seqkey):
 rr=sorted(rr,key=key)[:32];num=den=0
 for a,b in itertools.combinations(rr,2):
  x,y=a[seqkey],b[seqkey];num+=edlib.align(x,y,mode='NW')['editDistance'];den+=max(len(x),len(y))
 return len(rr),num/den if den else '',num,den
REG={r['array_id']:r for r in read(W/'array_catalog/regions.tsv') if r['species']=='TU' and r['qualifies_large']=='True'}
MEM={r['source_id']:r['array_id'] for r in read(W/'array_catalog/membership.tsv') if r['array_id'] in REG}
def process(ch):
 start=time.monotonic();out=D/'results'/ch;out.mkdir(exist_ok=True)
 if (out/'done.json').exists():return json.loads((out/'done.json').read_text())
 source=UNIT_ROOT/ch/'units.tsv.gz';groups=collections.defaultdict(list)
 for i,r in enumerate(read(source),2):
  if r['fragment'] in MEM:groups[r['fragment']].append(dict(source_line=i,**r))
 records=[];pools=collections.defaultdict(lambda:collections.defaultdict(list));fails=collections.Counter()
 for fid,rr in groups.items():
  rr.sort(key=lambda r:(int(r['start0']),int(r['end0'])));end=-1
  for i,r in enumerate(rr):
   a,b=int(r['start0']),int(r['end0']);ov=a<end or (i+1<len(rr) and int(rr[i+1]['start0'])<b);end=max(end,b)
   aid=MEM[fid];fam=r['family'];s=r['raw_sequence'];reason=''
   if ov:reason='overlap'
   elif r['within_fragment']!='True':reason='truncated'
   elif set(s)-set('ACGT'):reason='ambiguous_bases'
   elif not .72<=len(s)/len(F[fam])<=1.28:reason='length'
   v=('',0,0,0,'','','','');sc1=sc2=-1
   if not reason:
    sc1,sc2=scores(s);v=norm(s,fam)
   ns,nloc,nstart,nspan,eq,lo,hi,strand=v
   seqs={'strict':eq if ns=='unique' else '', 'equivalent':eq,'min_option':lo,'max_option':hi}
   # Normalization is audited independently of family assignment, never conflated.
   if not reason and r['pass_status']=='normalized':assert ns=='unique' and eq==r['oriented_circular_sequence'],(ch,a,ns)
   if not reason:
    old=('CL1' if sc1>sc2 else 'CL2') if max(sc1,sc2)>=.75 and abs(sc1-sc2)>=.10 else 'unresolved'
    assert old==r['source_family'],(a,old,r['source_family'])
   row=dict(array_id=aid,fragment=fid,seqid=ch,start0=a,end0=b,copy_coordinates=f'{ch}:{a}-{b}',family=fam,source_line=r['source_line'],old_status=r['pass_status'],old_normalization=r['normalization'],quality_exclusion=reason,CL1_score=sc1,CL2_score=sc2,phase_state=ns,optimal_options=nloc,distinct_offsets=nstart,distinct_spans=nspan,selected_direction=strand,sd_direction=r['sd_strand'],direction_disagrees_SD=(bool(strand) and strand!=r['sd_strand']),strict=seqs['strict'],equivalent=eq,min_option=lo,max_option=hi)
   records.append(row)
   own,other=(sc1,sc2) if fam=='CL1' else (sc2,sc1)
   for c in CONDS:
    fr=reason or ('family_conflict' if own<=other else 'family_score' if own<c[0] else 'family_margin' if own-other<c[1] else '')
    for ph in PHASES:
     name=sid(c,ph);why=fr or (ns if not seqs[ph] else 'retained')
     fails[aid,name,why]+=1
     if why=='retained':pools[aid][name].append(row)
 regional=[];common=[];sampling=[]
 for aid in sorted({MEM[f] for f in groups}):
  ar=REG[aid];base=pools[aid][BASE];basekeys={r['copy_coordinates'] for r in base}
  for c in CONDS:
   for ph in PHASES:
    name=sid(c,ph);rr=pools[aid][name];n,dist,num,den=metric(rr,ph);bp=sum(r['end0']-r['start0'] for r in rr)
    regional.append(dict(array_id=aid,family=ar['family'],chromosome=ar['chromosome'],scenario=name,score=c[0],margin=c[1],phase=ph,initial_units=sum(1 for r in records if r['array_id']==aid),eligible_units=len(rr),retained_bp=bp,satellite_bp=ar['supported_bp'],retained_fraction=bp/int(ar['supported_bp']),sampled_units=n,normalized_edit_distance=dist,edit_bases=num,denominator_bp=den))
    both=[r for r in rr if r['copy_coordinates'] in basekeys];n,di,nu,de=metric(both,ph);_,db,_,_=metric(both,'strict')
    common.append(dict(array_id=aid,family=ar['family'],scenario=name,common_copies=len(both),sampled_common_copies=n,baseline_distance=db,alternative_distance=di,common_bp=sum(r['end0']-r['start0'] for r in both)))
    for rank,r in enumerate(sorted(rr,key=key)[:32],1):sampling.append(dict(array_id=aid,scenario=name,rank=rank,copy_coordinates=r['copy_coordinates'],sequence=r[ph]))
 write(out/'unit_audit.tsv.gz',records);write(out/'regions.tsv',regional);write(out/'common_copies.tsv',common);write(out/'sampled_sequences.tsv.gz',sampling)
 write(out/'exclusions.tsv',[dict(array_id=a,scenario=b,reason=c,units=n) for (a,b,c),n in sorted(fails.items())])
 info=dict(chromosome=ch,records=len(records),seconds=time.monotonic()-start,peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
 (out/'done.json').write_text(json.dumps(info,indent=2));return info

def summarize():
 regional=[];commons=[];ex=[]
 for p in sorted((D/'results').glob('*/done.json')):
  regional+=read(p.parent/'regions.tsv');commons+=read(p.parent/'common_copies.tsv');ex+=read(p.parent/'exclusions.tsv')
 for r in regional:r['pass_fraction']=int(r['eligible_units'])/int(r['initial_units']) if int(r['initial_units']) else 0.
 write(D/'results/all_regions.tsv',regional);write(D/'results/common_copies.tsv',commons);write(D/'results/exclusions.tsv',ex)
 summ=[]
 for c in CONDS:
  for ph in PHASES:
   for fam in F:
    name=sid(c,ph);rr=[r for r in regional if r['family']==fam and r['scenario']==name];valid=[r for r in rr if r['normalized_edit_distance']!='']
    cc=[r for r in commons if r['family']==fam and r['scenario']==name and r['baseline_distance']!='']
    summ.append(dict(scenario=name,family=fam,regions=len(rr),comparable_regions=len(valid),eligible_units=sum(int(r['eligible_units']) for r in rr),initial_units=sum(int(r['initial_units']) for r in rr),retained_bp=sum(int(r['retained_bp']) for r in rr),satellite_bp=sum(int(r['satellite_bp']) for r in rr),median_retained_fraction=statistics.median(float(r['retained_fraction']) for r in rr),region_equal_distance=statistics.mean(float(r['normalized_edit_distance']) for r in valid) if valid else '',common_copy_regions=len(cc),common_baseline_distance=statistics.mean(float(r['baseline_distance']) for r in cc) if cc else '',common_alternative_distance=statistics.mean(float(r['alternative_distance']) for r in cc) if cc else ''))
 for r in summ:r['pass_fraction']=r['eligible_units']/r['initial_units'] if r['initial_units'] else 0.
 write(D/'results/family_summary.tsv',summ)
 baseline=[r for r in regional if r['scenario']==BASE]
 orig={r['array_id']:r for r in read(W/'array_metrics/within_array.tsv') if r['variant']=='gap25000_size50000' and r['category'].endswith('_all')}
 for r in baseline:
  x=orig[r['array_id']];assert int(r['eligible_units'])==int(x['eligible_copies']);assert abs(float(r['normalized_edit_distance'])-float(x['normalized_edit_distance']))<1e-14
 (D/'checks/baseline_reproduction.json').write_text(json.dumps(dict(passed=True,regions=len(baseline),scope='Exact original eligible counts and region estimates reproduced; existing files not changed'),indent=2))
 return summ
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--pilot',action='store_true');p.add_argument('--workers',type=int,default=12);args=p.parse_args()
 chs=sorted({r['seqid'] for r in REG.values()})
 if args.pilot:chs=[f'NC_{n}.1' for n in [133177,133186,133196]]
 t=time.monotonic()
 with ProcessPoolExecutor(args.workers) as pool:
  fs=[pool.submit(process,ch) for ch in chs]
  for f in as_completed(fs):print(json.dumps(f.result()),flush=True)
 summary=summarize();print(json.dumps({'wall_seconds':time.monotonic()-t,'processed_chromosomes':len(chs)}),flush=True)
 (D/'checks'/('pilot.json' if args.pilot else 'complete.json')).write_text(json.dumps(dict(wall_seconds=time.monotonic()-t,processed_chromosomes=len(chs),code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2))
