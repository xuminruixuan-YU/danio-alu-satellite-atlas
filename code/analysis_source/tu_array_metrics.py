"""Recompute copy divergence and localized methylation using large TU regions.

Reuse eligible raw copies and primary CpG observations, never old locus means.
"""
import csv,gzip,hashlib,itertools,json,math,sys,collections
from pathlib import Path
from statistics import mean,median
from datetime import datetime
from zoneinfo import ZoneInfo
W=Path(__file__).resolve().parents[1];ROOT=W.parents[1]
R=ROOT/'output/danio_five_part_followup_2026-11-04';E=ROOT/'output/danio_cl1_cl2_execution_2026-11-03'
sys.path.insert(0,str(R/'07_edit_distance_review_20261105/code'))
from metrics import Totals,one,FIELDS
O=W/'array_metrics';O.mkdir(exist_ok=True)
def rows(p):
    op=gzip.open if str(p).endswith('.gz') else open
    with op(p,'rt') as f:return list(csv.DictReader(f,delimiter='\t'))
def table(p,rr,fields=None):
    op=gzip.open if str(p).endswith('.gz') else open
    with op(p,'wt') as f:
        z=csv.DictWriter(f,fields or list(rr[0]),delimiter='\t');z.writeheader();z.writerows(rr)
def hashkey(u):return hashlib.sha256(('copy20261103|'+u['copy_coordinates']).encode()).digest()[:8]
def totals(pairs):
    t=Totals()
    for a,b in pairs:t.add(one(a['sequence'],b['sequence']))
    return t.row()

def main():
    sources={};units=collections.defaultdict(list);allcounts=collections.Counter()
    for p in sorted((E/'genomes/TU/catalog_units').glob('*/units.tsv.gz')):
        with p.open('rb') as f:sources[str(p)]=hashlib.file_digest(f,'sha256').hexdigest()
        grouped=collections.defaultdict(list)
        for x in rows(p):grouped[x['fragment']].append(x)
        for fid,uu in grouped.items():
            uu.sort(key=lambda x:(int(x['start0']),int(x['end0'])));end=-1
            for i,x in enumerate(uu):
                a,b=int(x['start0']),int(x['end0']);overlap=a<end or (i+1<len(uu) and int(uu[i+1]['start0'])<b);end=max(end,b)
                if overlap or x['pass_status']!='normalized':continue
                seq=x['oriented_circular_sequence'];assert seq and set(seq)<=set('ACGT')
                units[fid].append(dict(fragment=fid,seqid=x['seqid'],start0=a,end0=b,family=x['family'],copy_coordinates=f"{x['seqid']}:{a}-{b}",sequence=seq,subtype=x['source_subtype']))
                allcounts[x['family']]+=1
    assert dict(allcounts)=={'CL1':59774,'CL2':301082},dict(allcounts)
    variants=[(10000,50000),(25000,25000),(25000,50000),(25000,100000),(25000,200000),(50000,50000)]
    summaries=[];allwithin=[];sampled=[];between=[];primary_groups={};primary_arrays={};primary_members={}
    for gap,cutoff in variants:
        variant=f'gap{gap}_size{cutoff}';d=W/'array_catalog/variants'/f'TU_gap{gap}'
        regions={x['array_id']:x for x in rows(d/'regions.tsv') if int(x['supported_bp'])>=cutoff}
        membership=collections.defaultdict(list)
        for x in rows(d/'membership.tsv'):
            if x['array_id'] in regions:membership[x['array_id']].append(x['source_id'])
        for family in ['CL1','CL2']:
            cats=[family+'_all']+(['CL2_A','CL2_B'] if family=='CL2' else [])
            for cat in cats:
                wr=[];samples={}
                for aid,ar in sorted(regions.items()):
                    if ar['family']!=family:continue
                    pool=[u for fid in membership[aid] for u in units[fid] if cat.endswith('_all') or u['subtype']==cat[-1]]
                    pool.sort(key=hashkey);sample=pool[:32];samples[aid]=sample
                    metrics=totals(itertools.combinations(sample,2))
                    rec=dict(variant=variant,array_id=aid,category=cat,family=family,seqid=ar['seqid'],chromosome=ar['chromosome'],start0=ar['start0'],end0=ar['end0'],span_bp=ar['span_bp'],supported_bp=ar['supported_bp'],eligible_copies=len(pool),sampled_copies=len(sample),**metrics)
                    wr.append(rec);allwithin.append(rec)
                    sampled.extend(dict(variant=variant,array_id=aid,category=cat,rank=i+1,**u) for i,u in enumerate(sample))
                    if (gap,cutoff)==(25000,50000) and cat.endswith('_all'):
                        primary_groups[aid]=pool;primary_arrays[aid]=ar
                        primary_members.update({fid:aid for fid in membership[aid]})
                valid=[x for x in wr if x['normalized_edit_distance'] is not None]
                summaries.append(dict(variant=variant,family=family,category=cat,arrays=len(wr),arrays_with_comparable_pairs=len(valid),eligible_copies=sum(x['eligible_copies'] for x in valid),sampled_copies=sum(x['sampled_copies'] for x in valid),array_equal_edit_distance=mean(x['normalized_edit_distance'] for x in valid) if valid else None,array_equal_event_distance=mean(x['event_sensitivity'] for x in valid) if valid else None,copy_weighted_array_edit_distance=sum(x['eligible_copies']*x['normalized_edit_distance'] for x in valid)/sum(x['eligible_copies'] for x in valid) if valid else None))
                if (gap,cutoff)==(25000,50000) and cat.endswith('_all'):
                    for a,b in itertools.combinations(sorted(samples),2):
                        if min(len(samples[a]),len(samples[b]))<4:continue
                        ar,br=regions[a],regions[b]
                        between.append(dict(family=family,array1=a,array2=b,chromosome1=ar['chromosome'],chromosome2=br['chromosome'],stratum='same_chromosome' if ar['seqid']==br['seqid'] else 'different_chromosome',**totals(itertools.product(samples[a][:4],samples[b][:4]))))
        print(variant,'done',flush=True)
    table(O/'within_array.tsv',allwithin);table(O/'sampled_units.tsv.gz',sampled);table(O/'family_summary.tsv',summaries);table(O/'between_arrays.tsv',between)
    strata=[]
    for family in ['CL1','CL2']:
        for scope in ['same_chromosome','different_chromosome','all_between']:
            edges=[x for x in between if x['family']==family and (scope=='all_between' or x['stratum']==scope)]
            neighbors=collections.defaultdict(list)
            for x in edges:neighbors[x['array1']].append(x['normalized_edit_distance']);neighbors[x['array2']].append(x['normalized_edit_distance'])
            within=[totals(itertools.combinations(primary_groups[k][:4],2))['normalized_edit_distance'] for k in neighbors]
            strata.append(dict(family=family,stratum=scope,arrays=len(neighbors),array_pairs=len(edges),within_same_arrays_mean=mean(within) if within else None,between_equal_target_array_mean=mean(mean(v) for v in neighbors.values()) if neighbors else None,scope='shared copies and arrays; descriptive values, not independent pairs'))
    table(O/'within_between_summary.tsv',strata)
    # Relabel actual CpG rows through source fragments, then recalculate weights.
    cps=rows(R/'11_locus_methylation/cohort_review/20261106T0828/individual_CpG_primary.tsv')
    cp_groups=collections.defaultdict(list);cpout=[];runs=sorted({x['run'] for x in cps})
    for x in cps:
        aid=primary_members.get(x['fragment'])
        if aid and int(x['valid_calls'])>0:
            cp_groups[x['run'],aid].append(x);cpout.append(dict(array_id=aid,**x))
    methyl=[]
    for run in runs:
        for aid,ar in primary_arrays.items():
            pp=cp_groups[run,aid];calls=sum(int(x['valid_calls']) for x in pp)
            methyl.append(dict(run=run,array_id=aid,family=ar['family'],chromosome=ar['chromosome'],callable_CpGs=len(pp),valid_observations=calls,median_CpG_depth=median(int(x['valid_calls']) for x in pp) if pp else None,
                CpG_equal_5mC=mean(int(x['m'])/int(x['valid_calls']) for x in pp) if pp else None,observation_weighted_5mC=sum(int(x['m']) for x in pp)/calls if calls else None,scope='position-validated satellite read subset; unobserved sequence not imputed'))
    table(O/'array_methylation.tsv',methyl);table(O/'CpG_methylation.tsv.gz',cpout)
    family_meth=[]
    for run in runs:
        for fam in ['CL1','CL2']:
            vv=[x for x in methyl if x['run']==run and x['family']==fam and x['CpG_equal_5mC'] is not None]
            family_meth.append(dict(run=run,family=fam,arrays_with_CpGs=len(vv),callable_CpGs=sum(x['callable_CpGs'] for x in vv),array_equal_5mC=mean(x['CpG_equal_5mC'] for x in vv) if vv else None))
    table(O/'array_methylation_family_summary.tsv',family_meth)
    # Independently recompute a bounded subset of reported edit distances with DP.
    def dp(a,b):
        last=list(range(len(b)+1))
        for i,u in enumerate(a,1):
            cur=[i]
            for j,v in enumerate(b,1):cur.append(min(cur[-1]+1,last[j]+1,last[j-1]+(u!=v)))
            last=cur
        return last[-1]
    audit=[]
    for aid,pool in primary_groups.items():
        if len(pool)<2:continue
        a,b=pool[0]['sequence'],pool[min(len(pool),32)-1]['sequence'];distance=dp(a,b)
        assert one(a,b)[0]==distance
        audit.append(dict(array_id=aid,sequence1=a,sequence2=b,independent_ED=distance,denominator=max(len(a),len(b))))
    table(O/'independent_ED_checks.tsv',audit)
    (O/'computed.json').write_text(json.dumps(dict(at_hkt=datetime.now(ZoneInfo('Asia/Hong_Kong')).isoformat(),eligible_source_counts=dict(allcounts),unit_source_sha256=sources,primary_arrays=len(primary_arrays),independent_distance_checks=len(audit),primary_variant='gap25000_size50000',new_array_catalogue_final=False,old_locus_means_reused=False),indent=2)+'\n')

if __name__=='__main__':main()
