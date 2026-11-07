"""Recompute denominators, spatial summaries and sensitivity under array units."""
import os
import collections,csv,gzip,itertools,json,sys
from pathlib import Path
from statistics import mean,median
W=Path(__file__).resolve().parents[1]
ROOT=W.parents[1]
R=ROOT / os.environ.get('DANIO_ANNOTATION_ROOT', 'inputs/annotation')
E=ROOT / os.environ.get('DANIO_ASSEMBLY_RUN_ROOT', 'inputs/assembly_runs')
sys.path.insert(0,str(E/'code'))
from local_alternation import assess
sys.path.insert(0,str(R/'07_edit_distance_review_20261105/code'))
from metrics import Totals,one
def rows(p):
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:return list(csv.DictReader(f,delimiter='\t'))
def table(p,rr):
    with (gzip.open(p,'wt') if str(p).endswith('.gz') else open(p,'w')) as f:
        w=csv.DictWriter(f,list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)
def main():
    o=W/'array_summaries';o.mkdir(exist_ok=True)
    ar={x['array_id']:x for x in rows(W/'array_catalog/regions.tsv') if x['species']=='TU' and x['qualifies_large']=='True'}
    members={x['source_id']:x['array_id'] for x in rows(W/'array_catalog/membership.tsv') if x['array_id'] in ar}
    den={x['fragment']:int(x['reference_CpGs']) for x in rows(R/'11_locus_methylation/per_run/ERR16782314/reference_CpG_denominators.tsv')}
    arden=collections.Counter()
    for fid,aid in members.items():arden[aid]+=den[fid]
    results=[];coverage=[];allprimary=[];sens_arrays=[];runs=[f'ERR167823{i}' for i in range(14,20)]
    for run in runs:
        groups=collections.defaultdict(list)
        for x in rows(R/'11_locus_methylation/per_run'/run/'CpG_sensitivity.tsv.gz'):
            if x['fragment'] in members:groups[int(x['min_baseQ']),float(x['model_confidence']),members[x['fragment']]].append(x)
        gs=collections.defaultdict(list)
        for (q,t,aid),ss in groups.items():
            good=[x for x in ss if int(x['valid_calls'])>0]
            v=dict(run=run,array_id=aid,family=ar[aid]['family'],min_baseQ=q,model_confidence=t,reference_CpGs=arden[aid],observed_CpGs=len(ss),callable_CpGs=len(good),valid_calls=sum(int(x['valid_calls']) for x in good),CpG_equal_5mC=mean(int(x['m'])/int(x['valid_calls']) for x in good) if good else None,CpG_equal_5hmC=mean(int(x['h'])/int(x['valid_calls']) for x in good) if good else None)
            sens_arrays.append(v);gs[q,t,v['family']].append(v)
            if q==10 and t==.8:allprimary.extend(dict(array_id=aid,**x) for x in ss)
        for (q,t,fam),vv in gs.items():
            good=[x for x in vv if x['callable_CpGs']>0]
            results.append(dict(run=run,family=fam,min_baseQ=q,model_confidence=t,callable_CpGs=sum(x['callable_CpGs'] for x in good),arrays_with_calls=len(good),valid_calls=sum(x['valid_calls'] for x in good),array_equal_5mC=mean(x['CpG_equal_5mC'] for x in good) if good else None,array_equal_5hmC=mean(x['CpG_equal_5hmC'] for x in good) if good else None,encoding_limited_confidence095=t==.95))
        for fam in ['CL1','CL2']:
            ss=[x for x in allprimary if x['run']==run and x['family']==fam];good=[x for x in ss if int(x['valid_calls'])>0];depth=[int(x['valid_calls']) for x in good];reference=sum(arden[k] for k in ar if ar[k]['family']==fam)
            coverage.append(dict(run=run,family=fam,reference_CpGs=reference,observed_CpGs=len(ss),callable_CpGs=len(good),callable_fraction=len(good)/reference,valid_calls=sum(depth),single_observation_CpGs=depth.count(1),median_depth=median(depth),max_depth=max(depth),arrays_with_CpGs=len({x['array_id'] for x in good})))
    table(o/'methylation_sensitivity.tsv',results);table(o/'array_methylation_sensitivity.tsv.gz',sens_arrays);table(o/'methylation_coverage.tsv',coverage)
    table(o/'CpG_primary.tsv.gz',allprimary)
    bypos=collections.defaultdict(list)
    for x in allprimary:
        if int(x['valid_calls'])>0:bypos[x['family'],x['reference'],x['reference_C_position0']].append(x)
    common=[x for vv in bypos.values() if len(vv)==6 for x in vv];table(o/'common_six_CpGs.tsv',common)
    common_summary=[]
    for run in runs:
        for fam in ['CL1','CL2']:
            gg=collections.defaultdict(list)
            for x in common:
                if x['run']==run and x['family']==fam:gg[x['array_id']].append(x)
            common_summary.append(dict(run=run,family=fam,common_CpGs=sum(map(len,gg.values())),arrays=len(gg),array_equal_5mC=mean(mean(int(x['m'])/int(x['valid_calls']) for x in ll) for ll in gg.values())))
    table(o/'common_six_summary.tsv',common_summary)
    adj=[];patterns=[]
    for ch in sorted({x['seqid'] for x in ar.values()}):
        aa=sorted((x for x in ar.values() if x['seqid']==ch),key=lambda x:int(x['start0']))
        patterns.append(dict(seqid=ch,chromosome=aa[0]['chromosome'],AT_regions=sum(x['family']=='CL1' for x in aa),GC_regions=sum(x['family']=='CL2' for x in aa),order='-'.join('AT' if x['family']=='CL1' else 'GC' for x in aa),array_ids=';'.join(x['array_id'] for x in aa)))
        for a,b in zip(aa,aa[1:]):adj.append(dict(seqid=ch,array1=a['array_id'],family1=a['family'],array2=b['array_id'],family2=b['family'],gap_bp=int(b['start0'])-int(a['end0'])))
    table(o/'adjacent_arrays.tsv',adj);table(o/'chromosome_patterns.tsv',patterns)
    # Identical target arrays for within/same-/other-chromosome GC comparisons.
    pairs=rows(W/'array_metrics/between_arrays.tsv');samples=collections.defaultdict(list)
    for x in rows(W/'array_metrics/sampled_units.tsv.gz'):
        if x['variant']=='gap25000_size50000' and x['category'].endswith('_all') and int(x['rank'])<=4:samples[x['array_id']].append(x['sequence'])
    strata=[]
    for fam in ['CL1','CL2']:
        ss=collections.defaultdict(list);dd=collections.defaultdict(list)
        for x in pairs:
            if x['family']!=fam:continue
            for k in [x['array1'],x['array2']]:(ss if x['stratum']=='same_chromosome' else dd)[k].append(float(x['normalized_edit_distance']))
        shared=set(ss)&set(dd)
        ww=[]
        for aid in shared:
            t=Totals()
            for a,b in itertools.combinations(samples[aid],2):t.add(one(a,b))
            ww.append(t.row()['normalized_edit_distance'])
        strata.append(dict(family=fam,shared_arrays=len(shared),within=mean(ww) if ww else None,same_chromosome=mean(mean(ss[k]) for k in shared) if shared else None,different_chromosome=mean(mean(dd[k]) for k in shared) if shared else None))
    table(o/'shared_array_chromosome_comparison.tsv',strata)
    # Preserve contiguous chains; merging regions never bridges unknown units.
    chains=collections.defaultdict(list)
    for x in rows(E/'genomes/TU/local_alternation/contiguous_segments.tsv.gz'):
        if x['fragment'] in members:chains[members[x['fragment']]].append(x['order'])
    hor=[]
    for aid,x in ar.items():
        if x['family']=='CL2':
            for model in ['within_contiguous_segment','20_unit_local_blocks']:hor.append(dict(array_id=aid,**assess(chains[aid],aid,model)))
    table(o/'array_alternation.tsv',hor)
    (o/'computed.json').write_text(json.dumps(dict(large_arrays=len(ar),retained_fragments=len(members),CpG_rows=len(allprimary),shared_target_strata=strata,max_alternation=max(x['longest_alternating_units'] for x in hor)),indent=2)+'\n')
    print('all array summaries computed')
if __name__=='__main__':main()
