"""Re-extract raw copies for whole-nuclear supported congener regions.

Primary fixed panel and TU source-sequence/phase rules are matched. A local
monomer panel and omission of the copy-score screen are separate sensitivities.
Region membership is frozen before copy filtering; copy failure is not absence.
"""
import csv,gzip,json,hashlib,subprocess,os,sys,itertools,collections
from pathlib import Path
from statistics import mean
from concurrent.futures import ProcessPoolExecutor,as_completed
import pysam
W=Path(__file__).resolve().parents[1]
ROOT=W.parents[1]
E=ROOT / os.environ.get('DANIO_ASSEMBLY_RUN_ROOT', 'inputs/assembly_runs')
R=ROOT / os.environ.get('DANIO_ANNOTATION_ROOT', 'inputs/annotation')
sys.path.insert(0,str(E/'code'))
from catalog_units import normalize,FAM,rawfamily
from pilot import SD
sys.path.insert(0,str(R/'07_edit_distance_review_20261105/code'))
from metrics import total
from regional_catalog import rows,table,save,union,sha
O=W/'congener_metrics';O.mkdir(exist_ok=True)
def work(job):
    reg,tracts,local=job;d=O/reg['array_id'];d.mkdir(exist_ok=True)
    if (d/'checked.json').exists():return json.loads((d/'checked.json').read_text())
    meta=json.loads((E/'inputs/reference'/reg['assembly']/'ready.json').read_text())
    a,b=int(reg['start0']),int(reg['end0']);fam=reg['family'];L=len(FAM[fam])
    with pysam.FastaFile(meta['fasta']) as fa:
        lo=max(0,a-1000);hi=min(fa.get_reference_length(reg['seqid']),b+1000);seq=fa.fetch(reg['seqid'],lo,hi).upper()
    (d/'raw.fa').write_text('>region\n'+seq+'\n')
    templates={'fixed':W/'inputs/three.fa'}
    if local:
        (d/'local.fa').write_text('>local_unit\n'+local+'\n');templates['local']=d/'local.fa'
    obs=[];allunits=[]
    for panel,template in templates.items():
        out=d/panel;out.mkdir(exist_ok=True)
        cmd=[SD,str(d/'raw.fa'),str(template),'-t','2','-i','0','--second-best','-o',str(out)]
        if not (out/'final_decomposition.tsv').exists():
            with (out/'worker.log').open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env={**os.environ,'OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'1'})
            assert p.returncode==0,(reg['array_id'],panel,p.returncode)
        blocks=[]
        with (out/'final_decomposition.tsv').open() as f:
            for v in csv.reader(f,delimiter='\t'):
                s,e=lo+int(v[2]),lo+int(v[3])+1
                if s<b and e>a:blocks.append((s,e,v[1],float(v[4])))
        blocks.sort();units=[]
        for i,(s,e,tag,score) in enumerate(blocks):
            raw=seq[s-lo:e-lo];normalized='';phase='';sourcefam='unchecked'
            if not any(x<s<e<y for x,y in tracts):status='boundary_or_interruption'
            elif (i and blocks[i-1][1]>s) or (i+1<len(blocks) and blocks[i+1][0]<e):status='overlapping_block'
            elif set(raw)-set('ACGT'):status='ambiguous_sequence'
            elif not .72*L<=len(raw)<=1.28*L:status='not_single_unit_length'
            else:
                phase,normalized=normalize(raw,fam);sourcefam=rawfamily(raw)
                status='normalized' if normalized else 'phase_unresolved'
            u=dict(array_id=reg['array_id'],species=reg['species'],family=fam,panel=panel,seqid=reg['seqid'],start0=s,end0=e,tag=tag,native_identity=score,raw_sequence=raw,sequence=normalized,phase=phase,source_family=sourcefam,status=status)
            units.append(u);allunits.append(u)
        for screen in ['matched_TU','without_copy_family_screen']:
            pool=[x for x in units if x['status']=='normalized' and (screen!='matched_TU' or x['source_family']==fam)]
            pool.sort(key=lambda u:hashlib.sha256(('copy20261103|'+u['seqid']+':'+str(u['start0'])+'-'+str(u['end0'])).encode()).digest()[:8])
            sample=pool[:32]
            obs.append(dict(array_id=reg['array_id'],species=reg['species'],family=fam,panel=panel,copy_screen=screen,supported_bp=reg['supported_bp'],chromosome=reg['chromosome'],eligible_copies=len(pool),sampled_copies=len(sample),eligible_copy_bp=sum(int(x['end0'])-int(x['start0']) for x in pool),**total(itertools.combinations([x['sequence'] for x in sample],2))))
            table(d/(panel+'_'+screen+'_sample.tsv'),sample,fields=list(allunits[0]))
    with gzip.open(d/'units.tsv.gz','wt') as f:
        z=csv.DictWriter(f,list(allunits[0]),delimiter='\t');z.writeheader();z.writerows(allunits)
    result=dict(region=reg,local_template=local,source_sha256=hashlib.sha256(seq.encode()).hexdigest(),tracts=tracts,observations=obs,source_sequence_file=str(d/'raw.fa'),scope='actual complete source copies; fixed/native screen does not determine family region membership')
    save(d/'checked.json',result);return result
def main():
    rr=[x for x in rows(W/'congener_catalog/regions.tsv') if int(x['supported_bp'])>=5000]
    mm=rows(W/'congener_catalog/membership.tsv');src={x['source_id']:x for x in rows(W/'congener_identity/source_records.tsv')};jobs=[]
    for r in rr:
        ff=[x for x in mm if x['array_id']==r['array_id']];iv=union((int(x['start0']),int(x['end0'])) for x in ff)
        L=len(FAM[r['family']]);candidates=[src[x['source_id']] for x in ff if x['source_id'] in src and .85*L<=len(src[x['source_id']]['pattern'])<=1.15*L]
        candidates.sort(key=lambda x:(-(int(x['end0'])-int(x['start0'])),x['source_id']))
        local=candidates[0]['pattern'] if candidates else ''
        jobs.append((r,iv,local))
    results=[]
    with ProcessPoolExecutor(8) as pool:
        for future in as_completed([pool.submit(work,x) for x in jobs]):
            results.append(future.result());save(O/'progress.json',dict(completed=len(results),total=len(jobs)))
    obs=sorted([v for x in results for v in x['observations']],key=lambda x:(x['array_id'],x['panel'],x['copy_screen']));table(O/'within_arrays.tsv',obs)
    summary=[]
    for sp in ['aesculapii','kyathit']:
      for fam in ['CL1','CL2']:
       for cut in [5000,10000,15000,20000,25000,50000]:
        for panel in ['fixed','local']:
         for screen in ['matched_TU','without_copy_family_screen']:
            allregs=[x for x in rr if x['species']==sp and x['family']==fam and int(x['supported_bp'])>=cut]
            vv=[x for x in obs if x['species']==sp and x['family']==fam and int(x['supported_bp'])>=cut and x['panel']==panel and x['copy_screen']==screen]
            valid=[x for x in vv if x['pairs']>0]
            summary.append(dict(species=sp,family=fam,minimum_supported_bp=cut,panel=panel,copy_screen=screen,regions=len(allregs),regions_with_panel=len(vv),comparable_arrays=len(valid),eligible_copies=sum(x['eligible_copies'] for x in valid),sampled_copies=sum(x['sampled_copies'] for x in valid),array_equal_edit_distance=mean(x['normalized_edit_distance'] for x in valid) if valid else None,array_equal_event_distance=mean(x['event_sensitivity'] for x in valid) if valid else None,scope='recoverable assembly copies; species-specific size strata not quantitative abundance comparisons'))
    table(O/'family_summary.tsv',summary)
    save(O/'computed.json',dict(regions=len(results),primary=dict(size_bp=10000,panel='fixed',copy_screen='matched_TU'),size_sensitivity=[5000,15000,20000],old_window_means_reused=False,local_panel_missing=[x['region']['array_id'] for x in results if not x['local_template']],source_files={str(p.relative_to(W)):sha(p) for p in [W/'congener_catalog/regions.tsv',W/'congener_catalog/membership.tsv']}))
    print('computed congener regions',len(results),flush=True)
if __name__=='__main__':main()
