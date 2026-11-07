"""Complete congener two/three-template controls on every nuclear reference base.

This is the fixed-template arm, not template-independent family discovery or an
accepted array catalogue. Existing whole-nuclear TRF and unlabelled motif pools
are checked and kept as separate inputs. Never writes to previous result trees.
"""
import argparse,contextlib,csv,fcntl,hashlib,json,os,shutil,sys,time,traceback
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
from concurrent.futures import ThreadPoolExecutor,as_completed

W=Path(__file__).resolve().parents[1]
ROOT=W.parents[1]
E=ROOT/'output/danio_cl1_cl2_execution_2026-11-03'
OLD=ROOT/'output/danio_satellite_review_2026-11-04'
sys.path.insert(0,str(E/'code'))
import whole_fixed
from capacity import memory

O=W/'genome_sd';O.mkdir(parents=True,exist_ok=True)
now=lambda:datetime.now(ZoneInfo('Asia/Hong_Kong')).isoformat()
def save(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n');t.replace(p)
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
GIB=1024**3
@contextlib.contextmanager
def reserve(*args):
    avail,total=memory()
    if avail<max(256*GIB,.15*total)+6*GIB:raise RuntimeError('memory reserve would be exhausted')
    if shutil.disk_usage(W).free<302*GIB:raise RuntimeError('workspace reserve below 300 GiB + task allowance')
    yield
whole_fixed.reserve=reserve
whole_fixed.ROOT=W

def prepare():
    (W/'inputs').mkdir(exist_ok=True)
    for panel in ['two','three']:
        src=E/'inputs'/f'{panel}.fa';dst=W/'inputs'/f'{panel}.fa'
        if not dst.exists():shutil.copy2(src,dst)
        assert sha(src)==sha(dst)
    resources=[];alljobs=[]
    for species,acc in [('aesculapii','GCA_903798145.1'),('kyathit','GCA_903798195.1')]:
        ready=E/'inputs/reference'/acc/'ready.json';r=json.loads(ready.read_text())
        ref=Path(r['fasta']);assert ref.exists() and Path(str(ref)+'.fai').exists()
        lengths={x['seqid']:x['length'] for x in r['inventory'] if x['nuclear']}
        mpath=OLD/'results'/species/'trf_manifest.json';m=json.loads(mpath.read_text())
        assert m['reference']==str(ref) and dict(m['sequences'])==lengths
        checked=0;success=0
        for i,group in enumerate(m['groups']):
            d=OLD/'results'/species/'trf'/f'B{i:05d}';s=json.loads((d/'state.json').read_text())
            expected=sum(end-start for ch,start,end,left,right in group)
            assert s['status']=='success' and s['returncode']==0 and s['core_bp']==expected
            assert (d/'calls.tsv.gz').is_file() and (d/'native.dat.gz').is_file()
            checked+=expected;success+=1
        assert checked==sum(lengths.values())==m['core_bp']
        motif=OLD/'results'/species/'all_discovered_motifs.tsv.gz'
        assert motif.exists()
        res=dict(species=species,accession=acc,reference=str(ref),ready_sha256=sha(ready),
            compressed_source_MD5=r['md5'],reference_sha256=sha(ref),
            nuclear_records=len(lengths),nuclear_bp=checked,
            reused_TRF_batches=success,TRF_manifest=str(mpath),TRF_manifest_sha256=sha(mpath),
            unlabelled_motif_archive=str(motif),motif_sha256=sha(motif),
            scope='all nuclear assembled and unplaced records; organelle excluded')
        resources.append(res)
        for ch,n in lengths.items():
            last=0
            for start in range(0,n,whole_fixed.CORE):
                end=min(n,start+whole_fixed.CORE)
                assert start==last;last=end
                region=(ch,start,end,max(0,start-whole_fixed.PAD),min(n,end+whole_fixed.PAD))
                alljobs.append((species,ref,O/species/f'T{len([x for x in alljobs if x[0]==species]):05d}',region))
            assert last==n
    save(O/'inputs_checked.json',dict(at_hkt=now(),resources=resources,
        previous_analysis_read_only=True,full_SD_not_yet_verified=True))
    # Representative repeated and non-repeated cores from both species first.
    priorities=[]
    for species in ['aesculapii','kyathit']:
        tab=OLD/'results'/species/'selected_regions.tsv'
        rows=list(csv.DictReader(tab.open(),delimiter='\t'))
        selected=[next(x for x in rows if x['family']==fam) for fam in ['AT','GC']]
        for x in selected:
            hit=next(j for j in alljobs if j[0]==species and j[3][0]==x['seqid'] and j[3][1]<=int(x['start0'])<j[3][2])
            if hit not in priorities:priorities.append(hit)
        hit=next(j for j in alljobs if j[0]==species)
        if hit not in priorities:priorities.append(hit)
    ordered=priorities+[j for j in alljobs if j not in priorities]
    save(O/'plan.json',dict(at_hkt=now(),core_bp=whole_fixed.CORE,padding_bp=whole_fixed.PAD,
        panels=['two','three'],threads_per_job=4,print_min_identity=0,
        pilot_tiles=[str(j[2].relative_to(W)) for j in priorities],
        jobs=[dict(species=j[0],directory=str(j[2].relative_to(W)),region=j[3]) for j in ordered],
        expected_core_bp=sum(j[3][2]-j[3][1] for j in ordered),
        code_dependencies={str(p):sha(p) for p in [Path(__file__),Path(whole_fixed.__file__)]},
        fixed_panels={p:sha(W/'inputs'/f'{p}.fa') for p in ['two','three']},
        array_threshold=None,completion_meaning='full fixed-template decomposition only; candidate annotation and scientific review separate'))
    return ordered,len(priorities)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=24)
    ap.add_argument('--pilot-only',action='store_true');ap.add_argument('--prepare-only',action='store_true');a=ap.parse_args()
    assert 1<=a.workers<=36
    with (O/'worker.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        jobs,n=prepare()
        if a.prepare_only:print('prepared',len(jobs),'tiles;',n,'pilot tiles');return
        if a.pilot_only:jobs=jobs[:n]
        # A full run waits for bounded pilot validation before submitting others.
        if not a.pilot_only and not (O/'pilot_checked.json').exists():
            raise RuntimeError('run --pilot-only first; inspect failures before full submission')
        completed=bp=0;errors=[]
        with ThreadPoolExecutor(max_workers=a.workers) as pool:
            futures={pool.submit(whole_fixed.job,*j):j for j in jobs}
            for future in as_completed(futures):
                job=futures[future]
                try:
                    future.result();q=json.loads((job[2]/'complete.json').read_text())
                    assert q['core']==list(job[3]) and q['scanned_core_bp']==job[3][2]-job[3][1]
                    completed+=1;bp+=q['scanned_core_bp']
                except Exception as ex:errors.append(dict(species=job[0],tile=job[2].name,error=repr(ex)))
                save(O/'progress.json',dict(at_hkt=now(),stage='bounded_pilot' if a.pilot_only else 'whole_nuclear_fixed_SD',
                    successful_two_panel_tiles=completed,total_tiles=len(jobs),verified_core_bp=bp,
                    expected_core_bp=sum(j[3][2]-j[3][1] for j in jobs),workers=a.workers,threads_per_task=4,errors=errors))
        assert not errors,errors[:4]
        assert bp==sum(j[3][2]-j[3][1] for j in jobs)
        save(O/('pilot_checked.json' if a.pilot_only else 'fixed_SD_checked.json'),dict(at_hkt=now(),
            successful_tiles=completed,verified_core_bp=bp,panels=['two','three'],
            family_and_array_catalogue_accepted=False,scope='technical completion; no automatic family or array claims'))
        try:print('verified tiles',completed,'core bp',bp,flush=True)
        except BrokenPipeError:pass  # Completed receipts survive a closed reporting terminal.

if __name__=='__main__':
    try:main()
    except Exception as e:
        save(O/'failure.json',dict(at_hkt=now(),error=repr(e),traceback=traceback.format_exc()));raise
