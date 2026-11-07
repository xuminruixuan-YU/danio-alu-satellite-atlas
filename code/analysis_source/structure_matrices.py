"""ModDotPlot matrices on each large array plus 200 kb on both sides."""
import concurrent.futures,csv,fcntl,hashlib,json,os,shutil,subprocess,time
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import pysam
W=Path(__file__).resolve().parents[1];ROOT=W.parents[1]
E=ROOT / os.environ.get('DANIO_ASSEMBLY_RUN_ROOT', 'inputs/assembly_runs')
O=W/'structure';O.mkdir(exist_ok=True)
BIN=os.environ.get('MODDOTPLOT', 'moddotplot')
def now():return datetime.now(ZoneInfo('Asia/Hong_Kong')).isoformat()
def save(p,v):
    t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n');t.replace(p)
def sha(p):
    with open(p,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def job(x):
    d=O/x['array_id'];d.mkdir(exist_ok=True)
    if (d/'checked.json').exists():return json.loads((d/'checked.json').read_text())
    if x['species']=='TU':ref=json.loads((E/'external/reference/inputs.json').read_text())['source']
    else:ref=json.loads((E/'inputs/reference'/x['assembly']/'ready.json').read_text())['fasta']
    a,b=int(x['plot_start0']),int(x['plot_end0'])
    with pysam.FastaFile(ref) as f:seq=f.fetch(x['seqid'],a,b).upper()
    assert len(seq)==b-a
    fasta=d/'region.fa';fasta.write_text('>'+x['array_id']+'\n'+'\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')
    # 5 kb is fixed for all main large arrays; short-array views are separate.
    cmd=[BIN,'static','-f',str(fasta),'-k','17','-w','5000','-m','150','-id','80','-d','0.5','-o',str(d/'matrix'),'--no-plot']
    if shutil.disk_usage(W).free<300*1024**3:raise RuntimeError('disk reserve below 300GiB')
    # Existing successful native output may be checked again after a validator
    # fix, without rerunning matrices or confusing a check with computation.
    if not list((d/'matrix').rglob('*.bedpe')):
        with open(d/'program.log','w') as log:
            p=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'MPLCONFIGDIR':'/tmp/danio-array-revision-mpl','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
        assert p.returncode==0,(x['array_id'],p.returncode)
    beds=list((d/'matrix').rglob('*.bedpe'));assert len(beds)==1,(d,beds)
    n=0
    with open(beds[0]) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            qs,qe,rs,re=map(int,[r['query_start'],r['query_end'],r['reference_start'],r['reference_end']])
            assert 1<=qs<=qe and 1<=rs<=re and qs<=b-a and rs<=b-a
            assert qe-qs+1==5000 and re-rs+1==5000
            assert 80<=float(r['perID_by_events'])<=100
            n+=1
    z=dict(at_hkt=now(),array=x,command=cmd,version='0.9.8',input_fasta_SHA256=sha(fasta),BEDPE=str(beds[0]),BEDPE_SHA256=sha(beds[0]),matrix_cells=n,
           BEDPE_coordinate_semantics='exported BEDPE: one-based closed intervals and identity in percent; last bins clipped to input length when rendering',
           scope='windowed sequence self-similarity; 80% display cutoff is not a family admission or monomer identity threshold')
    save(d/'checked.json',z);return z
def main():
    lock=open(O/'worker.lock','a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    with open(W/'array_catalog/regions.tsv') as f:rr=[x for x in csv.DictReader(f,delimiter='\t') if x['qualifies_large']=='True']
    results=[]
    with concurrent.futures.ThreadPoolExecutor(8) as p:
        fs={p.submit(job,x):x for x in rr}
        for f in concurrent.futures.as_completed(fs):
            results.append(f.result());save(O/'progress.json',dict(at_hkt=now(),completed=len(results),total=len(rr),workers=8))
    save(O/'matrices_checked.json',dict(at_hkt=now(),arrays=len(results),matrices=results,figure_rendering_complete=False))
    print('checked',len(results),'matrices')
if __name__=='__main__':main()
