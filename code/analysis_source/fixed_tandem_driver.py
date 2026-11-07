"""One-shot dependent extraction of all new congener fixed-panel tandem runs."""
import os
import concurrent.futures,contextlib,csv,fcntl,gzip,hashlib,json,sys,time,collections,traceback
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
W=Path(__file__).resolve().parents[1];ROOT=W.parents[1]
E=ROOT / os.environ.get('DANIO_ASSEMBLY_RUN_ROOT', 'inputs/assembly_runs')
sys.path.insert(0,str(E/'code'))
import fixed_tandem
fixed_tandem.ROOT=W
O=W/'fixed_tandem';O.mkdir(parents=True,exist_ok=True)
def save(p,v):
    p.parent.mkdir(parents=True,exist_ok=True);t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(v,indent=2)+'\n');t.replace(p)
def now():return datetime.now(ZoneInfo('Asia/Hong_Kong')).isoformat()
def work(p):
    species,tile=p.parent.parent.name,p.parent.name;d=O/species/tile;d.mkdir(parents=True,exist_ok=True)
    receipt=d/'checked.json'
    if receipt.exists():return
    core=json.loads(p.read_text())['core'];info={}
    for panel in ['two','three']:
        info[panel]=fixed_tandem.extract(p.parent/panel/'final_decomposition.tsv.gz',panel,core,d)
        assert info[panel]['rows']==json.loads((p.parent/panel/'verified.json').read_text())['records']
    save(receipt,dict(at_hkt=now(),core=core,panels=info,biological_array_admission=False))
def main():
    lock=(O/'worker.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    expected=len(json.loads((W/'genome_sd/plan.json').read_text())['jobs'])
    done=set();completed=0
    with concurrent.futures.ThreadPoolExecutor(6) as pool:
        while True:
            paths=[p for p in (W/'genome_sd').glob('*/T*/complete.json') if str(p) not in done]
            for p in pool.map(work,paths):completed+=1
            done.update(str(p) for p in paths)
            save(O/'progress.json',dict(at_hkt=now(),verified_tiles=completed,expected_tiles=expected,waiting_upstream=not (W/'genome_sd/fixed_SD_checked.json').exists()))
            if (W/'genome_sd/failure.json').exists():raise RuntimeError('upstream failure requires diagnosis')
            if (W/'genome_sd/fixed_SD_checked.json').exists() and completed==expected:break
            time.sleep(15)
    # Exact union of core-clipped native tandem candidates, retaining both panels.
    intervals=collections.defaultdict(list)
    for p in O.glob('*/T*/*.runs.tsv.gz'):
        species=p.parent.parent.name;panel=p.name.split('.')[0]
        with gzip.open(p,'rt') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                intervals[species,panel,r['family'],r['seqid']].append((int(r['start0']),int(r['end0'])))
    dest=O/'diagnostic_segments.tsv'
    with dest.open('w') as f:
        writer=csv.writer(f,delimiter='\t');writer.writerow(['species','panel','family','seqid','start0','end0','span_bp','scope'])
        for key,iv in sorted(intervals.items()):
            merged=[]
            for a,b in sorted(iv):
                if merged and a<=merged[-1][1]:merged[-1][1]=max(b,merged[-1][1])
                else:merged.append([a,b])
            for a,b in merged:writer.writerow([*key,a,b,b-a,'native75_tandem_screen_not_family_or_array_confirmation'])
    save(O/'computed.json',dict(at_hkt=now(),tiles=completed,scope='full-reference fixed-panel candidate evidence; family and raw-sequence validation still required'))
if __name__=='__main__':
    try:main()
    except Exception as e:save(O/'failure.json',dict(at_hkt=now(),error=repr(e),traceback=traceback.format_exc()));raise
