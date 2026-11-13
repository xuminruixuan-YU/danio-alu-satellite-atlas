"""Recount fixed array-zone support from accepted reference-CpG projections.

Requires accepted projection artifacts (counts, reference CpG arrays, retained
primary modBAMs and returned candidate records). Does not rerun mapping or
relax filters. For deposited-summary reproduction use summarize.py instead.
"""
import argparse,json,gzip,itertools,sys,time
from bisect import bisect_right
from collections import Counter,defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo
import numpy as np
import pandas as pd
import pysam

REPO=Path(__file__).resolve().parents[2]
WIDTH=20000

def now():return datetime.now(ZoneInfo('Asia/Hong_Kong')).isoformat()
def save(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)

def site_statistics(positions, states, start, end):
    left, right = np.searchsorted(positions, [start, end])
    pos, calls = positions[left:right], states[left:right, :3]
    depth = calls.sum(axis=1, dtype=np.uint64)
    good = depth > 0
    pos, calls, depth = pos[good], calls[good], depth[good]
    observations = int(depth.sum())
    count = len(pos)
    if count:
        gap = int(max(pos[0] - start, end - pos[-1] - 1, np.diff(pos, prepend=pos[0]).max() - 1))
        center_count = int(np.searchsorted(pos, start + .75 * (end - start))
                           - np.searchsorted(pos, start + .25 * (end - start)))
        mean = float(np.mean(calls[:, 1] / depth) * 100)
        weighted = float(calls[:, 1].sum(dtype=np.uint64) / observations * 100)
    else:
        gap, center_count, mean, weighted = end - start, 0, np.nan, np.nan
    return dict(reference_CpGs=right-left, CpGs=count, observations=observations,
                mean=mean, weighted=weighted, max_gap=gap, central_CpGs=center_count,
                edge_flag=bool(gap > (end-start)/2 or center_count == 0))


def merge_intervals(keys):
    intervals = []
    for start, end in sorted((a, b) for _, a, b in keys):
        if intervals and start <= intervals[-1][1]:
            intervals[-1][1] = max(intervals[-1][1], end)
        else:
            intervals.append([start, end])
    return tuple(zip(*intervals))


def parse_parents(group, columns):
    parents = []
    for line in group:
        row = line.rstrip('\n').split('\t')
        if columns is None:
            score = next((int(tag[5:]) for tag in row[11:] if tag.startswith('AS:i:')), None)
            parents.append(dict(flag=int(row[1]), reference=row[2], start0=int(row[3])-1,
                                AS=score, CIGAR=row[5]))
        else:
            flag, ref, start, score, cigar = (row[index] for index in columns)
            parents.append(dict(flag=int(flag), reference=ref, start0=int(start),
                                AS=int(score) if score else None, CIGAR=cigar))
    return parents


def replay_part(part, targets, positions, progress):
    """Read paired saved streams and project only reads touching requested windows."""
    fmt = part['parents_format']
    opener = open if fmt == 'SAM' else gzip.open
    diagnostics = Counter()
    with opener(part['parents'], 'rt') as handle, pysam.AlignmentFile(part['bam'], 'rb', threads=2) as bam:
        columns = None
        if fmt == 'TSV':
            header = next(handle).rstrip('\n').split('\t')
            assert header[0] == 'read'
            columns = [header.index(field) for field in ('flag','reference','start0','AS','CIGAR')]
        records = (line for line in handle if not line.startswith('@')) if fmt == 'SAM' else handle
        parents = itertools.groupby(records, key=lambda line: line[:line.index('\t')])
        for seen, alignment in enumerate(bam, start=1):
            name, group = next(parents)
            assert alignment.query_name == name, (part['name'], seen, name, alignment.query_name)
            if seen % 100_000 == 0:
                guard()
                progress(seen)
            if alignment.is_unmapped or alignment.mapping_quality < 20:
                continue
            reference = alignment.reference_name
            target = targets.get(reference)
            if target is None:
                continue
            left, right = target['unions']
            interval = bisect_right(right, alignment.reference_start)
            if interval == len(right) or left[interval] >= alignment.reference_end:
                continue
            result = project(alignment, parse_parents(group, columns), positions, diagnostics)
            if result is None:
                continue
            _, _, states, _, coordinates = result
            valid = np.sort(coordinates[states < 3])
            if not len(valid):
                continue
            starts, ends = target['starts'], target['ends']
            a = np.searchsorted(starts, valid[0]-WIDTH, side='right')
            b = np.searchsorted(starts, valid[-1], side='right')
            if a < b:
                n = np.searchsorted(valid, ends[a:b]) - np.searchsorted(valid, starts[a:b])
                target['observations'][a:b] += n
                target['reads'][a:b] += n > 0
                for offset in np.flatnonzero(n > 0):
                    target['read_IDs'][a+int(offset)].add(alignment.query_name)
        assert next(parents, None) is None and seen == part['records'], part['name']
    progress(seen)
    return diagnostics['input_reads']


def new_read_counts(item, values, keys, progress_path):
    targets = {}
    by_ref = defaultdict(list)
    for key in keys:
        by_ref[key[0]].append(key)
    for reference, pending in by_ref.items():
        pending.sort(key=lambda x: (x[1], x[2]))
        targets[reference] = dict(keys=pending, unions=merge_intervals(pending),
                                  starts=np.array([x[1] for x in pending]),
                                  ends=np.array([x[2] for x in pending]),
                                  reads=np.zeros(len(pending), dtype=np.int64),
                                  read_IDs=[set() for _ in pending],
                                  observations=np.zeros(len(pending), dtype=np.int64))
    if not targets:
        return {}, 0
    run = item['run']
    positions = {ref: np.load(SURVEY / 'inputs/reference_CpGs' / f'{ref}.npy').astype(np.int64, copy=False) for ref in targets}
    parts = item['parts'] + [dict(name='missing_recovered_2.26', parents_format='TSV',
                                 bam=str(SURVEY / 'per_run' / run / 'new.primary.mod.bam'),
                                 parents=str(SURVEY / 'per_run' / run / 'new_alignment_records.tsv.gz'),
                                 records=item['missing_mapping_records'])]
    elapsed, processed, projected = time.monotonic(), 0, 0
    for part in parts:
        def progress(seen):
            save(progress_path, dict(at_HKT=now(), run=run, stage='recount_window_read_IDs', part=part['name'],
                                     source_records_checked=processed+seen, source_records=item['source_records'],
                                     pending_distinct_windows=len(keys), elapsed_seconds=time.monotonic()-elapsed))
        projected += replay_part(part, targets, positions, progress)
        processed += part['records']
    assert processed == item['source_records']
    reads = {}
    for target in targets.values():
        for key, count, observation_count, ids in zip(target['keys'], target['reads'], target['observations'], target['read_IDs']):
            assert int(observation_count) == values[key]['observations'], (run, key, observation_count, values[key])
            assert len(ids) == int(count), (run, key, 'duplicate primary read IDs')
            reads[key] = ids
    return reads, projected


def sample(item):
    run=item['run'];dest=D/'per_run'/run;dest.mkdir(exist_ok=True)
    started=time.monotonic()
    assert json.loads((SURVEY/'per_run'/run/'accepted.json').read_text())['passed']
    assert json.loads((SURVEY/'per_run'/run/'existing_region_reproduction.json').read_text())['passed']
    config=json.loads((SURVEY/'inputs/config.json').read_text())
    lengths={r['reference']:r['length'] for r in config['reference_records']}
    arrays=pd.read_csv(REPO/'data/array_study/array_catalog/variants/TU_gap25000/regions.tsv',sep='\t')
    arrays=arrays[arrays.supported_bp>=50000]
    old=pd.read_csv(REPO/'data/regional_display/per_run'/run/'valid_observations.tsv.gz',sep='\t',usecols=['read','reference','position0','call'])
    assert not old.duplicated(['read','reference','position0']).any()
    oldrefs={ref:g.sort_values('position0') for ref,g in old.groupby('reference',sort=False)}
    rows=[];values={};outer={};old_ids={};links={}
    current=None
    for ar in arrays.sort_values(['seqid','start0']).itertuples():
        ref=ar.seqid;s,e=int(ar.start0),int(ar.end0)
        assert s>=500000 and e+500000<=lengths[ref]
        if ref!=current:
            positions=np.load(SURVEY/'inputs/reference_CpGs'/f'{ref}.npy')
            with np.load(SURVEY/'per_run'/run/'counts'/f'{ref}.npz') as z:counts=z['counts']
            current=ref
        for flank in (200000,500000):
            for zone,a,b in [('left',s-flank,s),('inside',s,e),('right',e,e+flank)]:
                key=(ref,a,b)
                if key not in values:values[key]=site_statistics(positions,counts,a,b)
                rows.append(dict(run=run,array_id=ar.array_id,family=ar.family,chromosome=ar.chromosome,
                                 reference=ref,flank_bp=flank,zone=zone,start0=a,end0=b,**values[key]))
                if flank==200000:
                    g=oldrefs.get(ref)
                    if g is None:q=old.iloc[:0]
                    else:
                        lo,hi=np.searchsorted(g.position0.to_numpy(),[a,b]);q=g.iloc[lo:hi]
                    stats=values[key]
                    assert len(q)==stats['observations'] and q.position0.nunique()==stats['CpGs'],(run,key)
                    c=q.assign(modified=(q.call=='m')).groupby('position0').modified.agg(['sum','count'])
                    mean=float((c['sum']/c['count']).mean()*100) if len(c) else np.nan
                    np.testing.assert_allclose(mean,stats['mean'],rtol=0,atol=1e-9,equal_nan=True)
                    old_ids[key]=set(q.read)
        for zone,a,b,inner in [('left',s-500000,s-200000,(ref,s-200000,s)),('right',e+200000,e+500000,(ref,e,e+200000))]:
            outerkey=(ref,a,b);outer[outerkey]=site_statistics(positions,counts,a,b)
            links[(ar.array_id,zone)]=(inner,outerkey)
    del old,oldrefs
    global WIDTH
    WIDTH=300000
    pending=[k for k,v in outer.items() if v['observations']]
    print(run,'reused zones',len(old_ids),'outer 300 kb intervals',len(pending),flush=True)
    extra,projected=new_read_counts(item,outer,pending,dest/'progress.json')
    for row in rows:
        key=(row['reference'],row['start0'],row['end0'])
        if row['flank_bp']==200000 or row['zone']=='inside':
            row['reads']=len(old_ids[key]);row['read_count_source']='validated_200kb_observation_ledger'
        else:
            inner,outerkey=links[(row['array_id'],row['zone'])]
            row['reads']=len(old_ids[inner] | extra.get(outerkey,set()))
            row['read_count_source']='read_ID_union_of_200kb_ledger_and_outer_300kb_replay'
            assert row['observations']==values[inner]['observations']+outer[outerkey]['observations']
        assert 0<=row['reads']<=row['observations']
    frame=pd.DataFrame(rows)
    expected=pd.read_csv(REPO/'data/regional_display/results/background.tsv',sep='\t')
    expected=expected[expected.run==run].set_index('array_id')
    for r in frame[frame.flank_bp==200000].itertuples():
        src=expected.loc[r.array_id]
        for f in ['CpGs','observations','reads','max_gap','mean']:
            np.testing.assert_allclose(getattr(r,f),src[f'{r.zone}_{f}'],rtol=0,atol=1e-9,equal_nan=True)
    frame.to_csv(dest/'backgrounds.tsv',sep='\t',index=False)
    result=dict(passed=True,run=run,zones=len(rows),outer_intervals=len(outer),
                target_primary_reads_projected=projected,seconds=time.monotonic()-started,
                all_replayed_outer_observation_totals_match=True,
                all_200kb_counts_means_and_read_IDs_reproduced_from_ledger=True,
                exact_500kb_read_ID_union=True,duplicate_primary_read_IDs_in_replayed_intervals=0)
    save(dest/'checked.json',result);print(json.dumps(result),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=6)
    args=parser.parse_args();SURVEY=args.analysis_root.resolve();D=args.output.resolve()
    if D.exists() and any(D.iterdir()):parser.error('Choose a new or empty output directory')
    for name in ['per_run','data','checks']:(D/name).mkdir(parents=True,exist_ok=True)
    sys.path.insert(0,str(REPO/'code/genome_methylation'))
    import engine
    engine.D=SURVEY
    guard,project=engine.guard,engine.project
    config=json.loads((SURVEY/'inputs/config.json').read_text())
    assert config['assembly']=='GCF_049306965.1'
    assert config['MAPQ']==20 and config['score_fraction']==0.8 and config['baseQ']==10
    with ProcessPoolExecutor(max_workers=args.workers) as ex:receipts=list(ex.map(sample,config['items']))
    frames=[pd.read_csv(D/'per_run'/i['run']/'backgrounds.tsv',sep='\t') for i in config['items']]
    pd.concat(frames,ignore_index=True).to_csv(D/'data/zone_backgrounds.tsv',sep='\t',index=False)
    save(D/'checks/backgrounds.json',dict(passed=True,samples=receipts))
