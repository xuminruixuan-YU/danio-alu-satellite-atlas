"""Recover unthresholded fixed-panel units within each released family component."""
import bisect
import csv
import fcntl
import gzip
import json
import time
from collections import Counter,defaultdict
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache

import edlib
import pysam
from Bio import SeqIO

from capacity import reserve
from catalog_evidence import rawfamily,table
from genome_trf import obtain
from library import rc,write_tsv
from pilot import ROOT,atomic

G=ROOT/'genomes'
T={r.id:str(r.seq) for r in SeqIO.parse(ROOT/'inputs/three.fa','fasta')}
FAM={'CL1':T['CL1_186nt'],'CL2':T['CL2_half1']}
TICKET={'TU':'04','SAT':'05','Fish11':'05F'}
FIELDS='fragment seqid start0 end0 family sd_strand raw_tag native_identity source_subtype normalization raw_sequence oriented_circular_sequence within_fragment source_family pass_status'.split()

@lru_cache(maxsize=32768)
def normalize(s,family):
    t=FAM[family];L=len(t)
    opts=[]
    for strand,q in [('+',s),('-',rc(s))]:
        a=edlib.align(q,t*3,mode='HW',task='locations')
        loc={(l%L,h-l+1) for l,h in a['locations']}
        opts.append((a['editDistance'],strand,loc,q))
    cost=min(o[0] for o in opts);best=[o for o in opts if o[0]==cost]
    if len(best)!=1:return 'direction_tie',''
    _,strand,locations,q=best[0]
    if len(locations)!=1:return 'cycle_offset_or_span_tie',''
    start,span=next(iter(locations));ref=(t*3)[start:start+span]
    if start==0:return strand,q
    target=edlib.align(q,ref,mode='NW',task='path')
    aligned=edlib.getNiceAlignment(target,q,ref)
    ti=start;qi=0
    for base,tbase in zip(aligned['query_aligned'],aligned['target_aligned']):
        if ti==L and tbase!='-':
            return strand,q[qi:]+q[:qi]
        ti+=tbase!='-';qi+=base!='-'
    return 'origin_not_covered',''

def process(job):
    name,ch,grid=job;base=G/name;d=base/'catalog_units'/ch
    if (d/'checked.json').exists():return
    d.mkdir(exist_ok=True,parents=True)
    with reserve('annotated_units_'+name,1,2,.2):
        ps=defaultdict(list)
        for p in table(base/'catalog/fragments.tsv'):
            if p['seqid']==ch:ps[p['family']].append((int(p['start0']),int(p['end0']),p['fragment']))
        keys={fam:[v[0] for v in sorted(ps[fam])] for fam in FAM}
        for fa in ps:ps[fa].sort()
        totals=Counter()
        ref,_=obtain(name)
        with pysam.FastaFile(str(ref)) as f,gzip.open(d/'units.tsv.gz','wt',compresslevel=2) as out:
            w=csv.writer(out,delimiter='\t',lineterminator='\n');w.writerow(FIELDS)
            for tile,(_,a,b,lo,hi) in grid:
                if not any(l<b and h>a for vals in ps.values() for l,h,_ in vals):continue
                ph=base/'phase_comparison'/tile/'subtype_blocks.tsv.gz';sub={}
                if ph.exists():
                    sub={(int(r['start0']),int(r['end0']),r['strand']):r['call'] for r in table(ph)}
                seq=f.fetch(ch,lo,hi).upper();hits=defaultdict(set)
                with gzip.open(base/'fixed'/tile/'three/final_decomposition.tsv.gz','rt') as inp:
                    for v in csv.reader(inp,delimiter='\t'):
                        x,y=lo+int(v[2]),lo+int(v[3])+1
                        if a<=x<b:hits[x,y].add((v[1],float(v[4])))
                for (x,y),labels in sorted(hits.items()):
                    for family in FAM:
                        j=bisect.bisect_right(keys[family],x)-1
                        if j<0:continue
                        l,h,fid=ps[family][j]
                        if x>=h:continue
                        observed=seq[x-lo:y-lo]
                        rank=max(v for _,v in labels);tops=sorted(k for k,v in labels if v==rank)
                        strands={'-' if k.endswith("'") else '+' for k in tops}
                        strand=next(iter(strands)) if len(strands)==1 else '?'
                        size=len(observed)/len(FAM[family]);source='unchecked';norm='';aligned=''
                        if y>h:status='truncated_at_piece_end'
                        elif set(observed)-set('ACGT'):status='ambiguous_bases'
                        elif not .72<=size<=1.28:status='not_expected_single_unit_span'
                        else:
                            source=rawfamily(observed)
                            if source!=family:status='source_family_unresolved_or_conflicting'
                            else:
                                norm,aligned=normalize(observed,family)
                                status='normalized' if aligned else 'source_family_only'
                        label=sub.get((x,y,strand),'unknown') if family=='CL2' else 'not_applicable'
                        if source!=family or strand not in ['+','-'] or (aligned and norm!=strand):
                            if family=='CL2':label='unknown'
                        w.writerow([fid,ch,x,y,family,strand,'|'.join(tops),rank,label,norm,observed,aligned,y<=h,source,status])
                        totals[family,status]+=1
        atomic(d/'checked.json',dict(counts={f+':'+s:c for (f,s),c in totals.items()},scope='Includes low identity winning native blocks without a 90% cutoff. Owner is the reference start; clipping, ambiguous cycles, ties, and atypical unit spans remain reported. SD can miss/split units: evaluate recoverable subset only.'))

def main():
    lock=(G/'catalog_units.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    tasks=[]
    for n in TICKET:
        groups=defaultdict(list)
        for i,r in enumerate(json.loads((G/n/'fixed/plan.json').read_text())['jobs']):
            groups[r[0]].append((f'T{i:05d}',r))
        tasks.extend((n,ch,ts) for ch,ts in groups.items())
    with ProcessPoolExecutor(12) as pool:
        while tasks:
            ready=[v for v in tasks if (ROOT/'acceptance'/(TICKET[v[0]]+'.json')).is_file() and (G/v[0]/'catalog/freeze.json').is_file()]
            if ready:
                list(pool.map(process,ready));tasks=[x for x in tasks if x not in ready]
            else:time.sleep(15)
    atomic(G/'catalog_units_complete.json',dict(status='all released references processed'))

if __name__=='__main__':main()
