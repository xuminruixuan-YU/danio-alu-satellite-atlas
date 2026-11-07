"""Independent CIGAR coordinate audit using pysam aligned-pair enumeration.
Does not import the original custom coordinate helper.
"""
from pathlib import Path
import json,collections
import pandas as pd
import pysam
D=Path(__file__).resolve().parents[2];C=D/'data/regional_methylation/coordinate_example'
observed=pd.read_csv(C/'sensitivity_CpG_placements.tsv',sep='\t')
parents=pd.read_csv(C/'sensitivity_alignments.tsv',sep='\t')
sequence={}
with (C/'reads.fastq').open() as f:
    while True:
        name=f.readline().strip()
        if not name:break
        sequence[name[1:]]=f.readline().strip();f.readline();f.readline()
refs=sorted(parents.reference.unique());header=pysam.AlignmentHeader.from_dict({'SQ':[{'SN':r,'LN':200000000} for r in refs]})
maps={}
for (method,name),g in parents.groupby(['method','read']):
    pp=[]
    for r in g[g.score_ratio>=.8].itertuples():
        a=pysam.AlignedSegment(header);a.query_name=name;a.query_sequence=sequence[name];a.flag=int(r.flag);a.reference_id=refs.index(r.reference);a.reference_start=int(r.start0);a.cigarstring=r.CIGAR
        qlen=a.infer_query_length();assert qlen==len(sequence[name])
        pp.append((r.reference,a.is_reverse,qlen,dict(a.get_aligned_pairs(matches_only=True))))
    maps[method,name]=pp
counts=collections.Counter()
for x in observed.itertuples():
    matches=[]
    for ref,rev,length,qmap in maps[x.method,x.read]:
        qc=length-1-int(x.read_C0) if rev else int(x.read_C0);qg=qc-1 if rev else qc+1
        pc,pg=qmap.get(qc),qmap.get(qg)
        matches.append(pc is not None and pg is not None and abs(pc-pg)==1 and ref==x.original_reference and min(pc,pg)==x.original_reference_C0)
    assert bool(matches) and all(matches)==x.exact_original_CpG,(x.method,x.read,x.read_C0)
    counts[x.method]+=1
result=dict(status='all_coordinates_independently_reproduced',method='pysam aligned pairs, no custom coordinate helper',observations_checked=dict(counts))
print(json.dumps(result,indent=2))
