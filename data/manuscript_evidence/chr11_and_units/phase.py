"""Resolve CL2 subtype with raw diagnostic bases; compare rotated coordinates and order."""
from collections import Counter,defaultdict
import json,re,gzip,sys,time
from bisect import bisect_left
from pathlib import Path
import numpy as np,edlib
from Bio import SeqIO
from library import rc,write_tsv
from pilot import ROOT,PILOT,atomic
def rawcall(sequence,A,B,rev):
    s=rc(sequence) if rev else sequence
    if not 82<=len(s)<=102 or set(s)-set('ACGT'):return 'U','outside_length_or_nonACGT'
    d=[j for j in range(len(A)) if A[j]!=B[j]];votes=[]
    for t in [A,B]:
        aln=edlib.align(t,s,mode='NW',task='path')
        g=edlib.getNiceAlignment(aln,t,s)
        idx=-1;vals=[]
        for x,y in zip(g['query_aligned'],g['target_aligned']):
            if x!='-':
                idx+=1
                if idx in d:vals.append(y)
        votes.append('A' if vals==[A[i] for i in d] else 'B' if vals==[B[i] for i in d] else 'U')
    return (votes[0] if votes[0]==votes[1] else 'U','/'.join(votes))
def main():
    d=PILOT/'fixed_comparison';d.mkdir(exist_ok=True)
    r=list(SeqIO.parse(ROOT/'inputs/three.fa','fasta'));A=str(r[1].seq);B=str(r[2].seq)
    out=[];summary=[]
    manifest=json.loads((ROOT.parent/'CL1_CL2_local_handoff_2026-11-03/round3/seed_comparison/input/pilot_manifest.json').read_text())
    for obj in manifest:
        w=obj['id'];off=obj['start0'];s=str(next(SeqIO.parse(ROOT/'inputs/windows'/(w+'.fa'),'fasta')).seq).upper()
        maps={};starts={}
        for mode in ['two','rot92','three','reverse']:
            m=np.zeros(len(s),dtype='int8');score=np.full(len(s),-1,dtype=np.int16);ss=[]
            p=PILOT/'fixed_sd_input'/(w+'.'+mode)/'final_decomposition.tsv'
            for ln,line in enumerate(p.open(),1):
                v=line.split();tag=v[1].strip("'");a=int(v[2]);b=int(v[3])+1;ident=float(v[4])
                if not tag.startswith('CL2') or ident<75:continue
                rv=v[1].endswith("'");lab='U';why='dimer template'
                if mode in ('three','reverse'):
                    lab,why=rawcall(s[a:b],A,B,rv);sd='A' if tag.endswith('1') else 'B'
                    if lab!=sd:lab='U'
                    out.append(dict(window=w,mode=mode,raw_row=ln,start0=off+a,end0=off+b,direction='-' if rv else '+',native=tag,identity=ident,raw_votes=why,concordant_subtype=lab))
                val={'U':3,'A':1,'B':2}[lab];vscore=int(round(100*ident));best=vscore>score[a:b];tie=(vscore==score[a:b])&(m[a:b]!=val)
                m[a:b][tie]=3;m[a:b][best]=val;score[a:b][best]=vscore;ss.append((a,rv))
            maps[mode]=m;starts[mode]=ss
        a,b=maps['three'],maps['reverse'];common=(a>0)&(b>0);known=common&(a<3)&(b<3)
        nearest=[]
        for orient in [True,False]:
            sorted_b=sorted(x for x,rv in starts['rot92'] if rv==orient)
            for x,rv in starts['two']:
                if rv!=orient or not sorted_b:continue
                j=bisect_left(sorted_b,x);opts=sorted_b[max(j-1,0):j+1];dist=min(abs(q-x) for q in opts)
                nearest.append(dist)
        summary.append(dict(window=w,raw_confident_both_bp=int(known.sum()),raw_AB_discordant_bp=int((known&(a!=b)).sum()),unresolved_in_either_on_common_bp=int((common&~known).sum()),shared_CL2_coverage_within_this_diagnostic=int(common.sum()),rotation_nearest_start_classes=dict(Counter('<=5' if x<=5 else '6..70' if x<=70 else '71..105' if x<=105 else '>105' for x in nearest)),rotation_denominator=len(nearest)))
    write_tsv(d/'cl2_raw_diagnostics.tsv.gz',list(out[0]),out);write_tsv(d/'phase_summary.tsv',list(summary[0]),summary)
    atomic(d/'phase_checks.json',dict(diagnostic_positions0=[i for i in range(len(A)) if A[i]!=B[i]],diagnostic_scope='CL2-only blocks >=75 in diagnostic view; raw bases must agree using both template alignments and nominated orientation. Unknowns and phase changes retained. Shared denominator here omits competition with CL1 unlike whole-family summaries.'))
if __name__=='__main__':main()
