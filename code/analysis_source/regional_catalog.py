"""Explicit large-region grouping, retaining source fragments and interruptions.

The working rule is an analyst choice supported by a recorded sensitivity grid,
not an externally validated biological boundary or a fixed chromosome count.
Congener family-identity leads not already reviewed remain separate.
"""
import os
import csv,gzip,json,hashlib,collections
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import pysam

W=Path(__file__).resolve().parents[1];ROOT=W.parents[1]
R=ROOT / os.environ.get('DANIO_ANNOTATION_ROOT', 'inputs/annotation')
E=ROOT / os.environ.get('DANIO_ASSEMBLY_RUN_ROOT', 'inputs/assembly_runs')
OLD=ROOT / os.environ.get('DANIO_REPEAT_EVIDENCE_ROOT', 'inputs/repeat_evidence')
O=W/'array_catalog';O.mkdir(parents=True,exist_ok=True)
def rows(p):
    op=gzip.open if str(p).endswith('.gz') else open
    with op(p,'rt') as f:return list(csv.DictReader(f,delimiter='\t'))
def save(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n')
def table(p,rr,fields=None):
    with p.open('w') as f:
        z=csv.DictWriter(f,fields or list(rr[0]),delimiter='\t');z.writeheader();z.writerows(rr)
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def union(iv):
    out=[]
    for a,b in sorted(iv):
        if out and a<=out[-1][1]:out[-1][1]=max(out[-1][1],b)
        else:out.append([a,b])
    return out
def canon(s):
    s=s.upper();rc=s.translate(str.maketrans('ACGT','TGCA'))[::-1]
    return min(q[i:]+q[:i] for q in [s,rc] for i in range(len(q)))
def inputs():
    fr=rows(R/'05_catalog_review/fragments.tsv')
    meta=json.loads((E/'external/reference/inputs.json').read_text())
    # Assigned molecule labels are read from the actual accessioned report.
    report=R/'02_order_resources/TU_version_audit/GCF_049306965.1_GRCz12tu_assembly_report.txt'
    labels={}
    for line in report.read_text().splitlines():
        if not line.startswith('#') and line.strip():
            v=line.split('\t');labels[v[6]]=v[2] if v[1]=='assembled-molecule' else v[0]
    tu=[dict(source_id=x['fragment'],seqid=x['seqid'],start0=int(x['start0']),end0=int(x['end0']),family=x['family'],tier='reviewed_family_tandem',source='05_catalog_review/fragments.tsv') for x in fr]
    yield 'TU','GCF_049306965.1',Path(meta['source']),{k:int(v) for k,v in meta['available_sequences']},labels,tu,[]
    old_review=rows(R/'14_danio_scope_review/adjudication/local_region_adjudication.tsv')
    known=collections.defaultdict(set)
    for x in old_review:
        if x['supported_family'] in ['CL1','CL2']:known[canon(x['local_template'])].add(x['supported_family'])
    for name,acc in [('aesculapii','GCA_903798145.1'),('kyathit','GCA_903798195.1')]:
        meta=json.loads((E/'inputs/reference'/acc/'ready.json').read_text())
        lengths={x['seqid']:x['length'] for x in meta['inventory'] if x['nuclear']}
        labels={x['seqid']:(x['molecule'] if x['role']=='assembled-molecule' else 'unplaced') for x in meta['inventory'] if x['nuclear']}
        accepted=[];unresolved=[]
        for i,x in enumerate(rows(OLD/'results'/name/'family_TRF_candidates.tsv.gz')):
            fam='CL1' if x['family']=='AT' else 'CL2';pat=canon(x['pattern'])
            supported=known.get(pat,set())=={fam}
            y=dict(source_id=f'{name}.TRF.{i+1:05d}',seqid=x['seqid'],start0=int(x['start0']),end0=int(x['end0']),family=fam,
                tier='reviewed_motif_family_tandem' if supported else 'family_relation_review_pending',
                source=f"archive:results/{name}/family_TRF_candidates.tsv.gz",motif=pat,
                priority_score=x['score'],family_margin=x['margin'])
            (accepted if supported else unresolved).append(y)
        yield name,acc,Path(meta['fasta']),lengths,labels,accepted,unresolved

def partition(species,acc,ref,lengths,labels,fr,maxgap,minbp):
    gaps=[];out=[];membership=[]
    fa=pysam.FastaFile(str(ref))
    for family in ['CL1','CL2']:
        for seqid in sorted({x['seqid'] for x in fr}):
            # Overlapping source calls form elementary supported tracts first.
            tracts=union((x['start0'],x['end0']) for x in fr if x['family']==family and x['seqid']==seqid)
            opposite=union((x['start0'],x['end0']) for x in fr if x['family']!=family and x['seqid']==seqid)
            groups=[];current=[]
            for a,b in tracts:
                if current:
                    left=current[-1][1];seq=fa.fetch(seqid,left,a).upper();opp=sum(max(0,min(e,a)-max(s,left)) for s,e in opposite)
                    can_join=a-left<=maxgap and seq.count('N')==0 and opp==0
                    if a-left<=maxgap:
                        gaps.append(dict(species=species,family=family,seqid=seqid,start0=left,end0=a,span_bp=a-left,N_bp=seq.count('N'),non_ACGT_bp=len(seq)-sum(seq.count(c) for c in 'ACGT'),opposite_family_bp=opp,joined=can_join,sequence_sha256=hashlib.sha256(seq.encode()).hexdigest()))
                    if not can_join:groups.append(current);current=[]
                current.append((a,b))
            if current:groups.append(current)
            for i,iv in enumerate(groups,1):
                a,b=iv[0][0],iv[-1][1];covered=sum(e-s for s,e in iv);seq=fa.fetch(seqid,a,b).upper()
                ident=f'{species}.{seqid}.{family}.R{i:03d}'
                ff=[x for x in fr if x['family']==family and x['seqid']==seqid and x['start0']>=a and x['end0']<=b]
                selected=covered>=minbp
                out.append(dict(array_id=ident,species=species,assembly=acc,seqid=seqid,chromosome=labels.get(seqid,'unplaced'),family=family,
                    start0=a,end0=b,span_bp=b-a,supported_bp=covered,interruption_bp=b-a-covered,support_fraction=covered/(b-a),
                    source_records=len(ff),supported_tracts=len(iv),N_bp=seq.count('N'),non_ACGT_bp=len(seq)-sum(seq.count(c) for c in 'ACGT'),
                    qualifies_large=selected,minimum_supported_bp=minbp,maximum_interruption_bp=maxgap,
                    plot_start0=max(0,a-200000),plot_end0=min(lengths[seqid],b+200000),
                    boundary_scope='operational region; observed source evidence; not a validated biological endpoint',
                    census_scope='complete TU source catalogue' if species=='TU' else 'reviewed motif subset; full SD and broader motif adjudication pending'))
                membership.extend(dict(array_id=ident,source_id=x['source_id'],species=species,family=family,seqid=seqid,start0=x['start0'],end0=x['end0'],qualifies_large=selected) for x in ff)
    fa.close();return out,membership,gaps

def main():
    resources=list(inputs());primary=[];members=[];gaprows=[];sensitivity=[];summary=[];unresolved=[]
    for species,acc,ref,lengths,labels,fr,pending in resources:
        rr,mm,gg=partition(species,acc,ref,lengths,labels,fr,25000,50000)
        primary.extend(rr);members.extend(mm);gaprows.extend(gg)
        unresolved.extend(dict(species=species,**x) for x in pending)
        for gap in [10000,25000,50000]:
            if gap==25000:part,part_members=rr,mm
            else:part,part_members,_=partition(species,acc,ref,lengths,labels,fr,gap,50000)
            variant=O/'variants'/f'{species}_gap{gap}';variant.mkdir(parents=True,exist_ok=True)
            table(variant/'regions.tsv',part);table(variant/'membership.tsv',part_members)
            for cutoff in [10000,25000,50000,100000,200000]:
                for fam in ['CL1','CL2']:
                    vv=[x for x in part if x['family']==fam and x['supported_bp']>=cutoff]
                    sensitivity.append(dict(species=species,family=fam,maximum_gap_bp=gap,minimum_supported_bp=cutoff,regions=len(vv),placed_regions=sum(x['chromosome']!='unplaced' for x in vv),supported_bp=sum(x['supported_bp'] for x in vv),span_bp=sum(x['span_bp'] for x in vv),scope='TU complete; congener reviewed-motif subset pending whole-SD review'))
        for fam in ['CL1','CL2']:
            allfam=[x for x in rr if x['family']==fam];large=[x for x in allfam if x['qualifies_large']]
            summary.append(dict(species=species,assembly=acc,family=fam,source_records=sum(x['family']==fam for x in fr),all_supported_regions=len(allfam),large_regions=len(large),large_chromosomes=len({x['seqid'] for x in large if x['chromosome']!='unplaced'}),small_regions=len(allfam)-len(large),large_supported_bp=sum(x['supported_bp'] for x in large),large_span_bp=sum(x['span_bp'] for x in large),unresolved_source_records=sum(x['family']==fam for x in pending)))
    table(O/'regions.tsv',primary);table(O/'membership.tsv',members);table(O/'interruptions.tsv',gaprows)
    table(O/'summary.tsv',summary);table(O/'threshold_sensitivity.tsv',sensitivity)
    table(O/'unreviewed_family_leads.tsv',unresolved)
    save(O/'working_rule.json',dict(at_hkt=datetime.now(ZoneInfo('Asia/Hong_Kong')).isoformat(),
        stage='TU working regional rule; congener catalogue pending whole-SD and identity reconciliation',
        minimum_supported_bp=50000,maximum_interruption_bp=25000,unknown_N_or_opposite_family_interruptions='do not bridge',
        rationale='lowest tested50..200kb size range with stable TU large-region counts; gap25/50kb stability and source interruption quality checked; analyst operational rule, not a published universal array boundary',
        secondary_thresholds_bp=[10000,25000,100000,200000],secondary_max_gaps_bp=[10000,50000],
        family_coverage_excludes_interruptions=True,counts_not_forced=True,
        limitations=['congener short regions retained, not biological absence','a large region may contain multiple tandem tracts','family-sequence evidence is not centromere activity'],
        files={p.name:sha(p) for p in O.glob('*.tsv')},scientific_final_acceptance=False))
    print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
