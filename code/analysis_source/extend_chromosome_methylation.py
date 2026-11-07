"""Decode only independently anchored existing reads, including real flanks.

This is conditional coverage of the satellite-read subset, not an unbiased
whole-genome flank methylome. Each CpG must agree across all high-score parents.
The original in-fragment observations are reproduced exactly before acceptance.
"""
import bisect,collections,csv,gzip,json,sys,hashlib
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from zoneinfo import ZoneInfo
import pysam
W=Path(__file__).resolve().parents[1];ROOT=W.parents[1]
R=ROOT/'output/danio_five_part_followup_2026-11-04'
sys.path.insert(0,str(R/'code'))
from locus_cpg_coordinates import compile_placement,qualify_cpg
from mod_probability_calls import bounds,call
O=W/'chromosome_methylation_tracks';O.mkdir(exist_ok=True)
CLASSES=['C','m','h','unknown','quantization_tie','low_confidence','below_Q10']
def rows(p):
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:return list(csv.DictReader(f,delimiter='\t'))
def save(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n')
def table(p,rr,fields):
    with (gzip.open(p,'wt') if str(p).endswith('.gz') else open(p,'w')) as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(rr)
def worker(run):
    d=O/run;d.mkdir(exist_ok=True)
    if (d/'checked.json').exists():return json.loads((d/'checked.json').read_text())
    anchors={x['read']:x for x in rows(R/'11_locus_methylation/anchor_review'/run/'validated_main_anchors.tsv')}
    parents=collections.defaultdict(list)
    for x in rows(R/'11_locus_methylation/anchor_preparation'/run/'parent_high_score_placements.tsv'):
        if x['read'] in anchors:parents[x['read']].append(x)
    cp={};primary={}
    for name,pp in parents.items():
        assert len(pp)==int(anchors[name]['high_score_placements'])
        cp[name]=[compile_placement(x) for x in pp]
        ii=[i for i,x in enumerate(pp) if not int(x['flag'])&(256|2048|4)];assert len(ii)==1
        primary[name]=cp[name][ii[0]]
    large=rows(W/'chromosome_structure/windows.tsv')
    byref=collections.defaultdict(list)
    for x in large:byref[x['seqid']].append(x)
    fragments=collections.defaultdict(list)
    for x in rows(R/'05_catalog_review/fragments.tsv'):fragments[x['seqid']].append(x)
    starts={}
    for ref,ff in fragments.items():ff.sort(key=lambda x:int(x['start0']));starts[ref]=[int(x['start0']) for x in ff]
    # Exact comparison includes low quality/low confidence calls, not only 5mC.
    old={}
    for x in rows(R/'11_locus_methylation/per_run'/run/'position_qualified_observations.tsv.gz'):
        key=(x['read'],x['reference'],int(x['reference_C_position0']))
        old[key]=(x['primary_call'],x['m_byte'],x['h_byte'],int(x['read_C_position0']))
    reproduced={};seen=set();observations=[];sites=collections.defaultdict(collections.Counter)
    diagnostics=collections.Counter();info=json.loads((R/'10_family_methylation/per_run'/run/'observations_input.json').read_text())
    source=json.loads((ROOT/'output/danio_cl1_cl2_execution_2026-11-03/external/reference/inputs.json').read_text())['source']
    with pysam.FastaFile(source) as ref:
        for part in ['prefix','new']:
            with pysam.AlignmentFile(info[part+'_modBAM'],'rb') as bam:
                for a in bam:
                    name=a.query_name
                    if name not in anchors:continue
                    assert anchors[name]['part']==part and name not in seen;seen.add(name)
                    assert not a.is_unmapped and not a.is_secondary and not a.is_supplementary
                    pp=primary[name];assert pp.reverse==a.is_reverse and pp.reference==a.reference_name
                    seq=a.get_forward_sequence();q=a.get_forward_qualities();mods=a.modified_bases_forward or {}
                    forms=[s.split(',')[0] for s in a.get_tag('MM').split(';') if s]
                    supported='C+m.' in forms and 'C+h.' in forms and (not a.has_tag('dx') or a.get_tag('dx')!=1)
                    mm=dict(mods.get(('C',0,'m'),[]));hh=dict(mods.get(('C',0,'h'),[]))
                    qmap=dict(a.get_aligned_pairs(matches_only=True));context=ref.fetch(a.reference_name,a.reference_start,a.reference_end).upper()
                    pos=seq.find('CG')
                    while pos>=0:
                        qc=len(seq)-1-pos if a.is_reverse else pos;qg=qc-1 if a.is_reverse else qc+1
                        pc,pg=qmap.get(qc),qmap.get(qg)
                        if pc is not None and pg is not None and abs(pc-pg)==1:
                            c=min(pc,pg);reference=a.reference_name
                            if context[c-a.reference_start:c-a.reference_start+2]=='CG':
                                status=qualify_cpg(pos,reference,c,pp,cp[name]);diagnostics[status]+=1
                                if status=='position_concordant':
                                    mb,hb=mm.get(pos),hh.get(pos);b=bounds(mb,hb,mb is None and hb is None) if supported else None
                                    klass=call(b,.8) if min(q[pos],q[pos+1])>=10 else 'below_Q10'
                                    key=(name,reference,c);observed=(klass,'' if mb is None else str(mb),'' if hb is None else str(hb),pos)
                                    if key in old:assert old[key]==observed,(run,key,old[key],observed);reproduced[key]=observed
                                    i=bisect.bisect_right(starts.get(reference,[]),c)-1;frag=fragments[reference][i] if i>=0 else None
                                    if frag and not int(frag['start0'])<=c<c+1<int(frag['end0']):frag=None
                                    if frag and frag['family']==anchors[name]['family']:assert key in old,(run,key,'new unexpected original-fragment CpG')
                                    for region in byref[reference]:
                                        if int(region['plot_start0'])<=c<c+1<int(region['plot_end0']):
                                            zone='left_flank' if c<int(region['start0']) else 'right_flank' if c>=int(region['end0']) else ('satellite_'+frag['family']) if frag else 'inter_array_or_interruption'
                                            row=dict(run=run,array_id=region['array_id'],family=region['family'],read_family=anchors[name]['family'],reference=reference,reference_C_position0=c,read=name,read_C_position0=pos,zone=zone,call=klass,C_baseQ=q[pos],G_baseQ=q[pos+1],m_byte=observed[1],h_byte=observed[2],encoding_supported=supported)
                                            observations.append(row);sites[region['array_id'],region['family'],reference,c,zone][klass]+=1
                        pos=seq.find('CG',pos+1)
    assert seen==set(anchors),(run,len(seen),len(anchors))
    assert reproduced==old,(run,len(reproduced),len(old),list(set(old)-set(reproduced))[:3])
    out=[]
    for (aid,fam,ref,c,zone),counts in sorted(sites.items()):
        valid=sum(counts[k] for k in ['C','m','h'])
        out.append(dict(run=run,array_id=aid,family=fam,reference=ref,reference_C_position0=c,zone=zone,**{k:counts[k] for k in CLASSES},valid_calls=valid,observations=sum(counts.values()),five_mC_call_fraction=counts['m']/valid if valid else None))
    table(d/'observations.tsv.gz',observations,['run','array_id','family','read_family','reference','reference_C_position0','read','read_C_position0','zone','call','C_baseQ','G_baseQ','m_byte','h_byte','encoding_supported'])
    table(d/'CpGs.tsv.gz',out,['run','array_id','family','reference','reference_C_position0','zone',*CLASSES,'valid_calls','observations','five_mC_call_fraction'])
    receipt=dict(run=run,at_hkt=datetime.now(ZoneInfo('Asia/Hong_Kong')).isoformat(),anchored_reads=len(seen),exactly_reproduced_original_observations=len(reproduced),plotted_CpG_rows=len(out),observations=len(observations),zone_counts=dict(collections.Counter(x['zone'] for x in out)),coordinate_diagnostics=dict(diagnostics),scope='observed positions on anchored satellite-read subset; no claim of complete or unbiased flanking methylation coverage')
    save(d/'checked.json',receipt);return receipt
if __name__=='__main__':
    runs=[f'ERR167823{i}' for i in range(14,20)]
    with ProcessPoolExecutor(6) as p:
        results=list(p.map(worker,runs))
    save(O/'checked.json',dict(runs=results,baseline_fragment_calls_exactly_reproduced=True,scope='conditional localized tracks; gaps not imputed'))
    print(json.dumps(results,ensure_ascii=False,indent=2))
