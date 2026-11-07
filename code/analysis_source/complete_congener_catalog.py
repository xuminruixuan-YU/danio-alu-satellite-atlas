"""Reconcile complete scans with reviewed sequence evidence, retaining ambiguity.

Reuses TRF raw tandem evidence; fixed-template evidence remains an independent
method view, not automatic replacement of discovery or the family catalogue.
"""
import argparse,csv,gzip,json,collections,hashlib
from pathlib import Path
import pysam
from regional_catalog import W,ROOT,R,E,OLD,rows,save,table,sha,partition,canon,union
O=W/'congener_catalog';O.mkdir(exist_ok=True)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--catalogue-only',action='store_true');args=ap.parse_args()
    if not args.catalogue_only:
        assert (W/'genome_sd/fixed_SD_checked.json').exists()
        assert (W/'fixed_tandem/computed.json').exists()
    evidence={x['motif']:x for x in rows(W/'congener_identity/motif_evidence.tsv')}
    source=rows(W/'congener_identity/source_records.tsv')
    prior={}
    for x in rows(R/'14_danio_scope_review/adjudication/local_region_adjudication.tsv'):
        if x['supported_family']:prior[canon(x['local_template'])]=x['supported_family']
    # Explicit analyst-reviewed stratum: whole-period GC repeat-library matches
    # with corresponding full-query historical alignment. The minimum whole-
    # query score is 0.7784 in this set; no 90% identity admission rule is used.
    extra={k for k,x in evidence.items() if x['whole_query_best_family_lead']=='CL2' and x['background_top']=='SAT-1_DR' and float(x['background_one_period_span'])>=.8 and x['background_report_cap']=='False'}
    assert len(extra)==312
    save(O/'reviewed_motif_decision.json',dict(explicit_added_GC_motifs=sorted(extra),prior_supported_motifs=len(prior),basis=['TRF tandem organization in actual assembly','whole-comparable-unit historical GC relation','historical library alignment covers most or all of one period','full-query score and local identity have different denominators','divergent AT-related sequences retained unresolved, not called absent'],interpretation='analyst evidence stratum for this fixed collection; not a universal numerical homology threshold',minimum_observed_whole_query_score=min(float(evidence[k]['whole_query_best_score']) for k in extra)))
    regions=[];members=[];interruptions=[];review=[];summary=[];fixed=rows(W/'fixed_tandem/diagnostic_segments.tsv') if not args.catalogue_only else [];methodrows=[];variantrows=[]
    for spec,acc in [('aesculapii','GCA_903798145.1'),('kyathit','GCA_903798195.1')]:
        meta=json.loads((E/'inputs/reference'/acc/'ready.json').read_text());ref=Path(meta['fasta']);lengths={x['seqid']:x['length'] for x in meta['inventory'] if x['nuclear']};labels={x['seqid']:x['molecule'] if x['role']=='assembled-molecule' else 'unplaced' for x in meta['inventory'] if x['nuclear']}
        fr=[]
        with pysam.FastaFile(str(ref)) as fa:
            for x in [s for s in source if s['species']==spec]:
                family='CL1' if x['family']=='AT' else 'CL2';known=prior.get(x['canonical']);added=x['motif'] in extra and family=='CL2'
                supported=known==family or added
                a,b=int(x['start0']),int(x['end0']);seq=fa.fetch(x['seqid'],a,b).upper();assert len(seq)==b-a
                row=dict(species=spec,source_id=x['source_id'],motif=x['motif'],seqid=x['seqid'],start0=a,end0=b,family=family,accepted_family=family if supported else '',status='prior_full_unit_review' if known==family else 'full_unit_GC_historical_and_tandem' if added else 'family_relation_unresolved',period=x['period'],copies=x['TRF_copies'],raw_N_bp=seq.count('N'),raw_non_ACGT_bp=len(seq)-sum(seq.count(c) for c in 'ACGT'),raw_sha256=hashlib.sha256(seq.encode()).hexdigest(),raw_sequence_scope='reference-extracted; TRF boundaries retained, source calls can overlap')
                review.append(row)
                if supported:fr.append(dict(source_id=x['source_id'],seqid=x['seqid'],start0=a,end0=b,family=family))
        for x in rows(R/'14_danio_scope_review/adjudication/retrieval_control_adjudication.tsv'):
            if x['species']==spec:
                # Only the supported repeat component, not the padded window.
                fr.append(dict(source_id=x['region'],seqid=x['seqid'],start0=int(x['component_start0']),end0=int(x['component_end0']),family=x['supported_family']))
        rr,mm,gg=partition(spec,acc,ref,lengths,labels,fr,25000,10000)
        for x in rr:x['census_scope']='complete nuclear scan; supported subset and unresolved family candidates separately retained'
        regions+=rr;members+=mm;interruptions+=gg
        for gap in [10000,25000,50000]:
            pp,pm,pg=partition(spec,acc,ref,lengths,labels,fr,gap,10000)
            vd=O/'variants'/f'{spec}_gap{gap}';vd.mkdir(parents=True,exist_ok=True)
            table(vd/'regions.tsv',pp);table(vd/'membership.tsv',pm)
            for cutoff in [5000,10000,15000,20000,25000,50000,100000,200000]:
                for fam in ['CL1','CL2']:
                    ll=[x for x in pp if x['family']==fam and x['supported_bp']>=cutoff]
                    variantrows.append(dict(species=spec,family=fam,gap_bp=gap,minimum_supported_bp=cutoff,regions=len(ll),supported_bp=sum(x['supported_bp'] for x in ll)))
        for fam in ['CL1','CL2']:
            small=[x for x in rr if x['family']==fam];large=[x for x in small if x['qualifies_large']]
            summary.append(dict(species=spec,assembly=acc,family=fam,supported_regions=len(small),large_regions=len(large),supported_bp=sum(x['supported_bp'] for x in small),large_supported_bp=sum(x['supported_bp'] for x in large),maximum_supported_region_bp=max((x['supported_bp'] for x in small),default=0),unresolved_source_records=sum(x['species']==spec and x['family']==fam and not x['accepted_family'] for x in review),interpretation='observed assembly evidence; non-detection under size rule is not species-level absence'))
            for panel in ['two','three']:
                for ch in lengths:
                    ii=union((int(x['start0']),int(x['end0'])) for x in fixed if x['species']==spec and x['panel']==panel and x['family']==fam and x['seqid']==ch)
                    jj=union((x['start0'],x['end0']) for x in fr if x['seqid']==ch and x['family']==fam)
                    common=sum(max(0,min(b,d)-max(a,c)) for a,b in ii for c,d in jj)
                    if ii or jj:methodrows.append(dict(species=spec,panel=panel,family=fam,seqid=ch,SD_screen_bp=sum(b-a for a,b in ii),reviewed_TRF_bp=sum(b-a for a,b in jj),common_bp=common,SD_only_bp=sum(b-a for a,b in ii)-common,TRF_only_bp=sum(b-a for a,b in jj)-common,scope='native75 tandem screen is a diagnostic stratum; no required TRF-and-SD intersection; SD-only intervals retained for review'))
        # All diagnostic intervals are retained, including native-template-only
        # candidates. Neither mismatch nor unresolved identity is called absence.
    table(O/'regions.tsv',regions);table(O/'membership.tsv',members);table(O/'interruptions.tsv',interruptions)
    table(O/'source_family_review.tsv',review);table(O/'summary.tsv',summary);table(O/'method_coverage.tsv',methodrows);table(O/'threshold_sensitivity.tsv',variantrows)
    save(O/'computed.json',dict(input_scan_receipt=sha(W/'genome_sd/fixed_SD_checked.json') if not args.catalogue_only else None,native_tandem_receipt=sha(W/'fixed_tandem/computed.json') if not args.catalogue_only else None,full_nuclear_TRF=True,full_nuclear_SD=not args.catalogue_only,scientific_review_pending=True,scope='observed supported family catalogue plus explicit review candidates; no claim of exhaustive true satellite membership'))
    print(json.dumps(summary,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
