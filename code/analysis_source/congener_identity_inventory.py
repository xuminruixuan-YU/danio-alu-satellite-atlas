"""Full-genome TRF candidate evidence, beyond the earlier sampled windows.

All unlabelled discovery motifs remain upstream. Prioritization is not a
homology threshold; evidence annotation does not automatically admit arrays.
"""
import csv,gzip,json,hashlib,sys,collections,shutil
from pathlib import Path
from Bio import SeqIO
W=Path(__file__).resolve().parents[1];ROOT=W.parents[1]
R=ROOT/'output/danio_five_part_followup_2026-11-04';OLD=ROOT/'output/danio_satellite_review_2026-11-04'
sys.path.insert(0,str(R/'code'))
import order_candidate_identity as engine
from regional_catalog import canon
O=W/'congener_identity';O.mkdir(exist_ok=True);engine.O=O
def rows(p):
    with (gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)) as f:return list(csv.DictReader(f,delimiter='\t'))
def table(p,rr):
    with open(p,'w') as f:
        z=csv.DictWriter(f,list(rr[0]),delimiter='\t');z.writeheader();z.writerows(rr)
def main():
    seqs={};records=[]
    for species in ['aesculapii','kyathit']:
        for i,x in enumerate(rows(OLD/'results'/species/'family_TRF_candidates.tsv.gz')):
            s=canon(x['pattern']);mid=engine.mid(s);seqs[mid]=s
            records.append(dict(species=species,source_id=f'{species}.TRF.{i+1:05d}',motif=mid,canonical=s,**x))
    table(O/'source_records.tsv',records)
    with open(O/'candidate_motifs.fa','w') as f:
        for k,s in sorted(seqs.items()):f.write(f'>{k}\n{s}\n')
    # Reuse only identical motif and complete archived comparisons. Copy to new
    # evidence tree; never modify prior receipts or add new family labels there.
    reused=0;provenance=[]
    for base in [R/'14_order_survey/identity_evidence/motifs',R/'14_danio_scope_review/identity_evidence/motifs',R/'14_danio_scope_review/pilot_linear_evidence/motifs']:
        for p in base.glob('*/checked.json'):
            k=p.parent.name
            if k not in seqs or (O/'motifs'/k/'checked.json').exists():continue
            v=json.loads(p.read_text())
            if v.get('sequence')!=seqs[k]:continue
            names=['checked.json','whole_query_historical.tsv','repeat_background.tsv']
            if not all((p.parent/n).exists() for n in names):continue
            d=O/'motifs'/k;d.mkdir(parents=True,exist_ok=True)
            for n in names:shutil.copy2(p.parent/n,d/n)
            reused+=1;provenance.append(dict(motif=k,source=str(p.parent)))
    if provenance:table(O/'reused_evidence.tsv',provenance)
    engine.annotate(seqs,'whole_nuclear_TRF_candidates')
    evidence=[]
    for k,s in sorted(seqs.items()):
        v=json.loads((O/'motifs'/k/'checked.json').read_text())
        evidence.append(dict(motif=k,sequence=s,**{a:b for a,b in v.items() if a not in ['motif','sequence']}))
    table(O/'motif_evidence.tsv',evidence)
    (O/'computed.json').write_text(json.dumps(dict(unique_motifs=len(seqs),source_records=len(records),reused_motifs=reused,scope='whole-genome prioritized TRF candidates; homology interpretation and SD reconciliation pending; all lower-priority motifs retained in discovery archive'),indent=2)+'\n')
    print('evidence',len(seqs),'motifs',len(records),'records')
if __name__=='__main__':main()
