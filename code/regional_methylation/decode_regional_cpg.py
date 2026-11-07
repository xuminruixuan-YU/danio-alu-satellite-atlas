"""Decode regional CpGs from preserved primary modBAMs and complete parent records.

Provide all returned alignments in TSV(.gz) with columns read, flag, reference,
start0, AS and CIGAR. Primary modBAM inputs must be disjoint partitions of one
sample. No mapping, download or background job is launched by this decoder.
"""
from pathlib import Path
from collections import Counter,defaultdict
import argparse,csv,gzip,json,sqlite3,tempfile
import pysam
from locus_cpg_coordinates import compile_placement,qualify_cpg
from mod_probability_calls import bounds,call,selfcheck

FIELDS=['C','m','h','unknown','quantization_tie','low_confidence','below_Q10']
CFG={'MAPQ':20,'secondary_cap':100,'score_fraction':0.8,'minimum_baseQ':10,'mod_probability':0.8}
INTERVALS=defaultdict(list)
def inside(ref,c):return any(a<=c and c+1<b for a,b in INTERVALS.get(ref,[]))
def rows(path):
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as f:yield from csv.DictReader(f,delimiter='\t')

class Calls:
    def __init__(self,run):
        self.run=run;self.sites=defaultdict(Counter);self.counts=Counter()
        self.fasta=pysam.FastaFile(CFG['reference'])
    def add(self,a,parent_records,part):
        self.counts[part+'_target_primary']+=1
        if a.mapping_quality<CFG['MAPQ']:
            self.counts[part+'_below_MAPQ']+=1;return
        if any(int(x['flag'])&2048 for x in parent_records) or any('H' in x['CIGAR'] for x in parent_records):
            self.counts[part+'_supplementary_or_hardclip']+=1;return
        secondary=sum(bool(int(x['flag'])&256) and not int(x['flag'])&2048 for x in parent_records)
        if secondary>=CFG['secondary_cap']:
            self.counts[part+'_secondary_cap']+=1;return
        pri=[x for x in parent_records if not int(x['flag'])&(256|2048)]
        assert len(pri)==1 and int(pri[0]['start0'])==a.reference_start and pri[0]['CIGAR']==a.cigarstring
        score=pri[0]['AS']
        if score is None or int(score)<=0:
            self.counts[part+'_missing_score']+=1;return
        high=[x for x in parent_records if x['AS'] is not None and int(x['AS'])>=CFG['score_fraction']*int(score)]
        compiled=[compile_placement(x) for x in high]
        primary=next(p for p,x in zip(compiled,high) if not int(x['flag'])&(256|2048))
        seq=a.get_forward_sequence();qualities=a.get_forward_qualities()
        if not a.has_tag('MM') or not a.has_tag('ML') or qualities is None:
            self.counts[part+'_missing_mod_or_quality']+=1;return
        mods=a.modified_bases_forward
        if mods is None:raise ValueError('MM/ML failed to decode: '+a.query_name)
        forms=[s.split(',')[0] for s in a.get_tag('MM').split(';') if s]
        supported='C+m.' in forms and 'C+h.' in forms and (not a.has_tag('dx') or a.get_tag('dx')!=1)
        if not supported:self.counts[part+'_unsupported_encoding']+=1;return
        mm=dict(mods.get(('C',0,'m'),[]));hh=dict(mods.get(('C',0,'h'),[]))
        qmap=dict(a.get_aligned_pairs(matches_only=True))
        context=self.fasta.fetch(a.reference_name,a.reference_start,a.reference_end).upper()
        self.counts[part+'_eligible_reads']+=1
        pos=seq.find('CG')
        while pos>=0:
            qc=len(seq)-1-pos if a.is_reverse else pos;qg=qc-1 if a.is_reverse else qc+1
            pc,pg=qmap.get(qc),qmap.get(qg)
            if pc is not None and pg is not None and abs(pc-pg)==1:
                c=min(pc,pg);ref=a.reference_name
                if inside(ref,c) and context[c-a.reference_start:c-a.reference_start+2]=='CG':
                    status=qualify_cpg(pos,ref,c,primary,compiled);self.counts[part+'_'+status]+=1
                    if status=='position_concordant':
                        mb,hb=mm.get(pos),hh.get(pos)
                        b=bounds(mb,hb,mb is None and hb is None)
                        klass=call(b,CFG['mod_probability']) if min(qualities[pos],qualities[pos+1])>=CFG['minimum_baseQ'] else 'below_Q10'
                        self.sites[ref,c][klass]+=1
            pos=seq.find('CG',pos+1)
    def output(self,d):
        fields=['run','reference','reference_C_position0',*FIELDS,'valid_calls','observations','five_mC_call_fraction']
        with gzip.open(d/'CpGs.tsv.gz','wt') as f:
            writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t');writer.writeheader()
            for (ref,c),v in sorted(self.sites.items()):
                valid=sum(v[k] for k in ['C','m','h'])
                writer.writerow(dict(run=self.run,reference=ref,reference_C_position0=c,**{k:v[k] for k in FIELDS},
                                     valid_calls=valid,observations=sum(v.values()),five_mC_call_fraction=v['m']/valid if valid else None))
        self.fasta.close()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--regions',type=Path,required=True,help='JSON array with seqid, plot_start0 and plot_end0')
    p.add_argument('--run',required=True)
    p.add_argument('--primary-bam',type=Path,action='append',required=True)
    p.add_argument('--parents',type=Path,action='append',required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    if (args.output/'CpGs.tsv.gz').exists():raise FileExistsError('Refusing to overwrite prior counts')
    CFG['reference']=str(args.reference)
    for r in json.loads(args.regions.read_text()):INTERVALS[r['seqid']].append((int(r['plot_start0']),int(r['plot_end0'])))
    selfcheck();calls=Calls(args.run)
    with tempfile.TemporaryDirectory(dir=args.output,prefix='parents_') as tmp:
        con=sqlite3.connect(str(Path(tmp)/'parents.sqlite'))
        con.execute('create table alignments (read text, record text)')
        for path in args.parents:
            for r in rows(path):
                r['flag']=int(r['flag']);r['start0']=int(r['start0']);r['AS']=int(r['AS']) if r['AS'] else None
                con.execute('insert into alignments values (?,?)',(r['read'],json.dumps(r)))
        con.execute('create index read_index on alignments(read)');con.execute('create table used(read text primary key)');con.commit()
        for path in args.primary_bam:
            with pysam.AlignmentFile(path,'rb') as bam:
                for a in bam:
                    if a.is_secondary or a.is_supplementary or a.is_unmapped:continue
                    if not any(a.reference_start<end and a.reference_end>start for start,end in INTERVALS.get(a.reference_name,[])):continue
                    con.execute('insert into used values (?)',(a.query_name,))
                    parent=[json.loads(x[0]) for x in con.execute('select record from alignments where read=?',(a.query_name,))]
                    if not parent:raise ValueError('Missing complete parent records: '+a.query_name)
                    calls.add(a,parent,'input')
        calls.output(args.output)
        (args.output/'counts.json').write_text(json.dumps(dict(calls.counts),indent=2)+'\n')
        con.close()

if __name__=='__main__':main()
