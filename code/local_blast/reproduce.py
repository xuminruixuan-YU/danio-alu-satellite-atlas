"""Reproduce the fixed four-assembly BLAST search or verify its deposited summaries."""
from pathlib import Path
from collections import defaultdict
import argparse,csv,gzip,hashlib,json,re,subprocess,tempfile

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/local_assembly_blast'
COLUMNS='qseqid sseqid pident length mismatch gaps qstart qend sstart send evalue bitscore qlen slen sstrand'.split()

def summarize(path):
    groups=defaultdict(list)
    with path.open() as stream:
        for r in csv.DictReader(stream,fieldnames=COLUMNS,delimiter='\t'):
            r['query_covered_bp']=abs(int(r['qend'])-int(r['qstart']))+1
            groups[r['qseqid']].append(r)
    return {q:{'hsps':len(rs),'subject_sequences':len({r['sseqid']for r in rs}),
               'largest_query_coverage_HSP':max(rs,key=lambda r:(r['query_covered_bp'],float(r['bitscore']))),
               'best_evalue_HSP':min(rs,key=lambda r:float(r['evalue']))}for q,rs in groups.items()}

def verify():
    for record in json.loads((DATA/'assemblies.json').read_text()):
        acc=record['assembly'];folder=DATA/acc
        assert summarize(folder/'matches.tsv')==json.loads((folder/'summary.json').read_text())['summary'],acc
        with (folder/'sequence_inventory.tsv').open()as f:rows=list(csv.DictReader(f,delimiter='\t'))
        nuclear=[r for r in rows if r['nuclear']=='True']
        assert len(nuclear)==record['nuclear_records'] and sum(int(r['bp'])for r in nuclear)==record['nuclear_bp'],acc
        print(acc, 'deposited HSPs and sequence inventory verified')
    with gzip.open(DATA/'rasbora_candidate_distances.tsv.gz','rt')as f:rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==1722 and len({r['canonical_id']for r in rows})==861
    assert all(int(r['passes_d_le_0p30'])==0 and 10*int(r['edit_distance'])>3*int(r['max_length'])for r in rows)
    print('Rasbora: 861 candidate patterns; no full-unit match passes normalized distance <= 0.30')

def run_blast(accession,fasta,output,threads):
    record=next(r for r in json.loads((DATA/'assemblies.json').read_text())if r['assembly']==accession)
    with fasta.open('rb')as f:assert hashlib.file_digest(f,'md5').hexdigest()==record['md5'],'Exact source FASTA checksum differs'
    output.mkdir(parents=True,exist_ok=False)
    with tempfile.TemporaryDirectory(prefix='satellite-blast-')as tmp:
        temp=Path(tmp);nuclear=temp/'nuclear.fna';db=temp/'assembly'
        with gzip.open(fasta,'rt')as source,nuclear.open('w')as dest:
            keep=False
            for line in source:
                if line.startswith('>'):keep=not re.search(r'\bmitochondri(?:on|al)\b',line,re.I)
                if keep:dest.write(line)
        subprocess.run(['makeblastdb','-in',str(nuclear),'-dbtype','nucl','-parse_seqids','-blastdb_version','4','-out',str(db)],check=True)
        cmd=['blastn','-task','blastn','-query',str(DATA/'queries.fasta'),'-db',str(db),'-word_size','11','-evalue','10','-reward','2','-penalty','-3','-gapopen','5','-gapextend','2','-max_target_seqs','5000','-dust','no','-soft_masking','false','-num_threads',str(threads),'-outfmt','6 '+' '.join(COLUMNS),'-out',str(output/'matches.tsv'),'-max_hsps','3']
        subprocess.run(cmd,check=True)
    (output/'summary.json').write_text(json.dumps(summarize(output/'matches.tsv'),indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--accession');p.add_argument('--fasta',type=Path);p.add_argument('--output',type=Path);p.add_argument('--threads',type=int,default=4);a=p.parse_args()
    if a.accession:
        if not(a.fasta and a.output):p.error('--accession requires --fasta and a new --output directory')
        run_blast(a.accession,a.fasta,a.output,a.threads)
    else:verify()
