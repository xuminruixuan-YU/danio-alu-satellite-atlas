"""Separate reference-CpG availability from valid-observation availability."""
from pathlib import Path
import json,re,numpy as np,pandas as pd,pysam,gzip
D=Path(__import__('os').environ.get('DANIO_COMPARISON_DIR',str(Path(__file__).resolve().parents[2]/'data/regional_display')))
def main():
 cfg=json.loads((D/'inputs/new_config.json').read_text());arr=[x['array'] for x in json.loads((D/'inputs/array_windows.json').read_text())['matrices']];fa=pysam.FastaFile(__import__('os').environ.get('DANIO_REFERENCE_FASTA',cfg['reference']));refs={}
 with gzip.open(D/'inputs/reference_CpG_positions.tsv.gz','wt') as f:
  f.write('reference\tposition0\n')
  for ref in {a['seqid'] for a in arr}:
   rows=[a for a in arr if a['seqid']==ref];start=min(int(a['plot_start0']) for a in rows);end=max(int(a['plot_end0']) for a in rows);seq=fa.fetch(ref,start,end).upper();p=np.fromiter((start+m.start() for m in re.finditer('CG',seq)),dtype=np.int64);refs[ref]=p
   for x in p:f.write(f'{ref}\t{x}\n')
 fa.close();df=pd.read_csv(D/'results/windows.tsv.gz',sep='\t');rr=[];br=[]
 for a in arr:
  aid=a['array_id'];p=refs[a['seqid']];domain_end=int(a['plot_end0'])-1;g=df[df.array_id==aid];count=np.searchsorted(p,np.minimum(g.end0.to_numpy(),domain_end))-np.searchsorted(p,g.start0.to_numpy());z=g[['array_id','run','window_bp','step_bp','start0','end0','CpGs']].copy();z['reference_CpGs']=count;z['effective_fraction']=np.where(count>0,z.CpGs/count,np.nan);assert (z.CpGs<=count).all();rr.append(z)
  for name,x,y in [('left',int(a['plot_start0']),int(a['start0'])),('inside',int(a['start0']),int(a['end0'])),('right',int(a['end0']),domain_end)]:br.append(dict(array_id=aid,zone=name,reference_CpGs=int(np.searchsorted(p,y)-np.searchsorted(p,x))))
 pd.concat(rr).to_csv(D/'results/reference_CpG_coverage.tsv.gz',sep='\t',index=False);pd.DataFrame(br).to_csv(D/'results/reference_CpG_background.tsv',sep='\t',index=False)
if __name__=='__main__':main()
