from pathlib import Path
import json,numpy as np,pandas as pd
from concurrent.futures import ProcessPoolExecutor
D=Path(__import__('os').environ.get('DANIO_COMPARISON_DIR',str(Path(__file__).resolve().parents[2]/'data/regional_display')))
ARR={x['array']['array_id']:x['array'] for x in json.loads((D/'inputs/array_windows.json').read_text())['matrices']}
def summarize(p,a,b):
 i,j=np.searchsorted(p,[a,b]);p=p[i:j]
 if len(p):
  starts=np.r_[a,p+1];ends=np.r_[p,b];k=np.argmax(ends-starts);s,e=int(starts[k]),int(ends[k]);q=np.diff(np.searchsorted(p,np.linspace(a,b,5)))
 else:s,e=a,b;q=np.zeros(4,int)
 return dict(gap_start0=s,gap_end0=e,gap_bp=e-s,q1=int(q[0]),q2=int(q[1]),q3=int(q[2]),q4=int(q[3]),edge_flag=e-s>(b-a)/2 or q[1]+q[2]==0)
def work(run):
 dest=D/'per_run'/run;obs=pd.read_csv(dest/'valid_observations.tsv.gz',sep='\t',usecols=['read','reference','position0']);models=pd.read_csv(dest/'read_models.tsv.gz',sep='\t');opaque=set(models.loc[models.model_metadata.str.contains('modbase_model_version_id'),'read']);obs['known']=~obs.read.isin(opaque);refs={ref:(np.unique(g.position0),np.unique(g.loc[g.known,'position0'])) for ref,g in obs.groupby('reference')};w=pd.read_csv(dest/'windows.tsv.gz',sep='\t');rows=[]
 for z in w.itertuples():
  p,kp=refs.get(ARR[z.array_id]['seqid'],(np.array([],int),np.array([],int)));a,b=int(z.start0),int(z.end0);s=summarize(p,a,b);k=summarize(kp,a,b);assert s['gap_bp']==z.max_gap
  rows.append(dict(array_id=z.array_id,run=run,window_bp=z.window_bp,step_bp=z.step_bp,start0=a,end0=b,**s,**{'known_'+x:y for x,y in k.items()}))
 pd.DataFrame(rows).to_csv(dest/'coverage_detail.tsv.gz',sep='\t',index=False)
def main():
 with ProcessPoolExecutor(max_workers=6) as pool:list(pool.map(work,[f'ERR167823{i}' for i in range(14,20)]))
 pd.concat([pd.read_csv(p,sep='\t') for p in (D/'per_run').glob('ERR*/coverage_detail.tsv.gz')]).to_csv(D/'results/coverage_detail.tsv.gz',sep='\t',index=False)
if __name__=='__main__':main()
