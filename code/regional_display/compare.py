from pathlib import Path
import json,gzip,csv,collections,hashlib,sys
import os
BASELINE_ONLY=os.environ.get("BASELINE_ONLY")=="1"
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
D=Path(__import__('os').environ.get('DANIO_COMPARISON_DIR',str(Path(__file__).resolve().parents[2]/'data/regional_display')));ROOT=D.parents[1]
ARR=[x['array'] for x in json.loads((D/'inputs/array_windows.json').read_text())['matrices']]
REP=json.loads((D/'inputs/figure_manifest.json').read_text())['main_examples']
SCHEMES=[(10000,10000),(10000,5000),(20000,5000),(50000,10000),(100000,20000)]
def work(run):
 dest=D/'per_run'/run
 assert json.loads((dest/'checked.json').read_text())['all_site_class_counts_identical']
 df=pd.read_csv(dest/'valid_observations.tsv.gz',sep='\t');assert not df.duplicated(['read','reference','position0']).any()
 ids,unique=pd.factorize(df.read);df['id']=ids
 models=pd.read_csv(dest/'read_models.tsv.gz',sep='\t');unknown=set(models.loc[models.model_metadata.str.contains('modbase_model_version_id'),'read']);df['opaque']=df.read.isin(unknown)
 refs={};known_refs={}
 for ref,g in df.groupby('reference',sort=False):
  g=g.sort_values('position0');pos=g.position0.to_numpy();p,first,n=np.unique(pos,return_index=True,return_counts=True);m=np.add.reduceat((g.call=='m').to_numpy().astype(int),first)
  refs[ref]=(pos,g.id.to_numpy(),g.opaque.to_numpy(),p,n,m,first)
  kg=g[~g.opaque];kp,ki,kn=np.unique(kg.position0.to_numpy(),return_index=True,return_counts=True);km=np.add.reduceat((kg.call=='m').to_numpy().astype(int),ki) if len(ki) else np.array([])
  known_refs[ref]=(kg.position0.to_numpy(),kg.id.to_numpy(),kp,kn,km)
 def stats(ref,a,b):
  z=dict(CpGs=0,observations=0,reads=0,mean=np.nan,weighted=np.nan,span=0,max_gap=b-a,central_CpGs=0,opaque_fraction=np.nan,edge_flag=True,known_CpGs=0,known_reads=0,known_mean=np.nan,min_CpG_depth=0,median_CpG_depth=np.nan,p90_CpG_depth=np.nan,max_CpG_depth=0,max_read_observation_fraction=np.nan)
  if ref not in refs:return z
  pos,ids,opaque,p,n,m,first=refs[ref];i,j=np.searchsorted(p,[a,b]);lo,hi=np.searchsorted(pos,[a,b])
  if i==j:return z
  pp=p[i:j];nn=n[i:j];mm=m[i:j];gap=int(max([pp[0]-a,b-pp[-1]-1,*list(np.diff(pp)-1)]));central=int(np.searchsorted(pp,a+.75*(b-a))-np.searchsorted(pp,a+.25*(b-a)))
  read_ids,read_counts=np.unique(ids[lo:hi],return_counts=True)
  z.update(min_CpG_depth=int(nn.min()),median_CpG_depth=float(np.median(nn)),p90_CpG_depth=float(np.quantile(nn,.9)),max_CpG_depth=int(nn.max()),max_read_observation_fraction=float(read_counts.max()/nn.sum()))
  z.update(CpGs=j-i,observations=int(nn.sum()),reads=len(np.unique(ids[lo:hi])),mean=float(np.mean(mm/nn)*100),weighted=float(mm.sum()/nn.sum()*100),span=int(pp[-1]-pp[0]+1),max_gap=gap,central_CpGs=central,opaque_fraction=float(opaque[lo:hi].mean()),edge_flag=gap>(b-a)/2 or central==0)
  kpos,kids,kp,kn,km=known_refs[ref];kl,kh=np.searchsorted(kp,[a,b]);rl,rh=np.searchsorted(kpos,[a,b])
  if kh>kl:z.update(known_CpGs=int(kh-kl),known_reads=len(np.unique(kids[rl:rh])),known_mean=float(np.mean(km[kl:kh]/kn[kl:kh])*100))
  return z
 out=[];background=[]
 for ar in ARR:
  aid=ar['array_id'];ref=ar['seqid'];a,b=int(ar['plot_start0']),int(ar['plot_end0']);s,e=int(ar['start0']),int(ar['end0'])
  bg={k:stats(ref,x,y) for k,x,y in [('left',a,s),('inside',s,e),('right',e,b)]}
  background.append(dict(array_id=aid,run=run,**{k+'_'+t:v for k,d in bg.items() for t,v in d.items()}))
  for w,step in SCHEMES:
   if BASELINE_ONLY and (w,step)!=(10000,10000):continue
   if aid not in REP and (w,step) not in [(10000,10000),(20000,5000),(50000,10000)]:continue
   for x in range(a//step*step,b,step):
    lo,hi=max(a,x),min(b,x+w);q=stats(ref,lo,hi)
    r=dict(array_id=aid,family=ar['family'],chromosome=ar['chromosome'],run=run,window_bp=w,step_bp=step,start0=lo,end0=hi,partial=(hi-lo)!=w,center=(lo+hi)/2,inside=(lo>=s and hi<=e),**q)
    if aid in REP and lo-30000>=s and hi+30000<=e:
     l,rt=stats(ref,lo-30000,lo),stats(ref,hi,hi+30000)
     r.update(local_known_left=l['known_mean'],local_known_right=rt['known_mean'],local_known_depth=min(l['known_mean'],rt['known_mean'])-q['known_mean'],local_known_flank_CpGs=min(l['known_CpGs'],rt['known_CpGs']),local_known_flank_reads=min(l['known_reads'],rt['known_reads']),local_left=l['mean'],local_right=rt['mean'],local_depth=min(l['mean'],rt['mean'])-q['mean'],local_flank_CpGs=min(l['CpGs'],rt['CpGs']),local_flank_reads=min(l['reads'],rt['reads']),local_flank_edge=l['edge_flag'] or rt['edge_flag'])
    out.append(r)
 frame=pd.DataFrame(out);frame.to_csv(dest/('baseline_windows.tsv.gz' if BASELINE_ONLY else 'windows.tsv.gz'),sep='\t',index=False)
 pd.DataFrame(background).to_csv(dest/'background.tsv',sep='\t',index=False)
 # Exact baseline uses the same grid and actual CpG denominators.
 expected=pd.read_csv(D/'inputs/v43_windows.tsv',sep='\t');expected=expected[expected.run==run]
 base=frame[(frame.window_bp==10000)&(frame.step_bp==10000)].merge(expected,on=['array_id','run','start0','end0'],validate='one_to_one')
 assert len(base)==len(expected)
 for a,b in [('CpGs','observed_CpGs'),('observations','valid_observations')]:assert np.array_equal(base[a],base[b])
 assert np.allclose(base['mean'],base.mean_CpG_5mC_pct,equal_nan=True,atol=1e-9)
 old=pd.read_csv(dest/'window_read_counts.tsv',sep='\t');base['bin_index']=base.start0//10000
 assert np.array_equal(base.merge(old,on=['array_id','run','bin_index']).reads,base.merge(old,on=['array_id','run','bin_index']).new_unique_reads)
 (dest/'baseline_verified.json').write_text(json.dumps({'run':run,'windows':len(base),'CpG_counts_observations_unique_reads_and_means_match':True}))
 print(run,'baseline verified',len(base),'parameter windows',len(frame),flush=True)
 return len(frame)
def main():
 with ProcessPoolExecutor(max_workers=6) as ex:list(ex.map(work,[f'ERR167823{i}' for i in range(14,20)]))
 if BASELINE_ONLY:return
 pd.concat([pd.read_csv(p,sep='\t') for p in (D/'per_run').glob('ERR*/windows.tsv.gz')],ignore_index=True).to_csv(D/'results/windows.tsv.gz',sep='\t',index=False)
 pd.concat([pd.read_csv(p,sep='\t') for p in (D/'per_run').glob('ERR*/background.tsv')],ignore_index=True).to_csv(D/'results/background.tsv',sep='\t',index=False)
if __name__=='__main__':main()
