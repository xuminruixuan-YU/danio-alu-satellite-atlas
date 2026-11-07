from pathlib import Path
import os,json,collections,zipfile,hashlib,shutil,subprocess
os.environ.setdefault('MPLCONFIGDIR','/tmp/danio-param-mpl')
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages
from concurrent.futures import ProcessPoolExecutor
D=Path(__import__('os').environ.get('DANIO_COMPARISON_DIR',str(Path(__file__).resolve().parents[2]/'data/regional_display')));ROOT=D.parents[1]
ARR={x['array']['array_id']:x['array'] for x in json.loads((D/'inputs/array_windows.json').read_text())['matrices']}
REP=json.loads((D/'inputs/figure_manifest.json').read_text())['main_examples']
RUNS=[f'ERR167823{i}' for i in range(14,20)]
PARAMS=[(10,10,10),(10,10,20),(10,5,10),(10,5,20),(20,5,10),(20,5,20),(50,10,10),(50,10,20),(100,20,10),(100,20,20),(10,10,5),(20,5,5)]
font_manager.fontManager.addfont('/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf')
plt.rcParams.update({'font.family':['DejaVu Sans','Droid Sans Fallback'],'pdf.fonttype':42,'font.size':8,'axes.unicode_minus':False})
def runs_of(mask):
 mask=np.asarray(mask,bool);z=np.r_[False,mask,False];return list(zip(np.where(np.diff(z.astype(int))==1)[0],np.where(np.diff(z.astype(int))==-1)[0]))
def evaluate():
 df=pd.read_csv(D/'results/windows.tsv.gz',sep='\t');bg=pd.read_csv(D/'results/background.tsv',sep='\t');metrics=[];gaps=[];dips=[];shape=[]
 for (aid,run,w,step),g in df.groupby(['array_id','run','window_bp','step_bp']):
  g=g.sort_values('start0').reset_index(drop=True);b=bg[(bg.array_id==aid)&(bg.run==run)].iloc[0];ar=ARR[aid]
  thresholds=[10,20]+([5] if (w,step) in [(10000,10000),(20000,5000)] else []) if aid in REP else [10]
  for t in thresholds:
   for nr in [1,3,5]:
    ok=(g.CpGs>=t)&(g.reads>=nr);segs=runs_of(ok);no=runs_of(~ok)
    maxspan=max([g.iloc[j-1].end0-g.iloc[i].start0 for i,j in segs],default=0);centers=max([g.iloc[j-1].center-g.iloc[i].center for i,j in segs],default=0)
    metrics.append(dict(array_id=aid,run=run,window_bp=w,step_bp=step,min_CpGs=t,min_reads=nr,exploratory=t==5,total_windows=len(g),display_windows=int(ok.sum()),display_fraction=float(ok.mean()),longest_display_span=maxspan,longest_center_span=centers,longest_missing_span=max([g.iloc[j-1].end0-g.iloc[i].start0 for i,j in no],default=0),median_CpGs=float(g.loc[ok,'CpGs'].median()) if ok.any() else np.nan,median_reads=float(g.loc[ok,'reads'].median()) if ok.any() else np.nan,edge_concentrated_display=int((ok&g.edge_flag).sum()),median_observations=float(g.loc[ok,'observations'].median()) if ok.any() else np.nan,max_observation_gap=float(g.max_gap.max()),model_majority_windows=int((ok&(g.opaque_fraction>=.5)).sum())))
    for i,j in no:gaps.append(dict(array_id=aid,run=run,window_bp=w,step_bp=step,min_CpGs=t,min_reads=nr,start0=int(g.iloc[i].start0),end0=int(g.iloc[j-1].end0),window_count=j-i))
    if aid not in REP:continue
    candidate=ok&(g.reads>=max(3,nr))&g.inside&(g.local_depth>=20)&(g.local_flank_CpGs>=10)&(g.local_flank_reads>=max(3,nr))
    for i,r in g[candidate].iterrows():
     dips.append(dict(array_id=aid,run=run,window_bp=w,step_bp=step,min_CpGs=t,min_reads=nr,start0=int(r.start0),end0=int(r.end0),depth=r.local_depth,mean=r['mean'],left=r.local_left,right=r.local_right,reads=r.reads,CpGs=r.CpGs,edge_flag=bool(r.edge_flag or r.local_flank_edge),opaque_fraction=r.opaque_fraction))
    bgok=all(b[f'{z}_CpGs']>=10 and b[f'{z}_reads']>=max(3,nr) for z in ['left','right'])
    arrlow=ok&g.inside&(g['mean']<=min(b.left_mean,b.right_mean)-20) if bgok else pd.Series(False,index=g.index)
    aa=[(i,j) for i,j in runs_of(arrlow) if j-i>=2]
    for i,j in aa:shape.append(dict(array_id=aid,run=run,window_bp=w,step_bp=step,min_CpGs=t,min_reads=nr,start0=int(g.iloc[i].start0),end0=int(g.iloc[j-1].end0),windows=j-i,edge_flag=bool(g.iloc[i:j].edge_flag.any()),opaque_fraction=float(g.iloc[i:j].opaque_fraction.max()),type='array_fixed_flanks'))
 pd.DataFrame(metrics).to_csv(D/'results/continuity.tsv',sep='\t',index=False);pd.DataFrame(gaps).to_csv(D/'results/gaps.tsv.gz',sep='\t',index=False)
 dipcols=['array_id','run','window_bp','step_bp','min_CpGs','min_reads','start0','end0','depth','mean','left','right','reads','CpGs','edge_flag','opaque_fraction']
 dp=pd.DataFrame(dips,columns=dipcols);dp.to_csv(D/'results/all_local_candidates.tsv.gz',sep='\t',index=False)
 sp=pd.DataFrame(shape,columns=['array_id','run','window_bp','step_bp','min_CpGs','min_reads','start0','end0','windows','edge_flag','opaque_fraction','type']);sp.to_csv(D/'results/array_shape_intervals.tsv',sep='\t',index=False)
 # A single common coordinate is required for replication, never union across samples.
 consensus=[]
 for typ,tbl in [('local',dp),('array_fixed_flanks',sp)]:
  stable=tbl[(tbl.min_CpGs>=10)&(tbl.min_reads==5)&(~tbl.edge_flag)]
  stable=stable[stable.opaque_fraction<.5]
  for aid in REP:
   ar=ARR[aid];s,e=int(ar['start0']),int(ar['end0']);g=stable[stable.array_id==aid]
   for x in range((s+9999)//10000*10000,e,10000):
    hits=g[(g.start0<=x)&(g.end0>x)];by=hits.groupby('run').window_bp.nunique();n=int((by>=3).sum());byshort=hits[hits.window_bp<100000].groupby('run').window_bp.nunique()
    if n:consensus.append(dict(type=typ,array_id=aid,anchor0=x,samples_with_three_widths=n,runs=';'.join(by[by>=3].index),stable_candidate=n>=3,samples_three_widths_excluding100=int((byshort>=3).sum())))
 cc=pd.DataFrame(consensus,columns=['type','array_id','anchor0','samples_with_three_widths','runs','stable_candidate','samples_three_widths_excluding100']);cc.to_csv(D/'results/common_coordinate_support.tsv',sep='\t',index=False)
 summaries=[]
 for aid in REP:
  b=bg[bg.array_id==aid];eligible=np.ones(len(b),bool)
  for z in ['left','inside','right']:eligible&=(b[z+'_CpGs']>=10)&(b[z+'_reads']>=3)
  lower=eligible&(b.inside_mean<b.left_mean)&(b.inside_mean<b.right_mean)
  stable=cc[(cc.array_id==aid)&cc.stable_candidate];local=dp[(dp.array_id==aid)&(dp.min_CpGs>=10)&(dp.min_reads>=3)]
  cl='固定侧翼参照下重复的低值区段；局部谷形另判' if len(stable) else '仅局部或特定参数下出现下降' if len(local) else '内部均值较低，缺少清楚谷形' if lower.any() else '有效覆盖不足或未显示两侧支持的下降'
  summaries.append(dict(array_id=aid,classification=cl,mean_comparable_samples=int(eligible.sum()),inside_lower_than_both_samples=int(lower.sum()),stable_coordinate_anchors=len(stable),local_candidate_windows=len(local),interpretation='均为探索性描述，不能单独定位功能着丝粒'))
 pd.DataFrame(summaries).to_csv(D/'results/representative_judgments.tsv',sep='\t',index=False)
 return summaries

def draw_page(pdf,df,aid,w,step,t,zero=False):
 ar=ARR[aid];a,b=int(ar['plot_start0']),int(ar['plot_end0']);s,e=int(ar['start0']),int(ar['end0']);fig=plt.figure(figsize=(12,12));gs=fig.add_gridspec(7,3,height_ratios=[.4,1,1,1,1,1,1],left=.08,right=.98,bottom=.08,top=.9,hspace=.25,wspace=.23,width_ratios=[3,1,1]);band=fig.add_subplot(gs[0,0]);band.set_xlim(a/1e6,b/1e6);band.set_ylim(0,1);band.axis('off')
 members=pd.read_csv(D/'inputs/membership.tsv',sep='\t');col={'CL1':'#266b91','CL2':'#bd612b'}
 for r in members[members.seqid==ar['seqid']].itertuples():
  if r.start0<b and r.end0>a:band.plot([max(a,r.start0)/1e6,min(b,r.end0)/1e6],[.5,.5],lw=4,color=col[r.family],solid_capstyle='butt')
 for k,run in enumerate(RUNS):
  g=df[(df.array_id==aid)&(df.run==run)&(df.window_bp==w*1000)&(df.step_bp==step*1000)].sort_values('start0');x=g.center.to_numpy()/1e6
  for j in range(3):
   ax=fig.add_subplot(gs[k+1,j]);ax.set_xlim(a/1e6,b/1e6)
   for bound in [s,e]:ax.axvline(bound/1e6,color='#888',lw=.6)
   if j==0:
    for nr,c,lw in [(1,'#abc4d4',1.3),(3,'#c47d35',.9),(5,'#22605e',.7)]:
     y=np.where((g.CpGs>=t)&(g.reads>=nr),g['mean'],np.nan)
     if zero:
      if nr!=1:continue
      missing=~np.isfinite(y);y=np.nan_to_num(y);ax.plot(x[missing],y[missing],'x',color='#999',ms=2)
     ax.plot(x,y,color=c,lw=lw)
     isolated=[i for i,jj in runs_of(np.isfinite(y)) if jj-i==1]
     if isolated:ax.plot(x[isolated],y[isolated],linestyle='None',marker='_',markersize=3,color=c,clip_on=False,zorder=4)
    flag=(g.CpGs>=t)&g.edge_flag
    ax.scatter(x[flag],g.loc[flag,'mean'],marker='|',s=12,color='#a92323')
    ax.set_ylim(0,100);ax.set_yticks([0,50,100]);ax.set_ylabel(run[-2:]+'\n5mC (%)')
   else:
    val=g.CpGs if j==1 else g.reads;ax.plot(x,val,color='#555',lw=.6);ax.set_yscale('symlog',linthresh=1);ax.set_ylim(0,max(2,float(val.max())*1.1));ax.set_title(['','有效CpG数','不同读段ID数'][j],fontsize=8) if k==0 else None
   if k==5:ax.set_xlabel('TU参考坐标 (Mb)')
   else:ax.tick_params(labelbottom=False)
   ax.spines[['top','right']].set_visible(False)
 title=f"chr{ar['chromosome']} · {'富AT' if ar['family']=='CL1' else '富GC'} · {aid}\n窗口 {w} kb / 步长 {step} kb / 至少 {t} CpG"+(' · 较弱支持探索对照' if t==5 else '')+(' · 仅宽尺度趋势' if w==100 else '')+(' · 原零占位对照' if zero else '')
 fig.suptitle(title,fontsize=12,y=.98);fig.text(.08,.03,'蓝：≥1条；橙：≥3条；绿：≥5条读段。红短线：覆盖集中/中部空缺。灰竖线：阵列边界。\n主要图不足门槛处断线，实测0%保留；重叠窗口共享观测。覆盖轴为symlog；甲基化轴始终0—100%。独立单窗测量保留短横线。'+(' 灰×表示占位零。' if zero else ''),fontsize=8)
 pdf.savefig(fig);plt.close(fig)
def plot_job(job):
 kind,key=job;df=pd.read_csv(D/'results/windows.tsv.gz',sep='\t')
 if kind=='rep':
  out=D/'figures'/f'representative_{REP.index(key)+1}.pdf'
  with PdfPages(out) as pdf:
   for w,step,t in PARAMS:draw_page(pdf,df,key,w,step,t)
   draw_page(pdf,df,key,10,10,10,True)
 else:
  w,step=key;out=D/'figures'/f'all49_{w}kb_step{step}kb_CpG10.pdf'
  with PdfPages(out) as pdf:
   for aid in ARR:draw_page(pdf,df,aid,w,step,10)
 print(out.name,flush=True)
def main():
 summaries=evaluate()
 with ProcessPoolExecutor(max_workers=4) as ex:list(ex.map(plot_job,[('rep',a) for a in REP]+[('atlas',(w,s)) for w,s in [(10,10),(20,5),(50,10)]]))
if __name__=='__main__':main()
