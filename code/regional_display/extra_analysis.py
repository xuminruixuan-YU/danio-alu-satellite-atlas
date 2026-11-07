from pathlib import Path
import numpy as np,pandas as pd,json
D=Path(__import__('os').environ.get('DANIO_COMPARISON_DIR',str(Path(__file__).resolve().parents[2]/'data/regional_display')))
def main():
 df=pd.read_csv(D/'results/windows.tsv.gz',sep='\t');rows=[]
 for (aid,run),aa in df.groupby(['array_id','run']):
  base=aa[(aa.window_bp==10000)&(aa.step_bp==10000)].sort_values('center')
  for (w,st),g in aa.groupby(['window_bp','step_bp']):
   g=g.sort_values('center').reset_index(drop=True)
   # Compare display status at the original 10 kb bin centers, avoiding window-count inflation.
   x=base.center.to_numpy();y=g.center.to_numpy();ix=np.searchsorted(y,x).clip(0,len(y)-1);prev=(ix-1).clip(0,len(y)-1);ix=np.where(abs(y[prev]-x)<abs(y[ix]-x),prev,ix);q=g.iloc[ix].reset_index(drop=True);b=base.reset_index(drop=True)
   for nr in [1,3,5]:
    old=(b.CpGs>=10)&(b.reads>=nr);new=(q.CpGs>=10)&(q.reads>=nr);rescued=~old&new;common=old&new
    delta=q.loc[common,'mean'].to_numpy()-b.loc[common,'mean'].to_numpy()
    rows.append(dict(array_id=aid,run=run,window_bp=w,step_bp=st,min_reads=nr,baseline_centers=len(b),baseline_displayed=int(old.sum()),new_displayed_at_baseline_centers=int(new.sum()),newly_displayed_centers=int(rescued.sum()),newly_displayed_edge_flag=int((rescued&q.edge_flag).sum()),lost_displayed_centers=int((old&~new).sum()),common_displayed=int(common.sum()),median_absolute_change_pp=float(np.median(abs(delta))) if len(delta) else np.nan,p90_absolute_change_pp=float(np.quantile(abs(delta),.9)) if len(delta) else np.nan,baseline_p90_minus_p10=float(np.ptp(np.quantile(b.loc[old,'mean'],[.1,.9]))) if old.any() else np.nan,new_p90_minus_p10=float(np.ptp(np.quantile(q.loc[new,'mean'],[.1,.9]))) if new.any() else np.nan,max_newly_displayed_unobserved_gap=int(q.loc[rescued,'max_gap'].max()) if rescued.any() else 0))
 pd.DataFrame(rows).to_csv(D/'results/same_coordinate_parameter_effects.tsv',sep='\t',index=False)
 # Removing reads whose modification-model label is opaque is a metadata sensitivity,
 # not a newly calibrated modification filter.
 z=df[(df.CpGs>=10)&(df.reads>=5)].copy();z['known_model_sufficient']=(z.known_CpGs>=10)&(z.known_reads>=5);z['known_minus_all_pp']=z.known_mean-z['mean']
 z[['array_id','run','window_bp','step_bp','start0','end0','CpGs','reads','mean','opaque_fraction','known_CpGs','known_reads','known_mean','known_model_sufficient','known_minus_all_pp']].to_csv(D/'results/model_label_sensitivity.tsv.gz',sep='\t',index=False)
 dp=pd.read_csv(D/'results/all_local_candidates.tsv.gz',sep='\t');keys=['array_id','run','window_bp','step_bp','start0','end0']
 extra=['known_mean','known_CpGs','known_reads','local_known_depth','local_known_flank_CpGs','local_known_flank_reads']
 dp=dp.merge(df[keys+extra],on=keys,validate='many_to_one');dp['known_only_retains_condition']=(dp.known_CpGs>=dp.min_CpGs)&(dp.known_reads>=np.maximum(3,dp.min_reads))&(dp.local_known_depth>=20)&(dp.local_known_flank_CpGs>=10)&(dp.local_known_flank_reads>=np.maximum(3,dp.min_reads))
 dp.to_csv(D/'results/local_candidates_model_sensitivity.tsv.gz',sep='\t',index=False)
 # All tested local windows, including failures, remain in windows.tsv.gz.
 diag=[]
 for (aid,run,w,st),g in df[df.local_depth.notna()].groupby(['array_id','run','window_bp','step_bp']):
  for nr in [3,5]:
   support=(g.CpGs>=10)&(g.reads>=nr)&(g.local_flank_CpGs>=10)&(g.local_flank_reads>=nr)
   diag.append(dict(array_id=aid,run=run,window_bp=w,step_bp=st,min_reads=nr,tested_positions=len(g),supported_positions=int(support.sum()),drop20_positions=int((support&(g.local_depth>=20)).sum()),drop20_without_edge_flags=int((support&(g.local_depth>=20)&~g.edge_flag&~g.local_flank_edge.astype(bool)).sum())))
 pd.DataFrame(diag).to_csv(D/'results/local_search_denominators.tsv',sep='\t',index=False)
 q=dp[(dp.min_CpGs==10)&(dp.min_reads==5)&~dp.edge_flag]
 q.groupby(['array_id','run','window_bp','step_bp']).agg(candidate_windows=('start0','size'),first_start0=('start0','min'),last_end0=('end0','max'),minimum_depth_pp=('depth','min'),maximum_depth_pp=('depth','max'),minimum_reads=('reads','min'),minimum_CpGs=('CpGs','min')).to_csv(D/'results/local_candidate_parameter_ranges.tsv',sep='\t')
 rr=[]
 for (aid,run),g in df[(df.window_bp==20000)&(df.step_bp==5000)&df.inside&(df.CpGs>=10)&(df.reads>=5)&~df.edge_flag].groupby(['array_id','run']):
  r=g.loc[g['mean'].idxmin()];rr.append(r[['array_id','run','start0','end0','mean','CpGs','reads','max_gap','local_depth']])
 pd.DataFrame(rr).to_csv(D/'results/minimum_20kb_supported_windows.tsv',sep='\t',index=False)
 keys=['array_id','run','start0','end0']
 pd.DataFrame(rr).merge(df[(df.window_bp==20000)&(df.step_bp==5000)][keys+['min_CpG_depth','median_CpG_depth','max_CpG_depth','max_read_observation_fraction']],on=keys,validate='one_to_one').to_csv(D/'results/minimum_window_read_support.tsv',sep='\t',index=False)
if __name__=='__main__':main()
