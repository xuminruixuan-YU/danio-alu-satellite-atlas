"""Recompute deposited CpG-level windows in a separate output directory."""
from pathlib import Path
import argparse,os,shutil,subprocess,sys,json
import pandas as pd,numpy as np
R=Path(__file__).resolve().parents[2];T=R/'data/regional_display'
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--plots',action='store_true');args=p.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
assert out!=T.resolve(), 'Use a separate output directory'
shutil.copytree(T/'inputs',out/'inputs',dirs_exist_ok=True);(out/'results').mkdir(exist_ok=True);(out/'figures').mkdir(exist_ok=True)
for run in (T/'per_run').glob('ERR*'):
 dest=out/'per_run'/run.name;dest.mkdir(parents=True,exist_ok=True)
 for f in run.iterdir():
  if not (dest/f.name).exists():(dest/f.name).symlink_to(f.resolve())
env=os.environ.copy();env['DANIO_COMPARISON_DIR']=str(out);env['MPLCONFIGDIR']=str(out/'mplconfig')
subprocess.run([sys.executable,str(R/'code/regional_display/compare.py')],env=env,check=True)
keys=['array_id','run','window_bp','step_bp','start0','end0'];a=pd.read_csv(out/'results/windows.tsv.gz',sep='\t').sort_values(keys).reset_index(drop=True);b=pd.read_csv(T/'results/windows.tsv.gz',sep='\t').sort_values(keys).reset_index(drop=True)
assert list(a.columns)==list(b.columns) and a[keys].equals(b[keys])
for k in a:
 if pd.api.types.is_numeric_dtype(a[k]):assert np.allclose(a[k],b[k],equal_nan=True,atol=1e-8),k
 else:assert a[k].fillna('').equals(b[k].fillna('')),k
chosen=pd.read_csv(T/'regional_windows_20kb_step5kb.tsv.gz',sep='\t');meta=json.load(open(T/'inputs/array_windows.json'));ars={x['array']['array_id']:x['array'] for x in meta['matrices']}
for run,g in chosen.groupby('run'):
 obs=pd.read_csv(T/'per_run'/run/'valid_observations.tsv.gz',sep='\t');assert not obs.duplicated(['read','reference','position0']).any();refs={ref:np.unique(v.position0) for ref,v in obs.groupby('reference')}
 for r in g.itertuples():
  pp=refs[ars[r.array_id]['seqid']];edges=np.linspace(r.start0,r.end0,5);counts=np.diff(np.searchsorted(pp,edges))
  assert list(counts)==[r.q1,r.q2,r.q3,r.q4]
  ok=r.CpGs>=10 and r.reads>=3;assert ok==r.display3 and (ok and r.reads>=5)==r.support5
  assert (ok and r.plot_5mC_pct==r.mean) or (not ok and pd.isna(r.plot_5mC_pct))
chosen.to_csv(out/'regional_windows_20kb_step5kb.tsv.gz',sep='\t',index=False)
report={'passed':True,'all_parameter_windows':len(a),'chosen_windows':len(chosen),'CpG_means_counts_read_counts_and_spatial_fields_reproduced':True,'missing_and_zero_distinguished':True,'read_CpG_deduplication':True}
(out/'reproduction.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))
if args.plots:
 for part in ['tu','aesculapii','kyathit']:
  subprocess.run([sys.executable,str(R/'code/figures/render_regional.py'),part,'--assets',str(R/'data/array_study'),'--labels',str(R/'data/figure_source/figures'),'--windows',str(out/'regional_windows_20kb_step5kb.tsv.gz'),'--out',str(out/'figures')],env=env,check=True)
