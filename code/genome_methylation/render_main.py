"""Draw the pre-existing chr3/chr21 representatives from frozen genome windows."""
from pathlib import Path
import os,json,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data/genome_methylation'
OUT=Path(os.environ.get('DANIO_FIGURE_OUTPUT',ROOT/'figures/genome_methylation'));OUT.mkdir(parents=True,exist_ok=True)
font=Path('/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf')
if font.exists():font_manager.fontManager.addfont(str(font))
plt.rcParams.update({'font.family':['DejaVu Sans','Droid Sans Fallback'],'axes.unicode_minus':False,
  'pdf.fonttype':42,'svg.fonttype':'path','font.size':7.5,'axes.spines.top':False,'axes.spines.right':False})
cfg=json.loads((DATA/'parameters.json').read_text());runs=[i['run'] for i in cfg['items']]
fr=pd.read_csv(DATA/'fragments.tsv',sep='\t');fr=fr[fr.fragment.str.startswith('TU.')]
win=pd.concat([pd.read_csv(DATA/f'{r}.windows.tsv.gz',sep='\t',low_memory=False) for r in runs],ignore_index=True)
win=win[win.window_bp==200000]
colors={'CL1':'#287da7','CL2':'#bd622d'}

def merged(bounds,other):
    out=[]
    for a,b in sorted(bounds):
        if out and (a<=out[-1][1] or (a-out[-1][1]<=50000 and not any(x<a and y>out[-1][1] for x,y in other))):out[-1]=(out[-1][0],max(out[-1][1],b))
        else:out.append((a,b))
    return out

fig=plt.figure(figsize=(7.5,8.1))
grid=fig.add_gridspec(8,2,height_ratios=[.5,1,1,1,1,1,1,1.25],hspace=.26,wspace=.27)
audit=[]
for col,chrom in enumerate([3,21]):
    meta=next(r for r in cfg['reference_records'] if str(r['chromosome'])==str(chrom) and r['nuclear'])
    ref=meta['reference'];length=meta['length'];top=fig.add_subplot(grid[0,col]);f=fr[fr.seqid==ref]
    bounds={fam:list(zip(f.loc[f.family==fam,'start0'],f.loc[f.family==fam,'end0'])) for fam in colors}
    shade={fam:merged(bounds[fam],bounds['CL2' if fam=='CL1' else 'CL1']) for fam in colors}
    for fam,y in [('CL1',1),('CL2',0)]:
        top.broken_barh([(a/1e6,(b-a)/1e6) for a,b in bounds[fam]],(y,.6),facecolors=colors[fam],linewidths=0)
    top.set(ylim=(-.1,1.8),xlim=(0,length/1e6),yticks=[1.3,.3],yticklabels=['富AT','富GC'])
    top.tick_params(axis='x',bottom=False,labelbottom=False)
    top.set_title(f'{"AB"[col]}  第{chrom}号染色体',loc='left',fontsize=10,pad=7)
    coverage=[]
    for i,run in enumerate(runs):
        q=win[(win.run==run)&(win.reference==ref)].sort_values('start0')
        ax=fig.add_subplot(grid[i+1,col],sharex=top)
        for fam in colors:
            ax.broken_barh([(a/1e6,(b-a)/1e6) for a,b in shade[fam]],(0,1),transform=ax.get_xaxis_transform(),
                          facecolors=colors[fam],edgecolors='none',alpha=.75,linewidths=0,zorder=0)
        ok=q.display.to_numpy(bool);x=(q.start0.to_numpy()+q.end0.to_numpy())/2e6
        a=np.where(ok,q.CpG_equal_5mC.to_numpy()*100,np.nan)
        b=np.where(ok,q.observation_weighted_5mC.to_numpy()*100,np.nan)
        ax.plot(x,a,color='black',lw=.65);ax.plot(x,b,color='#d07a28',ls='--',lw=.65)
        lone=ok&~np.r_[False,ok[:-1]]&~np.r_[ok[1:],False]
        ax.scatter(x[lone],a[lone],color='black',s=6,zorder=4)
        ax.scatter(x[lone],b[lone],color='#d07a28',marker='x',s=8,zorder=5)
        red=ok&q.spatially_concentrated.to_numpy(bool)
        ax.plot(x[red],np.full(red.sum(),98),marker='|',ls='none',color='#b42b26',ms=3)
        ax.set(ylim=(-2,102),yticks=[0,50,100]);ax.grid(axis='y',alpha=.12)
        ax.set_ylabel(run[-2:],rotation=0,labelpad=9,va='center');ax.tick_params(axis='x',labelbottom=False)
        coverage.append(q.callable_reference_fraction.to_numpy()*100)
        audit.append(dict(chromosome=chrom,run=run,windows=len(q),displayed=int(ok.sum()),isolated=int(lone.sum()),red_marked=int(red.sum())))
    bottom=fig.add_subplot(grid[7,col],sharex=top)
    edges=np.r_[q.start0.to_numpy(),q.end0.iloc[-1]]/1e6
    im=bottom.pcolormesh(edges,np.arange(7)-.5,np.asarray(coverage),cmap='viridis',vmin=0,vmax=100)
    bottom.set(ylim=(5.5,-.5),yticks=range(6),yticklabels=[r[-2:] for r in runs],xlabel='TU参考坐标（Mb）')
    bottom.set_ylabel('可测CpG',fontsize=7.5)
fig.text(.019,.52,'CpG 5mC比例（%）',rotation=90,va='center',fontsize=9)
fig.legend(handles=[Line2D([0],[0],color='black',label='CpG等权'),Line2D([0],[0],color='#d07a28',ls='--',label='观测加权'),
    Patch(facecolor=colors['CL1'],alpha=.75,label='富AT'),Patch(facecolor=colors['CL2'],alpha=.75,label='富GC')],
    loc='lower center',bbox_to_anchor=(.5,.069),frameon=False,ncol=4,fontsize=8)
fig.text(.50,.048,'200 kb非重叠窗口；每窗≥10个有效CpG、≥3条独立读段；支持不足处断线。',ha='center',fontsize=7)
fig.text(.50,.027,'红标：窗内有效CpG分布不均。下方色带：有效CpG占参考CpG的比例。',ha='center',fontsize=7)
cax=fig.add_axes([.36,.124,.28,.013]);fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[0,50,100])
cax.tick_params(labelsize=6,pad=1)
fig.subplots_adjust(left=.105,right=.975,top=.95,bottom=.21)
for ext in ['pdf','png','svg']:fig.savefig(OUT/f'Fig5_chromosome_methylation.{ext}',dpi=300,facecolor='white')
plt.close(fig)
(OUT/'Fig5_chromosome_methylation.audit.json').write_text(json.dumps({'window_bp':200000,'heatmap_colorbar_center_x':0.5,'shade_alpha':.75,'chromosome_selection':'same chr3/chr21 as Fig4, not selected by dips','tracks':audit},ensure_ascii=False,indent=2))
print('Fig5',OUT)
