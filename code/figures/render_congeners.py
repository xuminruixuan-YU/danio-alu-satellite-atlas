import os,csv,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
P=Path(__file__).resolve().parents[2];W=P/'data/array_study';O=Path(os.environ.get('DANIO_FIGURE_OUT',str(P/'figures')));O.mkdir(parents=True,exist_ok=True)
font_manager.fontManager.addfont('/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf')
COL={'CL1':'#266b91','CL2':'#bd612b'};LAB={'CL1':'富AT卫星','CL2':'富GC卫星'}
def read(p):
 with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def output(fig,stem,dpi):
 for ext in ['pdf','png','svg']:fig.savefig(O/(stem+'.'+ext),dpi=dpi,bbox_inches='tight',facecolor='white')
 plt.close(fig)
def main():
 from matplotlib.patches import Patch
 plt.rcParams.update({'font.family':['DejaVu Sans','Droid Sans Fallback'],'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','font.size':9,'axes.spines.right':False,'axes.spines.top':False,'axes.unicode_minus':False})
 regions=read(W/'congener_catalog/regions.tsv');inventory=read(P/'data/figure_source/congener_chromosomes.tsv')
 for num,species in enumerate(['aesculapii','kyathit'],1):
  chrom=sorted([r for r in inventory if r['species']==species],key=lambda x:int(x['chromosome']));fig,a=plt.subplots(figsize=(8.3,8),layout='constrained')
  for i,ch in enumerate(chrom):
   a.plot([0,int(ch['length'])/1e6],[i,i],color='#ddd',lw=1.5)
   for x in [r for r in regions if r['species']==species and r['seqid']==ch['seqid']]:
    lo,hi=int(x['start0'])/1e6,int(x['end0'])/1e6;y=i+(-.12 if x['family']=='CL1' else .12);large=x['qualifies_large']=='True'
    a.plot([lo,hi],[y,y],lw=4 if large else 1.2,color=COL[x['family']],alpha=1 if large else .6);a.plot((lo+hi)/2,y,marker='|' if large else '.',ms=7 if large else 2.3,color=COL[x['family']],alpha=1 if large else .6)
  a.set_yticks(range(len(chrom)),[x['chromosome'] for x in chrom]);a.invert_yaxis();a.set_xlabel('参考坐标（Mb）');a.set_ylabel('染色体');a.set_title(f'Danio {species}',loc='left',pad=32,fontstyle='italic')
  a.legend(handles=[Patch(color=COL[f],label=LAB[f]) for f in COL],loc='lower right',bbox_to_anchor=(1,1.01),ncol=2,frameon=False)
  fig.text(.15,-.035,'竖线：卫星支持覆盖≥10 kb；浅色点：<10 kb。',fontsize=8)
  output(fig,f'Supplementary_FigS{num}_{species}_distribution',260)
 plt.rcParams.update({'font.size':10,'axes.spines.top':True,'axes.spines.right':True})
 fig,ax=plt.subplots(2,2,figsize=(8.4,6.4),layout='constrained')
 for j,sp in enumerate(['aesculapii','kyathit']):
  for i,fam in enumerate(['CL1','CL2']):
   v=sorted([int(x['supported_bp']) for x in regions if x['species']==sp and x['family']==fam],reverse=True);a=ax[j,i]
   a.plot(range(1,len(v)+1),[n/1000 for n in v],'.-',ms=3,color=COL[fam])
   if len(v)==1:
    a.lines[-1].set_marker('o');a.lines[-1].set_markersize(6);a.lines[-1].set_clip_on(False)
    a.set_xticks([1]);a.set_xlim(.5,1.5)
    a.annotate(f'{v[0]/1000:.3f} kb',(1,v[0]/1000),xytext=(8,8),textcoords='offset points',color=COL[fam],fontsize=9)
   a.axhline(10,color='black',ls='--',lw=.8)
   a.set_title(f'D. {sp} · '+['富AT','富GC'][i],fontsize=10);a.set_ylabel('卫星支持覆盖（kb）');a.set_xlabel('区域序号（按支持覆盖降序）');a.set_ylim(0,30)
 output(fig,'Supplementary_FigS3_region_sizes',240)
if __name__=='__main__':main()
