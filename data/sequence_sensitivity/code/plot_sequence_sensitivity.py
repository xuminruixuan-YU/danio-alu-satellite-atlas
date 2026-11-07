"""Render sequence-only supplementary figure S7 from frozen regional results."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
D=Path(__file__).resolve().parents[1]
O=Path(__import__('os').environ.get('DANIO_FIGURE_OUT',str(D.parents[1]/'figures')));O.mkdir(parents=True,exist_ok=True)
font_manager.fontManager.addfont('/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf')
plt.rcParams.update({'font.family':['DejaVu Sans','Droid Sans Fallback'],'pdf.fonttype':42,'svg.fonttype':'path'})
with (D/'results/all_regions.tsv').open() as f:rr=list(csv.DictReader(f,delimiter='\t'))
fig,axs=plt.subplots(1,2,figsize=(8,4.3),layout='constrained')
for fam,col in [('CL1','#266b91'),('CL2','#bd612b')]:
 ids=sorted({r['array_id'] for r in rr if r['family']==fam})
 for aid in ids:
  z=[next(r for r in rr if r['array_id']==aid and r['scenario']==n) for n in ['score75_margin10_strict','score75_margin10_equivalent']]
  for ax,k in zip(axs,['retained_fraction','normalized_edit_distance']):ax.plot([0,1],[100*float(r[k]) for r in z],color=col,alpha=.5,marker='o',ms=3,lw=.7)
 for ax in axs:ax.plot([],[],color=col,label='富AT' if fam=='CL1' else '富GC')
for ax,title,y in zip(axs,['A  可比较拷贝池覆盖','B  阵列内序列差异'],['占区域卫星覆盖（%）','归一化编辑距离（%）']):
 ax.set_xticks([0,1],['严格标准化','等价序列合并']);ax.set_ylabel(y);ax.set_title(title,fontsize=11);ax.legend(frameon=False);ax.spines[['top','right']].set_visible(False)
for ext in ['pdf','png','svg']:fig.savefig(O/f'Supplementary_FigS7_sequence_sensitivity.{ext}',dpi=220)
plt.close(fig)
