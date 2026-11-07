import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/danio-array-revision-mpl')
import csv,json
from pathlib import Path
from statistics import mean
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch,Rectangle
P=Path(__file__).resolve().parents[2];W=P/'data/array_study';R=P/'data/family_methylation';O=P/'figures';O.mkdir(exist_ok=True)
font_manager.fontManager.addfont('/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf')
plt.rcParams.update({'font.family':['DejaVu Sans','Droid Sans Fallback'],'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','font.size':9,'axes.spines.right':False,'axes.spines.top':False,'axes.unicode_minus':False})
COL={'CL1':'#266b91','CL2':'#bd612b'};LAB={'CL1':'富AT卫星','CL2':'富GC卫星'}
def rows(p):
    with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def output(fig,stem):
    for ext in ['pdf','png','svg']:fig.savefig(O/(stem+'.'+ext),dpi=260,bbox_inches='tight',facecolor='white')
    plt.close(fig)
def main():
    family=[]
    for run in [f'ERR167823{i}' for i in range(14,20)]:family+=rows(R/run/'family_call_primary.tsv')
    regional=rows(W/'array_metrics/array_methylation_family_summary.tsv');fig,ax=plt.subplots(1,2,figsize=(8,4.2),layout='constrained')
    for j,data in enumerate([family,regional]):
        field='five_mC_call_fraction' if j==0 else 'array_equal_5mC'
        for run in [f'ERR167823{i}' for i in range(14,20)]:
            ys=[100*float(next(x for x in data if x['run']==run and x['family']==fam)[field]) for fam in COL];ax[j].plot([0,1],ys,color='#999',lw=.9,zorder=0)
            for i,fam in enumerate(COL):ax[j].scatter(i,ys[i],color=COL[fam],s=22)
        ax[j].set_xticks([0,1],[LAB[f] for f in COL]);ax[j].set_xlim(-.35,1.35);ax[j].set_ylim(35,90);ax[j].set_ylabel('5mC比例（%）');ax[j].grid(axis='y',alpha=.15)
        ax[j].set_title(['A  家族归属读段','B  独立锚点支持的阵列等权均值'][j],loc='left',pad=15)
    output(fig,'Fig3_methylation')

if __name__=="__main__":main()
