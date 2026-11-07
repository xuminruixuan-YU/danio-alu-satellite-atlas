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
P=Path(__file__).resolve().parents[2];W=P/'data/array_study'
O=Path(os.environ.get('DANIO_FIGURE_OUT',str(P/'figures')));O.mkdir(parents=True,exist_ok=True)
font_manager.fontManager.addfont('/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf')
plt.rcParams.update({'font.family':['DejaVu Sans','Droid Sans Fallback'],'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','font.size':9,'axes.spines.right':False,'axes.spines.top':False,'axes.unicode_minus':False})
COL={'CL1':'#266b91','CL2':'#bd612b'};LAB={'CL1':'富AT卫星','CL2':'富GC卫星'}
def rows(p):
    with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def output(fig,stem):
    for ext in ['pdf','png','svg']:fig.savefig(O/(stem+'.'+ext),dpi=260,bbox_inches='tight',facecolor='white')
    plt.close(fig)
def main():
    fig,ax=plt.subplots(1,2,figsize=(8.5,4.2),layout='constrained');stats=rows(W/'array_metrics/within_between_summary.tsv');summ=rows(W/'array_metrics/family_summary.tsv')
    variants=['gap10000_size50000','gap25000_size50000','gap50000_size50000','gap25000_size25000','gap25000_size100000','gap25000_size200000']
    for fam in COL:
        s=next(x for x in stats if x['family']==fam and x['stratum']=='all_between');yy=[100*float(s[k]) for k in ['within_same_arrays_mean','between_equal_target_array_mean']]
        ax[0].plot([0,1],yy,marker='o',color=COL[fam],label=LAB[fam])
        for i,v in enumerate(yy):ax[0].annotate(f'{v:.2f}',(i,v),xytext=(0,8),textcoords='offset points',ha='center',fontsize=8)
        ys=[100*float(next(x for x in summ if x['variant']==v and x['category']==fam+'_all')['array_equal_edit_distance']) for v in variants];ax[1].plot(range(6),ys,marker='o',color=COL[fam])
    ax[0].set_xticks([0,1],['阵列内','阵列间']);ax[0].set_xlim(-.25,1.25);ax[0].set_title('A  阵列内与阵列间差异',loc='left',pad=22)
    ax[1].set_xticks(range(6),['10/50','25/50','50/50','25/25','25/100','25/200'],rotation=30);ax[1].set_xlabel('间断上限/卫星覆盖下限（kb）');ax[1].set_title('B  大型区域划分敏感性',loc='left',pad=22)
    for a in ax:a.set_ylim(0,19);a.set_ylabel('平均归一化编辑距离（%）');a.grid(axis='y',alpha=.18)
    fig.legend(handles=[Patch(color=COL[f],label=LAB[f]) for f in COL],loc='lower center',bbox_to_anchor=(.5,1.02),ncol=2,frameon=False);output(fig,'Fig2_divergence')

if __name__=="__main__":main()
