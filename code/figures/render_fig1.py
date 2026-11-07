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
P=Path(__file__).resolve().parents[2];W=P/'data/array_study';R=P/'data/five_part';O=Path(os.environ.get('DANIO_FIGURE_OUT',str(P/'figures')));O.mkdir(parents=True,exist_ok=True)
font_manager.fontManager.addfont('/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf')
plt.rcParams.update({'font.family':['DejaVu Sans','Droid Sans Fallback'],'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','font.size':9,'axes.spines.right':False,'axes.spines.top':False,'axes.unicode_minus':False})
COL={'CL1':'#266b91','CL2':'#bd612b'};LAB={'CL1':'富AT卫星','CL2':'富GC卫星'}
def rows(p):
    with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def output(fig,stem):
    for ext in ['pdf','png','svg']:fig.savefig(O/(stem+'.'+ext),dpi=260,bbox_inches='tight',facecolor='white')
    plt.close(fig)
def main():
    regions=[x for x in rows(W/'array_catalog/regions.tsv') if x['species']=='TU' and x['qualifies_large']=='True'];members=rows(W/'array_catalog/membership.tsv');large={x['array_id'] for x in regions}
    chrom=rows(P/'data/figure_source/chromosomes.tsv');fig=plt.figure(figsize=(8.3,9.8),layout='constrained');gs=fig.add_gridspec(2,1,height_ratios=[5.2,1.25]);a=fig.add_subplot(gs[0]);b=fig.add_subplot(gs[1])
    for i,ch in enumerate(chrom):
        a.plot([0,int(ch['length'])/1e6],[i,i],color='#ddd',lw=1.7)
        for r in [x for x in regions if x['seqid']==ch['seqid']]:
            y=i+(-.14 if r['family']=='CL1' else .14)
            for x in [v for v in members if v['array_id']==r['array_id']]:a.plot([int(x['start0'])/1e6,int(x['end0'])/1e6],[y,y],color=COL[x['family']],lw=2.7,solid_capstyle='butt')
    a.set_yticks(range(25),range(1,26));a.invert_yaxis();a.set_ylabel('染色体');a.set_xlabel('TU参考坐标（Mb）');a.set_xlim(-.5,82);a.set_title('A  大型阵列区域的分布',loc='left',pad=30)
    a.legend(handles=[Patch(color=COL[f],label=LAB[f]) for f in COL],loc='lower right',bbox_to_anchor=(1,1.015),ncol=2,frameon=False)
    for j,fam in enumerate(COL):
        b.axhline(j,color='#ddd',lw=1)
        rr=[x for x in regions if x['chromosome']=='3' and x['family']==fam]
        for k,x in enumerate(sorted(rr,key=lambda x:int(x['start0'])),1):
            lo,hi=int(x['start0'])/1e6,int(x['end0'])/1e6;b.plot([lo,hi],[j,j],color=COL[fam],lw=8,solid_capstyle='butt');b.text((lo+hi)/2,j+.18,('AT' if fam=='CL1' else 'GC')+str(k),ha='center',fontsize=8)
    b.set_yticks([0,1],[LAB[f] for f in COL]);b.set_ylim(-.55,1.6);b.set_xlim(45,60);b.set_xlabel('第3号染色体参考坐标（Mb）');b.set_title('B  第3号染色体的富GC—富AT—富GC排列',loc='left',pad=12);output(fig,'Fig1_distribution')

if __name__ == "__main__": main()
