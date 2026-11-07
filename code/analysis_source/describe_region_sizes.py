"""Report descriptive species-specific size strata and sensitivity, not a fitted biological breakpoint."""
import json,csv,collections
from pathlib import Path
from statistics import median
from regional_catalog import W,rows,table,save
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
font_manager.fontManager.addfont('/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf')
plt.rcParams.update({'font.family':['DejaVu Sans','Droid Sans Fallback'],'pdf.fonttype':42,'svg.fonttype':'none','axes.unicode_minus':False})
def main():
    rr=rows(W/'congener_catalog/regions.tsv');out=[];fig,ax=plt.subplots(2,2,figsize=(8.4,6.4),layout='constrained')
    for j,sp in enumerate(['aesculapii','kyathit']):
      for i,f in enumerate(['CL1','CL2']):
        v=sorted([int(x['supported_bp']) for x in rr if x['species']==sp and x['family']==f],reverse=True)
        out.append(dict(species=sp,family=f,all_regions=len(v),minimum=min(v),median=median(v),maximum=max(v),primary_cutoff_bp=10000,large=sum(n>=10000 for n in v),small=sum(n<10000 for n in v),interpretation='descriptive size stratum; continuous GC distribution, not an inferred natural class'))
        a=ax[j,i];a.plot(range(1,len(v)+1),[n/1000 for n in v],'.-',ms=3,color=['#266b91','#bd612b'][i]);a.axhline(10,color='black',ls='--',lw=.8)
        a.set_title(f'D. {sp} · '+['富AT','富GC'][i],fontsize=10);a.set_ylabel('卫星支持覆盖（kb）');a.set_xlabel('区域长度降序');a.set_ylim(0,30)
    dest=W/'publication/figures';dest.mkdir(parents=True,exist_ok=True)
    for ext in ['pdf','png','svg']:fig.savefig(dest/('Supplementary_FigS3_region_sizes.'+ext),bbox_inches='tight',dpi=240)
    plt.close(fig);table(W/'congener_catalog/size_distribution.tsv',out)
    save(W/'congener_catalog/size_rule.json',dict(user_decision='congener size threshold need not match zebrafish; separate large and small regions',analyst_rule=dict(TU_supported_bp=50000,congener_supported_bp=10000,gap_bp=25000),sensitivity_bp=[5000,15000,20000,25000,50000],basis=['D. aesculapii AT supported lengths 0.586,2.139,15.642,18.306,22.304 kb; 10 kb separates the observed long/short groups','both congener GC distributions form a continuum; 10 kb is a rounded descriptive stratum, not a natural breakpoint','same within-species cutoff for both families; exact covered length and gaps retained','do not force one AT or two GC per chromosome'],limits=['non-T2T assembly completeness constrains observed region size','different primary cutoffs disallow direct species abundance comparisons','D. kyathit AT has only one supported 424 bp component; kept as short evidence, not promoted to large']))
    print(out)
if __name__=='__main__':main()
