"""Sequence-only structural atlas; no methylation imputation for congeners."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/danio-array-revision-mpl')
import collections,csv,json,subprocess,concurrent.futures,hashlib
from pathlib import Path
import pysam,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager,colors
from matplotlib.collections import PolyCollection
from matplotlib.backends.backend_pdf import PdfPages
W=Path(__file__).resolve().parents[1]
ROOT=W.parents[1]
E=ROOT / os.environ.get('DANIO_ASSEMBLY_RUN_ROOT', 'inputs/assembly_runs')
O=W/'congener_structure';O.mkdir(exist_ok=True)
BIN=os.environ.get('MODDOTPLOT', 'moddotplot')
font_manager.fontManager.addfont('/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf')
plt.rcParams.update({'font.family':['DejaVu Sans','Droid Sans Fallback'],'pdf.fonttype':42,'svg.fonttype':'none','font.size':9,'axes.unicode_minus':False})
def rows(p):
    with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def save(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
def job(x):
    d=O/x['array_id'];d.mkdir(exist_ok=True)
    if (d/'checked.json').exists():
        z=json.loads((d/'checked.json').read_text());z['region']=x;save(d/'checked.json',z);return z
    info=json.loads((E/'inputs/reference'/x['assembly']/'ready.json').read_text());meta=[]
    with pysam.FastaFile(info['fasta']) as fa:
        n=fa.get_reference_length(x['seqid']);a,b=int(x['start0']),int(x['end0'])
        # Context maintains requested 200 kb flanks; a nested small-array view
        # resolves organization below the common 5 kb display scale.
        views=[('context',int(x['plot_start0']),int(x['plot_end0']),5000,150),('detail',max(0,a-2000),min(n,b+2000),250,50)]
        for label,left,right,win,sketch in views:
            seq=fa.fetch(x['seqid'],left,right).upper();name=x['array_id']+'_'+label;p=d/(label+'.fa');p.write_text('>'+name+'\n'+seq+'\n')
            out=d/label;cmd=[BIN,'static','-f',str(p),'-k','17','-w',str(win),'-m',str(sketch),'-id','80','-d','.5','-o',str(out),'--no-plot']
            if not list(out.rglob('*.bedpe')):
                with open(d/(label+'.log'),'w') as log:
                    q=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
                assert q.returncode==0,(name,q.returncode)
            beds=list(out.rglob('*.bedpe'));assert len(beds)==1
            rr=rows(beds[0]);assert all(80<=float(r['perID_by_events'])<=100 and 1<=int(r['query_start'])<=right-left for r in rr)
            meta.append(dict(view=label,start0=left,end0=right,window_bp=win,command=cmd,BEDPE=str(beds[0]),matrix_cells=len(rr),sequence_sha256=hashlib.sha256(seq.encode()).hexdigest()))
    z=dict(region=x,views=meta,scope='observed family-supported short or large assembly region; no methylation data for this species');save(d/'checked.json',z);return z
def draw(ax,v,reg):
    a,b=v['start0'],v['end0'];poly=[];ident=[]
    for x in rows(v['BEDPE']):
        lo=(a+int(x['query_start'])-1)/1e6;hi=(a+min(b-a,int(x['query_end'])))/1e6;lo2=(a+int(x['reference_start'])-1)/1e6;hi2=(a+min(b-a,int(x['reference_end'])))/1e6
        poly.append([((s+t)/2,abs(s-t)/2) for s,t in [(lo,lo2),(hi,lo2),(hi,hi2),(lo,hi2)]]);ident.append(float(x['perID_by_events']))
    c=PolyCollection(poly,array=np.array(ident),cmap='viridis',norm=colors.Normalize(80,100),edgecolors='none',rasterized=True);ax.add_collection(c)
    ax.set_xlim(a/1e6,b/1e6);ax.set_ylim(-(b-a)/1e6*.025,(b-a)/2e6);ax.set_aspect('equal',adjustable='box');ax.set_yticks([]);ax.spines[['left','top','right']].set_visible(False)
    ax.plot([int(reg['start0'])/1e6,int(reg['end0'])/1e6],[0,0],color='#266b91' if reg['family']=='CL1' else '#bd612b',lw=3)
    ax.set_xlabel('参考坐标（Mb）');ax.ticklabel_format(useOffset=False,axis='x',style='plain');ax.set_title('两侧各200 kb背景' if v['view']=='context' else '局部结构放大（250 bp窗口）',fontsize=10);return c
def main():
    rr=rows(W/'congener_catalog/regions.tsv');results=[]
    with concurrent.futures.ThreadPoolExecutor(8) as pool:
        for z in pool.map(job,rr):
            results.append(z);save(O/'progress.json',dict(complete=len(results),total=len(rr)))
    dest=W/'publication/figures';dest.mkdir(parents=True,exist_ok=True)
    for spec in ['aesculapii','kyathit']:
        ss=[x for x in results if x['region']['species']==spec]
        ss.sort(key=lambda x:(x['region']['seqid'],int(x['region']['start0'])))
        with PdfPages(dest/f'Supplementary_{spec}_structure.pdf') as pdf:
            for z in ss:
                reg=z['region'];fig,ax=plt.subplots(2,1,figsize=(8.3,8));fig.subplots_adjust(top=.88,bottom=.14,left=.11,right=.87,hspace=.5)
                for a,v in zip(ax,z['views']):c=draw(a,v,reg)
                fam='富AT' if reg['family']=='CL1' else '富GC';large='大型区域' if reg['qualifies_large']=='True' else '较短相关串联区段'
                fig.suptitle(f"Danio {spec} · {reg['seqid']} · {fam}\n{int(reg['start0'])+1}–{reg['end0']} bp，{large}",fontsize=11)
                cb=fig.colorbar(c,cax=fig.add_axes([.9,.72,.016,.12]));cb.set_label('窗口相似度（%）',fontsize=8)
                fig.text(.11,.05,'背景窗口5 kb，放大窗口250 bp；显示下限80%。彩线为有家族证据的区域。\n本页为组装序列结构；该物种缺少本研究可用的甲基化数据。',fontsize=8)
                pdf.savefig(fig,dpi=220,bbox_inches='tight');plt.close(fig)
    save(O/'checked.json',dict(regions=len(results),all_supported_regions_shown=True,large_regions=sum(x['qualifies_large']=='True' for x in rr),methylation_generated=False,views_per_region=2,window_scales=[5000,250],context_flanks_bp=200000,cells_rasterized_labels_vector=True))
    print('congener structure regions',len(results))
if __name__=='__main__':main()
