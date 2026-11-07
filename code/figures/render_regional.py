"""Render equal-scale BEDPE triangles and CpG-equal regional tracks.

No alignment, call filtering, interpolation or data-dependent region selection.
Coordinates u=(s+t)/2, v=abs(s-t)/2 use identical physical scale on both axes.
"""
from pathlib import Path
import argparse,json,math,sys,os
os.environ.setdefault('MPLCONFIGDIR','/tmp/danio-regional-plots')
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors,font_manager
from matplotlib.collections import PolyCollection
from matplotlib.ticker import MaxNLocator,ScalarFormatter
from pypdf import PdfReader,PdfWriter
sys.path.insert(0,str(Path(__file__).resolve().parent))
from pdf_image_compat import repair_pdf

font_manager.fontManager.addfont('/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf')
plt.rcParams.update({'font.family':['DejaVu Sans','Droid Sans Fallback'],'pdf.fonttype':42,'svg.fonttype':'none','font.size':9,'axes.unicode_minus':False})
COL={'CL1':'#266b91','CL2':'#bd612b'};LAB={'CL1':'富AT','CL2':'富GC'}
RUNS=[f'ERR167823{i}' for i in range(14,20)]
CPG_COL=['#ffffff','#d9d9d9','#c6dbef','#6baed6','#2171b5','#08306b']
READ_COL=['#ffffff','#d9d9d9','#c7e9c0','#74c476','#238b45','#005a32']
GEOMETRY=[];LINES=[]
def axinch(fig,x,y,w,h):
    fw,fh=fig.get_size_inches();return fig.add_axes([x/fw,y/fh,w/fw,h/fh])
def matrix_file(assets,value):
    p=Path(value)
    for key in ['congener_structure','structure']:
        if key in p.parts:return assets/Path(*p.parts[p.parts.index(key):])
    return assets/p
def triangle(fig,assets,v,x,y,cw,aid):
    a,b=int(v['start0']),int(v['end0']);ax=axinch(fig,x,y,cw,cw/2)
    rr=pd.read_csv(matrix_file(assets,v['BEDPE']),sep='\t')
    q0=(a+rr.query_start.to_numpy()-1)/1e6;q1=(a+np.minimum(b-a,rr.query_end.to_numpy()))/1e6
    r0=(a+rr.reference_start.to_numpy()-1)/1e6;r1=(a+np.minimum(b-a,rr.reference_end.to_numpy()))/1e6
    ss=np.stack([q0,q1,q1,q0],axis=1);tt=np.stack([r0,r0,r1,r1],axis=1)
    poly=np.stack([(ss+tt)/2,abs(ss-tt)/2],axis=2)
    heat=PolyCollection(poly,array=rr.perID_by_events.to_numpy(),cmap='viridis',norm=colors.Normalize(80,100),edgecolors='none',rasterized=True)
    ax.add_collection(heat);ax.set_xlim(a/1e6,b/1e6);ax.set_ylim(0,(b-a)/2e6);ax.set_aspect('equal',adjustable='box',anchor='SW')
    ax.set_yticks([]);ax.tick_params(axis='x',bottom=False,labelbottom=False);ax.spines[['left','right','top']].set_visible(False)
    fig.canvas.draw();box=ax.get_window_extent();scale_x=box.width/np.ptp(ax.get_xlim());scale_y=box.height/np.ptp(ax.get_ylim())
    assert abs(scale_x/scale_y-1)<1e-8
    GEOMETRY.append(dict(array_id=aid,view=v.get('view','context'),start0=a,end0=b,pixel_width=box.width,pixel_height=box.height,scale_x_over_y=scale_x/scale_y,matrix_cells=len(rr),matrix=str(matrix_file(assets,v['BEDPE']))))
    return ax,heat
def format_x(ax,label=True):
    ax.xaxis.set_major_locator(MaxNLocator(nbins=3));fmt=ScalarFormatter(useOffset=False);fmt.set_scientific(False);ax.xaxis.set_major_formatter(fmt)
    ax.tick_params(axis='x',labelsize=9,bottom=label,labelbottom=label)
    if label:ax.set_xlabel('参考坐标（Mb）',fontsize=9,labelpad=3)
def isolated(mask):
    mask=np.asarray(mask,bool);return np.flatnonzero(mask&~np.r_[False,mask[:-1]]&~np.r_[mask[1:],False])
def draw_curve(ax,g,threshold,color,lw):
    ok=(g.CpGs.to_numpy()>=10)&(g.reads.to_numpy()>=threshold);x=g.center.to_numpy()/1e6;y=np.where(ok,g['mean'],np.nan)
    ax.plot(x,y,color=color,lw=lw,solid_capstyle='round')
    ii=isolated(ok)
    # Short strokes identify measurements at a single window center; they do not span gaps.
    for i in ii:
        half=min(float(g.iloc[i].step_bp)*.18,float(g.iloc[i].end0-g.iloc[i].start0)*.18)/1e6
        ax.plot([x[i]-half,x[i]+half],[y[i],y[i]],color=color,lw=lw)
    LINES.append(dict(array_id=g.iloc[0].array_id,run=g.iloc[0].run,threshold=threshold,windows=int(ok.sum()),isolated_windows=len(ii),line_segments=int((ok[1:]&ok[:-1]).sum()),zero_windows=int((ok&(g['mean'].to_numpy()==0)).sum())))
    return ok
def strip(ax,g,key,palette,bounds):
    # Each center represents one sampling step, rather than filling the entire overlapping window.
    x=g.center.to_numpy()/1e6;step=g.step_bp.iloc[0]/1e6
    edges=np.r_[x[0]-step/2,(x[:-1]+x[1:])/2,x[-1]+step/2]
    cats=np.digitize(g[key].to_numpy(),bounds,right=False)
    ax.pcolormesh(edges,[0,1],cats[None,:],cmap=colors.ListedColormap(palette),vmin=0,vmax=5,shading='flat',rasterized=True)
    ax.set_yticks([]);ax.tick_params(bottom=False,labelbottom=False);ax.spines[:].set_visible(False)
def tu_panel(fig,assets,m,g,members,x,top,cw,label,large=False):
    ar=m['array'];a,b=int(ar['plot_start0']),int(ar['plot_end0']);font=10.5 if large else 9
    fig.text(x/fig.get_figwidth(),top/fig.get_figheight(),f"{label}  第{ar['chromosome']}号染色体 · {LAB[ar['family']]}阵列\n{int(ar['start0'])+1}–{ar['end0']} bp",fontsize=font,va='top',linespacing=1.35)
    y=top-(.44 if large else .50)-cw/2;tri,heat=triangle(fig,assets,dict(start0=a,end0=b,BEDPE=m['BEDPE']),x,y,cw,ar['array_id']);axes=[tri]
    band=axinch(fig,x,y-(.20 if large else .24),cw,.16 if large else .18);band.set_xlim(a/1e6,b/1e6);band.set_ylim(-.5,1.5);band.set_yticks([0,1],['GC','AT'],fontsize=7);band.spines[:].set_visible(False);band.tick_params(bottom=False,labelbottom=False)
    for r in members[(members.seqid==ar['seqid'])&(members.start0<b)&(members.end0>a)].itertuples():
        band.plot([max(a,r.start0)/1e6,min(b,r.end0)/1e6],[int(r.family=='CL1')]*2,lw=3,color=COL[r.family],solid_capstyle='butt')
    axes.append(band);t=y-(.27 if large else .31);track_h=.24 if large else .30;stride=.42 if large else .52
    for run in RUNS:
        q=g[(g.array_id==ar['array_id'])&(g.run==run)].sort_values('start0');assert len(q)
        ax=axinch(fig,x,t-track_h,cw,track_h);ax.set_xlim(a/1e6,b/1e6);ax.set_ylim(0,100);ax.set_yticks([0,100],['0','100'],fontsize=8.5 if large else 7)
        ax.set_ylabel(run[-2:],rotation=0,ha='right',va='center',fontsize=9,labelpad=5);ax.spines[['right','top','bottom']].set_visible(False);format_x(ax,False)
        ok=draw_curve(ax,q,3,'#303030',.85);draw_curve(ax,q,5,'#12699a',1.0)
        flag=ok&q.edge_flag.to_numpy()
        ax.plot(q.loc[flag,'center']/1e6,np.full(flag.sum(),1.035),linestyle='None',marker='|',ms=3,color='#b33136',transform=ax.get_xaxis_transform(),clip_on=False)
        for boundary in [int(ar['start0']),int(ar['end0'])]:ax.axvline(boundary/1e6,color='#aaa',lw=.45,zorder=0)
        axes.append(ax)
        for j,(key,palette,bounds,title) in enumerate([('CpGs',CPG_COL,[1,10,20,50,100],'C'),('reads',READ_COL,[1,3,5,10,20],'R')]):
            cov=axinch(fig,x,t-(.29 if large else .37)-j*(.055 if large else .065),cw,.045 if large else .055);cov.set_xlim(a/1e6,b/1e6);strip(cov,q,key,palette,bounds);axes.append(cov)
        t-=stride
    format_x(cov,True);cov.spines['bottom'].set_visible(True)
    fig.canvas.draw()
    for other in axes:
        assert np.allclose([other.get_position().x0,other.get_position().x1],[tri.get_position().x0,tri.get_position().x1])
    return heat
def legend(fig,y=.10,intro_y=None):
    fw,fh=fig.get_size_inches();fig.text(.06,y+.055 if intro_y is None else intro_y,'黑线：≥3条读段；蓝线：≥5条。红标：窗内有效CpG分布不均。',fontsize=9)
    for row,(title,palette,labels) in enumerate([('C 有效CpG数',CPG_COL,['0','1–9','10–19','20–49','50–99','≥100']),('R 独立读段数',READ_COL,['0','1–2','3–4','5–9','10–19','≥20'])]):
        yy=y-row*.028;fig.text(.06,yy,title,fontsize=8,va='center')
        for i,(col,txt) in enumerate(zip(palette,labels)):
            xx=.23+i*.115;ax=fig.add_axes([xx,yy-.006,.018,.012]);ax.set_facecolor(col);ax.set_xticks([]);ax.set_yticks([])
            for sp in ax.spines.values():sp.set_linewidth(.3)
            fig.text(xx+.022,yy,txt,fontsize=8,va='center')
def save(fig,out,name,svg=False):
    out.mkdir(parents=True,exist_ok=True)
    for ext in ['pdf','png']+(['svg'] if svg else []):
        p=out/(name+'.'+ext);fig.savefig(p,dpi=240,facecolor='white');
        if ext=='pdf':repair_pdf(p)
    plt.close(fig)
def join(paths,out):
    w=PdfWriter()
    for p in paths:w.append(PdfReader(p))
    with open(out,'wb') as f:w.write(f)
    repair_pdf(out)
def tu(args):
    metas=json.load(open(args.assets/'structure/matrices_checked.json'))['matrices'];metas.sort(key=lambda m:(int(m['array']['chromosome']),int(m['array']['start0'])))
    labels={x['array']['array_id']:x['figure'] for x in json.load(open(args.labels/'array_panel_manifest.json'))['regions']}
    g=pd.read_csv(args.windows,sep='\t');members=pd.read_csv(args.assets/'array_catalog/membership.tsv',sep='\t');reps=[m for m in metas if m['array']['chromosome'] in ['3','21']]
    if not args.supplement_only:
        fig=plt.figure(figsize=(9.0,11.15));slots=[(.80,10.85),(3.50,10.85),(6.20,10.85),(2.15,5.95),(4.85,5.95)]
        for i,(m,(x,t)) in enumerate(zip(reps,slots)):heat=tu_panel(fig,args.assets,m,g,members,x,t,2.05,chr(65+i),True)
        bar=fig.colorbar(heat,cax=fig.add_axes([.95,.88,.010,.075]));bar.ax.set_title('相似度\n(%)',fontsize=7,pad=5);bar.ax.tick_params(labelsize=7)
        legend(fig,.072,intro_y=.112);fig.text(.012,.62,'CpG 5mC比例（%）',rotation=90,fontsize=10);fig.text(.06,.015,'20 kb窗口 / 5 kb步长；≥10个有效CpG。留白为支持不足；真实0%保留。',fontsize=9)
        # Validate every footer line after actual font layout, rather than relying on spacing constants.
        fig.canvas.draw();footer=[t for t in fig.texts if t.get_position()[1]<.15 and t.get_rotation()==0];boxes=[t.get_window_extent(fig.canvas.get_renderer()) for t in footer]
        pairs=[(i,j) for i in range(len(boxes)) for j in range(i+1,len(boxes)) if boxes[i].overlaps(boxes[j])]
        assert not pairs, [(footer[i].get_text(),footer[j].get_text()) for i,j in pairs]
        (args.out/'main_footer_layout.json').write_text(json.dumps({'overlapping_text_pairs':pairs,'checked_texts':len(footer),'texts':[t.get_text() for t in footer]},ensure_ascii=False,indent=2))
        save(fig,args.out,'Fig4_structure_methylation',True)
        if args.main_only:return
    paths=[];groups=[]
    for chrom in sorted({m['array']['chromosome'] for m in metas},key=int):
        sub=[m for m in metas if m['array']['chromosome']==chrom];assert len(sub)<=4
        cw=2.95;fig=plt.figure(figsize=(11.7,14.2 if len(sub)==4 else 8.3))
        slots=[(2.30,13.70),(6.10,13.70),(2.30,7.60),(6.10,7.60)] if len(sub)==4 else [(.75+(3-len(sub))*1.75+i*3.50,7.80) for i in range(len(sub))]
        for (x,top),m in zip(slots,sub):heat=tu_panel(fig,args.assets,m,g,members,x,top,cw,labels[m['array']['array_id']])
        bar=fig.colorbar(heat,cax=fig.add_axes([.95,.79,.01,.12]));bar.ax.set_title('相似度\n(%)',fontsize=7,pad=5);bar.ax.tick_params(labelsize=7)
        legend(fig,.055)
        name=f'S4_chr{int(chrom):02d}';save(fig,args.out/'grouped_panels',name);paths.append(args.out/'grouped_panels'/(name+'.pdf'));groups.append(dict(chromosome=chrom,arrays=[m['array']['array_id'] for m in sub],figure=name))
    join(paths,args.out/'Supplementary_TU_array_structure_methylation.pdf')
    (args.out/'tu_group_manifest.json').write_text(json.dumps(dict(groups=groups,main_examples=[m['array']['array_id'] for m in reps],window_bp=20000,step_bp=5000,min_CpGs=10,min_reads=3,support_mark_reads=5,missing='NaN line breaks',means='equal effective reference CpG',triangles='equal physical scale'),ensure_ascii=False,indent=2))
def congeners(args):
    species=args.part;old=json.load(open(args.labels/'congener_panel_manifest.json'));rr=[r for r in old if r['region']['species']==species];group='S5' if species=='aesculapii' else 'S6';paths=[];groups=[]
    for ch in sorted({r['region']['chromosome'] for r in rr},key=lambda c:(not c.isdigit(),int(c) if c.isdigit() else 0)):
        subset=sorted([r for r in rr if r['region']['chromosome']==ch],key=lambda x:(x['region']['seqid'],int(x['region']['start0'])))
        for offset in range(0,len(subset),6):
            batch=subset[offset:offset+6];nc=min(3,len(batch));nr=math.ceil(len(batch)/nc);cw=3.15 if nc==3 else 3.8 if nc==2 else 5.0
            fw=.6+nc*cw+(nc-1)*.55+.85;rowheight=cw+1.90;fh=.65+nr*rowheight+.35;fig=plt.figure(figsize=(fw,fh))
            fig.suptitle(f'Danio {species} · '+(f'第{ch}号染色体' if ch.isdigit() else '未定位序列')+(f' · {offset//6+1}' if len(subset)>6 else ''),fontsize=12,y=1-.1/fh,va='top')
            for i,item in enumerate(batch):
                r=item['region'];data=json.load(open(args.assets/'congener_structure'/r['array_id']/'checked.json'));x=.6+(i%nc)*(cw+.55);top=fh-.57-(i//nc)*rowheight
                kind='较大阵列' if r['qualifies_large']=='True' else '较短区段';fig.text(x/fw,top/fh,f"{item['figure']}  {LAB[r['family']]} · {kind}\n{r['seqid']}: {int(r['start0'])+1}–{r['end0']} bp",fontsize=8.5,va='top',linespacing=1.5)
                for j,v in enumerate(data['views']):
                    bottom=top-.62-cw/2-j*(cw/2+.64);ax,heat=triangle(fig,args.assets,v,x,bottom,cw,r['array_id'])
                    ax.plot([int(r['start0'])/1e6,int(r['end0'])/1e6],[0,0],color=COL[r['family']],lw=2.5,clip_on=False)
                    ax.set_title('背景' if v['view']=='context' else '局部放大',fontsize=8,pad=2);format_x(ax,True);ax.tick_params(labelsize=7);ax.xaxis.label.set_fontsize(7)
            bar=fig.colorbar(heat,cax=fig.add_axes([(fw-.48)/fw,(fh-2.0)/fh,.10/fw,1.0/fh]));bar.ax.set_title('相似度\n(%)',fontsize=7,pad=5);bar.ax.tick_params(labelsize=7)
            name=f'{group}_chr{ch}_{offset//6+1}';save(fig,args.out/'grouped_panels',name);paths.append(args.out/'grouped_panels'/(name+'.pdf'));groups.append(dict(chromosome=ch,part=offset//6+1,figure=name,arrays=[x['region']['array_id'] for x in batch]))
    join(paths,args.out/f'Supplementary_{species}_structure.pdf')
    (args.out/f'{species}_group_manifest.json').write_text(json.dumps(dict(groups=groups,arrays=len(rr),triangles='equal physical scale'),ensure_ascii=False,indent=2))
def main():
    p=argparse.ArgumentParser();p.add_argument('part',choices=['tu','aesculapii','kyathit']);p.add_argument('--assets',type=Path,required=True);p.add_argument('--labels',type=Path,required=True);p.add_argument('--windows',type=Path);p.add_argument('--out',type=Path,required=True);p.add_argument('--supplement-only',action='store_true');p.add_argument('--main-only',action='store_true');args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    tu(args) if args.part=='tu' else congeners(args)
    (args.out/f"{args.part}{'_main' if args.main_only else ''}_geometry_audit.json").write_text(json.dumps(GEOMETRY,indent=2))
    if LINES:pd.DataFrame(LINES).to_csv(args.out/('track_line_main_audit.tsv' if args.main_only else 'track_line_audit.tsv'),sep='\t',index=False)
    print(args.part,'complete',len(GEOMETRY),'triangles',flush=True)
if __name__=='__main__':main()
