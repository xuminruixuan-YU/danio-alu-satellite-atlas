"""Reproduce Figure 4 and S4 from validated 500 kb flank tables and BEDPEs."""
from pathlib import Path
import os,sys,json
import pandas as pd
from types import SimpleNamespace
BASE=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(BASE/'code/figures'))
import render_regional as rr
from cover_pages import append_cover
T=BASE/'data/regional_500kb'
OUT=Path(os.environ.get('DANIO_REGIONAL_OUTPUT',BASE/'figures/regional_500kb'));OUT.mkdir(parents=True,exist_ok=True)

def main():
    args=SimpleNamespace(assets=T,labels=T,windows=T/'windows.tsv.gz',out=OUT,supplement_only=False,main_only=True)
    rr.tu(args)
    main_geometry=rr.GEOMETRY.copy();main_lines=rr.LINES.copy();rr.GEOMETRY.clear();rr.LINES.clear()
    metas=json.loads((T/'structure/matrices_checked.json').read_text())['matrices']
    metas.sort(key=lambda m:(int(m['array']['chromosome']),int(m['array']['start0'])))
    labels=json.loads((T/'labels.json').read_text());g=pd.read_csv(T/'windows.tsv.gz',sep='\t')
    members=pd.read_csv(T/'array_catalog/membership.tsv',sep='\t')
    from matplotlib.backends.backend_pdf import PdfPages
    body=OUT/'S4_body.pdf'
    with PdfPages(body) as pdf:
        for chrom in sorted({int(m['array']['chromosome']) for m in metas}):
            sub=[m for m in metas if int(m['array']['chromosome'])==chrom];count=len(sub)
            fig=rr.plt.figure(figsize=(11.7,14.2 if count==4 else 8.3))
            fig.text(.055,.986,'两侧各500 kb   20 kb窗口 / 5 kb步长   CpG 5mC（%）',va='top',fontsize=9)
            slots=[(2.30,13.70),(6.10,13.70),(2.30,7.60),(6.10,7.60)] if count==4 else [(.75+(3-count)*1.75+i*3.50,7.80) for i in range(count)]
            for (x,top),m in zip(slots,sub):heat=rr.tu_panel(fig,T,m,g,members,x,top,2.95,labels[m['array']['array_id']])
            bar=fig.colorbar(heat,cax=fig.add_axes([.95,.79,.01,.12]));bar.ax.set_title('相似度\n（%）',fontsize=7,pad=5);bar.ax.tick_params(labelsize=7)
            rr.legend(fig,.055);pdf.savefig(fig,dpi=240,facecolor='white');rr.plt.close(fig)
            print('S4 chr',chrom,flush=True)
    rr.repair_pdf(body)
    dest=OUT/'Supplementary_FigS4_TU_array_panels.pdf';append_cover(body,dest,4);rr.repair_pdf(dest)
    assert len(rr.GEOMETRY)==49 and len(main_geometry)==5
    for x in rr.LINES+main_lines:
        q=g[(g.array_id==x['array_id'])&(g.run==x['run'])]
        assert x['windows']==int(q['display3' if x['threshold']==3 else 'support5'].sum())
    (OUT/'regional_500kb_checked.json').write_text(json.dumps({'passed':True,'flank_bp':500000,'main_geometry':main_geometry,'supplement_geometry':rr.GEOMETRY,'main_lines':main_lines,'supplement_lines':rr.LINES},indent=2))
if __name__=='__main__':main()
