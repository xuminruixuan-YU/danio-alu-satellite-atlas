from pathlib import Path
import importlib.util,json
import pandas as pd,numpy as np
from matplotlib.backends.backend_pdf import PdfPages
spec=importlib.util.spec_from_file_location('rr',Path(__file__).with_name('render_regional.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
D=Path(__file__).resolve().parents[2];P=D/'data/regional_display';F=Path(__import__('os').environ.get('DANIO_FIGURE_OUT',str(D/'figures')));F.mkdir(parents=True,exist_ok=True)
w=pd.read_csv(P/'results/windows.tsv.gz',sep='\t');rep=json.load(open(P/'inputs/figure_manifest.json'))['main_examples'];ars={x['array']['array_id']:x['array'] for x in json.load(open(P/'inputs/array_windows.json'))['matrices']}
names=['第3号染色体富GC第一阵列','第3号染色体富AT阵列','第3号染色体富GC第二阵列','第21号染色体富GC阵列','第21号染色体富AT阵列']
with PdfPages(F/'Supplementary_FigS8_window_comparison.pdf') as pdf:
 for aid,name in zip(rep,names):
  ar=ars[aid];a,b=int(ar['plot_start0']),int(ar['plot_end0']);fig,axs=m.plt.subplots(6,3,figsize=(11.7,8.3));fig.subplots_adjust(left=.075,right=.97,top=.88,bottom=.22,hspace=.32,wspace=.20)
  fig.suptitle('图S8  '+name,fontsize=13,y=.98)
  for col,(ww,step) in enumerate([(10000,10000),(20000,5000),(50000,10000)]):
   for i,run in enumerate(m.RUNS):
    g=w[(w.array_id==aid)&(w.run==run)&(w.window_bp==ww)&(w.step_bp==step)].sort_values('start0');ax=axs[i,col];ax.set_xlim(a/1e6,b/1e6);ax.set_ylim(0,100);ax.set_yticks([0,100]);ax.tick_params(labelsize=7);m.format_x(ax,i==5)
    if i==0:ax.set_title(f'{ww//1000} kb窗口 / {step//1000} kb步长',fontsize=10)
    if col==0:ax.set_ylabel(run[-2:],rotation=0,labelpad=12,fontsize=9)
    m.draw_curve(ax,g,3,'#333333',.65);m.draw_curve(ax,g,5,'#12699a',.85)
    gg=g.copy();gg.loc[gg.CpGs<20,'CpGs']=0;m.draw_curve(ax,gg,5,'#a34891',.55)
    for boundary in [int(ar['start0']),int(ar['end0'])]:ax.axvline(boundary/1e6,color='#aaa',lw=.5,zorder=0)
    fl=(g.CpGs>=10)&(g.reads>=3)&g.edge_flag;ax.plot(g.loc[fl,'center']/1e6,np.full(fl.sum(),1.04),marker='|',ms=3,color='#b33136',ls='none',transform=ax.get_xaxis_transform(),clip_on=False)
  fig.text(.005,.62,'CpG 5mC比例（%）',rotation=90,fontsize=10)
  lines=['黑线：≥10个有效CpG、≥3条读段；蓝线：≥10个有效CpG、≥5条读段；紫线：≥20个有效CpG、≥5条读段。',
         '红标：窗内有效CpG分布不均；灰线：阵列边界。缺测处断线，实测0%保留，孤立窗口显示短横线。']
  if aid==rep[0]:
   lines+=['各页样本14—19对应ERR16782314—ERR16782319，均值为有效参考CpG等权平均。',
           '重叠窗口共享观测，宽窗口可合并短缺口并降低空间分辨率；定位及红标判据见补充方法15、16。',
           '逐窗覆盖和连续性见表S54—S56，区域与侧翼均值见表S57；完整图注见补充资料图S8。']
  for j,line in enumerate(lines):fig.text(.075,.15-j*.025,line,fontsize=8)
  pdf.savefig(fig,dpi=240);fig.savefig(F/f'S8_representative_{rep.index(aid)+1}.png',dpi=150);m.plt.close(fig)
m.repair_pdf(F/'Supplementary_FigS8_window_comparison.pdf')
print('S8: five representative arrays; 3 scales; 3 support conditions')
