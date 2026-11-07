from pathlib import Path
import json,html,math,re
import pandas as pd,numpy as np
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,KeepTogether
from reportlab.lib.pagesizes import A4
from openpyxl import Workbook,load_workbook
from openpyxl.styles import Font,PatternFill,Alignment
D=Path(__import__('os').environ.get('DANIO_COMPARISON_DIR',str(Path(__file__).resolve().parents[2]/'data/regional_display')));R=D/'results'
ARR={x['array']['array_id']:x['array'] for x in json.loads((D/'inputs/array_windows.json').read_text())['matrices']}
REP=json.loads((D/'inputs/figure_manifest.json').read_text())['main_examples']
NAMES=['chr3 富GC左区','chr3 富AT','chr3 富GC右区','chr21 富GC','chr21 富AT']
def table_md(headers,rows):
 return '\n'.join(['| '+' | '.join(map(str,headers))+' |','|'+'---|'*len(headers)]+['| '+' | '.join(map(str,row))+' |' for row in rows])
def main():
 w=pd.read_csv(R/'windows.tsv.gz',sep='\t');c=pd.read_csv(R/'continuity.tsv',sep='\t');e=pd.read_csv(R/'same_coordinate_parameter_effects.tsv',sep='\t');b=pd.read_csv(R/'background.tsv',sep='\t');ref=pd.read_csv(R/'reference_CpG_background.tsv',sep='\t');cc=pd.read_csv(R/'common_coordinate_support.tsv',sep='\t');md=pd.read_csv(R/'model_label_sensitivity.tsv.gz',sep='\t')
 report=[]
 def add(s):report.append(s.strip())
 add('# v4.3区域甲基化多参数独立比较报告')
 add('本轮独立重汇总有效CpG观测并生成比较图集；完整比较保留于独立目录；正式图采用的子方案见新版稿件。六个样本分别处理，图中的14—19为ERR16782314—ERR16782319。参考固定为TU GCF_049306965.1。')
 add('## 主要结论与建议')
 add('建议以20 kb窗口、5 kb步长、至少10个有效CpG作为可读性与分辨率的折中展示，并同时显示≥3和≥5条不同读段ID支持条件、缺口及覆盖集中标记；保留10 kb非重叠结果作为分辨率参照，50 kb作为宽尺度辅助。该方案已由用户确认用于新版图4和图S4，完整探索条件另存。100 kb仅适合观察更宽趋势；5 CpG条件的连续性增益有限，应保留为弱支持探索对照。')
 add('chr21富GC在六个样本中均有较宽尺度的内部相对下降，曲线证据与内部有效CpG均值相符，但不是一个跨样本固定位置、边界明确的窄低谷。chr21富AT在四个样本中出现明显局部下降，最低区段有约80 kb的位置偏移，另外两个样本的阵列内部覆盖稀疏；该结果不能概括成六个样本共有的同一局部谷。chr3两处富GC的有效CpG均值通常较低，但连续性、局部形态和样本覆盖存在差异；chr3富AT尤其受限于长缺口。')
 add('固定紧邻30 kb为左右背景的局部诊断中，没有同一10 kb参考锚点达到“至少3个样本、各自至少3种窗口宽度、≥5条读段且通过空间覆盖标记”的预设重复条件。该条件的失败不否定较宽的相对下降，也不证明不存在其他尺度的低值。现有数据及所测试尺度下，尚不足以用高—低—高模式支持着丝粒定位。即使可见低谷，其形状也不能单独确认功能核心。')
 add('## 输入、计量与基线核验')
 add('从v4.3实际使用的既有主比对BAM及父比对记录重新解码，复用原读段集合、参考版本、MAPQ≥20、已返回高分竞争比对的CpG坐标一致性、Q10和修饰概率0.8条件。未重新比对或放宽定位/调用条件。原审计仅有10 kb逐窗读段汇总，无法据此准确计算滑窗读段数；本轮因此补存逐读段—参考CpG有效观测账本，并断言每样本(read,reference,position0)唯一。')
 checked=json.loads((R/'replay_checked.json').read_text());add(f"六样本合计{sum(x['CpG_rows_reproduced'] for x in checked)}条样本×参考CpG分类记录成功重放；原52074个10 kb窗口的有效CpG数、有效观测次数、不同读段ID数及CpG等权比例全部一致。核验通过后才执行其他参数。证据：results/replay_checked.json、baseline_all_verified.json及per_run/ERR*/baseline_verified.json。")
 add('窗口值为每个有效参考CpG的m/(C+m+h)比例的等权平均×100；C、m、h分别为合格未修饰、5mC和5hmC分类。观测加权比例另列weighted，未用于主要曲线。不同读段ID数是技术支持数，不是独立生物学重复；也不自动保证不同ID绝无同源分子关系。重叠窗共享观测。所有滑窗重新集合原CpG及其计数，没有平均邻窗均值。')
 add('原10 kb表中有24054个窗口无有效CpG、1654个仅有1—9个；其占位零不进入均值或低谷诊断。≥10 CpG的原窗中有1个实测0%窗口，该窗位于样本17的chr3富GC左区46900000—46910000，仅1条读段支持，在≥1读段对照中保留实测零，在≥3/5读段条件下按支持不足留白。主要图不足门槛处断线，孤立合格窗口用短横线保留；每个代表区域第13页另保留原占位零对照并以灰×标记。')
 add('## 固定方案与覆盖诊断')
 add('方案在inputs/plan.md预先记录。五种窗宽/步长为10/10、10/5、20/5、50/10、100/20 kb，各用10与20 CpG门槛；10/10及20/5另加5 CpG。读段支持分别为1、3、5。坐标以染色体0点对齐步长，边缘窗截短于原绘图区间并标记partial；各区域的横轴、样本顺序和0%—100%纵轴固定。图上结构条是同一参考坐标中的AT/GC片段，灰竖线为操作性阵列边界。')
 add('阵列尺度背景固定为边界外左右各200 kb，内部为完整区域跨度；均值比较要求每侧与内部≥10 CpG、≥3读段。局部背景为候选窗紧邻的左右各30 kb，均须位于阵列内、≥10 CpG、≥3读段；描述性下降条件为低于两侧各至少20个百分点，再检查≥5读段。该门槛不是显著性检验或经验证生物学标准。所有符合范围的窗均扫描，含未通过者的完整结果在windows.tsv.gz和local_search_denominators.tsv。')
 add('最长无有效CpG区间包括两端；超过半窗或中央一半无有效CpG即加红色空间覆盖标记。另给每CpG观测深度的最小值、中位数、90分位和最大值，以及单条读段最高观测贡献。参考CpG数量单独计算，使“参考本来少CpG”与“有参考CpG却缺有效观测”可区分；无有效CpG区间不等同于BAM完全无覆盖。')
 q=c[c.array_id.isin(REP)&(c.min_reads==3)].groupby(['window_bp','step_bp','min_CpGs']).agg(total=('total_windows','sum'),show=('display_windows','sum'),edge=('edge_concentrated_display','sum'),span=('longest_center_span','median'))
 rows=[]
 for (win,st,t),r in q.iterrows():rows.append([f'{win//1000}/{st//1000}',t,f'{r.show/r.total*100:.2f}%',f'{int(r.edge)}/{int(r.show)}',f'{r.span/1000:.1f}'])
 add('以下为五代表区域×六样本，≥3读段条件；连续长度为每区域/样本最长连续可显示窗口中心跨度的中位数。重叠窗数量和跨度均不是独立数据量。')
 add(table_md(['窗/步长(kb)','CpG门槛','显示比例','集中标记/显示窗','连续中心跨度(kb)'],rows))
 add('读段支持也影响连续性。五代表区域的10/10方案在≥1、≥3、≥5读段时分别可显示3803、3567、3277个窗口；20/5为8013、7604、7083个；50/10为4270、4139、3985个。各数只在同一窗/步长内部比较。较小读段门槛下显示、提高门槛后消失的片段可直接在蓝/橙/绿曲线与读段计数轨道中核查，不应仅因合窗后连线就视为更可靠。')
 add('在原5424个10 kb中心位置上，20/5方案新增显示249个，其中86个有空间覆盖集中标记；50/10新增572个，其中192个有标记；100/20新增793个，其中298个有标记。原共同显示位置的甲基化变化绝对值，按区域/样本计算的中位数再取中位数，分别为2.54、4.37和4.76个百分点。10/5基本保留原中心对应的同一10 kb统计量，主要增加相邻展示位置。')
 add('对49阵列×六样本，在同一52074个原中心位置、≥10 CpG且≥3读段条件下，原方案可显示24616个；20/5为26398个，50/10为29330个。新增显示1782和4714个中，分别771和2320个带集中标记。说明合窗既能汇集短缺口两侧的有效测量，也会把长缺口边缘的一小部分观测扩展为一个可显示均值。没有新增原始覆盖。')
 add('扩大窗口降低局部振幅并改变低值位置：chr21富AT样本15、16、18在20 kb窗可见较深下降，50/100 kb会将其与邻近较高观测汇总，谷值、宽度及支持位置因而改变。低值宽度在本报告中是候选窗口或合格窗口段的宽度，不是精确功能边界。逐条件的分位差、共同位置变化和全部短/长缺口见same_coordinate_parameter_effects.tsv及gaps.tsv.gz。对于原始有效CpG空隙超过100 kb的区段，本轮最大窗口也没有获得内部观测；连续显示若跨入该空隙，应查看红色标记和max_gap。')
 add('## 各代表区域的独立判断')
 judgments=[
 '仅局部或特定条件下出现下降；有效CpG内部均值在5/6样本低于两侧。样本14内部仅175个有效CpG，最长有效CpG空隙约1.56 Mb；其他样本较连续。满足≥5读段、无集中标记的局部20个百分点候选仅见于样本18的20 kb条件，不能概括为稳定共同低谷。',
 '主要属于数据不足；仅样本16的内部有效CpG均值低于两侧。全区域参考有7996个CpG，各样本有效56—1046个，最大空隙约0.35—1.32 Mb。扩大窗口不能补足内部证据，也不能用占位零说明低甲基化。',
 '有效CpG均值在6/6样本均低于两侧，均值差较一致；清楚而重复的局部谷形证据有限。内部存在约0.26—0.60 Mb的最长空隙；满足较宽松读段支持的少数局部候选未通过≥5读段、无集中标记条件。宜归为均值差明确、谷形与连续范围受限。',
 '六样本均支持较宽尺度的内部相对下降，内部有效CpG均值约52.53%—55.97%。相对固定200 kb两侧背景，多尺度支持的低值锚点出现于若干分离区段；不能将其最小与最大坐标之间全部连为一个谷。局部30 kb背景下，各样本窄低值的位置与参数支持不同，没有达到预设跨样本局部重复条件。可描述宽尺度趋势，不能给出固定功能核心。',
 '仅局部或特定条件下出现明显下降。样本15、16、18、19可见阵列内部局部低值；样本14、17内部有效CpG仅244、245个，最长有效CpG空隙分别为533.2 kb和685.2 kb。整个内部有效CpG均值在0/6样本低于左右两侧，故不能写成整阵列低甲基化。固定远侧背景在36.64—36.67 Mb的部分锚点有三样本多尺度低值支持；改用预设紧邻30 kb背景后，同一锚点最多只在一个样本达到三种宽度，说明背景定义和位置漂移至关重要。']
 jr=[]
 for i,(aid,name,judgment) in enumerate(zip(REP,NAMES,judgments),1):
  ar=ARR[aid];g=b[b.array_id==aid];nc=int(ref[(ref.array_id==aid)&(ref.zone=='inside')].reference_CpGs.iloc[0]);low=int(((g.inside_mean<g.left_mean)&(g.inside_mean<g.right_mean)).sum())
  add(f'### {i}. {name}（{aid}）');add(f"阵列坐标0基半开区间[{ar['start0']}, {ar['end0']})；代表图集representative_{i}.pdf，第1页原方案、第5页20 kb/5 kb、第7页50 kb/10 kb、第9页100 kb/20 kb（均为10 CpG）。")
  add(judgment)
  rr=[]
  for z in g.sort_values('run').itertuples():rr.append([z.run[-2:],f'{z.left_mean:.2f}',f'{z.inside_mean:.2f}',f'{z.right_mean:.2f}',f'{z.inside_CpGs}/{nc}',f'{z.inside_max_gap/1000:.1f}'])
  add(table_md(['样本末位','左均值%','内部均值%','右均值%','内部有效/参考CpG','最长空隙kb'],rr))
  jr.append(dict(array_id=aid,name=name,inside_lower_samples=low,interpretation=judgment))
 pd.DataFrame(jr).to_csv(R/'reviewed_representative_judgments.tsv',sep='\t',index=False)
 add('### chr21富AT的最低支持窗口与位置漂移')
 mins=pd.read_csv(R/'minimum_20kb_supported_windows.tsv',sep='\t');mins=mins[mins.array_id==REP[4]]
 add('下面列出20/5方案、≥10 CpG、≥5读段、无空间集中标记时，各样本阵列内的最低窗口。它是可评估部分的极小值，不是显著性检验；样本14/17只有靠近边缘的少量区段可用，不能据此与其他样本的内部谷一一对应。')
 add(table_md(['样本末位','窗口0基范围','均值%','CpG数','读段数','相对近背景下降pp'],[[z.run[-2:],f'{z.start0}—{z.end0}',f'{z.mean:.2f}',z.CpGs,z.reads,'不适用' if pd.isna(z.local_depth) else f'{z.local_depth:.2f}'] for z in mins.itertuples()]))
 add('这四个样本（15、16、18、19）的最低支持20 kb窗口分别有11、17、26、24条读段，每CpG有效观测深度的中位数为4、6、7、10；单条读段最大观测贡献约9.3%—21.5%。因此这些低窗并非仅由一条读段构成；覆盖与定位的其他限制仍分别适用。逐窗读段贡献见minimum_window_read_support.tsv。')
 add('样本15、16、18的最低20 kb窗口中心分别约36.645、36.680、36.600 Mb；最大偏移约80 kb。样本19最低窗约36.675 Mb，但相对紧邻30 kb两侧仅下降8.40个百分点；它在较宽背景/窗宽下更容易满足下降条件。这些位置分别报告，没有拼成同一共有窄谷。全部局部候选的深度、窗口宽度和支持条件见all_local_candidates.tsv.gz及local_candidate_parameter_ranges.tsv；后者的首末坐标只表示候选分布范围，不能当连续谷宽。')
 add('## 读段、模型及定位可靠性限制')
 z=md[md.known_model_sufficient];add(f"排除修饰模型标签不明确的读段后，在所有合格（≥10 CpG、≥5读段）且仍满足该门槛的窗口中，有{int((z.known_minus_all_pp.abs()>5).sum())}个窗口变化超过5个百分点，最大变化{z.known_minus_all_pp.abs().max():.2f}个百分点（包括49区域的多种重叠窗，不是独立事件数）。原合格窗口中有{int((~md.known_model_sufficient).sum())}个失去该门槛。五代表区域也有少数位置敏感，不能概括为所有窗口均与模型标签无关。")
 add('在五代表区域主要CpG条件、≥5读段且无原空间集中标记的127个局部候选条件记录中，去除标签受限读段后，127个仍满足CpG/读段和相对两侧下降条件。记录包含不同门槛及重叠窗口，并非127个独立低谷；此项只复核计数与比例，没有校准模型准确性，也没有重新证明定位或去除后的全部空间覆盖条件。完整变化表包含失败记录。')
 add('本轮复用的是v4.3已通过筛选的区域观测。坐标一致性是相对于已经返回的高分竞争比对，不能证明穷尽所有重复阵列的竞争位置。独立非目录序列锚段证据与此区域集合的证据等级不同，连续性和读段数量不能提升到精确重复单元定位。区域统计提供的是条件性参考CpG趋势。两个样本某些阵列的大段空缺、不同样本低值位置偏移、局部模型标签敏感性及CpG有效比例，均需要在后续使用时保留。')
 add('本轮未按形态选择或改变参数，未放宽定位/修饰过滤，未跨长缺口插值。不同样本分别估计；没有将窗口或读段当独立生物学重复，也没有进行未经设计的低谷显著性检验。相对低甲基化可形成后续定位问题，功能着丝粒及精确边界仍需要独立定位证据。')
 add('## 交付与复算入口')
 add('图集共8份、212页：五代表区域各13页，以及49阵列三种方案各49页。figures/index.tsv逐页列出区域和参数。results/summary.xlsx提供摘要表；全部窗口、候选、缺口与模型敏感性表保存在results/；每个样本的有效观测账本与不同读段ID可在per_run/复核。inputs/plan.md、parameters.json、environment.json及source_provenance.json记录计划、实际参数和来源。详细运行方法见包根目录RUNNING.md。')
 add(f'本地结果目录：{D.resolve()}。独立ZIP在同级output目录。参考/BAM不重复打入包；包内有效观测账本足以重做窗口统计及出图。重放BAM解码需要RUNNING.md列出的原输入。')
 text='\n\n'.join(report)+'\n';(D/'report/comparison_report.md').write_text(text)
 # Simple embedded-font PDF. Markdown table rows are rendered as actual tables.
 pdfmetrics.registerFont(TTFont('CJK','/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf'))
 pdfmetrics.registerFont(TTFont('Latin','/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
 def fmt(s):
  chunks=re.split(r'([\u3400-\u9fff\u3000-\u303f\uff00-\uffef]+)',s)
  return ''.join('<font name="'+('CJK' if re.match(r'[\u3400-\u9fff\u3000-\u303f\uff00-\uffef]',c) else 'Latin')+'">'+html.escape(c)+'</font>' for c in chunks if c)
 normal=ParagraphStyle('normal',fontName='CJK',fontSize=9,leading=14,spaceAfter=7,wordWrap='CJK');styles={1:ParagraphStyle('h1',parent=normal,fontSize=16,leading=22,spaceBefore=8),2:ParagraphStyle('h2',parent=normal,fontSize=12,leading=18,spaceBefore=12),3:ParagraphStyle('h3',parent=normal,fontSize=10.5,leading=16,spaceBefore=10)}
 story=[]
 for block in report:
  if block.startswith('|'):
   lines=block.splitlines();data=[[Paragraph(fmt(x.strip()),ParagraphStyle('cell',parent=normal,fontSize=7.5,leading=11,spaceAfter=0)) for x in l.strip('|').split('|')] for l in [lines[0]]+lines[2:]];t=Table(data,colWidths=[480/len(data[0])]*len(data[0]),repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e6eef3')),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.5,colors.grey),('BOTTOMPADDING',(0,0),(-1,-1),5)]));story.extend([t,Spacer(1,8)])
  else:
   lev=len(block)-len(block.lstrip('#')) if block.startswith('#') else 0;story.append(Paragraph(fmt(block.lstrip('#').strip()).replace('\n','<br/>'),styles.get(lev,normal)))
 def footer(canvas,doc):canvas.setFont('Latin',8);canvas.drawRightString(540,25,str(doc.page))
 SimpleDocTemplate(str(D/'report/comparison_report.pdf'),pagesize=A4,leftMargin=55,rightMargin=55,topMargin=42,bottomMargin=42).build(story,onFirstPage=footer,onLaterPages=footer)
 # Human-readable summaries; detailed large ledgers remain in TSV.gz.
 wb=Workbook();wb.remove(wb.active)
 frames={'代表区域判读':pd.DataFrame(jr),'参数对照':c,'同坐标比较':e,'样本均值与覆盖':b,'局部搜索范围':pd.read_csv(R/'local_search_denominators.tsv',sep='\t'),'最低支持窗口':pd.read_csv(R/'minimum_20kb_supported_windows.tsv',sep='\t'),'参考CpG背景':ref,'文件索引':pd.read_csv(D/'figures/index.tsv',sep='\t')}
 for name,frame in frames.items():
  ws=wb.create_sheet(name);ws.append(list(frame.columns));ws.freeze_panes='B2';ws.auto_filter.ref=f'A1:{__import__("openpyxl").utils.get_column_letter(len(frame.columns))}{len(frame)+1}'
  for row in frame.itertuples(index=False,name=None):ws.append([None if pd.isna(x) else x for x in row])
  for row in ws:
   for cell in row:cell.font=Font(name='宋体' if isinstance(cell.value,str) else 'Times New Roman',size=10);cell.alignment=Alignment(vertical='top',wrap_text=True)
  for cell in ws[1]:cell.fill=PatternFill('solid',fgColor='DDEBF7');cell.font=Font(name='宋体',size=10,bold=True)
  for col in ws.columns:
   width=60 if col[0].value=='interpretation' else 28 if col[0].value in ['array_id','name'] else 18;ws.column_dimensions[col[0].column_letter].width=width
  ws.row_dimensions[1].height=45
  for row in ws.iter_rows(min_row=2):
   lines=max(math.ceil(sum(2 if ord(x)>127 else 1 for x in str(cell.value or ''))/(ws.column_dimensions[cell.column_letter].width-2)) for cell in row);ws.row_dimensions[row[0].row].height=max(30,lines*15+6)
 wb.save(R/'summary.xlsx');load_workbook(R/'summary.xlsx').close()
if __name__=='__main__':main()
