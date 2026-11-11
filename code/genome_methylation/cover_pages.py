"""Shared A4 caption-page format for supplementary atlases."""
from pathlib import Path
import re
from reportlab.platypus import SimpleDocTemplate, Paragraph
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.pagesizes import A4
from pypdf import PdfReader, PdfWriter

pdfmetrics.registerFont(TTFont('LegendCJK','/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf'))
pdfmetrics.registerFont(TTFont('LegendLatin','/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
def mixed(t):
    t=re.sub(r'(?<=\d) (?=(?:kb|bp)\b)','\u00a0',t)
    return re.sub(r'[\x20-\x7e\u2013\u2014]+',lambda m:'<font name="LegendLatin">'+m[0]+'</font>',t)

def append_cover(body,out,kind):
    if kind==4:
        title='图S4　斑马鱼逐阵列结构与甲基化'
        paras=[
         '49个大型阵列按染色体分组，各面板含一个阵列及两侧各500 kb；结构图、卫星注释及甲基化轨道共用TU GCF_049306965.1参考坐标。蓝色和橙色注释条分别为富AT和富GC，灰竖线为阵列边界。',
         '六条轨道来自TL成年斑马鱼前脑样本，14—19对应ERR16782314—ERR16782319。甲基化采用20 kb窗口、5 kb步长，纵轴为CpG 5mC比例（%）。每窗至少10个有效CpG且至少3条独立读段时显示黑线；至少5条支持时以蓝线标示同一均值。留白为支持不足，曲线在该处断开，实测0%保留；孤立窗口显示短横线。',
         '轨道下的C、R色带分别为有效CpG数和窗口独立读段数，分级见页内图例。红标表示窗内有效CpG分布不均：最长无有效CpG区间超过实际窗宽一半，或窗口中央一半缺少有效CpG观测。读段数按整个窗口去重，各CpG深度可不同。',
         '均值为有效参考CpG等权平均，重叠窗口共享观测。定位及调用条件见补充方法15、17，显示条件见补充方法16，500 kb侧翼的逐窗及覆盖统计见补充表S65、S66。'
        ]
    elif kind==9:
        title='图S9　TL前脑全染色体甲基化'
        paras=[
         '25条核染色体分别绘图，横轴为TU GCF_049306965.1参考坐标（Mb）。六条轨道来自TL成年斑马鱼前脑样本，14—19对应ERR16782314—ERR16782319。采用200 kb非重叠窗口，折线位于窗口中心，染色体末端按实际范围截取。纵轴为5mC比例（%）。',
         '顶部蓝色和橙色条带分别标示富AT和富GC卫星片段。折线中的同色底纹不透明度为75%；同家族片段间不超过50 kb且不跨另一家族的间隔在底纹中连接，片段坐标及统计范围保持原样。',
         '黑实线为有效参考CpG等权平均，橙虚线为有效观测加权比例。未修饰C、5mC和5hmC均计入有效观测的分母。每窗至少10个有效CpG且至少3条独立读段时显示；留白为支持不足，曲线在该处断开，实测0%保留。孤立有效窗分别用圆点和叉号表示。',
         '红标表示已显示窗口中，最长无有效CpG间隔超过实际窗宽一半，或窗口中央一半缺少有效CpG。下方色带的各行与样本对应，表示有效CpG占参考CpG的比例，深紫至黄色为低至高；缺少参考CpG的窗口留白。',
         '读段门槛按整个窗口判断。200 kb窗口的读段数界限和必要的精确复核见补充方法17；逐窗显示与覆盖统计见补充表S63，逐染色体统计见S64。'
        ]
    else:raise ValueError(kind)
    out=Path(out);cover=out.with_suffix('.cover.pdf')
    style=ParagraphStyle('body',fontName='LegendCJK',fontSize=11,leading=20,spaceAfter=16,wordWrap='CJK')
    head=ParagraphStyle('head',parent=style,fontSize=16,leading=25,spaceAfter=30)
    SimpleDocTemplate(str(cover),pagesize=A4,leftMargin=55,rightMargin=55,topMargin=70).build(
        [Paragraph(mixed(title),head)]+[Paragraph(mixed(p),style) for p in paras])
    assert len(PdfReader(cover).pages)==1
    writer=PdfWriter();writer.append(PdfReader(cover));writer.append(PdfReader(body))
    writer.add_metadata({'/Title':title})
    with out.open('wb') as f:writer.write(f)
    cover.unlink()
