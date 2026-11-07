"""Put shared figure explanations on one cover page for each atlas."""
from pathlib import Path
import json
import re
import shutil
from reportlab.platypus import SimpleDocTemplate, Paragraph
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.pagesizes import A4
from pypdf import PdfReader, PdfWriter
from pdf_image_compat import repair_pdf, inspect_predictors

P=Path(__file__).resolve().parents[2]
F=P/'figures'
pdfmetrics.registerFont(TTFont('LegendCJK','/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf'))
pdfmetrics.registerFont(TTFont('LegendLatin','/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
def mixed(t):
    t=re.sub(r'(?<=\d) (?=(?:kb|bp)\b)','\u00a0',t)
    return re.sub(r'[\x20-\x7e\u2013\u2014]+',lambda m:'<font name="LegendLatin">'+m[0]+'</font>',t)

def main():
    items=[(4,'斑马鱼逐阵列结构与甲基化','Supplementary_FigS4_TU_array_panels.pdf','Supplementary_TU_array_structure_methylation.pdf',25),
           (5,'D. aesculapii逐阵列结构','Supplementary_FigS5_aesculapii_structure.pdf','Supplementary_aesculapii_structure.pdf',34),
           (6,'D. kyathit逐阵列结构','Supplementary_FigS6_kyathit_structure.pdf','Supplementary_kyathit_structure.pdf',25)]
    (P/'checks').mkdir(exist_ok=True)
    checks=[]
    for n,title,name,alias,pages in items:
        import sys
        chosen={'tu':4,'aesculapii':5,'kyathit':6}
        if len(sys.argv)>1 and n!=chosen[sys.argv[1]]:continue
        if n==4:
            paras=[
                '49个大型阵列按染色体分组，各面板含一个阵列及两侧各200 kb；结构图、卫星注释及甲基化轨道共用TU参考坐标。蓝色和橙色注释条分别为富AT和富GC，灰竖线为阵列边界。',
                '甲基化轨道采用20 kb窗口、5 kb步长，纵轴为CpG 5mC比例（%）。每窗至少10个有效CpG且至少3条独立读段时显示黑线；至少5条支持时以蓝线标示同一均值。留白为支持不足，曲线在该处断开，实测0%保留；孤立窗口显示短横线。14—19对应ERR16782314—ERR16782319。',
                '轨道下的C、R色带分别为有效CpG数和窗口独立读段数，分级见页内图例。红标表示窗内有效CpG分布不均：最长无有效CpG区间超过实际窗宽一半，或窗口中央一半缺少有效CpG观测。读段数按整个窗口去重，各CpG深度可不同。',
                '均值为有效参考CpG等权平均，重叠窗口共享观测。定位及调用条件见补充方法15，显示条件见补充方法16，坐标和覆盖见补充表S3、S54、S55。'
            ]
        else:
            count=146 if n==5 else 64
            paras=[
                f'{count}个已注释区域按染色体分组，未定位序列另列。每个区域的背景图和局部放大图使用原始序列坐标；背景范围为两侧各200 kb，局部范围为两侧各2 kb，端部按实际序列长度截取。',
                '蓝线和橙线分别标示富AT和富GC卫星支持区域。面板标明区域编号、家族、序列记录、坐标及长度分层；横轴为参考坐标（Mb），色标为相似度（%）。',
                '完整三角形底长为高的两倍，两个坐标方向保持等比例。背景图采用5 kb窗口，局部图采用250 bp窗口，相似度显示下限均为80%。结构图参数见补充方法8、11，区域分层和支持覆盖见补充表S28。'
            ]
        body=P/'checks'/f'S{n}_body_layout_corrected.pdf'
        alias_pages=len(PdfReader(F/alias).pages)
        if alias_pages==pages:shutil.copy2(F/alias,body)
        else:assert alias_pages==pages+1 and body.exists()
        assert len(PdfReader(body).pages)==pages
        cover=P/'checks'/f'S{n}_legend_page.pdf'
        style=ParagraphStyle('body',fontName='LegendCJK',fontSize=11,leading=20,spaceAfter=16,wordWrap='CJK')
        head=ParagraphStyle('head',parent=style,fontSize=16,leading=25,spaceAfter=30)
        SimpleDocTemplate(str(cover),pagesize=A4,leftMargin=55,rightMargin=55,topMargin=70).build(
            [Paragraph(mixed(f'图S{n}　{title}'),head)]+[Paragraph(mixed(x),style) for x in paras])
        assert len(PdfReader(cover).pages)==1
        writer=PdfWriter();writer.append(PdfReader(cover));writer.append(PdfReader(body))
        writer.add_metadata({'/Title':f'图S{n} {title}','/Author':''})
        with (F/name).open('wb') as out:writer.write(out)
        repair_pdf(F/name)
        assert not inspect_predictors(PdfReader(F/name))
        shutil.copy2(F/name,F/alias)
        checks.append(dict(figure=n,body_pages=pages,legend_pages=1,shared_explanations_on_cover_only=True))
    (P/'checks/supplement_cover_layout.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2)+'\n')

if __name__=='__main__':main()
