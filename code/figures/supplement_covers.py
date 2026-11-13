"""Put shared explanations on congener atlas cover pages.

TU tracks and their cover use code/genome_methylation/render_regional_500kb.py.
"""
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
pdfmetrics.registerFont(TTFont('LegendLatinItalic','/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf'))
def mixed(t):
    t=re.sub(r'(?<=\d) (?=(?:kb|bp)\b)','\u00a0',t)
    t=re.sub(r'[\x20-\x7e\u2013\u2014]+',lambda m:'<font name="LegendLatin">'+m[0]+'</font>',t)
    return re.sub(r'Danio (?:aesculapii|kyathit)',lambda m:'<font name="LegendLatinItalic">'+m[0]+'</font>',t)

def main():
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='tu':
        raise SystemExit('TU: run python reproduce.py plot tu for the 500 kb panels and caption page.')
    items=[(5,'Danio aesculapii逐阵列结构','Supplementary_FigS4_aesculapii_structure.pdf','Supplementary_aesculapii_structure.pdf',34),
           (6,'Danio kyathit逐阵列结构','Supplementary_FigS5_kyathit_structure.pdf','Supplementary_kyathit_structure.pdf',25)]
    (P/'checks').mkdir(exist_ok=True)
    checks=[]
    for n,title,name,alias,pages in items:
        import sys
        chosen={'tu':4,'aesculapii':5,'kyathit':6}
        if len(sys.argv)>1 and n!=chosen[sys.argv[1]]:continue
        count=146 if n==5 else 64
        paras=[
            f'{count}个已注释区域按染色体分组，未定位序列另列。每个区域的背景图和局部放大图使用原始序列坐标；背景范围为两侧各200 kb，局部范围为两侧各2 kb，端部按实际序列长度截取。',
            '蓝线和橙线分别标示富AT和富GC卫星支持区域。面板标明区域编号、家族、序列记录、坐标及长度分层；横轴为参考坐标（Mb），色标为相似度（%）。',
            '完整三角形底长为高的两倍，两个坐标方向保持等比例。背景图采用5 kb窗口，局部图采用250 bp窗口，相似度显示下限均为80%。结构图参数见补充方法9、17，区域分层和支持覆盖见补充表S4。'
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
            [Paragraph(mixed(f'图S{({4:7,5:4,6:5}[n])}　{title}'),head)]+[Paragraph(mixed(x),style) for x in paras])
        assert len(PdfReader(cover).pages)==1
        writer=PdfWriter();writer.append(PdfReader(cover));writer.append(PdfReader(body))
        writer.add_metadata({'/Title':f'图S{({4:7,5:4,6:5}[n])} {title}','/Author':''})
        with (F/name).open('wb') as out:writer.write(out)
        repair_pdf(F/name)
        assert not inspect_predictors(PdfReader(F/name))
        shutil.copy2(F/name,F/alias)
        checks.append(dict(figure={4:7,5:4,6:5}[n],body_pages=pages,legend_pages=1,shared_explanations_on_cover_only=True))
    (P/'checks/supplement_cover_layout.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2)+'\n')

if __name__=='__main__':main()
