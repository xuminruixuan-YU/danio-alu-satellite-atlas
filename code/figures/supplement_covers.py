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
                '49个大型阵列按25条染色体分组，每个面板单独展示一个阵列及两侧各200 kb。上方为ModDotPlot三角自比对图，下方依次为卫星序列注释及六个样本的CpG 5mC比例轨道。上下图使用相同的TU参考坐标。',
                '结构图采用5 kb窗口，相似度显示下限为80%。蓝色和橙色注释条分别表示富AT和富GC卫星支持序列，灰色竖线标示阵列边界。',
                '轨道14—19依次对应ERR16782314—ERR16782319，纵轴为CpG 5mC比例（%）。20 kb窗口、5 kb步长，以有效参考CpG等权平均计算5mC比例。每窗至少含10个有效CpG且至少3条不同读段时以黑线显示，至少5条读段支持的相同均值另以蓝线标示。低于门槛时留白并断线，实测0%保留。孤立有效窗口用短横线表示。',
                'C、R色带分别为窗口有效CpG数和不同读段数，分级见页内图例。红短线提示最长无有效CpG区间超过实际窗宽一半，或中央一半缺少有效CpG。每窗读段数是整个窗口的去重数量，各CpG深度可不同；滑窗之间共享观测。',
                '轨道采用全读段比对中通过质量及已返回高分竞争位置坐标一致性筛选的CpG观测，描述参考坐标上的修饰分布。相邻阵列的侧翼发生重叠时，共同区间在各面板中分别显示。区域坐标沿用补充表S3。区域轨道与独立非目录序列锚段支持的阵列统计采用各自定义的读段集合。定位到阵列和定位到阵列内部具体CpG分别评价，竞争搜索的适用范围见补充方法15。'
            ]
        else:
            count=146 if n==5 else 64
            paras=[
                f'{count}个已注释区域按染色体分组，未定位序列另列。每个区域分别绘制背景自比对图和局部放大图，两幅图均使用该序列记录的原始坐标。',
                '背景图包括区域及两侧各200 kb，采用5 kb窗口；局部图包括两侧各2 kb，采用250 bp窗口。序列记录端部按实际长度截断，相似度显示下限为80%。',
                '三角图按相同参考距离等比绘制，完整三角形底与高之比为2∶1。背景图与局部图各自保留完整坐标范围。',
                '蓝线和橙线分别标示富AT和富GC卫星支持区域。各面板标明图号、家族、序列记录、区域坐标及长度分层；横轴为参考坐标（Mb），色标为相似度（%）。',
                '较大区域与较短区段的划分和支持覆盖量见补充表S28，分析方法见补充方法11、16。'
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
