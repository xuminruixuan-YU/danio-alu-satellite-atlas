from pathlib import Path
import os,json
from datetime import datetime
from zoneinfo import ZoneInfo
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import font_manager
BASE=Path(__file__).resolve().parents[2]
DATA=BASE/'data/genome_methylation'
def now():return datetime.now(ZoneInfo('Asia/Hong_Kong')).isoformat()
def save(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
def choose_font():
    font=Path('/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf')
    if font.exists():font_manager.fontManager.addfont(str(font))
    families={f.name for f in font_manager.fontManager.ttflist}
    for name in ['Noto Sans CJK SC','Source Han Sans CN','SimSun','Noto Sans CJK JP','WenQuanYi Zen Hei','Droid Sans Fallback']:
        if name in families:
            plt.rcParams['font.family']=[name,'DejaVu Sans'];break
    plt.rcParams.update({'axes.unicode_minus':False,'pdf.fonttype':42,'font.size':8})
def read_windows(config):
    return pd.concat([pd.read_csv(DATA/(i['run']+'.windows.tsv.gz'),sep='\t',low_memory=False) for i in config['items']],ignore_index=True)
