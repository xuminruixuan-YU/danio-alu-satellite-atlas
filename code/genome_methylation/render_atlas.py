"""Draw standalone chromosome atlases with a shared first-page reading guide."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[2]
SURVEY = BASE
DATA = BASE/'data/genome_methylation'
OUT = Path(os.environ.get('DANIO_ATLAS_OUTPUT', BASE/'figures/genome_methylation'))
OUT.mkdir(parents=True,exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(OUT/'tmp/matplotlib'))
sys.path.insert(0, str(SURVEY/'code'))

from plot_common import PdfPages, choose_font, plt, read_windows
from plot_common import now, save
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

FAMILIES = [('CL1', '富 AT', 'AT-rich', '#287da7'), ('CL2', '富 GC', 'GC-rich', '#bd622d')]
ORANGE = '#d07a28'
RED = '#b42b26'
ALPHA = .75
MAX_SHADE_GAP_BP = 50_000
SIZES = [200]
SIZE = (11.7, 12.7)
COLUMNS = ['reference', 'run', 'window_bp', 'start0', 'end0', 'display',
           'CpG_equal_5mC', 'observation_weighted_5mC', 'spatially_concentrated',
           'callable_reference_fraction']


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def swatch(fig, x, y, color, alpha):
    fig.add_artist(Rectangle((x,y-.009), .038,.018, transform=fig.transFigure,
                             facecolor=color, edgecolor='none', alpha=alpha))


def line_key(fig, y, color, style, marker):
    fig.add_artist(Line2D([.10,.174], [y,y], color=color, linestyle=style, lw=1.6,
                         transform=fig.transFigure))
    fig.add_artist(Line2D([.203], [y], color=color, linestyle='none', marker=marker,
                         markersize=5, transform=fig.transFigure))



def shade_intervals(intervals, opposite_intervals):
    """Join short interruptions for the background layer."""
    merged = []
    for start, end in sorted(intervals):
        if merged:
            left, right = merged[-1]
            gap = start - right
            if gap <= 0 or (gap <= MAX_SHADE_GAP_BP and
                            not any(a < start and b > right for a, b in opposite_intervals)):
                merged[-1] = (left, max(right, end))
                continue
        merged.append((start, end))
    return merged


def spans_mb(intervals):
    return [(start/1e6, (end-start)/1e6) for start, end in intervals]


def chromosome_page(meta, groups, runs, fragments):
    fig = plt.figure(figsize=SIZE)
    grid = fig.add_gridspec(8,1,height_ratios=[.45,1,1,1,1,1,1,1.25],hspace=.3)
    top = fig.add_subplot(grid[0])
    top.set_ylim(-.2,1.9)
    intervals = {}
    for y, (family, _, _, color) in zip([1,0], FAMILIES):
        f = fragments[fragments.family == family]
        assert (f.start0 >= 0).all() and (f.end0 <= meta['length']).all() and (f.end0 > f.start0).all()
        intervals[family] = list(zip(f.start0, f.end0))
        top.broken_barh(spans_mb(intervals[family]),(y,.6),facecolors=color,linewidths=0)
    shades = {family: shade_intervals(bounds, intervals['CL2' if family=='CL1' else 'CL1'])
              for family, bounds in intervals.items()}
    shade_spans = {family: spans_mb(bounds) for family, bounds in shades.items()}
    top.set_yticks([1.3,.3],['富 AT','富 GC'])
    top.set_xlim(0,meta['length']/1e6)
    top.tick_params(axis='x',which='both',bottom=False,labelbottom=False)
    top.set_title(f'第{meta["chromosome"]}号染色体',fontsize=14,pad=11)

    coverage, track_checks = [], []
    edges = None
    for i, run in enumerate(runs):
        q = groups.get_group(run).sort_values('start0')
        ax = fig.add_subplot(grid[i+1],sharex=top)
        for family, _, _, color in FAMILIES:
            overlay = ax.broken_barh(shade_spans[family],(0,1),transform=ax.get_xaxis_transform(),
                                     facecolors=color,edgecolors='none',alpha=ALPHA,linewidths=0,zorder=0)
            assert len(overlay.get_paths()) == len(shades[family])
        starts, ends = q.start0.to_numpy(), q.end0.to_numpy()
        assert starts[0] == 0 and ends[-1] == meta['length'] and np.array_equal(starts[1:],ends[:-1])
        boundaries = np.r_[starts,ends[-1]]/1e6
        if edges is None:
            edges = boundaries
        else:
            np.testing.assert_array_equal(edges,boundaries)
        x = (starts+ends)/2e6
        ok = q.display.to_numpy(dtype=bool)
        equal = np.where(ok,q.CpG_equal_5mC.to_numpy()*100,np.nan)
        weighted = np.where(ok,q.observation_weighted_5mC.to_numpy()*100,np.nan)
        assert np.all(np.isfinite(equal[ok])) and np.all(np.isfinite(weighted[ok]))
        ax.plot(x,equal,color='black',lw=1)
        ax.plot(x,weighted,color=ORANGE,lw=.8,ls='--')
        lone = ok & ~np.r_[False,ok[:-1]] & ~np.r_[ok[1:],False]
        ax.scatter(x[lone],equal[lone],s=8,color='black',zorder=4)
        ax.scatter(x[lone],weighted[lone],s=11,color=ORANGE,marker='x',zorder=5)
        concentrated = ok & q.spatially_concentrated.to_numpy(dtype=bool)
        ax.plot(x[concentrated],np.full(concentrated.sum(),98),marker='|',ls='none',color=RED,ms=4)
        ax.set_ylim(-2,102)
        ax.set_yticks([0,50,100])
        ax.set_ylabel(run[-2:]+'\n5mC (%)',rotation=0,labelpad=25,va='center')
        ax.tick_params(labelbottom=False)
        ax.grid(axis='y',alpha=.15)
        coverage.append(q.callable_reference_fraction.to_numpy()*100)
        track_checks.append(dict(run=run,windows=len(q),displayed=int(ok.sum()),
                                 isolated=int(lone.sum()),red_marked=int(concentrated.sum())))
    bottom = fig.add_subplot(grid[-1],sharex=top)
    im = bottom.pcolormesh(edges,np.arange(len(runs)+1)-.5,np.asarray(coverage),
                           vmin=0,vmax=100,cmap='viridis',shading='flat')
    bottom.set_ylim(len(runs)-.5,-.5)
    bottom.set_yticks(range(len(runs)),[r[-2:] for r in runs])
    bottom.set_ylabel('可测 CpG\n占参考（%）')
    bottom.set_xlabel('TU 参考坐标（Mb）')
    cax = bottom.inset_axes([1.025,0,.025,1])
    fig.colorbar(im,cax=cax).set_label('有效 CpG／参考 CpG（%）')
    fig.subplots_adjust(left=.12,right=.90,top=.94,bottom=.072)
    return fig, {k:len(v) for k,v in intervals.items()}, shades, track_checks


def render(window_kb, config, windows, fragments):
    original = SURVEY/'figures'/f'chromosome_atlas_{window_kb}kb.pdf'
    baseline = sha(original) if original.is_file() else None
    output = OUT/f'chromosome_atlas_{window_kb}kb_annotated.pdf'
    pending = output.with_name('.'+output.name)
    runs = [item['run'] for item in config['items']]
    pages = []
    with PdfPages(pending,metadata={'Title':f'TL 前脑染色体 5mC 图集 · {window_kb} kb'}) as pdf:
        for meta in config['reference_records']:
            if not meta['nuclear']:
                continue
            ref = meta['reference']
            fig, counts, shades, tracks = chromosome_page(meta, windows[windows.reference==ref].groupby('run'),
                                                           runs, fragments[fragments.seqid==ref])
            pdf.savefig(fig)
            plt.close(fig)
            pages.append(dict(chromosome=meta['chromosome'],fragments=counts,shade_intervals_bp=shades,tracks=tracks))
    assert len(pages) == 25
    assert sum(sum(p['fragments'].values()) for p in pages) == len(fragments) == 376
    assert baseline is None or sha(original) == baseline
    from cover_pages import append_cover
    append_cover(pending,output,9)
    pending.unlink()
    save(OUT/f'figure_{window_kb}kb_checked.json',dict(
        passed=True, at_HKT=now(), filename=output.name, sha256=sha(output),
        original_sha256=baseline, shade_alpha=ALPHA, max_shade_gap_bp=MAX_SHADE_GAP_BP,
        guide_page=1, chromosome_pages=pages))
    print(window_kb,'kb',output.name,output.stat().st_size,flush=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('kb',nargs='+',type=int,choices=SIZES)
    args=parser.parse_args()
    config=json.loads((DATA/'parameters.json').read_text())
    assert config['assembly']=='GCF_049306965.1' and len(config['items'])==6
    windows=read_windows(config)
    fragments=pd.read_csv(DATA/'fragments.tsv',sep='\t')
    fragments=fragments[fragments.fragment.str.startswith('TU.')]
    choose_font()
    for kb in args.kb:render(kb,config,windows[windows.window_bp==kb*1000][COLUMNS],fragments)

if __name__=='__main__':main()
