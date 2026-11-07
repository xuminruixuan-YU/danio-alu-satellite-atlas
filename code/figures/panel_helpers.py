"""Figure 4: unsmoothed per-sample CpG means in 10 kb bins.

The minimum 10 observed CpGs is a display-only choice, not an inferential
coverage standard. Scientific missing values are kept separate from their plotting placeholders.
"""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/danio-array-revision-mpl')
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors, font_manager
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.ticker import MaxNLocator, ScalarFormatter

P = Path(__file__).resolve().parents[2]
OUT = P / 'figures'
ROOT = P / 'data/array_study'
BIN = 10000
MIN_CPG = 10
RUNS = [f'ERR167823{i}' for i in range(14, 20)]
COL = {'CL1': '#266b91', 'CL2': '#bd612b'}
LAB = {'CL1': '富AT', 'CL2': '富GC'}
font_path = os.environ.get('DANIO_CJK_FONT', '/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf')
if not Path(font_path).is_file():
    raise FileNotFoundError('Set DANIO_CJK_FONT to a Chinese TrueType font file.')
font_manager.fontManager.addfont(font_path)
font_name = font_manager.FontProperties(fname=font_path).get_name()
plt.rcParams.update({'font.family': ['DejaVu Sans', font_name],
                     'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
                     'font.size': 9, 'axes.unicode_minus': False})


def read_rows(path):
    path=Path(path)
    if not path.is_absolute():path=ROOT/path
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def summarize(rows, a, b, size):
    groups = collections.defaultdict(list)
    for row in rows:
        if int(row['valid_calls']) > 0:
            groups[int(row['reference_C_position0']) // size].append(row)
    bins = []
    for k in range(a // size, (b - 1) // size + 1):
        rr = groups[k]
        start, end = max(a, k * size), min(b, (k + 1) * size)
        valid = sum(int(x['valid_calls']) for x in rr)
        n = len(rr)
        fractions = [int(x['m']) / int(x['valid_calls']) for x in rr]
        for x, fraction in zip(rr, fractions):
            assert int(x['valid_calls']) == sum(int(x[c]) for c in ('C', 'm', 'h'))
            assert abs(fraction - float(x['five_mC_call_fraction'])) < 1e-10
        positions = sorted(int(x['reference_C_position0']) for x in rr)
        assert len(set(positions)) == n
        bins.append(dict(bin_index=k, start0=start, end0=end,
                         midpoint0=(start + end) / 2, bin_size_bp=size,
                         observed_CpGs=n, valid_observations=valid,
                         five_mC_observations=sum(int(x['m']) for x in rr),
                         mean_CpG_5mC_pct=float(np.mean(fractions) * 100) if n else None,
                         call_weighted_5mC_pct=(100 * sum(int(x['m']) for x in rr) / valid) if valid else None,
                         observed_span_bp=positions[-1] - positions[0] + 1 if n else 0,
                         max_between_observed_CpG_gap_bp=max(np.diff(positions), default=0) if n else 0,
                         status='display' if n >= MIN_CPG else 'low_CpG_count' if n else 'no_observation'))
    assert sum(x['observed_CpGs'] for x in bins) == sum(int(r['valid_calls']) > 0 for r in rows)
    return bins


def panel(fig, spec, meta, tracks, members):
    arr = meta['array']
    a, b = int(arr['plot_start0']), int(arr['plot_end0'])
    gs = spec.subgridspec(8, 1, height_ratios=[3, .32, .55, .55, .55, .55, .55, .55], hspace=.14)
    ax = fig.add_subplot(gs[0])
    polygons, values = [], []
    for row in read_rows(meta['BEDPE']):
        x0 = (a + int(row['query_start']) - 1) / 1e6
        x1 = (a + min(int(row['query_end']), b - a)) / 1e6
        y0 = (a + int(row['reference_start']) - 1) / 1e6
        y1 = (a + min(int(row['reference_end']), b - a)) / 1e6
        polygons.append([((x+y)/2, abs(x-y)/2) for x,y in [(x0,y0),(x1,y0),(x1,y1),(x0,y1)]])
        values.append(float(row['perID_by_events']))
    heat = PolyCollection(polygons, array=np.asarray(values), cmap='viridis',
                          norm=colors.Normalize(80, 100), edgecolors='none', rasterized=True)
    ax.add_collection(heat)
    ax.set_xlim(a/1e6, b/1e6)
    ax.set_ylim(0, (b-a)/2e6)
    ax.set_aspect('auto')
    ax.set_yticks([])
    ax.tick_params(axis='x', labelbottom=False, bottom=False)
    ax.spines[['left','right','top']].set_visible(False)
    ax.set_title(f"{meta['panel_label']}  第{arr['chromosome']}号染色体 · {LAB[arr['family']]}阵列\n"
                 f"{int(arr['start0'])+1}–{arr['end0']} bp", fontsize=meta.get('title_fontsize',10), pad=8)
    band = fig.add_subplot(gs[1], sharex=ax)
    band.set_ylim(-.5,1.5)
    band.set_yticks([0,1], ['GC','AT'], fontsize=7)
    for row in members:
        if row['seqid'] == arr['seqid'] and int(row['start0']) < b and int(row['end0']) > a:
            band.plot([max(a,int(row['start0']))/1e6,min(b,int(row['end0']))/1e6],
                      [1 if row['family']=='CL1' else 0]*2, color=COL[row['family']], lw=3, solid_capstyle='butt')
    band.tick_params(axis='x', labelbottom=False, bottom=False)
    band.spines[['left','right','top','bottom']].set_visible(False)
    axes = [ax, band]
    panel_info = []
    total = None
    for j, run in enumerate(RUNS):
        track = fig.add_subplot(gs[j+2], sharex=ax)
        axes.append(track)
        bins = tracks[arr['array_id'], run]
        x = np.array([r['midpoint0']/1e6 for r in bins])
        y = np.array([r['mean_CpG_5mC_pct'] if r['status']=='display' else np.nan for r in bins])
        # Missing/insufficient bins remain NaN in scientific data. Their zero
        # drawing value is a solid-line zero placeholder, never a 0% call.
        display_y = np.nan_to_num(y, nan=0.0)
        xx = np.r_[a/1e6, x, b/1e6]
        yy = np.r_[display_y[0], display_y, display_y[-1]]
        ok = np.r_[np.isfinite(y[0]), np.isfinite(y), np.isfinite(y[-1])]
        track.plot(xx, yy, color='#202020', linewidth=.95, linestyle='-', marker=None)
        # A finite isolated bin gets a short horizontal stroke, not a point.
        for k,row in enumerate(bins):
            if np.isfinite(y[k]) and (k==0 or not np.isfinite(y[k-1])) and (k==len(y)-1 or not np.isfinite(y[k+1])):
                track.plot([row['start0']/1e6,row['end0']/1e6],[y[k],y[k]],color='#202020',lw=.95,marker=None)

        track.set_ylim(-7,107)
        track.set_yticks([0,100], ['0','100'], fontsize=6)
        track.set_ylabel(run[-2:], rotation=0, ha='right', va='center', labelpad=8, fontsize=8)
        track.tick_params(axis='x', labelbottom=False, bottom=False)
        track.spines[['right','top','bottom']].set_visible(False)
        track.axhline(100,color='#e9e9e9',lw=.4,zorder=0)
        for edge in [int(arr['start0']), int(arr['end0'])]:
            track.axvline(edge/1e6, color='#999999', lw=.5, ls='-')
        if not np.any(np.isfinite(y)):
            pass
        calls = np.array([r['valid_observations'] for r in bins])
        total = calls if total is None else total + calls
        # Matplotlib's finite path segments may only join adjacent grid bins.
        segments = [(bins[k]['bin_index'],bins[k+1]['bin_index'])
                    for k in range(len(bins)-1) if np.isfinite(y[k]) and np.isfinite(y[k+1])]
        assert all(right-left == 1 for left,right in segments)
        panel_info.append(dict(run=run, total_bins=len(bins), displayed_bins=int(np.isfinite(y).sum()),
                               line_segments=len(segments), observed_CpGs=sum(r['observed_CpGs'] for r in bins),
                               valid_observations=int(calls.sum())))
    track.set_xlabel('TU参考坐标（Mb）', fontsize=9)
    track.xaxis.set_major_locator(MaxNLocator(nbins=4))
    fmt = ScalarFormatter(useOffset=False); fmt.set_scientific(False)
    track.xaxis.set_major_formatter(fmt)
    track.tick_params(axis='x', labelbottom=True, bottom=True, labelsize=7)
    track.spines['bottom'].set_visible(True)
    fig.canvas.draw()
    for other in axes:
        assert other.get_xlim() == ax.get_xlim()
        assert np.allclose([other.get_position().x0,other.get_position().x1],
                           [ax.get_position().x0,ax.get_position().x1],atol=1e-6)
    return heat, dict(array_id=arr['array_id'], matrix_cells=len(values),
                     shared_limits0=[a,b], samples=panel_info)


