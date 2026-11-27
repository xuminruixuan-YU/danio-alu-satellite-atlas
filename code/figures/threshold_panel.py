"""Plot the chromosome counts for six TU array-merging settings."""

import csv
from collections import Counter

from matplotlib.ticker import MultipleLocator


SETTINGS = [(10, 50), (25, 25), (25, 50), (25, 100), (25, 200), (50, 50)]
COLORS = {0: '#B7BEC4', 1: '#286887', 2: '#8DBACB', 3: '#CDC3DD'}


def plot_threshold_panel(ax, source, family, primary, *, title_loc='left'):
    by_chromosome = {setting: {} for setting in SETTINGS}
    with source.open(newline='') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            setting = int(row['maximum_gap_kb']), int(row['minimum_supported_kb'])
            if row['family'] != family or setting not in by_chromosome:
                continue
            chromosome = int(row['chromosome'])
            if chromosome in by_chromosome[setting]:
                raise ValueError(f'Duplicate chromosome {chromosome} for setting {setting}')
            by_chromosome[setting][chromosome] = int(row['regions'])

    for setting, chromosomes in by_chromosome.items():
        if set(chromosomes) != set(range(1, 26)):
            raise ValueError(f'Expected 25 TU chromosomes for setting {setting}')

    distribution = [Counter(by_chromosome[setting].values()) for setting in SETTINGS]
    categories = sorted({n for counts in distribution for n in counts})
    width = .24 if len(categories) == 3 else .29
    for index, n in enumerate(categories):
        offset = (index - (len(categories) - 1) / 2) * width
        positions = [i + offset for i in range(len(SETTINGS))]
        values = [counts[n] for counts in distribution]
        ax.bar(positions, values, width * .9, label=f'{n}个阵列',
               color=COLORS.get(n), edgecolor='#324F61', linewidth=.65, zorder=3)
        for pos, value in zip(positions, values):
            if value:
                ax.annotate(str(value), (pos, value), xytext=(0, 2),
                            textcoords='offset points', ha='center', va='bottom', fontsize=11)

    labels = [f'{gap}/{support}' for gap, support in SETTINGS]
    primary_index = SETTINGS.index(primary)
    labels[primary_index] += '\n（主分析）'
    ax.set_xticks(range(len(SETTINGS)), labels, fontsize=10)
    ax.get_xticklabels()[primary_index].set_color('#234C68')
    ax.set_xlim(-.52, len(SETTINGS) - .48)
    ax.set_ylim(0, 29.5)
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.set_xlabel('最大间隔 / 最低卫星覆盖（kb）', fontsize=10)
    ax.set_ylabel('染色体数', fontsize=10)
    ax.set_title('(c)  不同划分条件下富AT大型阵列的染色体分布', loc=title_loc, pad=20)
    ax.tick_params(axis='x', length=0, pad=3)
    ax.tick_params(axis='y', labelsize=10, pad=2, length=3)
    ax.grid(axis='y', color='#DCE3E7', linewidth=.55, zorder=0)
    ax.legend(loc='lower right', bbox_to_anchor=(1, 1.015), ncol=len(categories),
              frameon=False, columnspacing=1.04, fontsize=9)
