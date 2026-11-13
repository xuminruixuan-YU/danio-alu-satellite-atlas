"""Summarize TU array thresholds and paired methylation differences from deposited tables."""

import argparse
from functools import lru_cache
from itertools import product
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
CATALOG = ROOT / 'data/array_study/array_catalog'
REGIONAL = ROOT / 'data/regional_display'
CONFIG = ROOT / 'metadata/array_comparisons.json'


def read_table(path):
    return pd.read_csv(path, sep='\t')


def summarize_thresholds(combinations, catalogs):
    counts_by_family = []
    counts_by_chromosome = []
    deposited = read_table(CATALOG / 'threshold_sensitivity.tsv')
    for maximum_gap, minimum_coverage in combinations:
        catalog = catalogs[maximum_gap]
        for family in ('CL1', 'CL2'):
            regions = catalog[(catalog.family == family) & (catalog.supported_bp >= minimum_coverage)]
            counts = regions.groupby('chromosome').size().reindex(range(1, 26), fill_value=0)
            parameters = {
                'maximum_gap_kb': maximum_gap // 1000,
                'minimum_supported_kb': minimum_coverage // 1000,
                'family': family,
            }
            supported_bp = int(regions.supported_bp.sum())
            expected = deposited.loc[
                (deposited.species == 'TU') & (deposited.family == family)
                & (deposited.maximum_gap_bp == maximum_gap)
                & (deposited.minimum_supported_bp == minimum_coverage)
            ].iloc[0]
            assert len(regions) == expected.regions and supported_bp == expected.supported_bp
            counts_by_family.append(dict(
                **parameters, regions=len(regions), chromosomes=int((counts > 0).sum()),
                chromosomes_one=int((counts == 1).sum()),
                chromosomes_multiple=int((counts > 1).sum()),
                chromosomes_zero=int((counts == 0).sum()), supported_bp=supported_bp,
            ))
            for chromosome, count in counts.items():
                counts_by_chromosome.append(dict(**parameters, chromosome=int(chromosome), regions=int(count)))
    return pd.DataFrame(counts_by_family), pd.DataFrame(counts_by_chromosome)


def summarize_individuals(config, arrays):
    width = config['window_bp']
    minimum_cpgs = config['minimum_effective_CpGs']
    minimum_reads = config['minimum_independent_reads']
    windows = read_table(ROOT / config['inputs']['window_estimates'])
    old_windows = read_table(ROOT / config['inputs']['baseline_window_check'])
    keys = ['array_id', 'run', 'start0', 'end0']
    shared = old_windows.merge(windows, on=keys, suffixes=('_200kb', '_500kb'), validate='one_to_one')
    for field in ('CpGs', 'observations', 'reads', 'mean'):
        np.testing.assert_allclose(shared[field + '_200kb'], shared[field + '_500kb'], equal_nan=True)
    assert set(windows.array_id) == set(arrays.index)
    assert set(windows.run) == set(config['runs'])
    backgrounds = read_table(ROOT / config['inputs']['cpg_equal_estimates'])
    backgrounds = backgrounds.set_index(['array_id','run','flank_bp','zone'],verify_integrity=True)
    baseline = read_table(REGIONAL / 'results/background.tsv').set_index(['array_id','run'])
    verified = 0
    for (aid,run,flank,zone),row in backgrounds.iterrows():
        if flank != 200000: continue
        source=baseline.loc[(aid,run)]
        for field in ('CpGs','observations','reads','max_gap','mean'):
            np.testing.assert_allclose(row[field],source[zone+'_'+field],rtol=0,atol=1e-9,equal_nan=True)
        verified += 1
    records = []
    flanks = [config['flank_bp']] + config['sensitivity_flank_bp']
    for flank in flanks:
        for (array_id, run), frame in windows.groupby(['array_id','run']):
            array = arrays.loc[array_id]
            start,end = int(array.start0),int(array.end0)
            zones = [('left',start-flank,start),('inside',start,end),('right',end,end+flank)]
            coordinates = dict(array_id=array_id,family=array.family,chromosome=int(array.chromosome),
                               run=run,start0=start,end0=end,flank_bp=flank)
            by_window = dict(**coordinates,estimator='window_equal')
            by_cpg = dict(**coordinates,estimator='reference_CpG_equal')
            for zone,lower,upper in zones:
                complete=frame.loc[(frame.start0>=lower)&(frame.end0<=upper)&((frame.end0-frame.start0)==width)]
                valid=complete.loc[(complete.CpGs>=minimum_cpgs)&(complete.reads>=minimum_reads)&complete['mean'].notna()]
                by_window[zone+'_mean_pct']=valid['mean'].mean()
                by_window[zone+'_total_windows']=len(complete)
                by_window[zone+'_valid_windows']=len(valid)
                by_window[zone+'_flagged_windows']=int(valid.edge_flag.sum())
                stats=backgrounds.loc[(array_id,run,flank,zone)]
                assert int(stats.start0)==lower and int(stats.end0)==upper
                qualified=stats.CpGs>=minimum_cpgs and stats.reads>=minimum_reads
                by_cpg[zone+'_mean_pct']=stats['mean'] if qualified else np.nan
                for source,field in [('CpGs','effective_CpGs'),('observations','observations'),
                                     ('reads','independent_reads'),('max_gap','max_unobserved_bp')]:
                    by_cpg[zone+'_'+field]=int(stats[source])
            records.extend([by_window,by_cpg])
    individuals=pd.DataFrame(records)
    assert len(individuals)==2*len(arrays)*len(config['runs'])*len(flanks)
    for side in ('left','right'):
        individuals[f'inside_minus_{side}_pp']=individuals.inside_mean_pct-individuals[side+'_mean_pct']
    individuals['inside_minus_both_pp']=individuals.inside_mean_pct-(individuals.left_mean_pct+individuals.right_mean_pct)/2
    return individuals,len(shared),verified


@lru_cache(maxsize=None)
def resampling_indices(n):
    return np.array(list(product(range(n), repeat=n)))


def summarize_contrasts(individuals, config):
    rows = []
    for (flank, array_id, estimator), group in individuals.groupby(['flank_bp','array_id', 'estimator'],sort=False):
        for side in ('left', 'right', 'both'):
            values = group[f'inside_minus_{side}_pp'].dropna().to_numpy()
            n = len(values)
            lower, upper = np.nan, np.nan
            draws = 0
            if n >= config['minimum_individuals_for_interval']:
                medians = np.median(values[resampling_indices(n)], axis=1)
                lower, upper = np.quantile(medians, config['interval_quantiles'])
                draws = n ** n
            rows.append(dict(
                array_id=array_id, family=group.family.iloc[0], chromosome=int(group.chromosome.iloc[0]),
                flank_bp=int(flank), estimator=estimator, comparison=side, evaluable_individuals=n,
                median_difference_pp=float(np.median(values)) if n else np.nan,
                ci95_low_pp=lower, ci95_high_pp=upper, individuals_below=int((values < 0).sum()),
                resampling_unit='individual', bootstrap_ordered_resamples=draws,
            ))
    return pd.DataFrame(rows)


def main():
    global ROOT,CATALOG,REGIONAL,CONFIG
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repository',type=Path,default=ROOT)
    parser.add_argument('--config',type=Path,default=None)
    parser.add_argument('--output', type=Path, required=True, help='New or empty directory for the four TSVs')
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        parser.error('Choose a new or empty output directory')

    ROOT=args.repository.resolve();CATALOG=ROOT/'data/array_study/array_catalog';REGIONAL=ROOT/'data/regional_display'
    CONFIG=args.config or ROOT/'metadata/array_comparisons.json'
    config = json.loads(CONFIG.read_text())
    combinations = config['threshold_combinations_bp']
    catalogs = {
        gap: read_table(CATALOG / f'variants/TU_gap{gap}/regions.tsv')
        for gap in sorted({gap for gap, _ in combinations})
    }
    thresholds, chromosomes = summarize_thresholds(combinations, catalogs)
    gap, coverage = config['main_threshold_combination_bp']
    main_catalog = catalogs[gap]
    arrays = main_catalog[main_catalog.supported_bp >= coverage].set_index('array_id', verify_integrity=True)
    assert len(arrays) == 49
    individuals, common_windows, baseline_zones = summarize_individuals(config, arrays)
    contrasts = summarize_contrasts(individuals, config)

    args.output.mkdir(parents=True, exist_ok=True)
    tables = {
        'threshold_summary.tsv': thresholds,
        'threshold_by_chromosome.tsv': chromosomes,
        'methylation_individual_contrasts.tsv': individuals,
        'methylation_contrast_summary.tsv': contrasts,
    }
    for name, table in tables.items():
        table.to_csv(args.output / name, sep='\t', index=False)
    print(json.dumps({'rows': {name: len(table) for name, table in tables.items()},
                      'verified_common_windows': common_windows,
                      'verified_200kb_zones': baseline_zones}))


if __name__ == '__main__':
    main()
