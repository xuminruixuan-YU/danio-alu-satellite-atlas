"""Verify deposited data, recalculate array statistics, and redraw array panels.

This entry point uses only files shipped in the repository. Full-genome and
raw-read stages additionally need the public inputs in metadata/public_inputs.tsv.
The analysis_source directory preserves the analysis algorithms; full raw-data
workflows require their external inputs and paths to be configured. This entry
point recalculates deposited array statistics and redraws structural panels.
"""
from pathlib import Path
import argparse
import collections
import csv
import gzip
import itertools
import json
import math
import statistics
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data/array_study'


def rows(path):
    path = Path(path)
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def close(actual, expected):
    assert math.isclose(float(actual), float(expected), rel_tol=1e-9, abs_tol=1e-12), (actual, expected)


def save(name, result):
    print(json.dumps(result, indent=2, ensure_ascii=False))


def verify():
    ids = {}
    tu = json.loads((DATA / 'structure/matrices_checked.json').read_text())['matrices']
    tu_ids = [r['array']['array_id'] for r in tu]
    assert len(tu_ids) == len(set(tu_ids)) == 49
    for row in tu:
        assert (DATA / row['BEDPE']).is_file(), row['BEDPE']
    ids['tu'] = len(tu_ids)
    near = json.loads((ROOT / 'data/figure_source/figures/congener_panel_manifest.json').read_text())
    for species in ('aesculapii', 'kyathit'):
        members = [r['region']['array_id'] for r in near if r['region']['species'] == species]
        assert len(members) == len(set(members))
        ids[species] = len(members)
        for member in members:
            config = json.loads((DATA / 'congener_structure' / member / 'checked.json').read_text())
            for view in config['views']:
                assert (DATA / view['BEDPE']).is_file(), view['BEDPE']
    assert ids == dict(tu=49, aesculapii=146, kyathit=64), ids
    windows = rows(ROOT / 'data/figure_source/methylation_windows.tsv')
    for r in windows:
        valid = r['status'] == 'display'
        assert valid == (int(r['observed_CpGs']) >= 10)
        assert (r['zero_placeholder'] == 'True') == (not valid)
        if valid:
            close(r['plot_5mC_pct'], r['mean_CpG_5mC_pct'])
        else:
            assert float(r['plot_5mC_pct']) == 0
        if not int(r['observed_CpGs']):
            assert r['mean_CpG_5mC_pct'] == ''
    save('integrity.json', dict(status='passed', arrays=ids,
                               methylation_windows=len(windows), missing_values_preserved=True))


def metrics():
    sys.path.insert(0, str(ROOT / 'code'))
    from metrics import total
    sampled = collections.defaultdict(list)
    for row in rows(DATA / 'array_metrics/sampled_units.tsv.gz'):
        if row['variant'] == 'gap25000_size50000':
            sampled[row['array_id'], row['category']].append(row['sequence'])
    tu = [r for r in rows(DATA / 'array_metrics/within_array.tsv')
          if r['variant'] == 'gap25000_size50000' and r['category'] in ('CL1_all', 'CL2_all')]
    assert len(tu) == 49
    recalculated = collections.defaultdict(list)
    for r in tu:
        out = total(itertools.combinations(sampled[r['array_id'], r['category']], 2))
        for key in ('edit_bases', 'max_length_bp', 'pairs', 'normalized_edit_distance', 'event_sensitivity'):
            close(out[key], r[key])
        recalculated[r['family']].append(out['normalized_edit_distance'])
    for r in rows(DATA / 'array_metrics/family_summary.tsv'):
        if r['variant'] == 'gap25000_size50000' and r['category'] in ('CL1_all', 'CL2_all'):
            close(statistics.mean(recalculated[r['family']]), r['array_equal_edit_distance'])
    congeners = [r for r in rows(DATA / 'congener_metrics/within_arrays.tsv')
                 if r['panel'] == 'fixed' and r['copy_screen'] == 'matched_TU']
    near = collections.defaultdict(list)
    checked_near = 0
    for r in congeners:
        seqs = [x['sequence'] for x in rows(DATA / 'congener_metrics' / r['array_id'] / 'fixed_matched_TU_sample.tsv')]
        out = total(itertools.combinations(seqs, 2))
        close(out['pairs'], r['pairs'])
        if out['pairs']:
            close(out['normalized_edit_distance'], r['normalized_edit_distance'])
            checked_near += 1
            if int(r['supported_bp']) >= 10000:
                near[r['species'], r['family']].append(out['normalized_edit_distance'])
    for r in rows(DATA / 'congener_metrics/family_summary.tsv'):
        if r['panel'] == 'fixed' and r['copy_screen'] == 'matched_TU' and r['minimum_supported_bp'] == '10000':
            vals = near[r['species'], r['family']]
            assert len(vals) == int(r['comparable_arrays'])
            if vals:
                close(statistics.mean(vals), r['array_equal_edit_distance'])
    methyl = collections.defaultdict(list)
    for r in rows(DATA / 'array_metrics/CpG_methylation.tsv.gz'):
        assert int(r['valid_calls']) == sum(int(r[k]) for k in ('C', 'm', 'h'))
        if int(r['valid_calls']) > 0:
            methyl[r['run'], r['array_id']].append(r)
    family = collections.defaultdict(list)
    methyl_checked = 0
    for r in rows(DATA / 'array_metrics/array_methylation.tsv'):
        rr = methyl[r['run'], r['array_id']]
        assert len(rr) == int(r['callable_CpGs'])
        assert len({(x['reference'], x['reference_C_position0']) for x in rr}) == len(rr)
        depth = sum(int(x['valid_calls']) for x in rr)
        assert depth == int(r['valid_observations'])
        if rr:
            mean = statistics.mean(int(x['m']) / int(x['valid_calls']) for x in rr)
            close(mean, r['CpG_equal_5mC'])
            close(sum(int(x['m']) for x in rr) / depth, r['observation_weighted_5mC'])
            family[r['run'], r['family']].append(mean)
            methyl_checked += 1
    for r in rows(DATA / 'array_metrics/array_methylation_family_summary.tsv'):
        vals = family[r['run'], r['family']]
        assert len(vals) == int(r['arrays_with_CpGs'])
        close(statistics.mean(vals), r['array_equal_5mC'])
    save('statistics.json', dict(status='passed', TU_primary_arrays=len(tu),
         TU_array_equal_edit_distance={k: statistics.mean(v) for k,v in recalculated.items()},
         congener_arrays_with_pairs=checked_near, methylation_array_sample_pairs=methyl_checked,
         scope='Recalculation from deposited sampled monomers and valid CpG calls; no raw-read realignment.'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('verify', help='Check matrix inputs, panel membership, and missing-value display flags.')
    sub.add_parser('metrics', help='Recalculate sampled-monomer edit distances and array methylation summaries.')
    plot = sub.add_parser('plot', help='Redraw grouped figures using deposited matrices and CpG calls; writes figures/.')
    plot.add_argument('part', choices=['tu', 'aesculapii', 'kyathit'])
    args = parser.parse_args()
    if args.command == 'verify':
        verify()
    elif args.command == 'metrics':
        metrics()
    else:
        subprocess.run([sys.executable, str(ROOT / 'code/figures/render_groups.py'), args.part], cwd=ROOT, check=True)
        subprocess.run([sys.executable, str(ROOT / 'code/figures/supplement_covers.py'), args.part], cwd=ROOT, check=True)


if __name__ == '__main__':
    main()
