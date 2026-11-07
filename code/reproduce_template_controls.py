"""Recompute S50 from saved per-core coverage counts (Python standard library)."""
from pathlib import Path
import csv,collections
R=Path(__file__).resolve().parents[1];D=R/'data/template_controls'
tot=collections.Counter()
for p in D.glob('*/*/family_bases.tsv'):
 for r in csv.DictReader(p.open(),delimiter='\t'):
  tot[(p.parts[-3],r['comparison'],r['family'],r['cutoff'],r['method'],r['coverage_class'])]+=int(r['covered_bp'])
for r in csv.DictReader((D/'coverage_summary.tsv').open(),delimiter='\t'):
 a,b=r['comparison'].split('_');key=(r['assembly'],r['comparison'],r['family'],r['cutoff'])
 assert tot[key+(a,'common')]==tot[key+(b,'common')]==int(r['common_bp'])
 assert tot[key+(a,'unique')]==int(r['baseline_only_bp'])
 assert tot[key+(b,'unique')]==int(r['control_only_bp'])
print('All S50 shared and scheme-specific coverage counts reproduced.')
