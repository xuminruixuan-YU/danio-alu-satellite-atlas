from pathlib import Path
import csv,gzip,collections
from permutation_model import null
p=Path(__file__).resolve().parent
chains=collections.defaultdict(list)
with gzip.open(p/'alternation_segments.tsv.gz','rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):chains[r['array_id']].append(r['order'])
with (p/'array_alternation.tsv').open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if not int(r['permutations']):continue
  c=chains[r['array_id']]
  if r['method']=='20_unit_local_blocks':c=[s[i:i+20] for s in c for i in range(0,len(s),20)]
  nt,nl=null(c,int(r['seed']))
  assert (1+int((nt>=int(r['alternating_pairs'])).sum()))/1000==float(r['p_ge_transitions'])
  assert (1+int((nl>=int(r['longest_alternating_units'])).sum()))/1000==float(r['p_ge_longest'])
print('All 46 testable models reproduce the archived empirical P values.')
