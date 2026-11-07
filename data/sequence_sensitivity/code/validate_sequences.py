"""Validate retained standardized sequences against exact raw unit records."""
import sequence_sensitivity as s
import json
count=0
for p in sorted((s.D/'results').glob('*/unit_audit.tsv.gz')):
 raw={str(i):r for i,r in enumerate(s.read(s.UNIT_ROOT/p.parent.name/'units.tsv.gz'),2)}
 for r in s.read(p):
  if not r['equivalent']:continue
  source=raw[r['source_line']]
  q=source['raw_sequence'];q=s.rc(q) if r['selected_direction']=='-' else q
  seq=r['equivalent']
  assert len(seq)==len(q)==int(r['end0'])-int(r['start0'])
  assert seq in q+q and not set(seq)-set('ACGT')
  count+=1
result={'all_equivalent_sequences_checked_against_source_rows':count,'length_and_exact_cyclic_sequence_preserved':True,'inserted_or_deleted_sequence_by_normalization':False}
(s.D/'checks/sequence_integrity.json').write_text(json.dumps(result,indent=2));print(count)
