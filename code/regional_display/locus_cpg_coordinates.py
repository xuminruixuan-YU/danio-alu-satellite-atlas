"""Position qualification for ticket 11; never infer CpG position from flank alone.

Input read positions use the original (forward) read sequence, as in the accepted
ticket 10 observations. CIGAR query positions use the SAM-oriented sequence.
This module only checks coordinates; it neither selects modification outcomes
nor assigns an anchor, family, or operational locus.
"""
from bisect import bisect_right
from dataclasses import dataclass
import re


@dataclass(frozen=True)
class Placement:
    reference: str
    flag: int
    query_length: int
    starts: tuple
    matches: tuple
    unusable_reason: str = ''

    @property
    def reverse(self):
        return bool(self.flag & 16)

    def reference_position(self, query_position):
        if self.unusable_reason or not 0 <= query_position < self.query_length:
            return None
        i = bisect_right(self.starts, query_position) - 1
        if i < 0:
            return None
        qs, qe, rs = self.matches[i]
        return rs + query_position - qs if query_position < qe else None


def compile_placement(record):
    cigar = record['CIGAR']
    terms = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    if ''.join(n + op for n, op in terms) != cigar or not terms:
        raise ValueError('missing or malformed mapped CIGAR')
    flag = int(record['flag'])
    reason = ('unmapped' if flag & 4 else
              'supplementary' if flag & 2048 else
              'hard_clipped_query_coordinates' if 'H' in cigar else '')
    q, r, intervals = 0, int(record['start0']), []
    for n, op in terms:
        n = int(n)
        if n <= 0:
            raise ValueError('nonpositive CIGAR operation')
        if op in 'M=X':
            intervals.append((q, q + n, r))
            q += n
            r += n
        elif op in 'IS':
            q += n
        elif op in 'DN':
            r += n
        # H/P consume neither stored query nor reference. H is excluded above.
    return Placement(record['reference'], flag, q,
                     tuple(x[0] for x in intervals), tuple(intervals), reason)


def qualify_cpg(original_read_C_position0, reference, reference_C_position0,
                primary, high_score_placements):
    """Require both bases at the observed reference CpG in every parent placement.

    This independent replay reuses only the coordinate predicate. It does not
    add a flank-anchor selection. The caller must provide every reported
    AS>=0.8*primary placement, including the primary. The caller
    also verifies actual read/reference CG sequence and confidence/base quality.
    A secondary alignment is allowed only when these exact two positions agree.
    """
    if not high_score_placements or primary not in high_score_placements:
        raise ValueError('primary and complete high-score placement set required')
    c = int(original_read_C_position0)
    target = int(reference_C_position0)
    if c < 0 or c + 1 >= primary.query_length:
        return 'query_CpG_outside_read'
    if primary.reverse:
        qc, qg = primary.query_length - 1 - c, primary.query_length - 2 - c
        expected = (target + 1, target)
    else:
        qc, qg = c, c + 1
        expected = (target, target + 1)
    for placement in high_score_placements:
        if placement.unusable_reason:
            return placement.unusable_reason
        if (placement.reference != reference or placement.reverse != primary.reverse
                or placement.query_length != primary.query_length):
            return 'parent_reference_orientation_or_length_conflict'
        observed = (placement.reference_position(qc), placement.reference_position(qg))
        if None in observed:
            return 'CpG_base_unaligned_in_parent'
        if observed != expected:
            return 'CpG_coordinate_conflict'
    return 'position_concordant'
