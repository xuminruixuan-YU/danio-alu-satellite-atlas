"""sat_alignment.py — Sequence phasing and alignment algorithms.

Contains pure functions for rotating/phasing sequences relative to a reference.
"""
import logging
from collections import Counter
from Bio.Seq import Seq
import numpy as np

logger = logging.getLogger(__name__)

try:
    import edlib
except ImportError:
    edlib = None


IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
_IUPAC_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)


def reverse_complement(seq):
    """Return the reverse complement for the full DNA IUPAC alphabet."""
    return seq.translate(_IUPAC_COMPLEMENT)[::-1]


def validate_dna_sequence(seq, record_id="sequence"):
    """Normalize and validate a non-empty, ungapped DNA IUPAC sequence."""
    normalized = str(seq).upper()
    if not normalized:
        raise ValueError(f"Empty sequence: {record_id}")
    invalid = sorted(set(normalized) - IUPAC_DNA)
    if invalid:
        raise ValueError(
            f"Unsupported DNA characters in {record_id}: {''.join(invalid)}"
        )
    return normalized


def _least_rotation(seq):
    """Lexicographically least rotation using Booth's linear-time algorithm."""
    doubled = seq + seq
    length = len(seq)
    left, right, offset = 0, 1, 0
    while left < length and right < length and offset < length:
        a = doubled[left + offset]
        b = doubled[right + offset]
        if a == b:
            offset += 1
            continue
        if a > b:
            left += offset + 1
            if left == right:
                left += 1
        else:
            right += offset + 1
            if left == right:
                right += 1
        offset = 0
    start = min(left, right)
    return doubled[start:start + length]


def canonical_circular_sequence(seq):
    """Canonicalize a sequence over all rotations of both orientations."""
    seq = validate_dna_sequence(seq)
    return min(_least_rotation(seq), _least_rotation(reverse_complement(seq)))


def _require_edlib():
    if edlib is None:
        raise RuntimeError(
            "satkit cluster requires edlib; install the project's declared dependencies"
        )


def circular_distance(a, b, max_edits=None, _canonical=False):
    """Exact circular/reverse-complement Levenshtein distance.

    When ``max_edits`` is provided, ``None`` means the exact distance exceeds
    that limit.  A returned integer is always the exact distance, never a
    clipped or heuristic value.
    """
    _require_edlib()
    if _canonical:
        if not a or not b:
            raise ValueError("Canonical sequences must be non-empty")
    else:
        a = validate_dna_sequence(a, "left sequence")
        b = validate_dna_sequence(b, "right sequence")
    if max_edits is not None:
        if not isinstance(max_edits, int) or max_edits < 0:
            raise ValueError("max_edits must be a non-negative integer or None")
        if abs(len(a) - len(b)) > max_edits:
            return None
    if a == b or (
        not _canonical
        and canonical_circular_sequence(a) == canonical_circular_sequence(b)
    ):
        return 0

    lower = abs(len(a) - len(b))
    if max_edits is None:
        best = edlib.align(a, b, mode="NW", task="distance")["editDistance"]
    else:
        linear = edlib.align(
            a, b, mode="NW", task="distance", k=max_edits
        )["editDistance"]
        best = linear if linear >= 0 else max_edits + 1

    # HW supplies a safe lower bound and likely phase locations. Every
    # candidate is verified as a one-period NW alignment; remaining rotations
    # are searched if the bound has not been met. Only ``a`` is rotated, which
    # exactly matches the public definition; swapping fixed/moving repeats the
    # same symmetric distance without adding a candidate phase.
    fixed = b
    for oriented in (a, reverse_complement(a)):
        doubled = oriented + oriented
        hw_args = {"mode": "HW", "task": "locations"}
        if max_edits is not None:
            hw_args["k"] = max_edits
        hw = edlib.align(fixed, doubled, **hw_args)
        if hw["editDistance"] < 0:
            continue
        bound = max(lower, hw["editDistance"])
        if best <= bound:
            continue
        checked = set()
        for begin, end in hw.get("locations", ())[:4]:
            for phase in (
                begin % len(oriented),
                (end - len(oriented) + 1) % len(oriented),
            ):
                checked.add(phase)
                candidate = doubled[phase:phase + len(oriented)]
                distance = edlib.align(
                    fixed, candidate, mode="NW", task="distance", k=best - 1
                )["editDistance"]
                if distance >= 0:
                    best = distance
        if best <= bound:
            continue
        for phase in range(len(oriented)):
            if best <= bound:
                break
            if phase in checked:
                continue
            candidate = doubled[phase:phase + len(oriented)]
            distance = edlib.align(
                fixed, candidate, mode="NW", task="distance", k=best - 1
            )["editDistance"]
            if distance >= 0:
                best = distance

    return None if max_edits is not None and best > max_edits else best


def phase_to_reference(query, reference):
    """Return the deterministic best circular/RC phase of ``query``.

    Ties prefer forward orientation, then the smallest rotation offset.
    """
    _require_edlib()
    query = validate_dna_sequence(query, "query sequence")
    reference = validate_dna_sequence(reference, "reference sequence")
    best_distance = circular_distance(query, reference)
    best = None
    for orientation, oriented in (
        ("forward", query),
        ("reverse_complement", reverse_complement(query)),
    ):
        doubled = oriented + oriented
        for rotation in range(len(oriented)):
            phased = doubled[rotation:rotation + len(oriented)]
            distance = edlib.align(
                phased,
                reference,
                mode="NW",
                task="distance",
                k=best_distance,
            )["editDistance"]
            if distance == best_distance:
                key = (0 if orientation == "forward" else 1, rotation)
                if best is None or key < best[0]:
                    best = (key, phased, orientation, rotation)
    if best is None:  # Defensive: circular_distance found a realizable phase.
        raise RuntimeError("Unable to recover an exact best circular phase")
    _, phased, orientation, rotation = best
    return phased, orientation, rotation, best_distance


def _get_best_phase_brute_force_dual(query_seq, ref_seq, aligner):
    """双向暴力穷举定相.

    Scoring backend: edlib edit distance (fast, C-backed) when available,
    falling back to the Biopython PairwiseAligner otherwise. The return
    contract is preserved: (best_seq, norm_score, matched) where norm_score
    is an identity in [0, 1] (higher = better) and matched is always True.
    """
    # Defensive: validate ref_seq
    if ref_seq is None:
        raise ValueError("Reference sequence is None - cannot perform phasing")

    q_len = len(query_seq)
    r_len = len(ref_seq)
    query_rc = str(Seq(query_seq).reverse_complement())

    if edlib is not None:
        # Minimize edit distance across all forward and reverse-complement
        # rotations. Identity = 1 - editDistance / max(len_query, len_ref).
        best_dist = None
        best_seq = None

        for i in range(q_len):
            rot = query_seq[i:] + query_seq[:i]
            dist = edlib.align(rot, ref_seq, mode="NW", task="distance")['editDistance']
            if best_dist is None or dist < best_dist:
                best_dist = dist
                best_seq = rot

        for i in range(q_len):
            rot = query_rc[i:] + query_rc[:i]
            dist = edlib.align(rot, ref_seq, mode="NW", task="distance")['editDistance']
            if best_dist is None or dist < best_dist:
                best_dist = dist
                best_seq = rot

        denom = max(q_len, r_len)
        norm_score = 1.0 - (best_dist / denom) if denom > 0 else 0.0
        if norm_score < 0.0:
            norm_score = 0.0
        if norm_score > 1.0:
            norm_score = 1.0
        return best_seq, norm_score, True

    # Fallback: original Biopython behavior (used only if edlib is missing).
    best_score = -1e9
    best_seq = None

    for i in range(q_len):
        rot = query_seq[i:] + query_seq[:i]
        score = aligner.score(ref_seq, rot)
        if score > best_score:
            best_score = score
            best_seq = rot

    for i in range(q_len):
        rot = query_rc[i:] + query_rc[:i]
        score = aligner.score(ref_seq, rot)
        if score > best_score:
            best_score = score
            best_seq = rot

    max_possible = min(q_len, r_len) * aligner.match_score
    norm_score = best_score / max_possible if max_possible > 0 else 0
    return best_seq, norm_score, True


def _get_best_phase_smart(query_seq, ref_seq, aligner, kmer_len):
    """
    Smart phasing using K-mer guidance to reduce search space.
    Instead of O(L^2), we find matching K-mers to propose offsets.
    Then we only verify those offsets.
    """
    k = kmer_len
    q_len = len(query_seq)

    # 1. Index Reference K-mers
    ref_kmers = {} # kmer -> [positions]
    r_len = len(ref_seq)

    for i in range(r_len - k + 1):
        sub = ref_seq[i:i+k]
        if sub not in ref_kmers:
            ref_kmers[sub] = []
        ref_kmers[sub].append(i)

    # 2. Scan Query K-mers and vote for offsets
    # We also need to check RC of Query

    offsets_fwd = Counter()
    offsets_rc = Counter()

    query_rc = str(Seq(query_seq).reverse_complement())

    # Scan Fwd
    for i in range(q_len - k + 1):
        sub = query_seq[i:i+k]
        if sub in ref_kmers:
            for ref_pos in ref_kmers[sub]:
                rot_idx = (i - ref_pos) % q_len
                offsets_fwd[rot_idx] += 1

    # Scan RC
    for i in range(q_len - k + 1):
        sub = query_rc[i:i+k]
        if sub in ref_kmers:
            for ref_pos in ref_kmers[sub]:
                rot_idx = (i - ref_pos) % q_len
                offsets_rc[rot_idx] += 1

    # 3. Select Top Candidates
    candidates = []

    # Top 3 Fwd
    for offset, count in offsets_fwd.most_common(3):
        candidates.append((offset, 1)) # (offset, orientation=1)

    # Top 3 RC
    for offset, count in offsets_rc.most_common(3):
        candidates.append((offset, -1)) # (offset, orientation=-1)

    # If no candidates found (no shared kmers), fallback to brute-force dual
    if not candidates:
        return _get_best_phase_brute_force_dual(query_seq, ref_seq, aligner)

    matched = True

    # 4. Verify candidates with edlib (fast) or Biopython (fallback)
    # Note: edlib path returns identity = 1 - editDistance/max(q_len, r_len);
    # Biopython fallback returns affine-ratio = score/(min_len * match_score).
    # These scales differ — edlib is the canonical path, Biopython is legacy.
    if edlib is not None:
        best_dist = None
        best_seq = None

        for rot_idx, orient in candidates:
            if orient == 1:
                rot = query_seq[rot_idx:] + query_seq[:rot_idx]
            else:
                rot = query_rc[rot_idx:] + query_rc[:rot_idx]

            dist = edlib.align(rot, ref_seq, mode="NW", task="distance")['editDistance']
            if best_dist is None or dist < best_dist:
                best_dist = dist
                best_seq = rot

        denom = max(q_len, r_len)
        if best_dist is None or best_seq is None:
            # edlib failed or returned no result — fall back to brute-force
            return _get_best_phase_brute_force_dual(query_seq, ref_seq, aligner)
        norm_score = 1.0 - (best_dist / denom) if denom > 0 else 0.0
        if norm_score < 0.0:
            norm_score = 0.0
        if norm_score > 1.0:
            norm_score = 1.0
    else:
        # Fallback: original Biopython verification
        best_score = -1e9
        best_seq = None

        for rot_idx, orient in candidates:
            if orient == 1:
                rot = query_seq[rot_idx:] + query_seq[:rot_idx]
            else:
                rot = query_rc[rot_idx:] + query_rc[:rot_idx]

            score = aligner.score(ref_seq, rot)
            if score > best_score:
                best_score = score
                best_seq = rot

        max_possible = min(q_len, r_len) * aligner.match_score
        norm_score = best_score / max_possible if max_possible > 0 else 0

    return best_seq, norm_score, matched
