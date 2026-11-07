import numpy as np
import math
from collections import Counter
from functools import lru_cache
from itertools import combinations
B=999

def score(s):
    if not len(s):
        return (0, 0)
    changed = 0
    longest = 1
    streak = 1
    for a, b in zip(s, s[1:]):
        streak = streak + 1 if a != b else 1
        changed += int(a != b)
        longest = max(longest, streak)
    return (changed, longest)

@lru_cache(maxsize=2048)
def exact_distribution(n, k):
    freq = Counter()
    for pos in combinations(range(n), k):
        s = bytearray(n)
        for x in pos:
            s[x] = 1
        freq[score(s)] += 1
    values = np.asarray(list(freq))
    p = np.asarray(list(freq.values()), dtype=float)
    return (values, p / p.sum())

def null(chains, seed):
    rng = np.random.default_rng(seed)
    total = np.zeros(B, dtype=np.int64)
    longest = np.zeros(B, dtype=np.int64)
    groups = Counter(((len(s), s.count('A')) for s in chains))
    for (n, k), multiplicity in groups.items():
        if k == 0 or k == n or n <= 2:
            tr, lg = (0, 1) if k == 0 or k == n else (1, 2)
            total += tr * multiplicity
            np.maximum(longest, lg, out=longest)
        elif math.comb(n, min(k, n - k)) <= 12000 and n < 150:
            vals, prob = exact_distribution(n, k)
            hist = rng.multinomial(multiplicity, prob, size=B)
            total += hist @ vals[:, 0]
            np.maximum(longest, np.where(hist > 0, vals[:, 1], 0).max(axis=1), out=longest)
        else:
            for _ in range(multiplicity):
                s = np.arange(n) < k
                for i in range(B):
                    vals = rng.permutation(s)
                    different = vals[:-1] != vals[1:]
                    total[i] += int(different.sum())
                    cut = np.r_[-1, np.flatnonzero(~different), n - 1]
                    longest[i] = max(longest[i], int(np.diff(cut).max()))
    return (total, longest)
