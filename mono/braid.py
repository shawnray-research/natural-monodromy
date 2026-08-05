"""
The braid carried by a closed vineyard.

Braiding Vineyards is about braids, not merely permutations: it proves that any
braid can be realized by a closed vineyard, using a constructed manifold. So the
right question to ask of a natural example is not only "what is the monodromy
order" but "which braid is it".

A closed vineyard with k vines that stay off the diagonal is k strands in
(birth, death) x S^1. Projecting to the birth coordinate turns it into a braid
diagram: strands cross when two births coincide, and the crossing sign is set by
which of the two has the larger death, that is, which strand passes in front.
Reading the crossings in order of t gives a word in the Artin generators.
"""

from __future__ import annotations

import numpy as np


def braid_word(births, deaths):
    """
    births, deaths : (T, k) arrays, the vines sampled around the loop, in a fixed
                     vine labeling (row t, column = vine index).

    Returns the braid word as a list of signed generator indices: +j means
    strand in position j crosses over strand j+1, -j means under. Positions are
    1-based, following the Artin convention.
    """
    T, k = births.shape
    order = list(np.argsort(births[0]))        # positions, left to right
    word = []
    crossings = []
    for t in range(1, T):
        changed = True
        guard = 0
        while changed and guard < 10 * k:
            changed = False
            guard += 1
            for j in range(k - 1):
                a, b = order[j], order[j + 1]
                if births[t, a] > births[t, b]:
                    # they crossed between t-1 and t; sign from the death coordinate
                    sgn = 1 if deaths[t, a] > deaths[t, b] else -1
                    word.append(sgn * (j + 1))
                    crossings.append((t, a, b, sgn))
                    order[j], order[j + 1] = order[j + 1], order[j]
                    changed = True
    return word, order, crossings


def word_to_string(word):
    out = []
    for w in word:
        j = abs(w)
        out.append(f"s{j}" if w > 0 else f"s{j}^-1")
    return " ".join(out) if out else "(empty)"


def permutation_of_word(word, k):
    """Underlying permutation of a braid word, as a check against the vineyard's."""
    p = list(range(k))
    for w in word:
        j = abs(w) - 1
        p[j], p[j + 1] = p[j + 1], p[j]
    return p


def reduce_word(word):
    """Free reduction: cancel adjacent s_j s_j^{-1}. Not a full braid
    normalization, just enough to report a tidier word."""
    out = []
    for w in word:
        if out and out[-1] == -w:
            out.pop()
        else:
            out.append(w)
    return out


def initial_order(births):
    """
    The strand order that braid_word starts from: positions left to right at the
    first slice.

    Any diagram of the braid MUST start from this order. Drawing from
    list(range(k)) instead silently produces non-adjacent crossings, in which two
    strands appear to swap straight through a third, which is not an Artin
    generator and not a braid diagram. That bug produced figure 3 panels (d), (e)
    and (f) in their first form, while the braid words themselves were correct.
    """
    return list(np.argsort(np.asarray(births)[0]))
