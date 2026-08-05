"""
Cross-check the braid, at both stages.

GUDHI verified the persistence but not what happens after it. Two things are
still mine: the extraction of a braid word from the tracked vines
(`mono.braid.braid_word`), and the identification of its closure.

  A  independent extraction. `braid_word` walks the strand order and records
     swaps as it bubbles them into place. The check here builds an EVENT LIST
     instead: every time a pair of birth values crosses is found first, the
     events are sorted by time, and the generator index is read off from the
     positions held at that moment. Different algorithm, same input.
  B  independent identification. The closure was identified by me through the
     reduced Burau matrix and by a conjugacy argument. SnapPy and spherogram
     are an outside implementation: hand them the word, ask what the closure is.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.nbody import FIGURE8, integrate
from mono.kde_exact import merge_tree_exact
from mono.braid import braid_word, word_to_string, reduce_word, initial_order

F8 = FIGURE8
MS, T, NS = F8["m"], F8["T"], 12000


def braid_word_events(B, D):
    """
    Independent extraction. Find every crossing of a pair of birth values,
    order the events in time, and read the generator index off the positions
    held at that instant. No incremental sorting.
    """
    K, k = B.shape
    events = []
    for i in range(k):
        for j in range(i + 1, k):
            g = B[:, i] - B[:, j]
            for tstep in range(K - 1):
                if np.sign(g[tstep]) != np.sign(g[tstep + 1]):
                    w = abs(g[tstep]) / (abs(g[tstep]) + abs(g[tstep + 1]))
                    events.append((tstep + w, i, j))
    events.sort()
    order = list(np.argsort(B[0]))
    word = []
    for (tt, i, j) in events:
        pi, pj = order.index(i), order.index(j)
        if abs(pi - pj) != 1:
            raise AssertionError(f"non-adjacent crossing at positions {pi},{pj}")
        p = min(pi, pj)
        a, b = order[p], order[p + 1]
        t0 = int(np.floor(tt))
        sgn = 1 if D[t0, a] > D[t0, b] else -1
        word.append(sgn * (p + 1))
        order[p], order[p + 1] = order[p + 1], order[p]
    return word


def main():
    print(__doc__)
    t, X, V = integrate(F8["x"], F8["v"], MS, T, NS)
    i0 = int(0.37 * NS)

    print("=" * 76)
    print("A. Two independent extractions of the braid word")
    print("=" * 76)
    for span, lab in ((NS // 3, "T/3"), (NS, "T")):
        for K in (2000, 4000):
            idx = np.linspace(i0, i0 + span, K).astype(int) % NS
            B = np.zeros((K, 3)); D = np.zeros((K, 3))
            for a, kk in enumerate(idx):
                d = merge_tree_exact(X[kk], MS, 0.30)
                for i in range(3):
                    B[a, i] = d[i]["birth"]
                    D[a, i] = 0.0 if d[i]["death"] is None else d[i]["death"]
            w1, _, _ = braid_word(B, D)
            w2 = braid_word_events(B, D)
            same = "IDENTICAL" if list(w1) == list(w2) else "DIFFER"
            print(f"  loop {lab:3s}  K={K:5d}  mono {word_to_string(reduce_word(w1)):34s}")
            print(f"  {'':13s}          events {word_to_string(reduce_word(w2)):34s}  {same}")

    print("\n" + "=" * 76)
    print("B. Independent identification of the closure")
    print("=" * 76)
    word = [-2, 1, -2, 1]
    print(f"  measured word: {word_to_string(word)}")
    try:
        import spherogram
        L = spherogram.Link(braid_closure=word)
        print(f"  spherogram: {L.exterior().identify() if hasattr(L,'exterior') else ''}")
        print(f"    components {len(L.link_components)}, crossings {len(L.crossings)}")
        print(f"    Alexander polynomial {L.alexander_polynomial()}")
        try:
            print(f"    Jones polynomial     {L.jones_polynomial()}")
        except Exception as e:
            print(f"    Jones unavailable: {e}")
    except Exception as e:
        print(f"  spherogram failed: {e}")
    try:
        import snappy
        M = snappy.Manifold("braid" + str(word))
        print(f"  snappy identification: {M.identify()}")
        print(f"    volume {M.volume():.6f}")
    except Exception as e:
        print(f"  snappy: {e}")

    print("\n  Reference: the figure-eight knot 4_1 has Alexander polynomial")
    print("  t - 3 + 1/t and hyperbolic volume 2.029883...")


if __name__ == "__main__":
    main()
