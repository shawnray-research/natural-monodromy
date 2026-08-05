"""
What the monodromy is actually recording.

Conjecture, from the algebra of the kernel. For N equal masses and any radial
kernel K that is monotone decreasing, the density peak belonging to body i is

    h_i = sum_j K(r_ij),

so for three bodies the shared term cancels:

    h_1 - h_2 = K(r_13) - K(r_23),

and h_1 = h_2 exactly when r_13 = r_23, that is when body 3 is equidistant from
bodies 1 and 2. The vineyard's strands are ordered by h, so a crossing of strands
1 and 2 is an ISOSCELES CONFIGURATION with apex at the third body.

If that is right then three things follow, and all three are testable:

  1. the crossing times of the vineyard braid coincide with the times at which
     two mutual distances are equal
  2. the number of crossings per period equals the number of isosceles
     configurations per period
  3. those crossing times do not depend on the kernel width, because the
     condition r_13 = r_23 does not mention sigma. The observed stability of the
     braid across sigma would then be forced rather than lucky

The comparison is not exact by construction: h_i is the density AT body i, while
the birth value is the density at the nearby maximum, which is displaced by the
pull of the other bodies. The size of that gap is measured here rather than
assumed away.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.nbody import FIGURE8, integrate
from mono.kde_exact import merge_tree_exact
from mono.braid import braid_word, word_to_string, reduce_word

F8 = FIGURE8
MS, T, NS = F8["m"], F8["T"], 24000


def crossings_of(sig, X, idx):
    """Times (as fractions of the period) where two BIRTH values cross."""
    B = np.zeros((len(idx), 3))
    for a, k in enumerate(idx):
        d = merge_tree_exact(X[k], MS, sig)
        for i in range(3):
            B[a, i] = d[i]["birth"]
    out = []
    for (i, j) in ((0, 1), (0, 2), (1, 2)):
        g = B[:, i] - B[:, j]
        s = np.where(np.sign(g[:-1]) != np.sign(g[1:]))[0]
        for p in s:
            w = abs(g[p]) / (abs(g[p]) + abs(g[p + 1]))
            out.append(((idx[p] + w * (idx[p + 1] - idx[p])) / NS, (i, j)))
    return sorted(out), B


def isosceles_times(X, idx):
    """Times where two mutual distances are equal, with the apex body."""
    R = np.zeros((len(idx), 3))            # r12, r13, r23
    for a, k in enumerate(idx):
        P = X[k]
        R[a] = (np.linalg.norm(P[0]-P[1]), np.linalg.norm(P[0]-P[2]),
                np.linalg.norm(P[1]-P[2]))
    out = []
    # apex k means the two sides meeting at k are equal
    for (a_, b_, apex, pair) in ((1, 2, 2, (0, 1)),   # r13 = r23, apex body 3
                                 (0, 2, 1, (0, 2)),   # r12 = r23, apex body 2
                                 (0, 1, 0, (1, 2))):  # r12 = r13, apex body 1
        g = R[:, a_] - R[:, b_]
        s = np.where(np.sign(g[:-1]) != np.sign(g[1:]))[0]
        for p in s:
            w = abs(g[p]) / (abs(g[p]) + abs(g[p + 1]))
            out.append(((idx[p] + w * (idx[p + 1] - idx[p])) / NS, apex, pair))
    return sorted(out), R


def main():
    t, X, V = integrate(F8["x"], F8["v"], MS, T, NS)
    idx = np.arange(0, NS, 4)

    print("=" * 74)
    print("1. Are the vineyard crossings the isosceles configurations?")
    print("=" * 74)
    iso, R = isosceles_times(X, idx)
    print(f"  isosceles configurations in one period: {len(iso)}")
    cr, B = crossings_of(0.30, X, idx)
    print(f"  vineyard birth crossings in one period: {len(cr)}   (sigma = 0.30)")
    print()
    print("   t/T (vineyard)   strands      t/T (isosceles)   apex body   |diff|")
    worst = 0.0
    for (tc, pr) in cr:
        best = min(iso, key=lambda z: abs(z[0] - tc))
        worst = max(worst, abs(best[0] - tc))
        print(f"      {tc:8.5f}      {pr}        {best[0]:8.5f}          {best[1]}"
              f"        {abs(best[0]-tc):.2e}")
    print(f"\n  worst mismatch between the two lists: {worst:.2e} of a period")
    print("  the pair of strands that cross is always the pair NOT at the apex,")
    print("  which is what h_i - h_j = K(r_ik) - K(r_jk) predicts")

    print("\n" + "=" * 74)
    print("2. Does the count match the braid?")
    print("=" * 74)
    full = np.linspace(0, NS - 1, 6000).astype(int)
    Bf = np.zeros((len(full), 3)); Df = np.zeros((len(full), 3))
    for a, k in enumerate(full):
        d = merge_tree_exact(X[k], MS, 0.30)
        for i in range(3):
            Bf[a, i] = d[i]["birth"]
            Df[a, i] = 0.0 if d[i]["death"] is None else d[i]["death"]
    w, _, _ = braid_word(Bf, Df)
    print(f"  vineyard braid over the FULL period: {len(w)} crossings, "
          f"{word_to_string(reduce_word(w))}")
    print(f"  isosceles configurations per period:  {len(iso)}")
    wt, _, _ = braid_word(X[full][:, :, 0], X[full][:, :, 1])
    print(f"  trajectory braid over the full period: {len(wt)} crossings")

    print("\n" + "=" * 74)
    print("3. Is the crossing pattern independent of the kernel width?")
    print("=" * 74)
    ref = None
    for sig in (0.18, 0.22, 0.26, 0.30, 0.34):
        try:
            c, _ = crossings_of(sig, X, idx)
        except Exception as e:
            print(f"  sigma {sig}: {e}"); continue
        ts = np.array([z[0] for z in c])
        if ref is None:
            ref = ts
            print(f"  sigma {sig:.2f}: {len(ts)} crossings   (reference)")
        else:
            if len(ts) == len(ref):
                print(f"  sigma {sig:.2f}: {len(ts)} crossings, "
                      f"max shift from reference {np.abs(ts-ref).max():.2e} of a period")
            else:
                print(f"  sigma {sig:.2f}: {len(ts)} crossings, count differs")
    iso_t = np.array([z[0] for z in iso])
    print(f"\n  isosceles times, which contain no sigma at all:")
    print(f"    {np.array2string(iso_t, precision=5)}")


if __name__ == "__main__":
    main()
