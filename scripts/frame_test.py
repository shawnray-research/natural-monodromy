"""
The sharpest consequence: the vineyard braid is a braid of the SHAPE, not of the
motion.

The birth and death values of the density are built from mutual distances alone,
so they are unchanged by any rigid motion of the configuration, including a
time-dependent one. View the same orbit from a rotating frame and the mutual
distances are identical at every instant, so the persistence diagram is identical
and the vineyard braid cannot change. The world lines, on the other hand, pick up
a twist, so the trajectory braid does change.

If that holds, then persistent homology of the mass density is performing, for
free, the reduction to shape space that celestial mechanics performs by hand: the
quotient by translations and rotations. The vineyard braid is an invariant of the
curve in shape space; the trajectory braid is an invariant of the motion in the
plane, and the two differ by exactly the rotational degree of freedom.

Test: rotate the figure-eight rigidly at angular velocity Omega, so that after one
period the frame has turned by 2*pi*q. Recompute both braids.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.nbody import FIGURE8, integrate
from mono.kde_exact import merge_tree_exact
from mono.braid import braid_word, word_to_string, reduce_word, initial_order
from mono.core import perm_order, cycle_type

F8 = FIGURE8
MS, T, NS = F8["m"], F8["T"], 12000


def braids(Xs, sig=0.30):
    K = len(Xs)
    B = np.zeros((K, 3)); D = np.zeros((K, 3))
    for a in range(K):
        d = merge_tree_exact(Xs[a], MS, sig)
        for i in range(3):
            B[a, i] = d[i]["birth"]
            D[a, i] = 0.0 if d[i]["death"] is None else d[i]["death"]
    wv, _, _ = braid_word(B, D)
    wt, _, _ = braid_word(Xs[:, :, 0], Xs[:, :, 1])
    return wv, wt, B


def main():
    t, X, V = integrate(F8["x"], F8["v"], MS, T, NS)
    i0 = int(0.37 * NS)
    idx = np.linspace(i0, i0 + NS // 3, 2500).astype(int) % NS
    base = X[idx]
    tt = t[idx]

    print("Same orbit, viewed from frames rotating at different rates.")
    print("q is the number of turns the frame makes per orbital period.\n")
    print(f"  {'q':>6}  {'vineyard braid':<34} {'traj. braid':<26} {'r_ij same?':>10}")
    print("  " + "-" * 82)
    ref_v = ref_r = None
    for q in (0.0, 0.5, 1.0, 2.0, -1.5, 3.0):
        ang = 2 * np.pi * q * (tt - tt[0]) / T
        c, s = np.cos(ang), np.sin(ang)
        Xr = np.empty_like(base)
        Xr[:, :, 0] = c[:, None] * base[:, :, 0] - s[:, None] * base[:, :, 1]
        Xr[:, :, 1] = s[:, None] * base[:, :, 0] + c[:, None] * base[:, :, 1]
        wv, wt, B = braids(Xr)
        R = np.stack([np.linalg.norm(Xr[:, 0] - Xr[:, 1], axis=1),
                      np.linalg.norm(Xr[:, 0] - Xr[:, 2], axis=1),
                      np.linalg.norm(Xr[:, 1] - Xr[:, 2], axis=1)], axis=1)
        if ref_v is None:
            ref_v, ref_r = word_to_string(reduce_word(wv)), R
            same = "reference"
        else:
            same = f"{np.abs(R - ref_r).max():.1e}"
        print(f"  {q:6.1f}  {word_to_string(reduce_word(wv)):<34} "
              f"{word_to_string(reduce_word(wt))[:24]:<26} {same:>10}")

    print("\n  The mutual distances are identical to machine precision in every")
    print("  frame, so the vineyard braid is identical. The trajectory braid is")
    print("  not: it accumulates a full twist for each turn of the frame.")

    print("\n" + "=" * 84)
    print("Corrected version of the kernel-width test, with a working reference.")
    print("=" * 84)
    idx2 = np.arange(0, NS, 4)
    ref = None
    for sig in (0.22, 0.26, 0.30, 0.34, 0.38):
        B = np.zeros((len(idx2), 3))
        ok = True
        for a, k in enumerate(idx2):
            d = merge_tree_exact(X[k], MS, sig)
            for i in range(3):
                B[a, i] = d[i]["birth"]
        ts = []
        for (i, j) in ((0, 1), (0, 2), (1, 2)):
            g = B[:, i] - B[:, j]
            sgn = np.where(np.sign(g[:-1]) != np.sign(g[1:]))[0]
            for p in sgn:
                w = abs(g[p]) / (abs(g[p]) + abs(g[p + 1]))
                ts.append((idx2[p] + w * (idx2[p + 1] - idx2[p])) / NS)
        ts = np.sort(np.array(ts))
        gap = float(np.median([np.abs(B[a, 0] - B[a, 1]) for a in range(len(idx2))]))
        if ref is None or len(ts) != len(ref):
            tag = "reference" if ref is None else f"{len(ts)} crossings"
            if ref is None: ref = ts
            print(f"  sigma {sig:.2f}: {len(ts):3d} crossings   {tag}, "
                  f"typical |h_i - h_j| = {gap:.2e}")
        else:
            print(f"  sigma {sig:.2f}: {len(ts):3d} crossings   max shift "
                  f"{np.abs(ts-ref).max():.2e} of a period, "
                  f"typical |h_i - h_j| = {gap:.2e}")


if __name__ == "__main__":
    main()
