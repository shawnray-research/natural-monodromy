"""
Are there knees in the vines, and is the monodromy carried by the critical points?

A knee is a corner in a vine, where a class changes which critical points it is
born at or dies at. Knees are the signature of PAIRING monodromy, where the
critical points stay put and the elder rule reassigns. Their absence points to
TRANSPORT monodromy, where the critical points move and trade places.

The answer here is that there are both, and they do different jobs.

Over T/3 the shape curve crosses four walls. At every wall two births coincide,
because a wall is an isosceles configuration and the reflection makes two peaks
equal. What differs is which two.

  The TOP two cross. One of them is the essential class, pinned at the infimum,
  and the other dies at a saddle near 1.03. Their births meet but the two diagram
  points stay 1.03 apart in death, so nothing is exchanged. The essential label
  passes from one body to the other, and the vine at death zero, whose birth is
  the running maximum of two crossing curves, turns a CORNER. This is a knee and
  it carries no transposition.

  The BOTTOM two cross. Neither is essential, both die at saddles, and at the
  wall their deaths agree as well as their births. The two diagram points collide
  to within 1e-07 and are EXCHANGED. This is a crossing and it carries a
  transposition.

Over T/3 there are two of each. The two transpositions compose to the 3-cycle, so
the monodromy comes entirely from the crossings, and the crossings are the maxima
themselves swapping places in the domain. Following the maxima by continuation
gives the same 3-cycle with a closing error of 4e-08, which is the direct answer:
yes, the monodromy is the monodromy of the individual critical points.

A wrong turn worth recording. The first version of this script identified which
saddle a class dies at by matching the merge tree's saddle position to the
nearest pair midpoint. In an elongated triangle two midpoints can be nearer the
same saddle, and the tie-break then flips every slice: it reported 175 "pairing
changes" that were an artifact of the matching. The check that exposed it was
that the two finite death values are never closer than 7.6e-04, so no genuine
reassignment could be happening without a visible jump in a vine.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.nbody import FIGURE8, integrate
from mono.kde_exact import merge_tree_exact, newton_critical

F8 = FIGURE8
MS, T, NS = F8["m"], F8["T"], 48000
SIG = 0.30


def main():
    print(__doc__)
    t, X, V = integrate(F8["x"], F8["v"], MS, T, NS)
    i0 = int(0.37 * NS)
    K = 24000
    idx = np.linspace(i0, i0 + NS // 3, K).astype(int) % NS

    B = np.zeros((K, 3)); D = np.zeros((K, 3)); E = np.zeros((K, 3), bool)
    for a, k in enumerate(idx):
        d = merge_tree_exact(X[k], MS, SIG)
        for i in range(3):
            B[a, i] = d[i]["birth"]
            D[a, i] = 0.0 if d[i]["death"] is None else d[i]["death"]
            E[a, i] = d[i]["death"] is None

    events = []
    for i, j in ((0, 1), (0, 2), (1, 2)):
        g = B[:, i] - B[:, j]
        for a in np.where(np.sign(g[:-1]) != np.sign(g[1:]))[0]:
            events.append((a, i, j))
    events.sort()

    print("=" * 78)
    print("The four walls over T/3")
    print("=" * 78)
    print(f"  {'s':>7} {'pair':>6} {'|db|':>10} {'|dd|':>10} {'rank':>15} "
          f"{'event':>28}")
    ncross = ncorner = 0
    for a, i, j in events:
        k = ({0, 1, 2} - {i, j}).pop()
        rank = "top two" if B[a, i] > B[a, k] else "bottom two"
        dd = abs(D[a, i] - D[a, j])
        if dd < 1e-2:
            kind = "collision, transposition"; ncross += 1
        else:
            kind = "corner, essential transfer"; ncorner += 1
        print(f"  {a/K:7.4f} {f'{i},{j}':>6} {abs(B[a,i]-B[a,j]):10.2e} "
              f"{dd:10.2e} {rank:>15} {kind:>28}")
    print(f"\n  over T/3: {ncross} crossings, {ncorner} corners")
    print(f"  over T:   {3*ncross} crossings, {3*ncorner} corners")

    print("\n" + "=" * 78)
    print("The corners, in detail")
    print("=" * 78)
    for a, i, j in events:
        if abs(D[a, i] - D[a, j]) < 1e-2:
            continue
        lo, hi = max(0, a - 400), min(K - 1, a + 400)
        bess = np.array([B[m, int(np.argmax(E[m]))] for m in range(lo, hi)])
        d1 = np.diff(bess)
        turn = abs(d1[len(d1)//2 + 5] - d1[len(d1)//2 - 5]) / (abs(d1).mean() + 1e-30)
        print(f"  s = {a/K:.4f}: essential passes from body "
              f"{int(np.argmax(E[lo]))} to {int(np.argmax(E[hi-1]))}, "
              f"deaths {D[a,i]:.5f} and {D[a,j]:.5f}")
        print(f"    slope change in the essential vine's birth, "
              f"relative to its own scale: {turn:.2f}")

    print("\n" + "=" * 78)
    print("The crossings, in detail")
    print("=" * 78)
    for a, i, j in events:
        if abs(D[a, i] - D[a, j]) >= 1e-2:
            continue
        lo, hi = max(0, a - 600), min(K - 1, a + 600)
        pi_b = np.array([B[lo, i], D[lo, i]]); pj_b = np.array([B[lo, j], D[lo, j]])
        pi_a = np.array([B[hi, i], D[hi, i]]); pj_a = np.array([B[hi, j], D[hi, j]])
        keep = np.linalg.norm(pi_b - pi_a) + np.linalg.norm(pj_b - pj_a)
        swap = np.linalg.norm(pi_b - pj_a) + np.linalg.norm(pj_b - pi_a)
        sep = float(np.linalg.norm([B[a, i] - B[a, j], D[a, i] - D[a, j]]))
        print(f"  s = {a/K:.4f}: bodies {i},{j} separate by {sep:.2e} at the wall")
        print(f"    cost if each keeps its point {keep:.5f}, if exchanged {swap:.5f}")

    print("\n" + "=" * 78)
    print("The maxima in the domain")
    print("=" * 78)
    fine = np.linspace(i0, i0 + NS // 3, 8000).astype(int) % NS
    M = np.array([newton_critical(X[fine[0]][i], X[fine[0]], MS, SIG)["p"]
                  for i in range(3)])
    start = M.copy(); travel = np.zeros(3)
    for k in fine[1:]:
        for i in range(3):
            c = newton_critical(M[i], X[k], MS, SIG)
            travel[i] += float(np.linalg.norm(c["p"] - M[i]))
            M[i] = c["p"]
    perm = [int(np.argmin(np.linalg.norm(start - M[i], axis=1))) for i in range(3)]
    err = max(np.linalg.norm(M[i] - start[perm[i]]) for i in range(3))
    print(f"  continuation permutation: {perm}")
    print(f"  arc length travelled: {np.array2string(travel, precision=4)}")
    print(f"  closing error: {err:.2e}")
    print("\n  The maxima travel about two units each and land on each other's")
    print("  starting points. The two transpositions in the diagram are those")
    print("  exchanges, so the monodromy is carried by the critical points.")


if __name__ == "__main__":
    main()
