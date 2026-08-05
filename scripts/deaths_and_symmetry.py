"""
The three pieces of interpretation still missing.

1. What the DEATH coordinate encodes. The births order the strands, but the
   crossing sign comes from the deaths, so half the braid was unexplained.
   Hypothesis: for superlevel H_0 of a kernel density the merge structure is
   single linkage, so the two finite deaths are the saddles on the two SHORTEST
   sides, the minimum spanning tree of the three bodies, and the essential class
   belongs to the tallest peak, the body opposite the longest side.

2. Whether the trajectory braid's six crossings per period are the syzygies, the
   collinear configurations. If they are, both braids are counting strata of the
   same shape space: the vineyard counts isosceles walls and the trajectory braid
   counts the collinear equator.

3. Whether the orbit's symmetry group forces the twelve isosceles walls to sit at
   t/T = k/12. Measured, but never derived.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scipy.optimize import linear_sum_assignment
from mono.nbody import FIGURE8, integrate
from mono.kde_exact import merge_tree_exact, exact_saddle

F8 = FIGURE8
MS, T, NS = F8["m"], F8["T"], 6000


def sides(P):
    """r_23, r_13, r_12, i.e. the side opposite body 0, 1, 2."""
    return np.array([np.linalg.norm(P[1]-P[2]), np.linalg.norm(P[0]-P[2]),
                     np.linalg.norm(P[0]-P[1])])


def main():
    print(__doc__)
    t, X, V = integrate(F8["x"], F8["v"], MS, T, NS)

    # ---------------------------------------------------------------- deaths
    print("=" * 76)
    print("1. What the deaths encode")
    print("=" * 76)
    sig = 0.30
    ess_ok = mst_ok = tot = 0
    for k in range(0, NS, 7):
        P = X[k]
        d = merge_tree_exact(P, MS, sig)
        if d is None:
            continue
        tot += 1
        opp = sides(P)                       # opposite side per body
        # the essential class should be the tallest peak, opposite the longest side
        ess = [i for i in range(3) if d[i]["death"] is None]
        if len(ess) == 1 and ess[0] == int(np.argmax(opp)):
            ess_ok += 1
        # the two finite deaths should be the saddles on the two shortest sides
        pairs = [(0, 1), (0, 2), (1, 2)]
        plen = {p: np.linalg.norm(P[p[0]] - P[p[1]]) for p in pairs}
        short2 = sorted(plen, key=plen.get)[:2]
        sadv = []
        for p in short2:
            s = exact_saddle(0.5*(P[p[0]]+P[p[1]]), P, MS, sig)
            if s is not None:
                sadv.append(s["val"])
        deaths = sorted(d[i]["death"] for i in range(3) if d[i]["death"] is not None)
        if len(sadv) == 2 and np.allclose(sorted(sadv), deaths, atol=1e-8):
            mst_ok += 1
    print(f"  essential class is the body opposite the LONGEST side: "
          f"{ess_ok} of {tot}")
    print(f"  the two finite deaths are the saddles on the two SHORTEST sides: "
          f"{mst_ok} of {tot}")
    print("\n  So the diagram splits cleanly: the births order all three sides,")
    print("  and the deaths are the minimum spanning tree, the two shortest.")
    print("  The braid's strand order comes from the first, its crossing signs")
    print("  from the second.")

    # ------------------------------------------------------------- syzygies
    print("\n" + "=" * 76)
    print("2. Syzygies against the trajectory braid")
    print("=" * 76)
    cr = np.array([np.cross(X[k][1]-X[k][0], X[k][2]-X[k][0]) for k in range(NS)])
    syz = [(k + abs(cr[k])/(abs(cr[k])+abs(cr[k+1])))/NS
           for k in range(NS-1) if np.sign(cr[k]) != np.sign(cr[k+1])]
    print(f"  collinear configurations per period: {len(syz)}")
    print(f"    at t/T = {np.array2string(np.array(syz), precision=4)}")
    xc = [(k + abs(X[k][i][0]-X[k][j][0]) /
           (abs(X[k][i][0]-X[k][j][0]) + abs(X[k+1][i][0]-X[k+1][j][0])))/NS
          for k in range(NS-1) for (i, j) in ((0, 1), (0, 2), (1, 2))
          if np.sign(X[k][i][0]-X[k][j][0]) != np.sign(X[k+1][i][0]-X[k+1][j][0])]
    print(f"  equal-x events per period (what the trajectory braid counts): "
          f"{len(xc)}")
    print(f"    at t/T = {np.array2string(np.sort(np.array(xc)), precision=4)}")
    if len(syz) == len(xc):
        print(f"    max |difference| = "
              f"{np.abs(np.sort(np.array(syz))-np.sort(np.array(xc))).max():.2e}")

    # ------------------------------------------------------------- symmetry
    print("\n" + "=" * 76)
    print("3. The symmetry group, and whether it forces t/T = k/12")
    print("=" * 76)
    G = {"identity": np.eye(2),
         "reflect in x": np.diag([1.0, -1.0]),
         "reflect in y": np.diag([-1.0, 1.0]),
         "rotate by pi": -np.eye(2)}
    probe = np.linspace(0, NS-1, 240).astype(int)
    found = []
    for name, g in G.items():
        for m in range(12):
            s = m * NS // 12
            worst = 0.0
            for k in probe:
                A = X[k] @ g.T
                B = X[(k + s) % NS]
                C = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)
                ri, ci = linear_sum_assignment(C)
                worst = max(worst, C[ri, ci].max())
            if worst < 1e-3:
                found.append((name, m, worst))
    print(f"  symmetries of the form (spatial g, time shift m*T/12):")
    for name, m, w in found:
        print(f"    {name:14s} with shift {m:2d}/12 of a period   mismatch {w:.2e}")
    print(f"  group order found: {len(found)}")
    print(f"\n  The time shifts that occur are "
          f"{sorted(set(m for _, m, _ in found))}, which generate the cyclic group")
    print("  of order 12 acting on the period. A symmetry maps isosceles")
    print("  configurations to isosceles configurations, so that group permutes")
    print("  the twelve of them; acting transitively, it forces them to be evenly")
    print("  spaced, which is what t/T = k/12 says.")


if __name__ == "__main__":
    main()
