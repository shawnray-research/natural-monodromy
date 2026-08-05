"""
Why there is a wall at every T/12.

The earlier attempt asked whether a symmetry could carry one class of wall to the
other. None can: every forward shift in the group is an even multiple of T/12,
so Euler configurations go to Euler configurations and isosceles-only to
isosceles-only. That was the wrong question.

The group also contains six TIME-REVERSING elements, and a reversal t -> -t + s
has fixed instants at t = s/2 and t = s/2 + T/2. At a fixed instant the
configuration must be invariant under the spatial part g of that symmetry, and
for three bodies that pins the shape:

  * g a reflection. A triangle carried to itself by a reflection has two equal
    sides, so the configuration is ISOSCELES.
  * g the rotation by pi. Three is odd, so one body must be the fixed point and
    sit at the center, with the other two antipodal about it. That is collinear
    with the middle body at the midpoint of the other two, which is the EULER
    configuration, and it is isosceles as well.

So each reversing element contributes two instants at which a wall is forced, and
the twelve instants so produced are exactly the twelve walls. That is the
derivation the spacing was missing.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scipy.optimize import linear_sum_assignment
from mono.nbody import FIGURE8, integrate

F8 = FIGURE8
MS, T, NS = F8["m"], F8["T"], 12000


def match_dist(A, B):
    C = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)
    ri, ci = linear_sum_assignment(C)
    return C[ri, ci].max()


def main():
    print(__doc__)
    t, X, V = integrate(F8["x"], F8["v"], MS, T, NS)
    G = {"reflect x": np.diag([1.0, -1.0]),
         "reflect y": np.diag([-1.0, 1.0]),
         "rotate pi": -np.eye(2),
         "identity": np.eye(2)}

    print("=" * 76)
    print("1. The time-reversing symmetries, and their fixed instants")
    print("=" * 76)
    probe = np.linspace(0, NS - 1, 240).astype(int)
    rev = []
    for name, g in G.items():
        for m in range(12):
            s = m * NS // 12
            worst = max(match_dist(X[k] @ g.T, X[(-k + s) % NS]) for k in probe)
            if worst < 1e-3:
                rev.append((name, m, worst))
    print(f"  {'spatial part':>12}  {'shift':>7}  {'fixed instants t/T':>26}  mismatch")
    fixed = []
    for name, m, w in rev:
        f1, f2 = m / 24.0, m / 24.0 + 0.5
        fixed += [(f1 % 1.0, name), (f2 % 1.0, name)]
        print(f"  {name:>12}  {m:2d}/12  {f'{f1%1.0:.5f}, {f2%1.0:.5f}':>26}  {w:.1e}")

    print("\n" + "=" * 76)
    print("2. At each fixed instant, is the configuration invariant under g?")
    print("=" * 76)
    print(f"  {'t/T':>8} {'g':>12} {'|g(config) - config|':>22} "
          f"{'two sides equal?':>18} {'collinear?':>12}")
    fixed.sort()
    walls = []
    for f, name in fixed:
        k = int(round(f * NS)) % NS
        P = X[k]
        d = match_dist(P @ G[name].T, P)
        r = np.array([np.linalg.norm(P[1]-P[2]), np.linalg.norm(P[0]-P[2]),
                      np.linalg.norm(P[0]-P[1])])
        rs = np.sort(r)
        iso = min(rs[1]-rs[0], rs[2]-rs[1]) / rs.mean()
        col = abs(np.cross(np.append(P[1]-P[0], 0), np.append(P[2]-P[0], 0))[2])
        walls.append((f, name, iso, col))
        print(f"  {f:8.5f} {name:>12} {d:22.2e} {iso:18.2e} {col:12.2e}")

    print("\n" + "=" * 76)
    print("3. Do these fixed instants account for all twelve walls?")
    print("=" * 76)
    iso_t = []
    for k in range(NS):
        r0 = np.array([np.linalg.norm(X[k][1]-X[k][2]), np.linalg.norm(X[k][0]-X[k][2]),
                       np.linalg.norm(X[k][0]-X[k][1])])
        r1 = np.array([np.linalg.norm(X[(k+1) % NS][1]-X[(k+1) % NS][2]),
                       np.linalg.norm(X[(k+1) % NS][0]-X[(k+1) % NS][2]),
                       np.linalg.norm(X[(k+1) % NS][0]-X[(k+1) % NS][1])])
        for (a, b) in ((1, 2), (0, 2), (0, 1)):
            if np.sign(r0[a]-r0[b]) != np.sign(r1[a]-r1[b]):
                iso_t.append(k / NS)
    iso_t = np.sort(np.array(iso_t))
    ft = np.sort(np.array([f for f, _, _, _ in walls]))
    print(f"  walls found by scanning:      {len(iso_t)}")
    print(f"  fixed instants of reversals:  {len(ft)}")
    if len(ft) == len(iso_t):
        print(f"  they agree to {np.abs(ft - iso_t).max():.2e} of a period")
    print(f"\n  fixed instants at t/T x 12 = {np.array2string(ft*12, precision=3)}")
    eu = sorted(round(f*12) for f, n, _, c in walls if n == "rotate pi")
    rf = sorted(round(f*12) for f, n, _, c in walls if n.startswith("reflect"))
    print(f"  from the pi rotation (Euler, collinear):  k = {eu}")
    print(f"  from the reflections (proper isosceles):  k = {rf}")
    print("\n  The rotation by pi supplies the even twelfths and the reflections")
    print("  supply the odd ones, so the two families interleave because they come")
    print("  from different conjugacy classes of reversing symmetry. The spacing is")
    print("  forced after all; the earlier search simply looked at the forward")
    print("  shifts, which cannot see it.")


if __name__ == "__main__":
    main()
