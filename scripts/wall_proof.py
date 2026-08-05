"""
Proving the two claims the note reported as measured.

Both are the same statement. Write L for the perpendicular bisector of q_i q_j,
sigma for the reflection in L, and H_i, H_j for the two closed half planes. Let
M_i be the maximum belonging to body i and h_i = rho(M_i) its birth value.

    PROPOSITION. Let K be strictly decreasing. Suppose sigma M_j lies in the
    basin of body i, so that h_i >= rho(sigma M_j), and that M_j is interior to
    H_j. If body k lies strictly inside H_i, that is r_ik < r_jk, then h_i > h_j.

    PROOF. sigma fixes the pair {q_i, q_j}, so the only term that moves is the
    one belonging to body k:

        rho o sigma (x) - rho(x) = K(|x - sigma q_k|) - K(|x - q_k|).

    For x interior to H_j we have |x - sigma q_k| < |x - q_k|, since sigma q_k
    lies in H_j and L is the perpendicular bisector of q_k and sigma q_k. As K is
    strictly decreasing, rho o sigma > rho on the interior of H_j. Hence

        h_i  >=  rho(sigma M_j)  =  rho o sigma (M_j)  >  rho(M_j)  =  h_j.   []

Two corollaries, which are exactly the two claims:

    the strand order is the order of the opposite sides, since h_i > h_j iff
    r_ik < r_jk, and r_jk is the side opposite body i;

    h_i = h_j forces r_ik = r_jk, so the walls are the isosceles configurations
    and nothing else.

No derivatives, no bandwidth, no Gaussian. The reflection argument the note
already gave is the equality case of this one.

A FIRST ATTEMPT took the hypothesis to be that h_i is the largest value of rho on
the whole half plane H_i, which is tidier and false. When r_ik < r_jk body k is
in H_i as well, and its peak is frequently the taller one: the margin runs to
-3.8e-03 at sigma = 0.30, and fails at every bandwidth tried. Only the single
value rho(sigma M_j) is ever needed, and sigma M_j sits near q_i by construction.

What this script checks is the surviving hypothesis, over the orbit, for all six
ordered pairs and a range of bandwidths. Its margin vanishes exactly at the
walls, which is the equality case: there sigma carries M_j onto M_i.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.nbody import FIGURE8, integrate
from mono.kde_exact import merge_tree_exact, newton_critical

F8 = FIGURE8
MS, T, NS = F8["m"], F8["T"], 3000


def rho(x, P, sigma):
    x = np.atleast_2d(x)
    d2 = ((x[:, None, :] - P[None, :, :]) ** 2).sum(2)
    return np.exp(-d2 / (2 * sigma ** 2)).sum(1)


def check_instant(P, sigma):
    """
    The hypothesis the proof actually needs, for every ordered pair (i, j):

        h_i >= rho(sigma M_j),   sigma = reflection in the bisector of q_i q_j,

    together with M_j lying strictly on its own side of that bisector.
    """
    d = merge_tree_exact(P, MS, sigma)
    if d is None:
        return None
    M = {i: newton_critical(P[i], P, MS, sigma)["p"] for i in range(3)}
    h = {i: d[i]["birth"] for i in range(3)}

    worst_hyp, worst_int = np.inf, np.inf
    for i in range(3):
        for j in range(3):
            if i == j:
                continue
            mid = 0.5 * (P[i] + P[j])
            nrm = (P[i] - P[j]) / np.linalg.norm(P[i] - P[j])
            # reflect M_j in the bisector
            sMj = M[j] - 2.0 * ((M[j] - mid) @ nrm) * nrm
            worst_hyp = min(worst_hyp, h[i] - float(rho(sMj, P, sigma)[0]))
            worst_int = min(worst_int, -((M[j] - mid) @ nrm))
    return worst_hyp, worst_int


def main():
    print(__doc__)
    t, X, V = integrate(F8["x"], F8["v"], MS, T, NS)

    print("=" * 78)
    print("Does the hypothesis hold along the figure-eight?")
    print("=" * 78)
    print("  margin   = h_i minus rho(sigma M_j), the one value the proof uses.")
    print("             Must be >= 0, i.e. sigma M_j is in body i's basin.")
    print("  interior = signed distance of M_j from the bisector, on its own side.")
    print("             Must be > 0.\n")
    print(f"  {'sigma':>7} {'instants':>9} {'min margin':>14} {'min interior':>14} "
          f"{'hypothesis':>12}")
    for sigma in (0.22, 0.24, 0.26, 0.28, 0.30, 0.32):
        m, q, n = np.inf, np.inf, 0
        for k in range(0, NS, 12):
            r = check_instant(X[k], sigma)
            if r is None:
                continue
            n += 1
            m = min(m, r[0]); q = min(q, r[1])
        ok = m >= -1e-9 and q > 0
        print(f"  {sigma:7.2f} {n:9d} {m:14.3e} {q:14.3e} "
              f"{'HOLDS' if ok else 'FAILS':>12}")

    print("\n" + "=" * 78)
    print("The two corollaries, checked directly as a sanity test")
    print("=" * 78)
    sigma = 0.30
    bad_order = bad_wall = tot = 0
    for k in range(0, NS, 4):
        P = X[k]
        d = merge_tree_exact(P, MS, sigma)
        if d is None:
            continue
        tot += 1
        h = np.array([d[i]["birth"] for i in range(3)])
        opp = np.array([np.linalg.norm(P[1]-P[2]), np.linalg.norm(P[0]-P[2]),
                        np.linalg.norm(P[0]-P[1])])
        if tuple(np.argsort(h)) != tuple(np.argsort(opp)):
            bad_order += 1
        for a, b, c in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
            gap = abs(h[a] - h[b])
            iso = abs(np.linalg.norm(P[a]-P[c]) - np.linalg.norm(P[b]-P[c]))
            if gap < 1e-9 and iso > 1e-6:
                bad_wall += 1
    print(f"  strand order equals opposite-side order:  "
          f"{tot - bad_order} of {tot}")
    print(f"  births coincide off the isosceles locus:  {bad_wall} times")
    print("\n  These now follow from the proposition rather than standing on their")
    print("  own, so the numbers are a check on the hypothesis, not the evidence.")


if __name__ == "__main__":
    main()
