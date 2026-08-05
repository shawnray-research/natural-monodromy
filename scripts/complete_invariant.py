"""
How much of the configuration does the persistence diagram actually keep?

I had been calling the birth values a measure of "crowding", which treats them as
a lossy summary. For three bodies they are not lossy at all.

Write u = K(r_23), v = K(r_13), w = K(r_12). The peak heights are

    h_1 = w + v,    h_2 = w + u,    h_3 = v + u,

a linear system in (u, v, w). Adding all three gives u + v + w = S with
S = (h_1 + h_2 + h_3)/2, hence

    u = S - h_1,    v = S - h_2,    w = S - h_3,

and since K is monotone it inverts, so

    r_23 = K^{-1}(S - g_1),   r_13 = K^{-1}(S - g_2),   r_12 = K^{-1}(S - g_3).

The three heights therefore determine the three mutual distances exactly, and a
triangle is determined by its three sides up to congruence.

CAUTION, and the reason this file no longer supports the strong claim. All of the
above is about h_i as written, which is the density read AT the bodies. The
actual birth value is the density at the MAXIMUM, and the maximum is displaced
towards the other two bodies. The second table below measures the gap: feeding
true birth values through this inversion leaves a 33 to 35 per cent relative
error at every bandwidth in the working window, and scripts/completeness_test.py
shows it does not converge anywhere the merge tree is still computable. So the
exact statement is about a model of the diagram, not about the diagram.

What survives is in part 2, and it is what the braid actually needs: the ORDER of
the birth values is the order of the opposite sides, at 858 of 858 instants. The
braid reads only that order, never the values.

Two consequences worth checking as well.

  * The strand order is then the ordering of the three side lengths, since
    h_1 > h_2 iff r_13 < r_23. The braid is literally a braid of the sides of the
    triangle, ordered by length.
  * The count is a coincidence of three. The heights give N numbers and the
    shape needs N(N-1)/2, equal only at N = 3. For four bodies and more the
    heights cannot determine the shape, so the completeness is special to the
    three-body problem.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.nbody import FIGURE8, integrate
from mono.kde_exact import merge_tree_exact

F8 = FIGURE8
MS, T, NS = F8["m"], F8["T"], 6000


def recover(h, sigma, k0=1.0):
    """
    Mutual distances from the three peak heights, by the linear solve.

    h_i includes the SELF term K(0), which must come off before the system is
    linear in the pairwise kernel values. Forgetting it was the first bug here.
    """
    g = np.asarray(h, dtype=float) - k0
    S = 0.5 * (g[0] + g[1] + g[2])
    u, v, w = S - g[0], S - g[1], S - g[2]          # K(r_23), K(r_13), K(r_12)
    out = []
    for k in (u, v, w):
        if k <= 0:
            return None
        out.append(sigma * np.sqrt(-2.0 * np.log(k)))   # invert the Gaussian
    return np.array(out)                                 # r_23, r_13, r_12


def main():
    print(__doc__)
    t, X, V = integrate(F8["x"], F8["v"], MS, T, NS)

    print("=" * 74)
    print("1. Recovering the triangle from the diagram, along the orbit.")
    print("=" * 74)
    print("  Using the density AT the bodies, where the algebra is exact:\n")
    print(f"  {'t/T':>7} {'true (r23, r13, r12)':>34} {'recovered':>34} {'error':>10}")
    worst = 0.0
    for frac in (0.05, 0.20, 0.37, 0.55, 0.80):
        k = int(frac * NS); P = X[k]
        R = np.array([np.linalg.norm(P[1]-P[2]), np.linalg.norm(P[0]-P[2]),
                      np.linalg.norm(P[0]-P[1])])
        sig = 0.30
        h = np.array([sum(np.exp(-np.linalg.norm(P[i]-P[j])**2/(2*sig**2))
                          for j in range(3)) for i in range(3)])
        rec = recover(h, sig)
        e = np.abs(rec - R).max(); worst = max(worst, e)
        print(f"  {frac:7.2f} {np.array2string(R, precision=5):>34} "
              f"{np.array2string(rec, precision=5):>34} {e:10.2e}")
    print(f"\n  worst error over these instants: {worst:.2e}")

    print("\n  Using the true BIRTH values from the merge tree, where the maxima")
    print("  are displaced from the bodies, so the recovery is approximate:\n")
    print(f"  {'sigma':>7} {'worst relative error in the recovered sides':>46}")
    for sig in (0.22, 0.26, 0.30):
        worst = 0.0
        for frac in np.linspace(0.02, 0.98, 40):
            k = int(frac * NS); P = X[k]
            R = np.array([np.linalg.norm(P[1]-P[2]), np.linalg.norm(P[0]-P[2]),
                          np.linalg.norm(P[0]-P[1])])
            d = merge_tree_exact(P, MS, sig)
            if d is None:
                continue
            h = np.array([d[i]["birth"] for i in range(3)])
            rec = recover(h, sig)
            if rec is None:
                continue
            worst = max(worst, float(np.abs(rec - R).max() / R.max()))
        print(f"  {sig:7.2f} {worst:46.3e}")

    print("\n" + "=" * 74)
    print("2. Is the strand order the ordering of the side lengths?")
    print("=" * 74)
    agree = tot = 0
    for k in range(0, NS, 7):
        P = X[k]
        opp = np.array([np.linalg.norm(P[1]-P[2]), np.linalg.norm(P[0]-P[2]),
                        np.linalg.norm(P[0]-P[1])])     # side opposite body i
        d = merge_tree_exact(P, MS, 0.30)
        if d is None:
            continue
        h = np.array([d[i]["birth"] for i in range(3)])
        tot += 1
        if tuple(np.argsort(h)) == tuple(np.argsort(opp)):
            agree += 1
    print(f"  order of the peak heights equals order of the opposite sides")
    print(f"  at {agree} of {tot} instants ({100*agree/max(tot,1):.1f} per cent)")

    print("\n" + "=" * 74)
    print("3. Why this is special to three bodies.")
    print("=" * 74)
    print(f"  {'N':>4} {'heights':>9} {'distances needed':>18}  {'determined?':>12}")
    for N in range(3, 8):
        print(f"  {N:4d} {N:9d} {N*(N-1)//2:18d}  "
              f"{'yes' if N == N*(N-1)//2 else 'no':>12}")
    print("\n  N = N(N-1)/2 only at N = 3, so the peak heights are a complete")
    print("  shape invariant for the three-body problem and for nothing else.")


if __name__ == "__main__":
    main()
