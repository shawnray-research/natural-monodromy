"""
Does the persistence diagram actually determine the configuration?

complete_invariant.py showed that the IDEALIZED heights h_i = sum_j K(r_ij),
which is the density read at the bodies, invert to give the three mutual
distances, exactly. It also showed, in its second table, that feeding the TRUE
birth values through the same inversion leaves a 33 per cent relative error at
every bandwidth in the working window, because the maxima are displaced from the
bodies. So the exact statement is about a model of the diagram, not the diagram.

That leaves the real question open, and it is the interesting one: never mind my
particular formula, is the map

    configuration modulo congruence  ->  the three birth values

injective on the orbit, and if so how well conditioned? Two tests that do not
depend on any inversion formula.

  A  convergence. Does the model inversion error fall to zero as the bandwidth
     falls, and does it do so anywhere near the window the note actually uses,
     sigma in [0.24, 0.30]?

  B  injectivity, measured directly. Take every pair of instants around the
     orbit. Compare the distance between their birth triples with the distance
     between their triangles, taken up to congruence as sorted side lengths. If
     the map is injective and well conditioned there is a positive lower bound
     on the ratio, a bi-Lipschitz constant. If some pair has nearly equal births
     and clearly different shapes, it is not injective and the completeness
     claim is dead on the real diagram.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.nbody import FIGURE8, integrate
from mono.kde_exact import merge_tree_exact

F8 = FIGURE8
MS, T, NS = F8["m"], F8["T"], 6000


def sides(P):
    """Sorted side lengths: a complete invariant of the triangle up to congruence."""
    return np.sort([np.linalg.norm(P[1] - P[2]),
                    np.linalg.norm(P[0] - P[2]),
                    np.linalg.norm(P[0] - P[1])])


def recover(h, sigma, k0=1.0):
    g = np.asarray(h, float) - k0
    S = 0.5 * g.sum()
    out = []
    for k in (S - g[0], S - g[1], S - g[2]):
        if k <= 0:
            return None
        out.append(sigma * np.sqrt(-2.0 * np.log(k)))
    return np.array(out)


def main():
    print(__doc__)
    t, X, V = integrate(F8["x"], F8["v"], MS, T, NS)

    print("=" * 78)
    print("A. Does the model inversion converge as the bandwidth falls?")
    print("=" * 78)
    print(f"  {'sigma':>7} {'worst rel. error':>18} {'in the note window?':>22}"
          f"   {'instants used':>14}")
    for sig in (0.40, 0.30, 0.24, 0.18, 0.12, 0.08, 0.05, 0.03):
        worst, n = 0.0, 0
        for frac in np.linspace(0.02, 0.98, 40):
            P = X[int(frac * NS)]
            d = merge_tree_exact(P, MS, sig)
            if d is None:
                continue
            h = np.array([d[i]["birth"] for i in range(3)])
            rec = recover(h, sig)
            if rec is None:
                continue
            R = np.array([np.linalg.norm(P[1]-P[2]), np.linalg.norm(P[0]-P[2]),
                          np.linalg.norm(P[0]-P[1])])
            worst = max(worst, float(np.abs(rec - R).max() / R.max())); n += 1
        inwin = "yes" if 0.24 <= sig <= 0.30 else "no"
        note = "" if n else "   <- NO VALID INSTANTS, the zero is empty"
        print(f"  {sig:7.2f} {worst:18.3e} {inwin:>22}   n={n:3d}{note}")
    print("\n  The note's window is [0.24, 0.30]. Read off whether the error is")
    print("  small anywhere inside it.")

    print("\n" + "=" * 78)
    print("B. Injectivity of the real diagram, with no inversion formula")
    print("=" * 78)
    for sig in (0.24, 0.30):
        B, S = [], []
        for k in range(0, NS, 6):
            d = merge_tree_exact(X[k], MS, sig)
            if d is None:
                continue
            B.append(np.sort([d[i]["birth"] for i in range(3)]))
            S.append(sides(X[k]))
        B, S = np.array(B), np.array(S)
        db = np.linalg.norm(B[:, None, :] - B[None, :, :], axis=2)
        ds = np.linalg.norm(S[:, None, :] - S[None, :, :], axis=2)
        iu = np.triu_indices(len(B), k=1)
        db, ds = db[iu], ds[iu]
        # scale-free: how far apart can two shapes be while their diagrams agree?
        keep = ds > 1e-9
        ratio = db[keep] / ds[keep]
        print(f"\n  sigma = {sig}, {len(B)} instants, {keep.sum()} distinct pairs")
        print(f"    range of shape separations   {ds.min():.3e} to {ds.max():.3e}")
        print(f"    range of diagram separations {db.min():.3e} to {db.max():.3e}")
        print(f"    lower Lipschitz ratio |dB|/|dS|: min {ratio.min():.3e}, "
              f"median {np.median(ratio):.3e}")
        # the damning case: shapes far apart whose diagrams nearly coincide
        far = ds > 0.10 * ds.max()
        if far.any():
            j = np.argmin(db[far] / 1.0)
            worst_db = db[far].min()
            print(f"    among pairs with shape separation above 10 per cent of the")
            print(f"    maximum, the SMALLEST diagram separation is {worst_db:.3e}")
            print(f"    -> {'COLLISION, not injective' if worst_db < 1e-6 else 'no collision, injective on this sample'}")

    print("\n" + "=" * 78)
    print("C. What survives exactly, regardless of the above")
    print("=" * 78)
    agree = tot = 0
    for k in range(0, NS, 7):
        P = X[k]
        opp = np.array([np.linalg.norm(P[1]-P[2]), np.linalg.norm(P[0]-P[2]),
                        np.linalg.norm(P[0]-P[1])])
        d = merge_tree_exact(P, MS, 0.30)
        if d is None:
            continue
        h = np.array([d[i]["birth"] for i in range(3)])
        tot += 1
        agree += tuple(np.argsort(h)) == tuple(np.argsort(opp))
    print(f"  ORDER of the births equals ORDER of the opposite sides: "
          f"{agree} of {tot}")
    print(f"  This is the fact the braid rests on, and it is exact. The braid")
    print(f"  needs only the ordering, never the values, so it is untouched by")
    print(f"  whatever A and B say about recovering the actual side lengths.")


if __name__ == "__main__":
    main()
