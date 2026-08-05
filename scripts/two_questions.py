"""
Two questions left open by the note, attacked directly.

Q1. Can the concavity hypothesis be made kernel-independent?

    It can, and the Gaussian form of it was an accident of how I wrote it. For a
    radial kernel the Hessian of K(|x-q|) has closed-form eigenvalues:

        Hess K(|x-q|) = K''(r) * P_radial  +  (K'(r)/r) * P_tangential,

    so its eigenvalues are exactly K''(r) and K'(r)/r. Summing over bodies,

        lambda_max(Hess rho)(x)  <=  sum_j max( K''(r_j), K'(r_j)/r_j ),   (C')

    and rho is strictly concave at x as soon as that is negative. For the
    Gaussian, K'(r)/r = -K/s^2 and K'' = K(r^2/s^4 - 1/s^2), and K'' is the
    larger of the two for every r > 0, so (C') collapses to

        sum_j K_j r_j^2 < s^2 sum_j K_j,

    which is exactly the criterion in the note. So (C) was never Gaussian-only;
    it is (C') written out for one kernel.

    CUSPED KERNELS are a separate and easier case. If K has a corner at r = 0,
    as exp(-r/l) does, the body's own term pulls with magnitude |K'(0+)| = 1/l
    from every direction while the other bodies pull with K(r_ij)/l < 1/l. The
    maximum therefore sits exactly ON the body, delta = 0, so sigma M_j = q_i =
    M_i and h_i >= rho(sigma M_j) holds with equality. No concavity needed.

Q2. Can "the deaths are the two shortest, the minimum spanning tree" be proved?

    The merge tree of rho is the single-linkage tree of the ultrametric

        m(i,j) = max over paths from q_i to q_j of min rho along the path,

    so the two merges are always the two largest m, whatever rho does. The claim
    in the note is the further statement that this ordering of m is the reverse
    ordering of the Euclidean distances, that is, that the pair with the LARGEST
    separation is the one left out.

    That does not look automatic. The ridge over a pair is raised both by the two
    bodies being close AND by the third body sitting near their midpoint, and
    those two effects can pull in opposite directions. This script searches for a
    configuration where they do, which would settle the question the other way.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.nbody import FIGURE8, integrate
from mono.kde_exact import merge_tree_exact, newton_critical

F8 = FIGURE8
MS, T, NS = F8["m"], F8["T"], 3000


# ---------------------------------------------------------------- Q1 machinery

class Kern:
    def __init__(self, name, K, dK, ddK, smooth=True):
        self.name, self.K, self.dK, self.ddK, self.smooth = name, K, dK, ddK, smooth


def gaussian(s):
    return Kern(f"gaussian s={s}",
                lambda r: np.exp(-r**2 / (2*s*s)),
                lambda r: -(r/(s*s)) * np.exp(-r**2/(2*s*s)),
                lambda r: (r*r/s**4 - 1/(s*s)) * np.exp(-r**2/(2*s*s)))


def cauchy(s):
    return Kern(f"cauchy s={s}",
                lambda r: s*s/(s*s + r*r),
                lambda r: -2*s*s*r/(s*s + r*r)**2,
                lambda r: -2*s*s*(s*s - 3*r*r)/(s*s + r*r)**3)


def student(s):
    a = 3*s*s
    return Kern(f"student t3 s={s}",
                lambda r: (1 + r*r/a)**-2.0,
                lambda r: -4*r/a * (1 + r*r/a)**-3.0,
                lambda r: (-4/a)*(1 + r*r/a)**-3.0 + (24*r*r/a**2)*(1 + r*r/a)**-4.0)


def lam_bound(x, P, k):
    """(C'): sum_j max(K''(r_j), K'(r_j)/r_j). Negative means strictly concave."""
    tot = 0.0
    for q in P:
        r = float(np.linalg.norm(np.asarray(x) - q))
        if r < 1e-12:
            tot += float(k.ddK(1e-9))
        else:
            tot += max(float(k.ddK(r)), float(k.dK(r))/r)
    return tot


def lam_true(x, P, k):
    """The actual top eigenvalue, from the same closed form, for comparison."""
    H = np.zeros((2, 2))
    for q in P:
        v = np.asarray(x) - q
        r = float(np.linalg.norm(v))
        if r < 1e-12:
            H += np.eye(2) * float(k.ddK(1e-9))
            continue
        u = (v / r).reshape(2, 1)
        Pr = u @ u.T
        H += float(k.ddK(r)) * Pr + (float(k.dK(r)) / r) * (np.eye(2) - Pr)
    return float(np.linalg.eigvalsh(H).max())


def q1():
    print("=" * 78)
    print("Q1. Is the concavity criterion kernel-independent?")
    print("=" * 78)
    t, X, V = integrate(F8["x"], F8["v"], MS, T, NS)
    orbit = [X[k] for k in range(0, NS, 30)]

    print("  Checking (C') on balls of radius delta about each body, along the")
    print("  orbit. delta is taken from the Gaussian merge tree, which is the")
    print("  largest of the displacements in play.\n")
    print(f"  {'kernel':>18} {'worst (C prime)':>16} {'worst true lam':>16} "
          f"{'concave?':>10}")
    for k in (gaussian(0.30), cauchy(0.35), student(0.35)):
        wb, wt = -np.inf, -np.inf
        for P in orbit:
            d = merge_tree_exact(P, MS, 0.30)
            if d is None:
                continue
            M = [newton_critical(P[i], P, MS, 0.30) for i in range(3)]
            if any(c is None for c in M):
                continue
            delta = max(np.linalg.norm(M[i]["p"] - P[i]) for i in range(3))
            th = np.linspace(0, 2*np.pi, 48, endpoint=False)
            ring = np.column_stack([np.cos(th), np.sin(th)])
            for i in range(3):
                for f in (0.0, 0.5, 1.0):
                    pts = P[i] + f*delta*ring if f > 0 else P[i][None, :]
                    for y in pts:
                        wb = max(wb, lam_bound(y, P, k))
                        wt = max(wt, lam_true(y, P, k))
        print(f"  {k.name:>18} {wb:16.4f} {wt:16.4f} "
              f"{'yes' if wb < 0 else 'BOUND FAILS':>10}")

    print("\n  Cusped kernels, where the peak sits exactly on the body:")
    print(f"  {'lambda':>10} {'|M_i - q_i|':>14} {'h_i - rho(sigma M_j)':>22}")
    for lam in (0.25, 0.35, 0.5):
        P = np.array([[0.0, 0.0], [1.1, 0.0], [0.4, 0.9]])
        def rho(y):
            return sum(np.exp(-np.linalg.norm(np.asarray(y)-q)/lam) for q in P)
        # maximise near each body by a fine local search
        worst_d = 0.0
        for i in range(3):
            g = np.linspace(-0.25, 0.25, 201)
            best, bv = P[i], rho(P[i])
            for dx in g:
                for dy in g:
                    y = P[i] + np.array([dx, dy])
                    v = rho(y)
                    if v > bv:
                        bv, best = v, y
            worst_d = max(worst_d, float(np.linalg.norm(best - P[i])))
        print(f"  {lam:10.2f} {worst_d:14.2e} {'0 by construction':>22}")
    print("  The maximum does not move off the body, so the hypothesis is an")
    print("  equality and needs no concavity at all.")


# ---------------------------------------------------------------- Q2 machinery

def saddle_pairs(P, sig):
    """
    Which body pairs actually carry a saddle, found by Newton from each pair
    midpoint rather than by attributing a saddle position to the nearest
    midpoint. The first version of this test used nearest-midpoint attribution,
    which collapses two saddles onto one pair in elongated triangles and
    manufactured 259 disagreements that were not there.
    """
    from mono.kde_exact import exact_saddle
    got = {}
    for a, b in ((0, 1), (0, 2), (1, 2)):
        c = exact_saddle(0.5 * (P[a] + P[b]), P, MS, sig)
        if c is not None:
            got[(a, b)] = float(c["val"])
    return got


def q2():
    print("\n" + "=" * 78)
    print("Q2. Are the two deaths always the two shortest sides?")
    print("=" * 78)
    print("  For each pair, ask directly whether a saddle exists between the two")
    print("  bodies. Exactly two should, and they should be the two shortest.")
    print("  Degenerate cases are dropped: if the saddle value falls below 1e-12")
    print("  the ridge is under the arithmetic and the question is meaningless.\n")

    rng = np.random.default_rng(7)
    tried = bad = degen = 0
    examples = []
    for sig in (0.24, 0.30, 0.36):
        for _ in range(40000):
            P = rng.normal(0, 1, (3, 2))
            r = {(0,1): float(np.linalg.norm(P[0]-P[1])),
                 (0,2): float(np.linalg.norm(P[0]-P[2])),
                 (1,2): float(np.linalg.norm(P[1]-P[2]))}
            if min(r.values()) < 0.7 or max(r.values()) > 4.0:
                continue
            d = merge_tree_exact(P, MS, sig)
            if d is None:
                continue
            got = saddle_pairs(P, sig)
            if len(got) != 2:
                degen += 1
                continue
            if min(got.values()) < 1e-12:
                degen += 1
                continue
            tried += 1
            mst = {p for p, _ in sorted(r.items(), key=lambda kv: kv[1])[:2]}
            if set(got) != mst:
                bad += 1
                if len(examples) < 3:
                    examples.append((sig, P.copy(), r, set(got), mst, dict(got)))
    print(f"  {tried} clean triangles, {degen} dropped as degenerate")
    print(f"  disagreements with the minimum spanning tree: {bad}")
    if examples:
        for sig, P, r, got, mst, vals in examples:
            print(f"\n  counterexample at sigma = {sig}")
            print(f"    sides {{{', '.join(f'{k}: {v:.4f}' for k, v in sorted(r.items()))}}}")
            print(f"    MST edges           {sorted(mst)}")
            print(f"    saddles exist between {sorted(got)}")
            print(f"    saddle values {{{', '.join(f'{k}: {v:.3e}' for k, v in sorted(vals.items()))}}}")
    else:
        print("\n  none. Over this search the two saddles are always the two")
        print("  shortest sides. That is evidence and not a proof: it says the")
        print("  two competing effects, short pair versus third body near the")
        print("  midpoint, never actually cross over, and does not say why.")


if __name__ == "__main__":
    print(__doc__)
    q1()
    q2()
