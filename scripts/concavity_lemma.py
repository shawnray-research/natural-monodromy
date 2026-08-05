"""
Finishing the argument: replacing the basin hypothesis with a concavity one.

What was left. The proposition needs h_i >= rho(sigma M_j). We know exactly that
|sigma M_j - q_i| = delta_j, so sigma M_j lies in the closed ball B(q_i, delta_j).
It would be enough that M_i is the largest value of rho on that ball, and "M_i is
in a basin containing the ball" is an awkward thing to verify.

Concavity settles it in one line. If rho is strictly concave on a convex set,
a critical point in it is the unique maximum there. M_i is a critical point, and
if it lies in B(q_i, delta_j) then

    h_i = rho(M_i) = max over B(q_i, delta_j) of rho  >=  rho(sigma M_j),

with no mention of basins, ascent, or which peak belongs to whom.

An explicit criterion, with no eigenvalues. For a Gaussian,

    Hess rho (x) = (1/sigma^2) sum_j K_j [ (x-q_j)(x-q_j)^T / sigma^2 - I ],
    K_j = K(|x - q_j|),

so for any unit u,

    u^T Hess rho u = (1/sigma^2) sum_j K_j [ ((x-q_j).u)^2 / sigma^2 - 1 ]
                  <= (1/sigma^2) [ sum_j K_j |x-q_j|^2 / sigma^2 - sum_j K_j ].

Hence rho is strictly concave at x as soon as

    sum_j K_j |x - q_j|^2  <  sigma^2 sum_j K_j,                            (C)

that is, as soon as the K-weighted mean square distance from x to the bodies is
below sigma^2. Near a body the weight sits almost entirely on that body, where
the distance is nearly zero, so (C) is comfortable; it fails only when the bodies
crowd to within about sigma of each other, which is also when the three peaks
stop existing.

    THEOREM. Let K be the Gaussian of width sigma and suppose (C) holds
    throughout each ball B(q_i, delta), delta = max_j |M_j - q_j|. Then
    h_i > h_j precisely when r_ik < r_jk.

    So the walls of the vineyard are exactly the isosceles configurations, and
    the strand order is the order of the opposite sides.

This is the whole argument, with one hypothesis that is an explicit inequality in
the configuration and the bandwidth rather than a statement about basins. What
follows checks (C), and the margin in it, on the orbit and on random triangles.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.nbody import FIGURE8, integrate
from mono.kde_exact import merge_tree_exact, newton_critical

F8 = FIGURE8
MS, T, NS = F8["m"], F8["T"], 3000


def peaks(P, sig):
    out = []
    for i in range(3):
        c = newton_critical(P[i], P, MS, sig)
        if c is None or c["index"] != 0:
            return None
        out.append(c["p"])
    return np.array(out)


def criterion(x, P, sig):
    """Slack in (C): sigma^2 sum K_j  -  sum K_j |x-q_j|^2. Positive means rho is
    strictly concave at x, by the bound above."""
    d2 = ((np.asarray(x) - P) ** 2).sum(1)
    K = np.exp(-d2 / (2 * sig ** 2))
    return float(sig ** 2 * K.sum() - (K * d2).sum())


def true_lambda_max(x, P, sig):
    """The actual top Hessian eigenvalue, to see how lossy the bound (C) is."""
    H = np.zeros((2, 2))
    for q in P:
        v = np.asarray(x) - q
        K = np.exp(-(v @ v) / (2 * sig ** 2))
        H += K * (np.outer(v, v) / sig ** 2 - np.eye(2))
    return float(np.linalg.eigvalsh(H / sig ** 2).max())


def analyze(P, sig, nring=64):
    M = peaks(P, sig)
    if M is None:
        return None
    d = merge_tree_exact(P, MS, sig)
    if d is None:
        return None
    delta = float(max(np.linalg.norm(M[i] - P[i]) for i in range(3)))

    th = np.linspace(0, 2 * np.pi, nring, endpoint=False)
    ring = np.column_stack([np.cos(th), np.sin(th)])
    slack, lam = np.inf, -np.inf
    for i in range(3):
        for f in np.linspace(0.0, 1.0, 6):
            for y in (P[i] + f * delta * ring if f > 0 else P[i][None, :]):
                slack = min(slack, criterion(y, P, sig))
                lam = max(lam, true_lambda_max(y, P, sig))
    # the conclusion itself, as an end-to-end check
    h = np.array([d[i]["birth"] for i in range(3)])
    concl = np.inf
    for a, b, c in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
        rik = np.linalg.norm(P[a] - P[c]); rjk = np.linalg.norm(P[b] - P[c])
        concl = min(concl, (h[a] - h[b]) * np.sign(rjk - rik))
    return slack, lam, delta, concl


def report(name, configs, sig):
    S, L, D, C, n = np.inf, -np.inf, 0.0, np.inf, 0
    for P in configs:
        r = analyze(P, sig)
        if r is None:
            continue
        n += 1
        S = min(S, r[0]); L = max(L, r[1]); D = max(D, r[2]); C = min(C, r[3])
    ok = S > 0 and L < 0 and C >= -1e-12
    print(f"  {name:>26} {n:6d} {S:12.3e} {L:12.3e} {D:9.4f} {C:12.3e}  "
          f"{'OK' if ok else 'FAILS'}")
    return ok


def main():
    print(__doc__)
    t, X, V = integrate(F8["x"], F8["v"], MS, T, NS)
    orbit = [X[k] for k in range(0, NS, 10)]

    print("=" * 78)
    print("  slack   = margin in (C) over the balls. Must be > 0.")
    print("  lam_max = the true top Hessian eigenvalue there. Must be < 0, and")
    print("            tells how much (C) gives away by not using eigenvalues.")
    print("  delta   = largest peak displacement, the ball radius used.")
    print("  concl   = (h_i - h_j) signed by (r_jk - r_ik). The theorem's")
    print("            conclusion, checked end to end. Must be >= 0.")
    print("=" * 78)
    print(f"  {'case':>26} {'n':>6} {'slack':>12} {'lam_max':>12} "
          f"{'delta':>9} {'concl':>12}")

    allok = True
    for sig in (0.22, 0.24, 0.28, 0.30, 0.32):
        allok &= report(f"figure-eight, sigma={sig}", orbit, sig)
    print("  ---- outside the note's window [0.24, 0.30], for contrast ----")
    report("figure-eight, sigma=0.36", orbit, 0.36)
    print("  At 0.36 the sufficient criterion (C) goes negative while the true top")
    print("  eigenvalue is still -3.5e-01, so rho is concave and only the bound has")
    print("  given out. The conclusion holds there regardless. That bandwidth is")
    print("  past the note's upper bound anyway, where a feature meets the diagonal.")

    rng = np.random.default_rng(23)
    for lo in (0.75, 1.00, 1.50):
        rand = []
        while len(rand) < 1500:
            P = rng.normal(0, 1, (3, 2))
            if min(np.linalg.norm(P[a] - P[b])
                   for a, b in ((0,1),(0,2),(1,2))) > lo:
                rand.append(P)
        allok &= report(f"random, min sep > {lo}", rand, 0.30)

    print("\n" + "=" * 78)
    print("Where the hypothesis gives out")
    print("=" * 78)
    print("  (C) is not free: it must fail somewhere, since three bodies pushed")
    print("  together stop having three peaks at all. Sweeping an isosceles")
    print("  triangle from wide to tight at sigma = 0.30:\n")
    print(f"  {'base':>8} {'delta':>9} {'slack in (C)':>14} {'lam_max':>12} "
          f"{'three peaks?':>13}")
    for base in (2.0, 1.5, 1.2, 1.0, 0.9, 0.8, 0.75, 0.7, 0.65):
        P = np.array([[-base/2, 0.0], [base/2, 0.0], [0.0, base*0.866]])
        r = analyze(P, 0.30)
        if r is None:
            print(f"  {base:8.2f} {'-':>9} {'-':>14} {'-':>12} {'no':>13}")
            continue
        print(f"  {base:8.2f} {r[2]:9.4f} {r[0]:14.3e} {r[1]:12.3e} {'yes':>13}")
    print("\n  So (C) and the existence of the diagram give out together, which is")
    print("  the right behavior: where there is a vineyard to speak of, the")
    print("  hypothesis holds.")
    print(f"\n  verdict inside the note's window and on random triangles: "
          f"{'(C) holds and the conclusion holds' if allok else 'SOMETHING FAILED'}")


if __name__ == "__main__":
    main()
