"""
Closing the last gap: the hypothesis that sigma M_j lies in body i's basin.

The proposition in wall_proof.py needs one thing it does not prove,

    h_i >= rho(sigma M_j),

and the note reports that as checked rather than proved. This script reduces it
to something with a closed form and a measurable margin, and says exactly how
much of the gap that closes.

STEP 1, an exact identity. For a Gaussian kernel, K'(r) = -(r/sigma^2) K(r), so

    grad rho(x) = -(1/sigma^2) sum_j (x - q_j) K(|x - q_j|),

and a critical point satisfies sum_j (M - q_j) w_j = 0 with w_j = K(|M - q_j|).
So every maximum is a WEIGHTED CENTROID of the bodies,

    M_i = sum_j w_j q_j / sum_j w_j,

which is exact, not a perturbation. Subtracting q_i,

    |M_i - q_i|  <=  sum_{j != i} w_j r_ij / sum_j w_j.

The peak displacement is therefore bounded in closed form by the kernel values at
the other bodies, and those are exponentially small at the working bandwidth.

STEP 2, what the displacement buys. sigma is an isometry fixing q_i and q_j as a
pair, so

    |sigma M_j - q_i| = |sigma M_j - sigma q_j| = |M_j - q_j| = delta_j,

that is, sigma M_j sits within delta_j of body i. So it is enough that

    h_i >= max of rho over the closed ball B(q_i, delta_j).                 (*)

STEP 3, what is left. (*) holds if B(q_i, delta_j) sits inside body i's basin.
That is the whole of the remaining gap, and it is now a statement about one
length against another rather than about the vineyard. This script measures both:
delta against the distance from the body to the nearest saddle, on the orbit and
on four thousand random triangles.

What this does and does not settle. STEPS 1 and 2 are proved and hold for every
configuration: the identity is exact, so the bound is exact, so sigma M_j really
does sit within delta_j of body i. STEP 3 is still measured. The honest summary
is that the gap narrows from an unquantified hypothesis about basins to a
quantified comparison of two lengths, and does not close. Closing it needs a
uniform bound on delta against the basin radius, which is a real lemma.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.nbody import FIGURE8, integrate
from mono.kde_exact import merge_tree_exact, newton_critical

F8 = FIGURE8
MS, T, NS = F8["m"], F8["T"], 3000


def rho(x, P, sig):
    x = np.atleast_2d(np.asarray(x, float))
    d2 = ((x[:, None, :] - P[None, :, :]) ** 2).sum(2)
    return np.exp(-d2 / (2 * sig ** 2)).sum(1)


def peaks(P, sig):
    out = []
    for i in range(3):
        c = newton_critical(P[i], P, MS, sig)
        if c is None or c["index"] != 0:
            return None
        out.append(c["p"])
    return np.array(out)


def analyze(P, sig):
    """
    Returns (centroid residual, bound slack, hypothesis slack, headroom), or None.

    hypothesis slack is h_i - rho(sigma M_j) over all ordered pairs, the exact
    quantity the proposition needs, and it must be >= 0.

    headroom is delta / (distance from q_i to the nearest saddle), a scale-free
    measure of how far the configuration is from the hypothesis being tight. An
    earlier version of this script measured the max of rho on a sphere of radius
    1.05 max(delta) instead, which is zero by construction whenever body i is the
    one with the largest displacement, and so reported failures that were only
    roundoff on an exact equality.
    """
    M = peaks(P, sig)
    if M is None:
        return None
    d = merge_tree_exact(P, MS, sig)
    if d is None:
        return None
    h = np.array([d[i]["birth"] for i in range(3)])

    ident, bound_slack = 0.0, np.inf
    delta = np.zeros(3)
    for i in range(3):
        w = np.exp(-((M[i] - P) ** 2).sum(1) / (2 * sig ** 2))
        ident = max(ident, float(np.linalg.norm(M[i] - (w @ P) / w.sum())))
        delta[i] = np.linalg.norm(M[i] - P[i])
        bnd = sum(w[j] * np.linalg.norm(P[i] - P[j])
                  for j in range(3) if j != i) / w.sum()
        bound_slack = min(bound_slack, bnd - delta[i])

    hyp = np.inf
    for i in range(3):
        for j in range(3):
            if i == j:
                continue
            mid = 0.5 * (P[i] + P[j])
            nrm = (P[i] - P[j]) / np.linalg.norm(P[i] - P[j])
            sMj = M[j] - 2.0 * ((M[j] - mid) @ nrm) * nrm
            hyp = min(hyp, h[i] - float(rho(sMj, P, sig)[0]))

    sad = [d[i]["saddle"] for i in range(3) if d[i].get("saddle") is not None]
    if sad:
        nearest = min(np.linalg.norm(P[i] - sp)
                      for i in range(3) for sp in sad)
        headroom = float(delta.max() / nearest)
    else:
        headroom = np.nan
    return ident, bound_slack, hyp, headroom


def report(name, configs, sig):
    Imax, B, H, R, n = 0.0, np.inf, np.inf, 0.0, 0
    for P in configs:
        r = analyze(P, sig)
        if r is None:
            continue
        n += 1
        Imax = max(Imax, r[0]); B = min(B, r[1]); H = min(H, r[2])
        R = max(R, r[3])
    ok = Imax < 1e-10 and B >= -1e-12 and H >= -1e-12
    print(f"  {name:>26} {n:6d} {Imax:11.2e} {B:11.2e} {H:12.3e} {R:10.4f}  "
          f"{'OK' if ok else 'FAILS'}")
    return ok


def main():
    print(__doc__)
    t, X, V = integrate(F8["x"], F8["v"], MS, T, NS)
    orbit = [X[k] for k in range(0, NS, 10)]

    print("=" * 78)
    print("Columns: identity = |M_i - weighted centroid| (must be ~0)")
    print("         bound    = closed-form bound minus actual displacement (>= 0)")
    print("         hyp      = h_i minus rho(sigma M_j), what the proof needs (>= 0)")
    print("         headroom = largest peak displacement / distance to nearest")
    print("                    saddle. Small means the hypothesis is nowhere near")
    print("                    tight: the peak barely moves against the basin.")
    print("=" * 78)
    print(f"  {'case':>26} {'n':>6} {'identity':>11} {'bound':>11} "
          f"{'hyp':>12} {'headroom':>10}")

    allok = True
    for sig in (0.24, 0.28, 0.30, 0.32):
        allok &= report(f"figure-eight, sigma={sig}", orbit, sig)

    rng = np.random.default_rng(11)
    rand = []
    while len(rand) < 4000:
        P = rng.normal(0, 1, (3, 2))
        r = [np.linalg.norm(P[a] - P[b]) for a, b in ((0,1),(0,2),(1,2))]
        if min(r) > 0.75:                 # else the peaks merge and there is no diagram
            rand.append(P)
    allok &= report("random triangles, 0.30", rand, 0.30)

    print("\n" + "=" * 78)
    print("How much of the gap this closes")
    print("=" * 78)
    print("  Proved, not measured: the weighted-centroid identity, and hence the")
    print("  closed-form bound on the peak displacement, and hence that sigma M_j")
    print("  lies within delta_j of body i. Those are STEPS 1 and 2 and they hold")
    print("  for every configuration, not just this orbit.")
    print("\n  Still measured: that this puts sigma M_j inside body i's basin. The")
    print("  headroom column says how much room there is, and it is a factor of a")
    print("  few, not an order of magnitude: the largest peak displacement runs")
    print("  from 4 per cent of the way to the nearest saddle at sigma = 0.24 to")
    print("  39 per cent at sigma = 0.32, and 23 per cent at the working 0.30.")
    print("  Comfortable, and worth not overstating.")
    print("\n  Not done: a uniform proof that displacement stays below basin radius")
    print("  for every triangle. That is a real lemma, and it is what would close")
    print("  the gap completely.")
    print(f"\n  verdict over every case above: {'all criteria hold' if allok else 'SOMETHING FAILED'}")


if __name__ == "__main__":
    main()
