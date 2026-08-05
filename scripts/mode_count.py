"""
How many modes does the density actually have?

The number of modes of an isotropic Gaussian mixture is NOT the number of
components. Carreira-Perpinan and Williams conjectured that a homoscedastic
mixture of M components in more than one dimension has at most M modes, and
Duistermaat gave a counterexample: three equal Gaussians at the vertices of an
equilateral triangle in the plane have, for a range of variances, three modes
near the vertices AND a fourth at the center.

That is precisely this setup, and the figure-eight passes near the equilateral
Lagrange configuration. Worse, `exact_maxima` in mono/kde_exact.py seeds Newton
from the body positions and returns one maximum per body, so it CANNOT find a
mode at the centroid. If a fourth mode exists anywhere on the orbit inside the
bandwidth window, the merge tree is missing a maximum, the sublevel-set
connectivity is wrong, and the elder-rule pairing among the three that were found
can be wrong too. The braid word would inherit that.

There is a second reason to check. Earlier counting found either 2 or 3 saddles
depending on the instant, and the 3-saddle case was interpreted as three saddles
plus a central local MINIMUM. Four maxima and three saddles satisfies the same
Euler count, so that interpretation was an assumption rather than a measurement.

This script does not seed from the bodies. It runs Newton from every cell of a
dense grid covering the configuration, deduplicates the limit points, and
classifies each by the sign of the Hessian eigenvalues.
"""
import sys, os, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.nbody import FIGURE8, integrate
from mono.kde_exact import newton_critical, rho_and_derivs

F8 = FIGURE8
MS, T = F8["m"], F8["T"]


def all_critical(pts, masses, sigma, pad=2.5, n=26, tol=1e-6):
    """
    Every critical point, found without any knowledge of where the bodies are.
    Newton from each cell of a grid covering the configuration plus a margin of
    `pad` bandwidths, then deduplicated.
    """
    lo = pts.min(axis=0) - pad * sigma
    hi = pts.max(axis=0) + pad * sigma
    xs = np.linspace(lo[0], hi[0], n)
    ys = np.linspace(lo[1], hi[1], n)
    found = []
    for x in xs:
        for y in ys:
            c = newton_critical(np.array([x, y]), pts, masses, sigma)
            if c is None or c["grad"] > 1e-8:
                continue
            p = c["p"]
            if not (lo[0] - sigma <= p[0] <= hi[0] + sigma
                    and lo[1] - sigma <= p[1] <= hi[1] + sigma):
                continue
            if any(np.linalg.norm(p - q["p"]) < tol for q in found):
                continue
            found.append(c)
    nmax = sum(1 for c in found if c["index"] == 0)
    nsad = sum(1 for c in found if c["index"] == 1)
    nmin = sum(1 for c in found if c["index"] == 2)
    return found, nmax, nsad, nmin


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", type=int, default=120)
    ap.add_argument("--grid", type=int, default=26)
    a = ap.parse_args()

    print(__doc__)
    NS = 3000
    t, X, V = integrate(F8["x"], F8["v"], MS, T, NS)

    print("=" * 76)
    print("Control: Duistermaat's configuration, to show the finder can see a")
    print("fourth mode when there is one.")
    print("=" * 76)
    eq = np.array([[np.cos(a_), np.sin(a_)] for a_ in
                   (np.pi/2, np.pi/2 + 2*np.pi/3, np.pi/2 + 4*np.pi/3)])
    print(f"  {'sigma':>7} {'maxima':>7} {'saddles':>8} {'minima':>7}   "
          f"{'center is':>12}   Morse")
    for sig in (0.40, 0.50, 0.60, 0.70, 0.80, 0.90):
        found, nmx, nsd, nmn = all_critical(eq, np.ones(3), sig, n=a.grid)
        ctr = [c for c in found if np.linalg.norm(c["p"]) < 1e-4]
        lab = {0: "a MAXIMUM", 1: "a saddle", 2: "a minimum"}.get(
            ctr[0]["index"], "-") if ctr else "nothing"
        print(f"  {sig:7.2f} {nmx:7d} {nsd:8d} {nmn:7d}   {lab:>12}   "
              f"{nmx - nsd + nmn}")

    print("\n" + "=" * 76)
    print("The figure-eight, along the whole orbit, at each bandwidth.")
    print("=" * 76)
    ks = np.linspace(0, NS - 1, a.frames).astype(int)
    print(f"  {'sigma':>7} {'max':>16} {'saddles':>12} {'minima':>10}  "
          f"{'Morse':>7}  {'worst frame':>12}")
    for sig in (0.22, 0.26, 0.30, 0.34, 0.38, 0.45):
        mx, sd, mn, bad = set(), set(), set(), None
        morse = set()
        for k in ks:
            found, nmx, nsd, nmn = all_critical(X[k], MS, sig, n=a.grid)
            mx.add(nmx); sd.add(nsd); mn.add(nmn)
            # per-FRAME Morse count. Taking the Cartesian product of the sets
            # collected across frames, as an earlier version did, invents
            # combinations that never occur at any single instant and made the
            # finder look broken at large sigma when it was not.
            morse.add(nmx - nsd + nmn)
            if nmx != 3 and bad is None:
                bad = (k / NS, nmx, nsd, nmn)
        print(f"  {sig:7.2f} {str(sorted(mx)):>16} {str(sorted(sd)):>12} "
              f"{str(sorted(mn)):>10}  {str(sorted(morse)):>7}  "
              f"{'-' if bad is None else f't/T={bad[0]:.3f} {bad[1]}max'}")

    print("\n  How close does the figure-eight come to the equilateral shape?")
    worst = 1e9
    for k in range(NS):
        P = X[k]
        r = np.array([np.linalg.norm(P[0]-P[1]), np.linalg.norm(P[0]-P[2]),
                      np.linalg.norm(P[1]-P[2])])
        worst = min(worst, (r.max() - r.min()) / r.mean())
    print(f"    minimum relative spread of the three mutual distances: {worst:.4f}")
    print("    (0 would be exactly equilateral)")


if __name__ == "__main__":
    main()
