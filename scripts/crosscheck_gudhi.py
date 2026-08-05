"""
Cross-check the two core persistence routines against GUDHI.

Every number in the write-up runs through mono/, so a bug in
`merge_tree_exact` or `extended_persistence_circle` would propagate everywhere
and no amount of internal consistency would catch it. GUDHI is an independent
implementation using a different algorithm: a cubical complex on a sampled grid,
against my analytic Newton refinement plus an elder-rule sweep.

Two comparisons.

  A  the planar density used for the figure-eight. Superlevel H_0 of rho is
     sublevel H_0 of -rho, so GUDHI's cubical complex on -rho, negated back,
     should converge to the exact merge tree as the grid refines.
  B  the circle machinery used for the coastline, the X-ray scan and the
     rotating waves. A function on S^1, with GUDHI's periodic cubical complex.

Convergence with grid refinement is the test. The exact routine returns machine
precision values; the cubical complex can only resolve them to the grid, so the
right signature is an error falling like the grid spacing.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import gudhi
from mono.nbody import FIGURE8, integrate
from mono.kde_exact import merge_tree_exact
from mono.core import extended_persistence_circle

F8 = FIGURE8
MS, T = F8["m"], F8["T"]


def density_grid(P, masses, sigma, lim, n):
    g = np.linspace(-lim, lim, n)
    X, Y = np.meshgrid(g, g, indexing="ij")
    F = np.zeros_like(X)
    for p, m in zip(P, masses):
        F += m * np.exp(-((X - p[0])**2 + (Y - p[1])**2) / (2 * sigma**2))
    return F


def gudhi_superlevel_h0(F):
    """Superlevel H_0 pairs (birth, death) via GUDHI on -F."""
    cc = gudhi.CubicalComplex(top_dimensional_cells=(-F).flatten(),
                              dimensions=list(F.shape))
    cc.compute_persistence(homology_coeff_field=2)
    out = []
    for dim, (b, d) in cc.persistence():
        if dim != 0:
            continue
        out.append((-b, -d))          # back to superlevel convention
    return out


def main():
    print(__doc__)
    t, X, V = integrate(F8["x"], F8["v"], MS, T, 3000)

    print("=" * 76)
    print("A. Planar density, GUDHI cubical complex against the exact merge tree")
    print("=" * 76)
    for frac in (0.10, 0.37, 0.62):
        P = X[int(frac * 3000)]
        sig = 0.30
        d = merge_tree_exact(P, MS, sig)
        mine_b = np.sort([d[i]["birth"] for i in range(3)])
        mine_d = np.sort([d[i]["death"] for i in range(3)
                          if d[i]["death"] is not None])
        print(f"\n  t/T = {frac}")
        print(f"    exact births {np.array2string(mine_b, precision=8)}")
        print(f"    exact deaths {np.array2string(mine_d, precision=8)}")
        print(f"    {'grid':>6} {'max |birth diff|':>18} {'max |death diff|':>18} "
              f"{'spacing':>10}")
        for n in (400, 800, 1600):
            F = density_grid(P, MS, sig, 3.0, n)
            pr = gudhi_superlevel_h0(F)
            fin = sorted([b for b, dd in pr if np.isfinite(dd)], reverse=True)
            dea = sorted([dd for b, dd in pr if np.isfinite(dd)])
            ess = [b for b, dd in pr if not np.isfinite(dd)]
            gb = np.sort(np.array(fin + ess))
            gd = np.sort(np.array(dea))
            eb = np.abs(gb - mine_b).max() if len(gb) == 3 else np.nan
            ed = np.abs(gd - mine_d).max() if len(gd) == 2 else np.nan
            print(f"    {n:6d} {eb:18.2e} {ed:18.2e} {6.0/n:10.4f}")

    print("\n" + "=" * 76)
    print("B. Function on the circle, GUDHI periodic cubical against mono.core")
    print("=" * 76)
    rng = np.random.default_rng(4)
    for trial in range(3):
        M = 2000
        th = np.linspace(0, 2*np.pi, M, endpoint=False)
        f = np.zeros(M)
        for k in (1, 2, 3, 5):
            f += rng.normal(0, 1) * np.cos(k*th + rng.uniform(0, 2*np.pi))
        r = extended_persistence_circle(f)
        allp = np.sort([f[r["maxs"][j]] - f[r["mins"][i]] for (i, j, _) in r["pairs"]])
        # mono reports the ESSENTIAL class as a pair, global max to global min,
        # while gudhi reports it separately with infinite death. Drop it before
        # comparing, or the two lists are offset by one and look like they
        # disagree when they do not.
        ess_val = float(f.max() - f.min())
        mine = np.sort([v for v in allp if abs(v - ess_val) > 1e-9])
        cc = gudhi.PeriodicCubicalComplex(top_dimensional_cells=(-f).flatten(),
                                          dimensions=[M],
                                          periodic_dimensions=[True])
        cc.compute_persistence(homology_coeff_field=2)
        pers = [(-b, -d) for dim, (b, d) in cc.persistence() if dim == 0]
        theirs = np.sort([b - d for b, d in pers if np.isfinite(d)])
        n_ess = sum(1 for b, d in pers if not np.isfinite(d))
        print(f"\n  trial {trial}: mono finds {len(mine)} pairs, "
              f"gudhi finds {len(theirs)} finite + {n_ess} essential")
        print(f"    mono essential class {ess_val:.6f} = global range, "
              f"gudhi reports {n_ess} essential")
        if len(mine) == len(theirs):
            print(f"    {len(mine)} finite pairs, agree to "
                  f"{np.abs(mine - theirs).max():.2e}")
            print(f"    mono   {np.array2string(mine, precision=8)}")
            print(f"    gudhi  {np.array2string(theirs, precision=8)}")
        else:
            print(f"    COUNT MISMATCH: mono {len(mine)} finite, "
                  f"gudhi {len(theirs)} finite")


if __name__ == "__main__":
    main()
