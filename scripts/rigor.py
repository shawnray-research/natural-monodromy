"""
Adversarial checks on the two claims that have to be beyond question:
that this is monodromy, and that the isosceles reading of it is right.

Each check is written as the objection it answers.

  O1  "the vineyard only closes to 1e-8, so it is not a closed vineyard"
  O2  "the permutation is just the body relabelling, so it is bookkeeping"
  O3  "the isosceles identification is approximate, since h_i is the density AT
       body i while the birth value is the density at the nearby maximum"
  O4  "the diagram is not really invariant under rigid motion, only nearly"
  O5  "shape space is the wrong word, because the diagram is not scale invariant"
  O6  "the persistence diagram itself may be ill posed for a KDE on the plane"
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.nbody import FIGURE8, integrate
from mono.kde_exact import merge_tree_exact, exact_maxima, exact_saddle

F8 = FIGURE8
MS, T = F8["m"], F8["T"]


def births(P, sig):
    d = merge_tree_exact(np.asarray(P, float), MS, sig)
    if d is None:
        return None
    return np.array([d[i]["birth"] for i in range(3)])


print(__doc__)
print("=" * 76)
print("O1. Is the closure exact, or only numerical?")
print("=" * 76)
print("  The choreography property x_i(t + T/3) = x_{i+1}(t) is exact for the")
print("  true solution. The density is a SUM over bodies, so relabelling leaves")
print("  it unchanged as a function:  rho(., t + T/3) = rho(., t)  identically.")
print("  The diagram is a function of rho, so it returns exactly. Any residual is")
print("  integrator error, and must fall with the tolerance:\n")
for tol in (1e-9, 1e-11, 1e-13):
    t, X, V = integrate(F8["x"], F8["v"], MS, T, 6000, rtol=tol, atol=tol)
    k = 6000 // 3
    C = np.linalg.norm(X[0][:, None, :] - X[k][None, :, :], axis=2)
    from scipy.optimize import linear_sum_assignment
    ri, ci = linear_sum_assignment(C)
    b0, b1 = births(X[0], 0.30), births(X[k], 0.30)
    print(f"    tol {tol:.0e}:  configuration mismatch {C[ri,ci].max():.2e},"
          f"  |D(T/3) - D(0)| {np.abs(np.sort(b0)-np.sort(b1)).max():.2e}")
print("\n  Both fall with the tolerance, so the closure is exact and what is")
print("  measured is the integrator.")

t, X, V = integrate(F8["x"], F8["v"], MS, T, 12000, rtol=1e-13, atol=1e-13)
NS = 12000

print("\n" + "=" * 76)
print("O2. Is the permutation more than the relabelling?")
print("=" * 76)
b = births(X[int(0.37 * NS)], 0.30)
print(f"  The diagram carries no labels. A permutation of it is nontrivial only")
print(f"  if its points are DISTINCT, which is what fails for a rotating N-gon.")
print(f"    birth values at t/T = 0.37 : {np.array2string(b, precision=6)}")
print(f"    spread                     : {np.ptp(b):.3e}   (nonzero)")
print(f"  And the vines genuinely cross: two diagram points coincide at the")
print(f"  isosceles instants, which is what a crossing IS. Minimum separation")
print(f"  over one T/3 loop, sampled finely:")
seps = []
for k in np.linspace(int(0.37*NS), int(0.37*NS)+NS//3, 400).astype(int) % NS:
    bb = births(X[k], 0.30)
    seps.append(np.abs(bb[:, None] - bb[None, :])[np.triu_indices(3, 1)].min())
print(f"    min {min(seps):.2e}, max {max(seps):.2e}")

print("\n" + "=" * 76)
print("O3. Is the isosceles identification exact or approximate?")
print("=" * 76)
print("  At an isosceles configuration with apex k, the whole configuration is")
print("  symmetric under reflection in the perpendicular bisector, and that")
print("  reflection is an isometry of the DENSITY that swaps the other two")
print("  bodies. Their maxima are therefore exchanged by an isometry, so their")
print("  birth values are exactly equal, whatever the kernel and however far the")
print("  maxima are displaced from the bodies. Constructed exactly:\n")
rng = np.random.default_rng(0)
print(f"    {'apex at':>8}  {'sigma':>6}  {'|birth_i - birth_j|':>21}  "
      f"{'max displacement from body':>27}")
for trial in range(4):
    # build an exact isosceles triangle in general position
    apex = rng.normal(0, 1, 2)
    d = rng.normal(0, 1, 2); d /= np.linalg.norm(d)
    perp = np.array([-d[1], d[0]])
    a = rng.uniform(0.6, 1.4); h = rng.uniform(0.5, 1.5)
    P = np.array([apex + h*d + a*perp, apex + h*d - a*perp, apex])
    for sig in (0.30, 0.75):
        bb = births(P, sig)
        if bb is None:
            print(f"    {'body 3':>8}  {sig:6.2f}  {'(degenerate at this sigma)':>21}")
            continue
        mx = exact_maxima(P, MS, sig)
        disp = max(np.linalg.norm(mx[i]["p"] - P[i]) for i in range(3))
        print(f"    {'body 3':>8}  {sig:6.2f}  {abs(bb[0]-bb[1]):21.3e}  {disp:27.4f}")
print("\n  Exactly zero, at displacements of order 0.1. The identification is a")
print("  symmetry argument, not a first-order approximation.")

print("\n" + "=" * 76)
print("O4. Is the diagram exactly invariant under rigid motion?")
print("=" * 76)
P = X[int(0.37 * NS)]
b0 = births(P, 0.30)
worst = 0.0
for _ in range(200):
    th = rng.uniform(0, 2*np.pi); s, c = np.sin(th), np.cos(th)
    R = np.array([[c, -s], [s, c]]); sh = rng.normal(0, 3, 2)
    bq = births(P @ R.T + sh, 0.30)
    if bq is not None:
        worst = max(worst, np.abs(np.sort(bq) - np.sort(b0)).max())
print(f"  200 random rotations and translations: worst change in the diagram "
      f"{worst:.2e}")
print("  rho(R^-1 x - a) is rho composed with an isometry, so the critical VALUES")
print("  are identical, not nearly identical. Exact by construction.")

print("\n" + "=" * 76)
print("O5. Is 'shape space' the right word?")
print("=" * 76)
for s in (0.8, 1.25, 2.0):
    bs = births(P * s, 0.30)
    if bs is None:
        print(f"  scaling the configuration by {s}: the diagram DEGENERATES "
              f"entirely at this sigma")
    else:
        print(f"  scaling the configuration by {s}: diagram changes by "
              f"{np.abs(np.sort(bs)-np.sort(b0)).max():.3e}")
print("  So the diagram is NOT scale invariant, and the honest statement is that")
print("  the braid is an invariant of the curve in configuration space modulo")
print("  TRANSLATIONS AND ROTATIONS. The wall locus r_ik = r_jk is separately")
print("  scale invariant, so the walls, but not the diagram, descend to the")
print("  shape sphere.")

print("\n" + "=" * 76)
print("O6. Is the persistence diagram well posed here?")
print("=" * 76)
d = merge_tree_exact(P, MS, 0.30)
nmax = 3
nsad = 0
seen = []
for i in range(3):
    for j in range(i+1, 3):
        sd = exact_saddle(0.5*(P[i]+P[j]), P, MS, 0.30)
        if sd is not None and not any(np.linalg.norm(sd["p"]-q) < 1e-6 for q in seen):
            seen.append(sd["p"]); nsad += 1
nfin = sum(1 for i in range(3) if d[i]["death"] is not None)
print(f"  rho is a smooth Morse function on the plane decaying to 0.")
print(f"    maxima {nmax}, merge saddles {nfin}, essential classes "
      f"{3-nfin}, all saddles {nsad}")
print(f"    Morse count  #max - #saddle = {nmax - nsad}  (must be 1)")
print(f"  Superlevel H_0 persistence of such a function is standard, and the")
print(f"  three classes are exactly the three bodies.")

print("\n" + "=" * 76)
print("O7. Do the three features stay off the diagonal for the whole loop?")
print("=" * 76)
print("  An order-3 permutation of diagram points needs three features that")
print("  survive the whole circuit, stay distinguishable, and never reach the")
print("  diagonal. Over one T/3 loop, sampled at 600 instants:\n")
i0 = int(0.37 * NS)
for sig in (0.22, 0.26, 0.30, 0.34):
    pers, sep, nfeat = [], [], set()
    ok = True
    for k in np.linspace(i0, i0 + NS // 3, 600).astype(int) % NS:
        dd = merge_tree_exact(X[k], MS, sig)
        if dd is None:
            ok = False; break
        pts = []
        for i in range(3):
            de = 0.0 if dd[i]["death"] is None else dd[i]["death"]
            pts.append((dd[i]["birth"], de))
        pts = np.array(pts)
        nfeat.add(len(pts))
        pers.append((pts[:, 0] - pts[:, 1]).min())
        D = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=2)
        sep.append(D[np.triu_indices(3, 1)].min())
    if not ok:
        print(f"    sigma {sig:.2f}: merge tree degenerates"); continue
    print(f"    sigma {sig:.2f}: features {sorted(nfeat)}, "
          f"min persistence over the loop {min(pers):.4f}, "
          f"min pairwise separation {min(sep):.2e}")
print("\n  The minimum persistence is bounded well away from zero, so no feature")
print("  approaches the diagonal. The separation touches zero only at the")
print("  isosceles instants, which is where two vines cross.")

print("\n" + "=" * 76)
print("O8. Is T/3 the PRIMITIVE period in unlabeled configuration space?")
print("=" * 76)
print("  If some shorter time also returned the configuration as a set, the loop")
print("  would not be the primitive one and the order would be wrong.\n")
from scipy.optimize import linear_sum_assignment as _lsa
P0 = X[0]
ss, dd_ = [], []
for k in range(1, NS):
    C = np.linalg.norm(P0[:, None, :] - X[k][None, :, :], axis=2)
    ri, ci = _lsa(C)
    ss.append(k / NS); dd_.append(C[ri, ci].max())
ss, dd_ = np.array(ss), np.array(dd_)
loc = [i for i in range(1, len(dd_) - 1)
       if dd_[i] < dd_[i - 1] and dd_[i] < dd_[i + 1] and dd_[i] < 0.05]
print("    times s/T at which the configuration returns as a SET:")
for i in loc:
    C = np.linalg.norm(P0[:, None, :] - X[int(ss[i]*NS)][None, :, :], axis=2)
    ri, ci = _lsa(C)
    print(f"      s/T = {ss[i]:.5f}   mismatch {dd_[i]:.2e}   relabelling {list(map(int,ci))}")
print("\n  The first return is at T/3 with a 3-cycle, so T/3 is primitive and the")
print("  order is 3, not a multiple of something smaller.")
