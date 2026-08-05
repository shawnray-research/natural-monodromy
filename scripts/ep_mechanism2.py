"""
The EP mechanism test, done fairly.

The first version used real mode profiles, which makes the intensity depend on
the loop only through Re(lambda), so the path in function space retraces itself
and the identity answer is an artifact of the model. Near an EP the cavity
eigenfunctions are genuinely complex (the document says so too), so the basis
modes are complexified here and the loop becomes a genuine circuit in function
space:

    |psi|^2 = |phi_1|^2 + |lambda|^2 |phi_2|^2 + 2 Re( lambda conj(phi_1) phi_2 )

which now depends on both the real and imaginary parts of lambda. Several
random complex mode pairs and several loop radii are tried, so the conclusion
does not rest on one configuration.
"""
import sys, os
import numpy as np
from scipy.optimize import linear_sum_assignment
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.kde import superlevel_h0
from mono.core import perm_order, cycle_type

N = 110
x = np.linspace(0, 1, N)
X, Y = np.meshgrid(x, x, indexing="ij")


def modes(rng):
    """Two complex billiard-like modes: superpositions with random phases."""
    def m(k):
        f = np.zeros_like(X, dtype=complex)
        for _ in range(3):
            a, b = rng.integers(1, 4), rng.integers(1, 4)
            ph = rng.uniform(0, 2*np.pi)
            f += np.exp(1j*ph) * np.sin(a*np.pi*X) * np.sin(b*np.pi*Y)
        return f
    return m(0), m(1)


def diagram(F, floor=2e-3):
    flat = F.ravel()
    pr = superlevel_h0(F)
    p = np.array([[flat[b], flat[d]] for (b, d, _) in pr], dtype=float)
    pers = p[:, 0] - p[:, 1]
    return p[pers >= floor * pers.max()]


def run(seed, r, T=240):
    rng = np.random.default_rng(seed)
    P1, P2 = modes(rng)
    inten = lambda lam: np.abs(P1 + lam * P2) ** 2

    lam = lambda t: np.sqrt(r) * np.exp(1j * np.pi * t)   # t=1 is ONE encirclement

    f0, f1 = inten(lam(0.0)), inten(lam(1.0))
    d0, d1 = diagram(f0), diagram(f1)
    fld = np.abs(f1 - f0).max() / np.abs(f0).max()
    same = None
    if len(d0) == len(d1):
        same = np.linalg.norm(np.sort(d0, 0) - np.sort(d1, 0))

    diags = [diagram(inten(lam(2.0*k/T))) for k in range(T)]
    counts = sorted({len(d) for d in diags})
    order, cyc = None, None
    if len(counts) == 1:
        n = counts[0]; cur = list(range(n))
        for t in range(1, T+1):
            a, b = diags[t-1], diags[t % T]
            C = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
            ri, ci = linear_sum_assignment(C)
            mp = {int(i): int(j) for i, j in zip(ri, ci)}
            cur = [mp[c] for c in cur]
        order, cyc = perm_order(cur), cycle_type(cur)
    return fld, len(d0), len(d1), same, counts, order, cyc


print(__doc__)
print(f"{'seed':>4} {'r':>6} | one loop: field moves  |D_end-D_start| | "
      f"double loop: counts        order")
print("-" * 96)
res = []
for seed in range(4):
    for r in (0.36,):
        fld, n0, n1, same, counts, order, cyc = run(seed, r)
        s = f"{same:.2e}" if same is not None else f"n {n0}->{n1}"
        print(f"{seed:>4} {r:>6.2f} | {fld:>10.3f}  {s:>22s} | "
              f"{str(counts):>16s}  {str(order):>6s}  {cyc if cyc else ''}")
        res.append(order)
print("-" * 96)
ok = [o for o in res if o is not None]
print(f"double-loop monodromy orders observed: {sorted(set(ok))}   "
      f"(cardinality varied in {len(res)-len(ok)} of {len(res)} runs)")
