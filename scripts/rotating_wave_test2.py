"""
Symmetry versus transport, computed with the exact analytic merge tree.

Tracking is by OBJECT LABEL (each bump/body carries its own vine), so the
permutation is not something a matcher guessed. The question is then whether
that permutation is observable in the DIAGRAM, which is what vineyard monodromy
is about. It is observable only if the permuted features have distinct diagram
points.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.kde_exact import merge_tree_exact
from mono.nbody import FIGURE8, integrate
from mono.core import perm_order


def pts_of(X, m, sig):
    d = merge_tree_exact(np.asarray(X, dtype=float), np.asarray(m, dtype=float), sig)
    out = []
    for i in range(len(X)):
        b = d[i]["birth"]
        de = d[i]["death"]
        out.append((b, np.nan if de is None else de))
    return np.array(out, dtype=float)


def min_sep(P):
    """Closest pair among the FINITE diagram points."""
    Q = P[np.isfinite(P).all(axis=1)]
    if len(Q) < 2:
        return np.inf
    D = np.linalg.norm(Q[:, None, :] - Q[None, :, :], axis=2)
    return float(D[np.triu_indices(len(Q), 1)].min())


print(__doc__)
print("=" * 78)
print("SYMMETRY: k identical bumps exchanged by a rotation of the instantaneous field")
print("=" * 78)
for k in (2, 3, 4, 5):
    seps = []; spans = []
    for s in range(60):
        th = 2*np.pi*s/(k*60)
        X = np.array([[np.cos(th + 2*np.pi*i/k), np.sin(th + 2*np.pi*i/k)] for i in range(k)])
        P = pts_of(X, [1.0]*k, 0.42)
        seps.append(min_sep(P))
        fin = P[np.isfinite(P).all(axis=1)]
        spans.append(np.ptp(fin[:, 0]) if len(fin) else 0.0)
    print(f"  k={k}: closest pair of diagram points over the loop = {min(seps):.3e}")
    print(f"        spread of birth values among the k bumps      = {max(spans):.3e}")
    print(f"        -> the rotation is an isometry of the instantaneous field, so the")
    print(f"           k bumps carry EQUAL values. The permutation is real on labels")
    print(f"           and invisible in the diagram.")

print()
print("=" * 78)
print("TRANSPORT: figure-eight choreography, bodies exchanged by a TIME SHIFT")
print("=" * 78)
F8 = FIGURE8
NS = 6000
t, Xtr, V = integrate(F8["x"], F8["v"], F8["m"], F8["T"], NS)
i0 = 0; i1 = NS // 3
idx = np.linspace(i0, i1, 240).astype(int)
seps = []; P0 = None; Pend = None
for a, kk in enumerate(idx):
    P = pts_of(Xtr[kk], F8["m"], 0.30)
    seps.append(min_sep(P))
    if a == 0: P0 = P
    Pend = P
print(f"  closest pair of diagram points over T/3 = {min(seps):.3e}")
fin = P0[np.isfinite(P0).all(axis=1)]
print(f"  spread of birth values among the 3 bodies at t=0 = {np.ptp(P0[:,0]):.3e}")
print(f"  body i's diagram point at t=0   : {np.array2string(P0[:,0], precision=5)}")
print(f"  body i's diagram point at t=T/3 : {np.array2string(Pend[:,0], precision=5)}")
print(f"  -> body i lands on body i+1's starting values, and those values DIFFER,")
print(f"     so the permutation is a permutation of distinct diagram points.")
