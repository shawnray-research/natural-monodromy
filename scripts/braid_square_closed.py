"""
Settling: is the vineyard braid the square of the trajectory braid?

The claim came from the figure-eight, where the vineyard word (s2^-1 s1)^2 is the
square of the trajectory word s1 s2^-1 up to conjugacy. It was withdrawn when the
supporting N=5 comparison turned out to compare permutations written in two
different strand labellings. It is closed here, structurally.

PROPOSITION. For any choreography, both braids have the same underlying
permutation over T/N, namely a single N-cycle.

  trajectory: x_i(t + T/N) = x_{i+1}(t) is the definition of a choreography, so
              body i ends where body i+1 began.
  vineyard:   diagram point i belongs to body i. At t = T/N the configuration is
              identical as a set, so the field is identical, and body i sits at
              body i+1's former position carrying body i+1's former birth and
              death. So vine i ends on vine i+1's starting diagram point.

COROLLARY. If the vineyard braid were the square of the trajectory braid up to
conjugacy then, applying the homomorphism B_N -> S_N,

        N-cycle  =  (N-cycle)^2.

The square of an N-cycle is an N-cycle when N is odd, and splits into two cycles
of length N/2 when N is even. So the relation is FALSE for every even N, with no
computation and no numerics.

AND THE N=3 EVIDENCE WAS VACUOUS. At N=3 the square of a 3-cycle is a 3-cycle, so
the permutation invariant cannot distinguish anything; and both exponent sums are
0, so e(V) = 2 e(T) reads 0 = 0. Neither conjugation invariant available at N=3
is capable of detecting a difference. The agreement was never evidence.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.nbody import FIGURE8, integrate
from mono.kde_exact import merge_tree_exact
from mono.braid import braid_word, initial_order
from mono.core import cycle_type

def lens(p):
    return sorted(len(c) for c in cycle_type(p))

def sq(p):
    return [p[p[i]] for i in range(len(p))]

def perm_of(word, k, init):
    order = list(init)
    for x in word:
        j = abs(x) - 1
        order[j], order[j + 1] = order[j + 1], order[j]
    p = [0] * k
    for pos, vine in enumerate(order):
        p[list(init).index(vine)] = pos
    return p

print(__doc__)
print("=" * 72)
print("Step 1. The corollary, as pure group theory.")
print("=" * 72)
for N in range(3, 11):
    c = [(i + 1) % N for i in range(N)]          # an N-cycle
    print(f"  N={N:2d}: N-cycle {lens(c)}   its square {lens(sq(c))}   "
          f"{'compatible' if lens(c)==lens(sq(c)) else 'INCOMPATIBLE -> relation false'}")

print("\n" + "=" * 72)
print("Step 2. Check the proposition on the real orbit (N=3).")
print("=" * 72)
F8 = FIGURE8; ms, T, NS = F8["m"], F8["T"], 12000
_, X, _ = integrate(F8["x"], F8["v"], ms, T, NS)
i0 = int(0.37 * NS); idx = np.linspace(i0, i0 + NS // 3, 2000).astype(int)
B = np.zeros((len(idx), 3)); D = np.zeros((len(idx), 3))
for a, kk in enumerate(idx):
    d = merge_tree_exact(X[kk], ms, 0.30)
    for i in range(3):
        B[a, i] = d[i]["birth"]
        D[a, i] = 0.0 if d[i]["death"] is None else d[i]["death"]
wv, _, _ = braid_word(B, D)
wt, _, _ = braid_word(X[idx][:, :, 0], X[idx][:, :, 1])
pv = perm_of(wv, 3, initial_order(B))
pt = perm_of(wt, 3, initial_order(X[idx][:, :, 0]))
ev = sum(1 if x > 0 else -1 for x in wv)
et = sum(1 if x > 0 else -1 for x in wt)
print(f"  vineyard   permutation cycle lengths {lens(pv)}   exponent sum {ev:+d}")
print(f"  trajectory permutation cycle lengths {lens(pt)}   exponent sum {et:+d}")
print(f"  both are 3-cycles, as the proposition requires: "
      f"{lens(pv) == [3] and lens(pt) == [3]}")
print(f"  square of trajectory permutation: {lens(sq(pt))}")
print(f"  permutation test at N=3: {lens(pv)} vs {lens(sq(pt))} -> "
      f"{'passes, and passes automatically' if lens(pv)==lens(sq(pt)) else 'fails'}")
print(f"  exponent test at N=3: e(V)={ev:+d} vs 2e(T)={2*et:+d} -> "
      f"{'passes, and both sides are zero' if ev==2*et else 'fails'}")

print("\n" + "=" * 72)
print("VERDICT")
print("=" * 72)
print("  The relation is FALSE for every even N, by the corollary.")
print("  At N=3, the only case where it was ever observed, both available")
print("  conjugation invariants are incapable of detecting a difference.")
print("  So the claim is refuted in general and was never supported at N=3.")
print("  No statement is made for odd N >= 5, and none is needed: a relation")
print("  that fails for every even N is not a relation.")
