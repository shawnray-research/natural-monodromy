"""
Independent re-derivation of every number that will be shown to Edelsbrunner.

Nothing here is recalled from notes. The orbit is integrated, the density field
is built, the merge tree is computed analytically (no grid), the braid is read
off the measured crossings, and the Alexander polynomial is computed from the
reduced Burau representation. Discipline adopted after the July 28 retraction:
every figure quoted must come from a run that is saved.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.nbody import FIGURE8, integrate
from mono.kde_exact import merge_tree_exact
from mono.braid import braid_word, word_to_string, permutation_of_word, reduce_word
from mono.core import perm_order, cycle_type

def burau_alexander(word):
    """Alexander polynomial of the closure of a 3-braid, via reduced Burau."""
    import sympy as sp
    t = sp.symbols('t')
    def burau(i, inv):
        M = sp.Matrix([[-t, 0], [1, 1]]) if i == 1 else sp.Matrix([[1, t], [0, -t]])
        return M.inv() if inv else M
    B = sp.eye(2)
    for w in word:
        B = B * burau(abs(w), w < 0)
    num = sp.expand(sp.det(sp.simplify(B) - sp.eye(2)))
    alex = sp.cancel(sp.simplify(num / sp.expand(1 + t + t**2)))
    norm = sp.expand(sp.simplify(-t * alex))
    coeffs = [norm.coeff(t, k) for k in (1, 0, -1)]
    return norm, coeffs


def main():
    F8 = FIGURE8
    m, T, N, NS = F8["m"], F8["T"], 3, 12000
    print("FIGURE-EIGHT THREE-BODY CHOREOGRAPHY")
    print("=" * 70)
    t, X, V = integrate(F8["x"], F8["v"], m, T, NS)
    print(f"  integrator DOP853, rtol=atol=1e-13, {NS} steps")
    print(f"  orbit closure |x(T)-x(0)|            = "
          f"{np.abs(X[-1]-X[0]).max():.3e}")
    # choreography relabelling after T/3
    k3 = NS // 3
    A = X[0]; Bc = X[k3]
    cost = np.linalg.norm(A[:,None,:]-Bc[None,:,:],axis=2)
    from scipy.optimize import linear_sum_assignment
    ri, ci = linear_sum_assignment(cost)
    print(f"  cloud identical after T/3, max mismatch = "
          f"{cost[ri,ci].max():.3e},  relabelling {list(ci)}")

    for SIG in (0.20, 0.25, 0.30):
        i0 = int(0.37*NS); i1 = i0 + NS//N
        K = 3000
        idx = np.linspace(i0, i1, K).astype(int)
        B = np.zeros((K,N)); D = np.zeros((K,N))
        for a,kk in enumerate(idx):
            d = merge_tree_exact(X[kk], m, SIG)
            for i in range(N):
                B[a,i] = d[i]["birth"]
                D[a,i] = 0.0 if d[i]["death"] is None else d[i]["death"]
        w,_,cr = braid_word(B,D)
        wt,_,crt = braid_word(X[idx][:,:,0], X[idx][:,:,1])
        pv = permutation_of_word(w,N); pt = permutation_of_word(wt,N)
        print(f"\n  sigma = {SIG}")
        print(f"    vineyard braid   : {word_to_string(reduce_word(w))}"
              f"   ({len(w)} crossings)  perm {pv} order {perm_order(pv)}")
        print(f"    trajectory braid : {word_to_string(reduce_word(wt))}"
              f"   ({len(wt)} crossings)  perm {pt} order {perm_order(pt)}")
        if SIG == 0.30:
            keep = reduce_word(w)
    norm, coeffs = burau_alexander(keep)
    print(f"\n  Alexander polynomial of the closure : {norm}")
    print(f"  coefficients [t, 1, 1/t]           : {coeffs}")
    print(f"  figure-eight knot 4_1 has t - 3 + 1/t, coefficients [1, -3, 1]: "
          f"{'MATCH' if coeffs == [1, -3, 1] else 'NO MATCH'}")
    print(f"\n  exponent sums (conjugacy invariants, labeling independent):")
    print(f"    vineyard   {sum(1 if x>0 else -1 for x in keep)}")

if __name__ == "__main__":
    main()
