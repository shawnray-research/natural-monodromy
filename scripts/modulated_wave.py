"""
The corrected theorem, tested.

The proposed theorem said: symmetry exchanges the extrema, therefore monodromy.
The computation says the opposite: if an isometry of the INSTANTANEOUS field
exchanges two extrema, they carry equal values and the diagram cannot see the
exchange. Symmetry is what kills it.

What survives is a rotating wave carrying a STATIONARY amplitude modulation:

    f(theta, t) = A(theta) * h(theta - omega t),     h of period 2*pi/n

  - h has n crests, so f(theta, t + T/n) = f(theta, t) EXACTLY: the loop closes,
    and its length T/n is the system's own period, not a choice
  - crest j moves to where crest j+1 was, so the features cyclically permute
  - the envelope A is fixed in the lab frame, so at any instant the n crests have
    DIFFERENT heights: the configuration is not C_n symmetric even though the
    crest positions are equally spaced
  - no crest is created or destroyed, so D7 holds by construction

That is transport rather than symmetry, and it is the same mechanism as the
choreography: the exchanged objects are at different phases of one process.
A(theta) constant recovers the degenerate case as a control.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.core import extended_persistence_circle, perm_order, cycle_type
from scipy.optimize import linear_sum_assignment

M = 3000
TH = np.linspace(0, 2*np.pi, M, endpoint=False)


def field(n, t, mod, phase=0.0):
    """n crests rotating through a fixed envelope 1 + mod*cos(theta)."""
    A = 1.0 + mod * np.cos(TH + phase)
    h = np.cos(n * (TH - 2*np.pi*t/n))
    return A * h


def diagram(f):
    r = extended_persistence_circle(f)
    if r is None:
        return None
    return np.array([[f[r["mins"][i]], f[r["maxs"][j]]] for (i, j, _) in r["pairs"]])


def run(n, mod, steps=400):
    diags = []
    for s in range(steps):
        d = diagram(field(n, s/steps, mod))
        if d is None:
            return None
        diags.append(d)
    counts = {len(d) for d in diags}
    if len(counts) != 1:
        return {"ok": False, "counts": sorted(counts)}
    N = counts.pop()
    seps = []
    for d in diags:
        D = np.linalg.norm(d[:, None, :] - d[None, :, :], axis=2)
        seps.append(D[np.triu_indices(len(d), 1)].min())
    cur = list(range(N))
    for t in range(1, steps+1):
        a, b = diags[t-1], diags[t % steps]
        C = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
        ri, ci = linear_sum_assignment(C)
        mp = {int(i): int(j) for i, j in zip(ri, ci)}
        cur = [mp[c] for c in cur]
    closure = np.abs(field(n, 0.0, mod) - field(n, 1.0, mod)).max()
    return {"ok": True, "n_pts": N, "order": perm_order(cur),
            "cycles": cycle_type(cur), "min_sep": min(seps), "closure": closure}


print(__doc__)
print(f"{'crests':>7} {'modulation':>11} {'closure':>11} {'pts':>4} "
      f"{'min sep':>10} {'order':>6}   cycles")
print("-" * 78)
for n in (3, 4, 5, 6):
    for mod in (0.0, 0.05, 0.20, 0.45):
        r = run(n, mod)
        if r is None or not r.get("ok"):
            print(f"{n:>7} {mod:>11.2f}   D7 fails, counts {r['counts'] if r else '-'}")
            continue
        tag = "  <- degenerate control" if mod == 0.0 else ""
        print(f"{n:>7} {mod:>11.2f} {r['closure']:>11.2e} {r['n_pts']:>4} "
              f"{r['min_sep']:>10.2e} {r['order']:>6}   {r['cycles']}{tag}")
