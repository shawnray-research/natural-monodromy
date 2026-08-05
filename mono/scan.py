"""
Plane-wide monodromy scanner.

For a closed planar curve M, scan the plane of observation points; at each
grid vertex run a small loop and record the permutation the vineyard induces
on the diagram points.  The result is a `monodromy map' of the plane for M:
the (isolated) points where a nontrivial permutation appears are exactly the
monodromy-critical A_1^2/A_1^2 singularities of the generalized symmetry set.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import linear_sum_assignment

from .core import radial_diagram, circle_loop, perm_order, cycle_type


def _persistent_points(d, tau):
    """Diagram points whose persistence exceeds tau * (largest persistence)."""
    pts = d["points"]
    if len(pts) == 0:
        return pts, np.array([], dtype=int)
    pers = pts[:, 1] - pts[:, 0]
    thr = tau * pers.max()
    keep = np.where(pers >= thr)[0]
    return pts[keep], keep


def monodromy_on_loop(M, loop, tau=0.05, return_stack=False):
    """
    Robust monodromy permutation over a closed loop.

    Points with persistence below `tau` times the maximum persistence in their
    own slice are discarded: by Observation 3.3 of arXiv:2607.01046 the vines
    created/annihilated at a focal-set crossing have persistence far below all
    others and are unlinked from the rest, so they cannot carry monodromy.

    Returns dict(perm, order, cycles, ok, reason, n).
    """
    stack = [radial_diagram(M, p) for p in loop]
    if any(s is None for s in stack):
        return {"ok": False, "reason": "degenerate slice", "perm": None,
                "order": None, "cycles": None, "n": 0}

    filt = [_persistent_points(s, tau) for s in stack]
    counts = {len(f[0]) for f in filt}
    if len(counts) != 1:
        return {"ok": False, "reason": f"cardinality varies {sorted(counts)}",
                "perm": None, "order": None, "cycles": None, "n": 0}
    n = counts.pop()
    if n < 2:
        return {"ok": True, "reason": "fewer than 2 points", "perm": list(range(n)),
                "order": 1, "cycles": [], "n": n}

    cur = list(range(n))
    T = len(stack)
    for t in range(1, T + 1):
        a = filt[t - 1][0]
        b = filt[t % T][0]
        C = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
        ri, ci = linear_sum_assignment(C)
        mp = {int(x): int(y) for x, y in zip(ri, ci)}
        cur = [mp[c] for c in cur]

    out = {"ok": True, "reason": "", "perm": cur, "order": perm_order(cur),
           "cycles": cycle_type(cur), "n": n}
    if return_stack:
        out["stack"] = stack
        out["filt"] = filt
    return out


def scan_plane(M, xlim, ylim, nx, ny, loop_steps=64, loop_frac=0.75, tau=0.05,
               progress=False):
    """
    Scan a rectangle of observation points.  At each grid vertex, run a circular
    loop of radius loop_frac * (grid spacing) and record the monodromy order.

    Returns dict with X, Y (grid), order (ny x nx int array, 0 = undetermined),
    npts, and a list of `hits' (x, y, order, cycles).
    """
    xs = np.linspace(*xlim, nx)
    ys = np.linspace(*ylim, ny)
    hx = (xlim[1] - xlim[0]) / max(nx - 1, 1)
    hy = (ylim[1] - ylim[0]) / max(ny - 1, 1)
    rad = loop_frac * min(hx, hy)

    order = np.zeros((ny, nx), dtype=int)
    npts = np.zeros((ny, nx), dtype=int)
    hits = []
    base = circle_loop([0.0, 0.0], rad, loop_steps)
    for iy, y in enumerate(ys):
        for ix, x in enumerate(xs):
            loop = base + np.array([x, y])
            r = monodromy_on_loop(M, loop, tau=tau)
            npts[iy, ix] = r["n"]
            if r["ok"]:
                order[iy, ix] = r["order"]
                if r["order"] and r["order"] > 1:
                    hits.append({"x": float(x), "y": float(y),
                                 "order": r["order"], "cycles": r["cycles"],
                                 "n": r["n"], "radius": float(rad)})
        if progress:
            print(f"    row {iy+1}/{ny}  hits so far: {len(hits)}", flush=True)

    return {"xs": xs, "ys": ys, "order": order, "npts": npts, "hits": hits,
            "radius": rad}


def refine_hit(M, x, y, radii=(0.5, 0.25, 0.12, 0.06), base_radius=1.0,
               loop_steps=256, tau=0.05):
    """Re-test a hit at several shrinking loop radii; genuine local monodromy
    must persist as the loop shrinks onto the singularity."""
    out = []
    for fr in radii:
        loop = circle_loop([x, y], fr * base_radius, loop_steps)
        r = monodromy_on_loop(M, loop, tau=tau)
        out.append({"radius": fr * base_radius, "order": r["order"],
                    "ok": r["ok"], "reason": r["reason"], "n": r["n"],
                    "cycles": r["cycles"]})
    return out
