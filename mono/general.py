"""
Monodromy detection for an ARBITRARY two-parameter family of scalar fields on
the circle.

The A_1^2/A_1^2 theorem of arXiv:2607.01046 is stated for the radial transform,
but its mechanism is not special to distance functions.  Stripped of the
distance-function geometry it says: in a two-parameter family f_(u,v) of Morse
functions, monodromy is generated at an isolated point of the control plane
where a **birth-birth wall** (two minima exchange value) crosses a
**death-death wall** (two maxima exchange value), *and* the elder-rule pairing
changes exactly twice around it.

The birth-birth wall is the analogue of the symmetry-set branch bitangent at
two minima; the death-death wall the analogue of the branch bitangent at two
maxima; their transverse crossing is the analogue of A_1^2/A_1^2.  This module
detects exactly that for any field-generating callable, so the same certified
machinery can be pointed at physical fields (quantum densities, potentials,
measured profiles) rather than only at shapes.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import linear_sum_assignment

from .core import extended_persistence_circle, perm_order, cycle_type


def field_diagram(f):
    """Extended persistence diagram of a generic field on S^1, with pairing."""
    r = extended_persistence_circle(np.asarray(f, dtype=float))
    if r is None:
        return None
    pts = np.array([[f[r["mins"][i]], f[r["maxs"][j]]] for (i, j, _) in r["pairs"]])
    return {"points": pts, "pairs": r["pairs"], "mins": r["mins"],
            "maxs": r["maxs"], "f": np.asarray(f, dtype=float)}


def pairing_key(d):
    return tuple(sorted((i, j) for (i, j, _) in d["pairs"]))


def rank_pairing_key(points):
    """
    The elder-rule pairing of a SET of diagram points, keyed by rank rather than
    by critical-point index.

    pairing_key above identifies a pairing by the indices of the critical points
    on the curve. That is exact when every critical point is kept, which is the
    case for a clean synthetic field, but it is wrong for thresholded measured
    data: sub-threshold noise features appear and vanish around the loop and
    shift the indices of the significant ones, so the pairing appears to change
    at almost every step even when the significant features are pairing up in
    exactly the same way.

    Keying on (rank of the birth among births, rank of the death among deaths)
    is invariant to that relabelling and still detects the two events that
    matter: crossing a birth-birth wall swaps two birth ranks, and crossing a
    death-death wall swaps two death ranks.
    """
    if len(points) == 0:
        return ()
    b = np.asarray(points)[:, 0]
    d = np.asarray(points)[:, 1]
    rb = np.empty(len(b), dtype=int); rb[np.argsort(b, kind="stable")] = np.arange(len(b))
    rd = np.empty(len(d), dtype=int); rd[np.argsort(d, kind="stable")] = np.arange(len(d))
    return tuple(sorted(zip(rb.tolist(), rd.tolist())))


def loop_points(cx, cy, radius, steps):
    t = np.linspace(0, 2 * np.pi, steps, endpoint=False)
    return np.column_stack([cx + radius * np.cos(t), cy + radius * np.sin(t)])


def monodromy_on_param_loop(field_fn, pts, tau=0.05):
    """
    field_fn(u, v) -> 1-D array on S^1.
    pts            -> (T x 2) array of control-plane points forming a closed loop.
    """
    diags = []
    for u, v in pts:
        d = field_diagram(field_fn(u, v))
        if d is None:
            return {"ok": False, "reason": "degenerate field", "order": None,
                    "n": 0, "pairing_changes": None}
        diags.append(d)

    keys = [pairing_key(d) for d in diags]
    changes = sum(1 for k in range(len(keys)) if keys[k] != keys[k - 1])
    uniq = []
    for k in keys:
        if not uniq or k != uniq[-1]:
            uniq.append(k)

    filt = []
    for d in diags:
        p = d["points"]
        if len(p) == 0:
            filt.append(p)
            continue
        pers = p[:, 1] - p[:, 0]
        filt.append(p[pers >= tau * pers.max()])
    counts = {len(f) for f in filt}
    if len(counts) != 1:
        return {"ok": False, "reason": f"cardinality varies {sorted(counts)}",
                "order": None, "n": 0, "pairing_changes": changes,
                "pairing_sequence": uniq}
    n = counts.pop()
    if n < 2:
        return {"ok": True, "reason": "", "order": 1, "perm": list(range(n)),
                "cycles": [], "n": n, "pairing_changes": changes,
                "pairing_sequence": uniq}
    cur = list(range(n))
    T = len(filt)
    for t in range(1, T + 1):
        a, b = filt[t - 1], filt[t % T]
        C = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
        ri, ci = linear_sum_assignment(C)
        mp = {int(x): int(y) for x, y in zip(ri, ci)}
        cur = [mp[c] for c in cur]
    return {"ok": True, "reason": "", "order": perm_order(cur), "perm": cur,
            "cycles": cycle_type(cur), "n": n, "pairing_changes": changes,
            "pairing_sequence": uniq}


def reorganization(field_fn, pts):
    """Normalized maximum L1 excursion of the field around the loop.
    ~0 = field barely moves (no wall crossed); ~2 = field completely replaced
    (features annihilate rather than permute).  Monodromy lives in between."""
    f0 = np.asarray(field_fn(*pts[0]), dtype=float)
    m = 0.0
    for u, v in pts:
        f = np.asarray(field_fn(u, v), dtype=float)
        m = max(m, float(np.abs(f - f0).sum() / np.abs(f0).sum()))
    return m


def scan_param_plane(field_fn, ulim, vlim, nu, nv, loop_steps=32, loop_frac=0.8,
                     tau=0.05, progress=False):
    """Loop-scan a rectangle of the control plane; return grid of monodromy orders."""
    us = np.linspace(*ulim, nu)
    vs = np.linspace(*vlim, nv)
    hu = (ulim[1] - ulim[0]) / max(nu - 1, 1)
    hv = (vlim[1] - vlim[0]) / max(nv - 1, 1)
    rad = loop_frac * min(hu, hv)
    order = np.zeros((nv, nu), dtype=int)
    changes = np.full((nv, nu), -1, dtype=int)
    hits = []
    for j, v in enumerate(vs):
        for i, u in enumerate(us):
            r = monodromy_on_param_loop(field_fn, loop_points(u, v, rad, loop_steps),
                                        tau=tau)
            if r["pairing_changes"] is not None:
                changes[j, i] = r["pairing_changes"]
            if r["ok"] and r["order"]:
                order[j, i] = r["order"]
                if r["order"] > 1:
                    hits.append({"u": float(u), "v": float(v), "order": r["order"],
                                 "pairing_changes": r["pairing_changes"],
                                 "n": r["n"], "radius": float(rad)})
        if progress:
            print(f"    row {j+1}/{nv}, hits {len(hits)}", flush=True)
    return {"us": us, "vs": vs, "order": order, "pairing_changes": changes,
            "hits": hits, "radius": rad}


def localize_param(field_fn, u0, u1, v0, v1, depth=24, steps=32, tau=0.05):
    """Quadrisection localization in the control plane."""
    def rect(a0, a1, b0, b1, s):
        t = np.linspace(0, 1, s, endpoint=False)
        bot = np.column_stack([a0 + (a1 - a0) * t, np.full_like(t, b0)])
        rgt = np.column_stack([np.full_like(t, a1), b0 + (b1 - b0) * t])
        top = np.column_stack([a1 - (a1 - a0) * t, np.full_like(t, b1)])
        lft = np.column_stack([np.full_like(t, a0), b1 - (b1 - b0) * t])
        return np.vstack([bot, rgt, top, lft])

    r = monodromy_on_param_loop(field_fn, rect(u0, u1, v0, v1, steps), tau=tau)
    if not r["ok"] or not r["order"] or r["order"] == 1:
        return None
    for _ in range(depth):
        um, vm = 0.5 * (u0 + u1), 0.5 * (v0 + v1)
        nxt = None
        for q in [(u0, um, v0, vm), (um, u1, v0, vm),
                  (u0, um, vm, v1), (um, u1, vm, v1)]:
            rr = monodromy_on_param_loop(field_fn, rect(*q, steps), tau=tau)
            if rr["ok"] and rr["order"] and rr["order"] > 1:
                nxt = (q, rr)
                break
        if nxt is None:
            break
        (u0, u1, v0, v1), r = nxt
    return (0.5 * (u0 + u1), 0.5 * (v0 + v1),
            0.5 * max(u1 - u0, v1 - v0), r["order"], r["pairing_changes"])


def wall_certificate(field_fn, u, v):
    """
    Analogue of C3 for a general field: at a birth-birth x death-death crossing,
    two minima of the field must coincide in value AND two maxima must coincide.
    """
    d = field_diagram(field_fn(u, v))
    if d is None:
        return None
    f = d["f"]
    mn = np.sort(f[d["mins"]])
    mx = np.sort(f[d["maxs"]])
    rng = float(f.max() - f.min())
    gm = np.min(np.diff(mn)) / rng if len(mn) > 1 else None
    gM = np.min(np.diff(mx)) / rng if len(mx) > 1 else None
    return {"n_min": len(mn), "n_max": len(mx),
            "min_gap_rel": None if gm is None else float(gm),
            "max_gap_rel": None if gM is None else float(gM),
            "min_values": mn.tolist(), "max_values": mx.tolist()}


def monodromy_abs(field_fn, pts, floor):
    """
    As monodromy_on_param_loop, but with an ABSOLUTE persistence floor rather
    than a fraction of the largest persistence in each diagram.

    For measured data the right floor is set by the noise of the instrument, not
    by the brightest feature present: a relative threshold silently changes the
    feature inventory whenever the strongest feature grows or shrinks along the
    loop, which is precisely the failure that D7 is meant to exclude.
    """
    diags = []
    for u, v in pts:
        d = field_diagram(field_fn(u, v))
        if d is None:
            return {"ok": False, "reason": "degenerate field", "order": None,
                    "n": 0, "pairing_changes": None}
        diags.append(d)

    filt = []
    for d in diags:
        p = d["points"]
        if len(p) == 0:
            filt.append(p)
            continue
        filt.append(p[(p[:, 1] - p[:, 0]) >= floor])

    keys = [rank_pairing_key(f) for f in filt]
    changes = sum(1 for k in range(len(keys)) if keys[k] != keys[k - 1])

    counts = {len(f) for f in filt}
    if len(counts) != 1:
        return {"ok": False, "reason": f"cardinality varies {sorted(counts)}",
                "order": None, "n": 0, "pairing_changes": changes}
    n = counts.pop()
    if n < 2:
        return {"ok": True, "reason": "", "order": 1, "perm": list(range(n)),
                "cycles": [], "n": n, "pairing_changes": changes}
    cur = list(range(n))
    T = len(filt)
    for t in range(1, T + 1):
        a, b = filt[t - 1], filt[t % T]
        C = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
        ri, ci = linear_sum_assignment(C)
        mp = {int(x): int(y) for x, y in zip(ri, ci)}
        cur = [mp[c] for c in cur]
    return {"ok": True, "reason": "", "order": perm_order(cur), "perm": cur,
            "cycles": cycle_type(cur), "n": n, "pairing_changes": changes,
            "min_gap": float(min(
                np.linalg.norm(f[:, None, :] - f[None, :, :], axis=2)[
                    np.triu_indices(len(f), 1)].min() for f in filt if len(f) > 1))}


def localize_abs(field_fn, u0, u1, v0, v1, floor, depth=26, steps=64):
    """
    Quadrisection localization of a monodromy generator, with an absolute
    persistence floor.

    Certification must be applied at the generator, not at the scan grid point
    that first noticed it. Shrinking a loop about a center that is offset from
    the singularity correctly returns the identity as soon as the loop no longer
    encloses it, which looks exactly like a spurious detection and is not one.
    """
    def rect(a0, a1, b0, b1, s):
        t = np.linspace(0, 1, s, endpoint=False)
        bot = np.column_stack([a0 + (a1 - a0) * t, np.full_like(t, b0)])
        rgt = np.column_stack([np.full_like(t, a1), b0 + (b1 - b0) * t])
        top = np.column_stack([a1 - (a1 - a0) * t, np.full_like(t, b1)])
        lft = np.column_stack([np.full_like(t, a0), b1 - (b1 - b0) * t])
        return np.vstack([bot, rgt, top, lft])

    r = monodromy_abs(field_fn, rect(u0, u1, v0, v1, steps), floor)
    if not r["ok"] or not r["order"] or r["order"] == 1:
        return None
    order0 = r["order"]
    for _ in range(depth):
        um, vm = 0.5 * (u0 + u1), 0.5 * (v0 + v1)
        nxt = None
        for q in [(u0, um, v0, vm), (um, u1, v0, vm),
                  (u0, um, vm, v1), (um, u1, vm, v1)]:
            rr = monodromy_abs(field_fn, rect(*q, steps), floor)
            if rr["ok"] and rr["order"] and rr["order"] > 1:
                nxt = (q, rr); break
        if nxt is None:
            break
        (u0, u1, v0, v1), r = nxt
    return {"u": 0.5 * (u0 + u1), "v": 0.5 * (v0 + v1),
            "box": 0.5 * max(u1 - u0, v1 - v0),
            "order": r["order"], "order_outer": order0,
            "changes": r["pairing_changes"], "n": r["n"]}


def monodromy_transport(field_fn, pts, floor, tol_frac=0.12):
    """
    Monodromy of a parameter loop, tracked by TRANSPORT rather than by assignment.

    Every measured-data scan in this project used monodromy_abs, which follows
    vines by scipy.optimize.linear_sum_assignment on the (birth, death)
    coordinates. That tracker was then shown to return the identity on rotating
    waves whose true monodromy is an n-cycle, in 11 of 11 cases, with the closest
    pair of diagram points never below 9.6e-03, so the failure is not a
    near-collision artifact: vines exchange in the death coordinate while staying
    separated in birth, and distance-based matching follows the wrong branch.
    Every negative obtained with that tracker is therefore suspect.

    Here a vine is followed by the identity of its MAXIMUM, which is a critical
    point of the field and moves continuously with the parameters, so it is
    carried from step to step by nearest position along the domain circle. This
    is the same principle as mono.exact, applied to a two-parameter loop.

    Returns the permutation of vines, where vine i is the one whose maximum is
    i-th in position order at the first loop point.
    """
    slices = []
    for u, v in pts:
        f = np.asarray(field_fn(u, v), dtype=float)
        r = extended_persistence_circle(f)
        if r is None:
            return {"ok": False, "reason": "degenerate field", "order": None, "n": 0}
        keep = [(i, j) for (i, j, _) in r["pairs"]
                if f[r["maxs"][j]] - f[r["mins"][i]] >= floor]
        if not keep:
            return {"ok": False, "reason": "no features above floor",
                    "order": None, "n": 0}
        mx = np.array([int(r["maxs"][j]) for (_, j) in keep])
        bd = np.array([[f[r["mins"][i]], f[r["maxs"][j]]] for (i, j) in keep])
        o = np.argsort(mx)
        slices.append({"mx": mx[o], "bd": bd[o], "n": len(f)})

    counts = {len(s["mx"]) for s in slices}
    if len(counts) != 1:
        return {"ok": False, "reason": f"cardinality varies {sorted(counts)}",
                "order": None, "n": 0}
    k = counts.pop()
    if k < 2:
        return {"ok": True, "order": 1, "perm": list(range(k)), "cycles": [],
                "n": k, "max_jump": 0.0}

    ncur = slices[0]["n"]
    tol = tol_frac * ncur
    cur = list(range(k))
    worst = 0.0
    T = len(slices)
    for t in range(1, T + 1):
        a, b = slices[t - 1]["mx"], slices[t % T]["mx"]
        D = np.abs(a[:, None] - b[None, :])
        D = np.minimum(D, ncur - D)
        ri, ci = linear_sum_assignment(D)      # on POSITIONS, not on (birth, death)
        worst = max(worst, float(D[ri, ci].max()))
        if D[ri, ci].max() > tol:
            return {"ok": False, "reason": f"transport jump {D[ri,ci].max():.1f} "
                    f"exceeds tol {tol:.1f}", "order": None, "n": k}
        mp = {int(i): int(j) for i, j in zip(ri, ci)}
        cur = [mp[c] for c in cur]
    return {"ok": True, "order": perm_order(cur), "perm": cur,
            "cycles": cycle_type(cur), "n": k, "max_jump": worst}
