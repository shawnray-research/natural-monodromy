"""
Localization and certification of monodromy-generating singularities.

Three *independent* certificates are computed for a candidate point p*:

  (C1) geometric matching certificate -- Hungarian transport of diagram points
       around a small loop about p* yields a nontrivial permutation;

  (C2) combinatorial certificate, the elder-rule pairing of critical points
       changes exactly twice around the loop (Definition 3.5 of
       arXiv:2607.01046: a `monodromy critical' A_1^2/A_1^2 has two
       consecutive changes of pairing; 0 or 4 changes give the identity);

  (C3) singularity certificate: at p* the squared distance to M has two
       coincident local minima (birth-birth) *and* two coincident local
       maxima (death-death), i.e. p* is the common center of two bitangent
       circles: an A_1^2/A_1^2 point of the symmetry set whose smaller circle
       is bitangent at minima and larger circle bitangent at maxima.

Agreement of C1, C2, C3 makes a detection essentially unimpeachable.
"""

from __future__ import annotations

import numpy as np

from .core import radial_diagram, perm_order, cycle_type
from .scan import monodromy_on_loop


def rect_loop(x0, x1, y0, y1, steps_per_side=64):
    """Counter-clockwise boundary of a rectangle, as a closed loop of points."""
    s = np.linspace(0, 1, steps_per_side, endpoint=False)
    bot = np.column_stack([x0 + (x1 - x0) * s, np.full_like(s, y0)])
    rgt = np.column_stack([np.full_like(s, x1), y0 + (y1 - y0) * s])
    top = np.column_stack([x1 - (x1 - x0) * s, np.full_like(s, y1)])
    lft = np.column_stack([np.full_like(s, x0), y1 - (y1 - y0) * s])
    return np.vstack([bot, rgt, top, lft])


def localize(M, x0, x1, y0, y1, depth=28, steps=64, tau=0.05):
    """
    Quadrisection search for a monodromy-generating singularity inside the
    given rectangle.  Returns (cx, cy, half_width, order) or None.
    """
    r = monodromy_on_loop(M, rect_loop(x0, x1, y0, y1, steps), tau=tau)
    if not r["ok"] or not r["order"] or r["order"] == 1:
        return None
    for _ in range(depth):
        xm = 0.5 * (x0 + x1)
        ym = 0.5 * (y0 + y1)
        quads = [(x0, xm, y0, ym), (xm, x1, y0, ym),
                 (x0, xm, ym, y1), (xm, x1, ym, y1)]
        nxt = None
        for q in quads:
            rr = monodromy_on_loop(M, rect_loop(*q, steps), tau=tau)
            if rr["ok"] and rr["order"] and rr["order"] > 1:
                nxt = (q, rr)
                break
        if nxt is None:
            break
        (x0, x1, y0, y1), r = nxt
    return (0.5 * (x0 + x1), 0.5 * (y0 + y1), 0.5 * max(x1 - x0, y1 - y0),
            r["order"])


def pairing_changes(M, loop):
    """
    Certificate C2: number of elder-rule pairing changes around a loop, and
    the sequence of pairings.  A monodromy-critical A_1^2/A_1^2 shows exactly
    two changes; a non-critical one shows zero or four.
    """
    seq = []
    for p in loop:
        d = radial_diagram(M, p)
        if d is None:
            return None, None
        seq.append(tuple(sorted((i, j) for (i, j, _) in d["pairs"])))
    changes = sum(1 for k in range(len(seq)) if seq[k] != seq[k - 1])
    uniq = []
    for s in seq:
        if not uniq or s != uniq[-1]:
            uniq.append(s)
    return changes, uniq


def singularity_certificate(M, p, n_report=6):
    """
    Certificate C3: report the critical values of d(.,p)^2|_M and how close the
    two smallest minima and the two largest maxima are to coinciding.

    At a true A_1^2/A_1^2 monodromy generator the *pair of minima* involved
    and the *pair of maxima* involved each become equal simultaneously.
    Returns a dict with the critical values and the two relevant gaps,
    normalized by the diameter of the critical-value range.
    """
    d = radial_diagram(M, np.asarray(p, dtype=float))
    if d is None:
        return None
    f = d["f"]
    mn = np.sort(f[d["mins"]])
    mx = np.sort(f[d["maxs"]])[::-1]
    rng = float(f.max() - f.min())
    out = {
        "n_min": len(mn), "n_max": len(mx),
        "min_values": mn[:n_report].tolist(),
        "max_values": mx[:n_report].tolist(),
        "range": rng,
    }
    # all pairwise gaps between minima and between maxima, normalized
    gm = np.abs(mn[:, None] - mn[None, :]) + np.eye(len(mn)) * 1e18
    gM = np.abs(mx[:, None] - mx[None, :]) + np.eye(len(mx)) * 1e18
    out["min_gap_rel"] = float(gm.min() / rng) if len(mn) > 1 else None
    out["max_gap_rel"] = float(gM.min() / rng) if len(mx) > 1 else None
    # which pair of minima / maxima are closest (radii of the two circles)
    if len(mn) > 1:
        i, j = np.unravel_index(np.argmin(gm), gm.shape)
        out["inner_radius"] = float(np.sqrt(0.5 * (mn[i] + mn[j])))
    if len(mx) > 1:
        i, j = np.unravel_index(np.argmin(gM), gM.shape)
        out["outer_radius"] = float(np.sqrt(0.5 * (mx[i] + mx[j])))
    return out


def full_certificate(M, cx, cy, radius, steps=512, tau=0.05):
    """Run C1, C2, C3 at a candidate point and return a combined report."""
    from .core import circle_loop
    loop = circle_loop([cx, cy], radius, steps)
    c1 = monodromy_on_loop(M, loop, tau=tau)
    n_changes, seq = pairing_changes(M, loop)
    c3 = singularity_certificate(M, [cx, cy])
    return {
        "point": (float(cx), float(cy)),
        "loop_radius": float(radius),
        "C1_perm": c1["perm"], "C1_order": c1["order"],
        "C1_cycles": c1["cycles"], "C1_ok": c1["ok"], "C1_reason": c1["reason"],
        "C1_n_diagram_points": c1["n"],
        "C2_pairing_changes": n_changes,
        "C2_pairing_sequence": seq,
        "C3": c3,
    }
