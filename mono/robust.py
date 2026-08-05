"""
Correct robustness testing for a located monodromy generator.

Naive test (wrong): re-run `localize` on a big rectangle around the original
point.  This fails for a reason the theory predicts -- Figure 4.1 of
arXiv:2607.01046 shows that a *large* loop enclosing an A_1^2/A_1^2 singularity
that locally generates monodromy may itself have trivial monodromy, because the
permutations contributed by several enclosed singularities can cancel.  A big
rectangle returning "no monodromy" therefore says nothing about whether the
singularity survived.

Correct test: re-scan a small neighborhood with *small* loops, find the nearest
surviving monodromy generator, and report its drift.
"""

from __future__ import annotations

import numpy as np

from .core import circle_loop
from .scan import monodromy_on_loop
from .localize import localize, full_certificate

# C3 is an *identification* certificate, not the primary evidence; its numerical
# gate must be loose enough not to be resolution-sensitive.  C1 (stable nontrivial
# permutation), C2 (exactly two pairing changes) and the causal toggle carry the proof.
GAP_GATE = 1e-3


def find_nearest_monodromy(M, center, search_radius, ngrid=13, tau=0.05,
                           certify=True):
    """
    Scan a square of side 2*search_radius about `center` with small loops and
    return the certified monodromy generator closest to `center`, or None.
    """
    cx, cy = center
    xs = np.linspace(cx - search_radius, cx + search_radius, ngrid)
    ys = np.linspace(cy - search_radius, cy + search_radius, ngrid)
    h = xs[1] - xs[0]
    rad = 0.8 * h
    cands = []
    base = circle_loop([0.0, 0.0], rad, 40)
    for y in ys:
        for x in xs:
            r = monodromy_on_loop(M, base + np.array([x, y]), tau=tau)
            if r["ok"] and r["order"] and r["order"] > 1:
                cands.append((x, y))
    best = None
    for (x, y) in sorted(cands, key=lambda p: (p[0] - cx) ** 2 + (p[1] - cy) ** 2):
        loc = localize(M, x - rad, x + rad, y - rad, y + rad, depth=28, steps=40)
        if loc is None:
            continue
        lx, ly, hw, order = loc
        if not certify:
            best = {"x": lx, "y": ly, "order": order,
                    "drift": float(np.hypot(lx - cx, ly - cy))}
            break
        cert = full_certificate(M, lx, ly, rad / 30, steps=256)
        c3 = cert["C3"]
        if (cert["C1_order"] and cert["C1_order"] > 1
                and cert["C2_pairing_changes"] == 2
                and c3["min_gap_rel"] is not None
                and c3["min_gap_rel"] < GAP_GATE and c3["max_gap_rel"] < GAP_GATE):
            best = {"x": lx, "y": ly, "order": cert["C1_order"],
                    "pairing_changes": cert["C2_pairing_changes"],
                    "min_gap_rel": c3["min_gap_rel"],
                    "max_gap_rel": c3["max_gap_rel"],
                    "drift": float(np.hypot(lx - cx, ly - cy))}
            break
    return best
