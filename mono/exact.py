"""
Exact, purely combinatorial vine tracking. No geometric matching anywhere.

The earlier tracker matched diagram points between consecutive slices with
`scipy.optimize.linear_sum_assignment` on their (birth, death) coordinates.
That is a heuristic: it can mis-associate when slices are coarse, and it makes
the answer depend on the number of slices.

Here a vine is followed by identity. Persistence pairs are recorded as pairs of
*critical points of the curve*, and critical points move continuously with the
parameter, so they are carried from slice to slice by nearest position along the
curve, cyclically. Between consecutive slices the elder-rule pairing can change
in only two ways: at a birth-birth wall two births exchange order while the
deaths are untouched, so each vine continues as the pair carrying its own DEATH;
at a death-death wall the mirror statement holds. Which occurred is decided by
asking which side's labels were preserved, not by measuring anything.

Two bookkeeping points matter and were got wrong on the first attempt.

* A vine is followed by its POSITION in the slice's pair list, while labels are
  used only to decide the successor. The permutation is then read off directly
  at the end, because the last slice is the first slice. Comparing labels
  instead of positions is wrong precisely when there is monodromy, since a
  critical point that travels once around the curve keeps its label but arrives
  where a different one started.

* Critical points are created and destroyed in min-max pairs whenever the
  parameter crosses the locus where the field stops being Morse. Newly created
  and newly destroyed labels are excluded when deciding the swap kind, and a
  vine whose birth or death critical point is destroyed reaches the diagonal and
  ends there. By Definition 2.6 of arXiv:2607.01046 such a vine forms its own
  closed loop and cannot permute with the vines that stay off the diagonal, so
  the permutation is reported on the survivors with the number that died beside
  it.
"""

from __future__ import annotations

import numpy as np

from .core import perm_order, cycle_type


def match_cyclic(a_idx, b_idx, n, tol):
    """Mutual nearest-neighbor matching of two sets of curve indices, cyclically."""
    if len(a_idx) == 0 or len(b_idx) == 0:
        return {}, 0.0
    D = np.abs(np.asarray(a_idx)[:, None] - np.asarray(b_idx)[None, :])
    D = np.minimum(D, n - D)
    ja = np.argmin(D, axis=1)
    ib = np.argmin(D, axis=0)
    out, worst = {}, 0.0
    for i, j in enumerate(ja):
        if ib[j] == i and D[i, j] <= tol:
            out[i] = int(j)
            worst = max(worst, float(D[i, j]))
    return out, worst


def _swap_kind(old_pairs, new_pairs):
    """
    Which side re-paired between two slices, using only labels common to both so
    that creation and destruction of critical points does not confuse the test.
    """
    if set(old_pairs) == set(new_pairs):
        return "none"
    ob = {p[0] for p in old_pairs} & {p[0] for p in new_pairs}
    od = {p[1] for p in old_pairs} & {p[1] for p in new_pairs}
    o = {p for p in old_pairs if p[0] in ob and p[1] in od}
    nw = {p for p in new_pairs if p[0] in ob and p[1] in od}
    if o == nw:
        return "none"
    changed_old = o - nw
    changed_new = nw - o
    births_preserved = sorted(p[0] for p in changed_old) == sorted(p[0] for p in changed_new)
    deaths_preserved = sorted(p[1] for p in changed_old) == sorted(p[1] for p in changed_new)
    if deaths_preserved and not births_preserved:
        return "births"
    if births_preserved and not deaths_preserved:
        return "deaths"
    if deaths_preserved and births_preserved:
        return "either"
    return "both"


def track_exact(slices, n_curve, tol_frac=0.05):
    """
    slices : list of dicts with 'mins', 'maxs' (curve indices) and 'pairs'
             (list of (i_min_pos, j_max_pos, kind)); slice[0] is revisited at
             the end of the circuit.

    Returns the permutation of the vines that survive one full circuit.
    """
    T = len(slices)
    tol = max(4.0, tol_frac * n_curve)
    nextlab = [0]

    def fresh(k):
        out = list(range(nextlab[0], nextlab[0] + k))
        nextlab[0] += k
        return out

    lab_min = fresh(len(slices[0]["mins"]))
    lab_max = fresh(len(slices[0]["maxs"]))

    def lpairs(sl, lm, lM):
        return [(lm[i], lM[j]) for (i, j, _) in sl["pairs"]]

    k0 = len(slices[0]["pairs"])
    pos = {v: v for v in range(k0)}                 # vine -> index in current slice
    lab = {v: lpairs(slices[0], lab_min, lab_max)[v] for v in range(k0)}
    alive = set(range(k0))
    worst_move = 0.0
    ambiguous = 0
    both_steps = 0

    for t in range(1, T + 1):
        a, b = slices[t - 1], slices[t % T]
        mm, w1 = match_cyclic(a["mins"], b["mins"], n_curve, tol)
        mM, w2 = match_cyclic(a["maxs"], b["maxs"], n_curve, tol)
        worst_move = max(worst_move, w1, w2)

        nlm = [None] * len(b["mins"])
        nlM = [None] * len(b["maxs"])
        for i, j in mm.items():
            nlm[j] = lab_min[i]
        for i, j in mM.items():
            nlM[j] = lab_max[i]
        for j in range(len(nlm)):
            if nlm[j] is None:
                nlm[j] = fresh(1)[0]
        for j in range(len(nlM)):
            if nlM[j] is None:
                nlM[j] = fresh(1)[0]

        old_lp = lpairs(a, lab_min, lab_max)
        new_lp = lpairs(b, nlm, nlM)
        kind = _swap_kind(old_lp, new_lp)
        if kind == "both":
            both_steps += 1
        smin, smax = set(nlm), set(nlM)

        newalive = set()
        for v in list(alive):
            m, M = lab[v]
            if m not in smin or M not in smax:
                continue                              # reached the diagonal
            if (m, M) in new_lp:
                pos[v] = new_lp.index((m, M)); lab[v] = (m, M); newalive.add(v); continue
            if kind == "births":
                c = [q for q in range(len(new_lp)) if new_lp[q][1] == M]
            elif kind == "deaths":
                c = [q for q in range(len(new_lp)) if new_lp[q][0] == m]
            else:
                c = [q for q in range(len(new_lp))
                     if new_lp[q][0] == m or new_lp[q][1] == M]
            if len(c) == 1:
                pos[v] = c[0]; lab[v] = new_lp[c[0]]; newalive.add(v)
            elif len(c) > 1:
                ambiguous += 1
                pos[v] = c[0]; lab[v] = new_lp[c[0]]; newalive.add(v)
        alive = newalive
        lab_min, lab_max = nlm, nlM

    # sigma: the one-period successor map on slice-0 diagram positions.
    # It is a partial injection; `None` means the vine reached the diagonal.
    sigma = {v: (pos[v] if v in alive else None) for v in range(k0)}

    # The permuting core is the largest subset on which sigma is a bijection and
    # no member ever reaches the diagonal. A vine that survives one circuit but
    # lands where a dying vine started will itself die on the next circuit, so it
    # is not part of a closed loop either; peeling those away repeatedly is what
    # this fixed point computes.
    core = {v for v in range(k0) if sigma[v] is not None}
    while True:
        drop = {v for v in core if sigma[v] not in core}
        if not drop:
            break
        core -= drop
    core = sorted(core)

    base = {"n_started": k0, "n_survivors": len(alive), "n_core": len(core),
            "max_move": worst_move, "ambiguous_steps": ambiguous,
            "both_side_steps": both_steps,
            "sigma": {v: sigma[v] for v in range(k0)}}
    if not core:
        return {"ok": True, "perm": [], "order": 1, "cycles": [], "n": 0,
                "note": "no vine survives repeated circuits off the diagonal", **base}
    idx = {v: i for i, v in enumerate(core)}
    perm = [idx[sigma[v]] for v in core]
    return {"ok": True, "perm": perm, "order": perm_order(perm),
            "cycles": cycle_type(perm), "n": len(core), **base}
