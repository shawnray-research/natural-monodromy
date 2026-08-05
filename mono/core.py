"""
Core machinery for detecting vineyard monodromy of the radial persistence
transform of a closed planar curve.

Setting (following Chambers-Fillmore-Stephenson-Wintraecken, "Braiding Vineyards",
SODA 2026, and Chambers-Fillmore-Mukherjee-Roy-Stephenson-Wintraecken,
"The Singular Source of Vineyard Monodromy", arXiv:2607.01046):

    M      a smooth closed curve in R^2, sampled as a cyclic polyline
    gamma  a loop of observation points in R^2
    f_t    = d_E(., gamma(t))^2 restricted to M
    V(M,gamma) = the closed vineyard of extended persistence diagrams of f_t

Monodromy = the permutation of off-diagonal diagram points induced by
transporting them once around gamma.

The implementation is exact-combinatorial: persistence pairs are recorded as
pairs of *critical point indices on the curve*, not merely as (birth, death)
values, so vines can be followed without any geometric matching heuristic.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import linear_sum_assignment


# --------------------------------------------------------------------------
# curves
# --------------------------------------------------------------------------

def resample_closed(curve: np.ndarray, n: int) -> np.ndarray:
    """Arclength-resample a closed polyline (n x 2, first != last) to n points."""
    P = np.asarray(curve, dtype=float)
    Q = np.vstack([P, P[:1]])
    seg = np.linalg.norm(np.diff(Q, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    total = s[-1]
    t = np.linspace(0.0, total, n, endpoint=False)
    x = np.interp(t, s, Q[:, 0])
    y = np.interp(t, s, Q[:, 1])
    return np.column_stack([x, y])


def smooth_closed(curve: np.ndarray, sigma_pts: float) -> np.ndarray:
    """Periodic Gaussian smoothing of a closed polyline (keeps it closed)."""
    if sigma_pts <= 0:
        return curve.copy()
    n = len(curve)
    k = int(np.ceil(4 * sigma_pts))
    x = np.arange(-k, k + 1)
    g = np.exp(-0.5 * (x / sigma_pts) ** 2)
    g /= g.sum()
    out = np.empty_like(curve)
    for d in range(curve.shape[1]):
        ext = np.concatenate([curve[-k:, d], curve[:, d], curve[:k, d]])
        out[:, d] = np.convolve(ext, g, mode="same")[k:k + n]
    return out


# --------------------------------------------------------------------------
# critical points of the squared-distance function on the curve
# --------------------------------------------------------------------------

def critical_points(f: np.ndarray):
    """
    Local minima and maxima of a cyclic sequence f (assumed generic: no ties
    between neighbors).  Returns (mins, maxs) as arrays of indices, in
    increasing index order.  They necessarily alternate cyclically.
    """
    n = len(f)
    prev = np.roll(f, 1)
    nxt = np.roll(f, -1)
    mins = np.where((f < prev) & (f < nxt))[0]
    maxs = np.where((f > prev) & (f > nxt))[0]
    return mins, maxs


class _UF:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        self.p[self.find(a)] = self.find(b)


def extended_persistence_circle(f: np.ndarray):
    """
    Extended persistence of a generic function f on S^1 (given as a cyclic
    array).  Uses the elder rule sweep.

    Returns a dict with
        mins, maxs : index arrays of the critical points on the curve
        pairs      : list of (i_min, i_max, kind) with kind in {'ord','ext'};
                     i_min / i_max are *positions in the mins / maxs arrays*
        n_crit     : number of minima ( == number of maxima )

    Every minimum is matched to exactly one maximum: (m-1) ordinary H0 pairs
    plus 1 extended pair (global min, global max).  This is the complete
    pairing of critical points, i.e. the full above-diagonal part of the
    extended diagram; the below-diagonal part mirrors it for a curve without
    boundary (Remark 3.7 of arXiv:2607.01046).
    """
    mins, maxs = critical_points(f)
    m = len(mins)
    if m == 0 or m != len(maxs):
        return None

    # Establish the cyclic alternation: maxs_sorted[j] lies between
    # mins_sorted[j] and mins_sorted[j+1] (cyclically).
    # mins and maxs are both in increasing index order and alternate.
    # Find the offset so that maxs[j] is the first maximum after mins[j].
    off = int(np.searchsorted(maxs, mins[0]))
    maxs_al = np.roll(maxs, -off)          # maxs_al[j] follows mins[j]

    uf = _UF(m)
    birth_of_root = {j: f[mins[j]] for j in range(m)}
    rep_of_root = {j: j for j in range(m)}   # root -> index of its oldest min

    order = np.argsort(f[maxs_al], kind="stable")
    pairs = []
    for j in order:
        a = uf.find(j)
        b = uf.find((j + 1) % m)
        if a == b:
            # closes the circle: H1 born here; in extended persistence this
            # maximum pairs with the global minimum.
            pairs.append((rep_of_root[a], int(j), "ext"))
            continue
        # elder rule: the younger (larger birth) component dies
        if birth_of_root[a] > birth_of_root[b]:
            young, old = a, b
        else:
            young, old = b, a
        pairs.append((rep_of_root[young], int(j), "ord"))
        uf.union(young, old)
        r = uf.find(old)
        birth_of_root[r] = birth_of_root[old]
        rep_of_root[r] = rep_of_root[old]

    return {"mins": mins, "maxs": maxs_al, "pairs": pairs, "n_crit": m}


def radial_diagram(M: np.ndarray, p: np.ndarray):
    """Extended persistence of d(., p)^2 restricted to the closed curve M."""
    f = ((M - p) ** 2).sum(axis=1)
    res = extended_persistence_circle(f)
    if res is None:
        return None
    res["f"] = f
    res["points"] = np.array(
        [[f[res["mins"][i]], f[res["maxs"][j]]] for (i, j, _) in res["pairs"]]
    )
    res["kinds"] = [k for (_, _, k) in res["pairs"]]
    return res


# --------------------------------------------------------------------------
# combinatorial signature (constant on components of the complement of the
# generalized symmetry set -- Prop. 2.7 / Cor. 2.10 of arXiv:2607.01046)
# --------------------------------------------------------------------------

def signature(M: np.ndarray, p: np.ndarray):
    """
    A hashable combinatorial descriptor of the diagram at observation point p:
    the number of critical points, the rank-order of all critical values, and
    the elder-rule pairing.  Two points of the plane in the same connected
    component of the complement of the generalized symmetry set have equal
    signatures.
    """
    r = radial_diagram(M, p)
    if r is None:
        return None
    m = r["n_crit"]
    vals = np.concatenate([r["f"][r["mins"]], r["f"][r["maxs"]]])
    ranks = tuple(np.argsort(np.argsort(vals)).tolist())
    pr = tuple(sorted((i, j) for (i, j, _) in r["pairs"]))
    return (m, ranks, pr)


# --------------------------------------------------------------------------
# vineyard over a loop, and the monodromy permutation
# --------------------------------------------------------------------------

def circle_loop(center, radius, n_steps):
    t = np.linspace(0.0, 2.0 * np.pi, n_steps, endpoint=False)
    c = np.asarray(center, dtype=float)
    return c + radius * np.column_stack([np.cos(t), np.sin(t)])


def vineyard(M: np.ndarray, loop: np.ndarray):
    """Compute the stack of diagrams over the (closed) loop of observation points."""
    return [radial_diagram(M, p) for p in loop]


def _match(prev_pts: np.ndarray, cur_pts: np.ndarray):
    """Hungarian matching of diagram points between consecutive time slices."""
    if len(prev_pts) == 0 or len(cur_pts) == 0:
        return {}
    C = np.linalg.norm(prev_pts[:, None, :] - cur_pts[None, :, :], axis=2)
    ri, ci = linear_sum_assignment(C)
    return {int(a): int(b) for a, b in zip(ri, ci)}


def monodromy_permutation(M: np.ndarray, loop: np.ndarray, return_stack=False):
    """
    Track every off-diagonal diagram point once around `loop` (which must be a
    closed loop of observation points; loop[0] is *not* repeated at the end).

    Returns (perm, info).  perm[i] = j means the vine starting at diagram
    point i of slice 0 ends at diagram point j of slice 0.  perm is None if
    the number of diagram points is not constant around the loop (i.e. the
    loop crosses the focal set), see `monodromy_permutation_robust`.
    """
    stack = vineyard(M, loop)
    if any(s is None for s in stack):
        return None, {"reason": "degenerate slice"}
    counts = {len(s["pairs"]) for s in stack}
    info = {"counts": sorted(counts), "n_slices": len(stack)}
    if len(counts) != 1:
        info["reason"] = "diagram cardinality changes (focal set crossing)"
        if return_stack:
            info["stack"] = stack
        return None, info

    n = len(stack[0]["pairs"])
    cur = list(range(n))          # cur[i] = index in current slice of vine i
    for t in range(1, len(stack) + 1):
        a = stack[t - 1]["points"]
        b = stack[t % len(stack)]["points"]
        mp = _match(a, b)
        cur = [mp[c] for c in cur]
    perm = cur
    if return_stack:
        info["stack"] = stack
    return perm, info


def perm_order(perm):
    """Order of a permutation given as a list."""
    if perm is None:
        return None
    n = len(perm)
    seen = [False] * n
    from math import gcd
    order = 1
    for i in range(n):
        if seen[i]:
            continue
        L = 0
        j = i
        while not seen[j]:
            seen[j] = True
            j = perm[j]
            L += 1
        order = order * L // gcd(order, L)
    return order


def cycle_type(perm):
    if perm is None:
        return None
    n = len(perm)
    seen = [False] * n
    cycles = []
    for i in range(n):
        if seen[i]:
            continue
        c = []
        j = i
        while not seen[j]:
            seen[j] = True
            c.append(j)
            j = perm[j]
        cycles.append(tuple(c))
    return sorted(cycles, key=len, reverse=True)
