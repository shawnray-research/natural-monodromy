"""
Analytically refined critical points of a Gaussian mixture.

Reading birth and death values off a grid is not good enough here. The three
peak heights of the figure-eight density span about 0.09 while a grid samples
them with error that does not shrink fast enough, and the resulting braid word
was unstable: 110, 134, 152, 90 crossings at resolutions 161, 241, 401, 601.
A braid word that changes with the mesh is not a braid word.

For a Gaussian mixture the gradient and Hessian are analytic, so every critical
point can be refined by Newton to machine precision. The grid is then used only
for the COMBINATORICS, which pairs with which, a discrete and robust question,
while every VALUE that enters the braid is computed exactly.
"""

from __future__ import annotations

import numpy as np


def rho_and_derivs(p, pts, masses, sigma):
    """Value, gradient and Hessian of the Gaussian mixture at a point p."""
    d = p[None, :] - pts
    r2 = (d ** 2).sum(axis=1)
    w = masses * np.exp(-r2 / (2 * sigma ** 2))
    val = w.sum()
    g = -(w[:, None] * d).sum(axis=0) / sigma ** 2
    H = np.zeros((2, 2))
    for k in range(len(pts)):
        outer = np.outer(d[k], d[k]) / sigma ** 2 - np.eye(2)
        H += w[k] * outer / sigma ** 2
    return val, g, H


def newton_critical(p0, pts, masses, sigma, iters=60, tol=1e-14):
    """Newton on grad rho = 0. Converges to a critical point of any index."""
    p = np.array(p0, dtype=float)
    for _ in range(iters):
        val, g, H = rho_and_derivs(p, pts, masses, sigma)
        if np.linalg.norm(g) < tol:
            break
        try:
            step = np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            return None
        p = p - step
        if not np.all(np.isfinite(p)):
            return None
    val, g, H = rho_and_derivs(p, pts, masses, sigma)
    ev = np.linalg.eigvalsh(H)
    idx = int((ev > 0).sum())          # 0 = max, 1 = saddle, 2 = min
    return {"p": p, "val": float(val), "grad": float(np.linalg.norm(g)), "index": idx}


def exact_maxima(pts, masses, sigma):
    """One maximum per body, refined from the body position."""
    out = {}
    for i, x in enumerate(pts):
        c = newton_critical(x, pts, masses, sigma)
        if c is None or c["index"] != 0:
            return None
        out[i] = c
    # they must be distinct
    P = np.array([out[i]["p"] for i in range(len(pts))])
    for i in range(len(P)):
        for j in range(i+1, len(P)):
            if np.linalg.norm(P[i]-P[j]) < 1e-8:
                return None
    return out


def exact_saddle(seed, pts, masses, sigma):
    c = newton_critical(seed, pts, masses, sigma)
    if c is None or c["index"] != 1:
        return None
    return c


def diagram_exact(pts, masses, sigma, grid_diagram):
    """
    Combine the grid's combinatorics with analytic values.

    grid_diagram : output of mono.kde.diagram_with_bodies, used only to learn
                   which body pairs with which saddle, and roughly where that
                   saddle is.
    """
    mx = exact_maxima(pts, masses, sigma)
    if mx is None:
        return None
    out = {}
    for e in grid_diagram:
        b = e["body"]
        if b in out:
            continue
        if e["kind"] == "ext":
            out[b] = {"birth": mx[b]["val"], "death": None, "kind": "ext"}
            continue
        s = exact_saddle(np.array([e["dx"], e["dy"]]), pts, masses, sigma)
        if s is None:
            return None
        out[b] = {"birth": mx[b]["val"], "death": s["val"], "kind": "ord",
                  "saddle": s["p"]}
    return out


def merge_tree_exact(pts, masses, sigma):
    """
    Persistence pairing of the Gaussian mixture computed entirely analytically,
    with no grid anywhere.

    The grid was previously used to decide the elder rule, and that is exactly
    where it broke: two peak heights are nearly equal, so which peak counts as
    the oldest flipped with mesh resolution, which flipped which body carried
    the essential class, which flipped the braid crossing signs. The saddle
    values themselves were identical across meshes, so the fault was the
    pairing, not the refinement.

    Maxima are found by Newton from each body; saddles by Newton from each pair
    midpoint, keeping those of Morse index 1. The merge tree is then built by
    sweeping saddles downward with union-find and the elder rule applied to the
    exact peak heights.
    """
    n = len(pts)
    mx = exact_maxima(pts, masses, sigma)
    if mx is None:
        return None

    saddles = []
    for i in range(n):
        for j in range(i + 1, n):
            s = exact_saddle(0.5 * (pts[i] + pts[j]), pts, masses, sigma)
            if s is not None:
                saddles.append({"val": s["val"], "p": s["p"], "pair": (i, j)})
    if not saddles:
        return None
    saddles.sort(key=lambda s: -s["val"])

    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    birth_of = {i: i for i in range(n)}       # root -> body owning the oldest peak
    out = {}
    for s in saddles:
        i, j = s["pair"]
        ri, rj = find(i), find(j)
        if ri == rj:
            continue
        bi, bj = birth_of[ri], birth_of[rj]
        # elder rule on EXACT peak heights
        if mx[bi]["val"] >= mx[bj]["val"]:
            elder, young = ri, rj
        else:
            elder, young = rj, ri
        out[birth_of[young]] = {"birth": mx[birth_of[young]]["val"],
                                "death": s["val"], "kind": "ord",
                                "saddle": s["p"]}
        parent[young] = elder
        birth_of[find(elder)] = birth_of[elder] if \
            mx[birth_of[elder]]["val"] >= mx[birth_of[young]]["val"] else birth_of[young]
    root = find(0)
    ess = birth_of[root]
    out[ess] = {"birth": mx[ess]["val"], "death": None, "kind": "ext"}
    if len(out) != n:
        return None
    return out
