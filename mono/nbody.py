"""
N-body choreographies and the vineyards of their point clouds.

A choreography is a periodic solution of the Newtonian N-body problem in which
all bodies traverse one closed curve equally spaced in time,

    x_i(t) = q(t + (i-1) T/N),

so the UNLABELLED configuration has period T/N while the labelled one has period
T. Over the short period the bodies cyclically permute. That permutation is the
defining property of the solution, not something imposed.

For H_0 of a point cloud under the Vietoris-Rips (equivalently, for H_0, the
Cech or the Euclidean MST) filtration, every point is born at radius 0 and the
deaths are the edge lengths of the minimum spanning tree. Diagram points are
therefore distinguished only by their death value, so vines must be followed by
**MST edge identity**, that is, by which pair of bodies the edge joins. Sorting
the death values instead would make crossings invisible by construction, since
a sorted sequence cannot cross itself.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from .core import perm_order, cycle_type


# Chenciner-Montgomery figure eight, equal masses, G = 1.
FIGURE8 = {
    "x": np.array([[0.97000436, -0.24308753],
                   [-0.97000436, 0.24308753],
                   [0.0, 0.0]]),
    "v": np.array([[0.46620369, 0.43236573],
                   [0.46620369, 0.43236573],
                   [-0.93240737, -0.86473146]]),
    "T": 6.32591398292621,
    "m": np.ones(3),
}


def accel(x, m, soft=0.0):
    n = len(x)
    d = x[None, :, :] - x[:, None, :]
    r2 = (d ** 2).sum(axis=2) + soft ** 2
    np.fill_diagonal(r2, np.inf)
    inv = r2 ** -1.5
    return (d * (inv * m[None, :])[:, :, None]).sum(axis=1)


def integrate(x0, v0, m, T, n_steps, rtol=1e-13, atol=1e-13):
    n = len(m)

    def rhs(t, y):
        x = y[: 2 * n].reshape(n, 2)
        v = y[2 * n:].reshape(n, 2)
        return np.concatenate([v.ravel(), accel(x, m).ravel()])

    y0 = np.concatenate([x0.ravel(), v0.ravel()])
    ts = np.linspace(0.0, T, n_steps + 1)
    sol = solve_ivp(rhs, (0.0, T), y0, t_eval=ts, method="DOP853",
                    rtol=rtol, atol=atol, dense_output=False)
    X = sol.y[: 2 * n].T.reshape(-1, n, 2)
    V = sol.y[2 * n:].T.reshape(-1, n, 2)
    return sol.t, X, V


def mst_edges(P):
    """
    Euclidean minimum spanning tree of a small point set, returned as a list of
    (i, j, length) sorted by length. These lengths are exactly the finite deaths
    in H_0 of the Rips filtration, and (i, j) is the identity used to track the
    vine.
    """
    n = len(P)
    D = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)
    inT = [0]
    out = set(range(1, n))
    edges = []
    while out:
        best = None
        for a in inT:
            for b in out:
                if best is None or D[a, b] < best[2]:
                    best = (a, b, D[a, b])
        edges.append((min(best[0], best[1]), max(best[0], best[1]), float(best[2])))
        inT.append(best[1])
        out.discard(best[1])
    return sorted(edges, key=lambda e: e[2])


def vineyard_h0(X):
    """MST edges at every time slice: the H_0 vineyard of the moving point cloud."""
    return [mst_edges(P) for P in X]


def monodromy_by_edge_identity(X, sigma):
    """
    Track each MST edge by the pair of bodies it joins, around the loop, and
    apply the body permutation `sigma` at the end.

    `sigma` is the cyclic relabelling the choreography performs over one short
    period: after T/N the body that was at slot i sits where sigma(i) was. The
    unlabelled cloud is identical, so the diagram is identical, and the question
    is which diagram point each vine has landed on.
    """
    T = len(X) - 1
    slices = [mst_edges(P) for P in X]
    start = [(e[0], e[1]) for e in slices[0]]
    cur = list(range(len(start)))
    tracks = [[start[i]] for i in range(len(start))]

    labels = [tuple(e[:2]) for e in slices[0]]
    for t in range(1, T + 1):
        nxt = [tuple(e[:2]) for e in slices[t]]
        newlab = []
        for lab in labels:
            if lab in nxt:
                newlab.append(lab)
            else:
                # the MST changed which pair it uses: the edge that vanished is
                # replaced by the one that appeared
                gone = [l for l in labels if l not in nxt]
                came = [l for l in nxt if l not in labels]
                newlab.append(came[gone.index(lab)] if lab in gone and came else lab)
        labels = newlab
        for k, lab in enumerate(labels):
            tracks[k].append(lab)

    end = labels
    mapped = [tuple(sorted((sigma[a], sigma[b]))) for (a, b) in end]
    pos = {e: i for i, e in enumerate(start)}
    if any(e not in pos for e in mapped):
        return {"ok": False, "reason": "end edges are not a relabelling of the start",
                "start": start, "end": end, "mapped": mapped}
    perm = [pos[e] for e in mapped]
    return {"ok": True, "perm": perm, "order": perm_order(perm),
            "cycles": cycle_type(perm), "start": start, "end": end,
            "mapped": mapped, "tracks": tracks}
