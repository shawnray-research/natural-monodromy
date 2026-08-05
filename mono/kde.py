"""
Kernel density of a moving point cloud, and its vineyard.

Why not H_0 of the point cloud itself: under Rips or Cech every point is born at
radius 0, so the whole diagram sits on the line birth = 0 and the deaths are the
minimum spanning tree edge lengths. Order statistics of continuous functions are
continuous AND ordered, so those vines can touch but never cross. Vines confined
to a line have nowhere to go, and H_0 of a point cloud can never exhibit
monodromy. Verified on the figure-eight: order 1.

Smoothing the cloud into a density fixes exactly that. The field

    rho(x; t) = sum_i m_i exp(-|x - x_i(t)|^2 / 2 sigma^2)

is a genuine scalar field on the plane whose maxima sit at the bodies and whose
saddles sit between them, so births and deaths both vary and the diagram is
honestly two-dimensional. This is also the kernel density estimate whose
vineyard Turner's paper suggests studying.

Superlevel-set H_0 is computed by an elder-rule union-find sweep in decreasing
density, which records the grid cell of every birth and every death, so vines
can be followed by identity: each maximum is labelled by the body it belongs to.
"""

from __future__ import annotations

import numpy as np

from .core import perm_order, cycle_type


def kde_grid(pts, masses, sigma, lim, res):
    xs = np.linspace(-lim, lim, res)
    ys = np.linspace(-lim, lim, res)
    X, Y = np.meshgrid(xs, ys)
    F = np.zeros_like(X)
    for p, m in zip(pts, masses):
        F += m * np.exp(-((X - p[0]) ** 2 + (Y - p[1]) ** 2) / (2 * sigma ** 2))
    return F, X, Y


class _UF:
    def __init__(self):
        self.p = {}

    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a


def superlevel_h0(F):
    """
    Elder-rule superlevel sweep. Returns pairs (birth_flat, death_flat, kind)
    where birth_flat is the grid cell of a local maximum and death_flat is the
    saddle at which its component merges into an older one. The global maximum
    is the essential class.
    """
    res = F.shape[0]
    flat = F.ravel()
    order = np.argsort(-flat, kind="stable")
    uf = _UF()
    birth_of = {}
    added = np.zeros(flat.size, dtype=bool)
    pairs = []
    for p in order:
        p = int(p)
        added[p] = True
        r, c = divmod(p, res)
        roots = []
        for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1),
                       (-1, -1), (-1, 1), (1, -1), (1, 1)):
            rr, cc = r + dr, c + dc
            if 0 <= rr < res and 0 <= cc < res:
                q = rr * res + cc
                if added[q]:
                    roots.append(uf.find(q))
        roots = list(dict.fromkeys(roots))
        uf.p[p] = p
        if not roots:
            birth_of[p] = p
            continue
        elder = max(roots, key=lambda q: flat[birth_of[q]])
        for q in roots:
            if q != elder:
                pairs.append((birth_of[q], p, "ord"))
                uf.p[q] = elder
        uf.p[p] = elder
        birth_of[uf.find(elder)] = birth_of[elder]
    root = uf.find(int(order[0]))
    pairs.append((birth_of[root], int(order[-1]), "ext"))
    return pairs


def diagram_with_bodies(F, X, Y, pts):
    """
    Persistence pairs, with each birth labelled by the body whose maximum it is.
    Labelling by body is what makes the tracking exact: the bodies are the
    physical identities, and a choreography permutes them by construction.
    """
    pairs = superlevel_h0(F)
    xr, yr = X.ravel(), Y.ravel()
    fr = F.ravel()
    out = []
    for (b, d, k) in pairs:
        bx, by = xr[b], yr[b]
        j = int(np.argmin(((pts[:, 0] - bx) ** 2 + (pts[:, 1] - by) ** 2)))
        out.append({"body": j, "birth": float(fr[b]), "death": float(fr[d]),
                    "kind": k, "bcell": int(b), "dcell": int(d),
                    "bx": float(bx), "by": float(by),
                    "dx": float(xr[d]), "dy": float(yr[d])})
    return out


def monodromy_over_loop(traj, masses, sigma, lim, res, sigma_perm):
    """
    traj        : (T+1, N, 2) body positions over the loop, traj[-1] the cloud
                  identical to traj[0] as a SET
    sigma_perm  : the relabelling, body i at the start sits where body
                  sigma_perm[i] is at the end

    Vines are labelled by body. After the loop the labels are pushed through
    sigma_perm and compared with the start.
    """
    d0 = diagram_with_bodies(*kde_grid(traj[0], masses, sigma, lim, res)[:1],
                             *kde_grid(traj[0], masses, sigma, lim, res)[1:], traj[0])
    F, X, Y = kde_grid(traj[0], masses, sigma, lim, res)
    d0 = diagram_with_bodies(F, X, Y, traj[0])
    F1, X1, Y1 = kde_grid(traj[-1], masses, sigma, lim, res)
    d1 = diagram_with_bodies(F1, X1, Y1, traj[-1])
    return d0, d1
