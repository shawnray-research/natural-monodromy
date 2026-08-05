"""Loading and conditioning real measured planar outlines."""

from __future__ import annotations

import json
import numpy as np

from .core import resample_closed, smooth_closed


def load_natural_earth_rings(path):
    """Every exterior ring of a Natural Earth land/country geojson, as lon/lat arrays."""
    with open(path) as f:
        d = json.load(f)
    rings = []
    for ft in d["features"]:
        g = ft["geometry"]
        name = (ft.get("properties") or {}).get("NAME") or \
               (ft.get("properties") or {}).get("name") or "?"
        if g["type"] == "Polygon":
            polys = [g["coordinates"]]
        elif g["type"] == "MultiPolygon":
            polys = g["coordinates"]
        else:
            continue
        for poly in polys:
            ext = np.asarray(poly[0], dtype=float)
            if len(ext) > 3:
                if np.allclose(ext[0], ext[-1]):
                    ext = ext[:-1]
                rings.append({"name": name, "lonlat": ext})
    return rings


def project_local(lonlat):
    """Equirectangular projection about the ring's centroid, scaled by cos(lat)
    so that the result is a faithful local planar shape (km-like units)."""
    lon = lonlat[:, 0]
    lat = lonlat[:, 1]
    lon0 = np.mean(lon)
    lat0 = np.mean(lat)
    x = (lon - lon0) * np.cos(np.deg2rad(lat0)) * 111.32
    y = (lat - lat0) * 110.57
    return np.column_stack([x, y])


def normalize(P):
    """Center and scale a closed polyline to unit RMS radius."""
    Q = P - P.mean(axis=0)
    s = np.sqrt((Q ** 2).sum(axis=1).mean())
    return Q / s


def condition(P, n=1500, smooth_frac=0.012):
    """
    Resample to n points and apply periodic Gaussian smoothing at a scale given
    as a fraction of the number of samples.  Real coastlines are fractal; the
    smoothing selects a scale, exactly as in scale-space shape analysis.  The
    monodromy findings are reported across a range of scales.
    """
    Q = resample_closed(P, n)
    if smooth_frac > 0:
        Q = smooth_closed(Q, smooth_frac * n)
        Q = resample_closed(Q, n)
    return normalize(Q)


def is_simple(P, step=1):
    """Self-intersection test on a closed polyline (vectorized, O(n^2/step))."""
    n = len(P)
    Q = np.vstack([P, P[:1]])
    A = Q[:-1]
    B = Q[1:]
    D = B - A
    for i in range(0, n, step):
        p, r = A[i], D[i]
        idx = np.arange(n)
        idx = idx[(idx != i) & (idx != (i - 1) % n) & (idx != (i + 1) % n)]
        q = A[idx]
        s = D[idx]
        rxs = r[0] * s[:, 1] - r[1] * s[:, 0]
        qp = q - p
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (qp[:, 0] * s[:, 1] - qp[:, 1] * s[:, 0]) / rxs
            u = (qp[:, 0] * r[1] - qp[:, 1] * r[0]) / rxs
        good = np.isfinite(t) & (t > 1e-9) & (t < 1 - 1e-9) & (u > 1e-9) & (u < 1 - 1e-9)
        if good.any():
            return False
    return True
