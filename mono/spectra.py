"""
Adapter that lets the validated S^1 machinery act on measured spectra and
profiles, which live on an interval rather than on a circle.

One-point compactification.  Prepend a single value LOW strictly below the
whole profile and read the result cyclically:

    LOW -> f_0 -> f_1 -> ... -> f_{n-1} -> LOW

LOW is then a strict local minimum, the two ends of the profile are joined
through it, and every other critical point of f is unchanged.  On the resulting
circle the numbers of minima and maxima agree, as they must, and the extra pair
(LOW, global max) is exactly the essential class of the profile.  Nothing about
the interior structure of the spectrum is altered, so the pairing that the
elder-rule sweep produces for the interior critical points is the ordinary
sublevel-set pairing of the profile.

The convention does treat a monotone run into either end of the recorded window
as terminating at a peak.  That is the honest reading of a truncated
measurement, but to keep it from mattering the analysis windows below are cut
at spectral minima wherever possible.
"""

from __future__ import annotations

import numpy as np

from .core import extended_persistence_circle


def compactify(f, pad_frac=0.25):
    """Profile on an interval -> function on S^1, by one added low value."""
    f = np.asarray(f, dtype=float)
    rng = float(f.max() - f.min())
    low = float(f.min()) - max(pad_frac * rng, 1e-12)
    return np.concatenate([[low], f])


def spectrum_diagram(f, pad_frac=0.25):
    """
    Extended persistence of a measured profile on an interval.

    Returns the same dict shape as mono.general.field_diagram, with indices
    referred back to the ORIGINAL profile (index -1 denotes the added low
    value, i.e. the essential class).
    """
    g = compactify(f, pad_frac)
    r = extended_persistence_circle(g)
    if r is None:
        return None
    pts = np.array([[g[r["mins"][i]], g[r["maxs"][j]]] for (i, j, _) in r["pairs"]])
    return {"points": pts, "pairs": r["pairs"],
            "mins": r["mins"] - 1, "maxs": r["maxs"] - 1,
            "f": np.asarray(f, dtype=float), "g": g}


def smooth(f, sigma):
    """Gaussian smoothing of a profile with reflecting ends (no wraparound)."""
    f = np.asarray(f, dtype=float)
    if sigma <= 0:
        return f
    k = int(np.ceil(4 * sigma))
    x = np.arange(-k, k + 1)
    w = np.exp(-0.5 * (x / sigma) ** 2)
    w /= w.sum()
    ext = np.concatenate([f[k:0:-1], f, f[-2:-k - 2:-1]])
    if len(ext) != len(f) + 2 * k:
        ext = np.pad(f, k, mode="reflect")
    return np.convolve(ext, w, mode="same")[k:k + len(f)]


def n_extrema(f):
    d = spectrum_diagram(f)
    if d is None:
        return 0, 0
    return len(d["mins"]), len(d["maxs"])


def field_fn_from_cube(cube, wl_slice, sigma=0.0, normalize=True):
    """
    Build a field_fn(u, v) -> profile, suitable for mono.general, from a
    spectral cube of shape (n_wavelength, ny, nx).  u, v are pixel coordinates
    and are interpolated bilinearly so that loops can be traversed continuously
    rather than only on the pixel lattice.
    """
    nz, ny, nx = cube.shape
    sub = cube[wl_slice]

    def field(u, v):
        u = float(np.clip(u, 0, nx - 1.001))
        v = float(np.clip(v, 0, ny - 1.001))
        i0, j0 = int(np.floor(u)), int(np.floor(v))
        a, b = u - i0, v - j0
        s = ((1 - a) * (1 - b) * sub[:, j0, i0] + a * (1 - b) * sub[:, j0, i0 + 1]
             + (1 - a) * b * sub[:, j0 + 1, i0] + a * b * sub[:, j0 + 1, i0 + 1])
        if sigma > 0:
            s = smooth(s, sigma)
        if normalize:
            m, M = s.min(), s.max()
            if M > m:
                s = (s - m) / (M - m)
        return s

    return field
