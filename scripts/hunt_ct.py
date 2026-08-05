"""
Monodromy hunt in directly measured X-ray projections.

Why this setting, after eight negatives. The recurring obstruction, recorded in
notes/RECORD.md, is that for measured FIELDS the stability of the feature inventory
and its richness pull against each other: the threshold that stabilises the
inventory strips it below the complexity the mechanism needs. The one escape is
a field whose features are tied to PERSISTENT OBJECTS rather than to level sets,
because a fixed set of objects cannot change its inventory.

An X-ray projection is exactly that. Each dense structure inside the specimen
contributes a peak to the profile, at position

    s_i(theta) = x_i cos(theta) + y_i sin(theta),

so the peaks sweep sinusoids and cross one another as the gantry turns. The
structures are physical and permanent, so D7 is nearly free, and a projection of
a real specimen has many strong peaks, so D9 is free too. Both failure modes
that killed the galaxy, solar and atmospheric attempts are absent by
construction.

  domain      the detector coordinate u, an interval, one-point compactified.
              This is exact rather than a device: the attenuation really is zero
              in air on both sides of the specimen, so the two ends genuinely sit
              at the same value and joining them adds no structure.
  field       the MEASURED attenuation -log(I / I_air). Not a reconstruction:
              these are the raw detector counts.
  parameters  (gantry angle theta, detector row v). Both physical, both densely
              sampled by the instrument, neither chosen by me.
  loop        small, a few tenths of a degree by a few hundred micrometres.

A one-parameter loop in theta alone cannot work, and it is worth saying why,
because the gantry circle is tempting: after a full rotation every structure
returns to itself, so every vine closes on its own diagram point, and by the
same codimension count as N1 the order is 1. The forced loop is real but the
permutation is not. Two parameters are needed, and the cone-beam geometry
supplies the second for free.

Noise floor. Features are kept only if their persistence exceeds a multiple of
the measured photon noise, estimated from an air region of the same panel after
the same smoothing. This is an absolute floor, set by the instrument, not a
fraction of the strongest feature.

Smoothing scale. Set by the imaging system, not chosen: the detector and focal
spot blur the projection over a couple of detector pixels, so structure below
that scale is not resolved by the measurement in the first place.
"""
import sys, os, json, argparse
import numpy as np
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.general import (field_diagram, monodromy_abs, loop_points)
from mono.spectra import compactify

DEG_PER_PROJ = 0.5      # AngleInterval from the scanner metadata
UM_PER_ROW = 50.0       # detector PixelSize 0.050 mm
VSTEP = 2               # rows kept when the patch was fetched


class Walnut:
    """Measured attenuation profile along u, as a function of (theta, v)."""

    def __init__(self, sigma_pix=3.0, path="data/ct"):
        self.A = np.load(f"{path}/walnut_patch.npy")        # (theta, v, u)
        self.vs = np.load(f"{path}/walnut_patch_v.npy")
        self.sigma_pix = sigma_pix
        self.S = gaussian_filter1d(self.A.astype(np.float64), sigma_pix,
                                   axis=2, mode="nearest")
        self.nth, self.nv, self.nu = self.S.shape

    def noise(self):
        """Scatter of the smoothed attenuation in air, i.e. the detection floor."""
        air = np.concatenate([self.S[:, :, :40], self.S[:, :, -40:]], axis=2)
        return float(np.std(air))

    def profile(self, k, j):
        """Bilinear in (theta index k, row index j)."""
        k = float(np.clip(k, 0, self.nth - 1 - 1e-9))
        j = float(np.clip(j, 0, self.nv - 1 - 1e-9))
        k0, j0 = int(k), int(j)
        a, b = k - k0, j - j0
        k1, j1 = min(k0 + 1, self.nth - 1), min(j0 + 1, self.nv - 1)
        return ((1 - a) * (1 - b) * self.S[k0, j0]
                + a * (1 - b) * self.S[k1, j0]
                + (1 - a) * b * self.S[k0, j1]
                + a * b * self.S[k1, j1])

    def field(self, k, j):
        return compactify(self.profile(k, j))


def counts_on_loop(W, k, j, r, floor, steps=32):
    th = np.linspace(0, 2 * np.pi, steps, endpoint=False)
    cs = set()
    for t in th:
        d = field_diagram(W.field(k + r * np.cos(t), j + r * np.sin(t)))
        if d is None:
            return None
        p = d["points"]
        cs.add(int((p[:, 1] - p[:, 0] >= floor).sum()))
    return sorted(cs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sigma", type=float, default=3.0)
    ap.add_argument("--nsig", type=float, default=6.0, help="noise multiples")
    ap.add_argument("--rad", type=float, default=1.5, help="loop radius, index units")
    ap.add_argument("--steps", type=int, default=40)
    ap.add_argument("--grid", type=float, default=4.0)
    a = ap.parse_args()

    W = Walnut(a.sigma)
    nz = W.noise()
    floor = a.nsig * nz
    print(f"patch {W.S.shape}  (theta x row x detector)")
    print(f"  theta span {(W.nth-1)*DEG_PER_PROJ:.1f} deg, "
          f"row span {(W.nv-1)*VSTEP*UM_PER_ROW/1000:.2f} mm")
    print(f"  smoothing {a.sigma} detector pixels "
          f"({a.sigma*UM_PER_ROW:.0f} um at the detector)")
    print(f"  air noise sigma = {nz:.5f} in attenuation; "
          f"floor = {a.nsig} sigma = {floor:.5f}")
    print(f"  loop radius {a.rad*DEG_PER_PROJ:.2f} deg x "
          f"{a.rad*VSTEP*UM_PER_ROW:.0f} um")

    d0 = field_diagram(W.field(W.nth // 2, W.nv // 2))
    pers = d0["points"][:, 1] - d0["points"][:, 0]
    print(f"\n  center diagram: {len(pers)} points, "
          f"{int((pers>=floor).sum())} above the noise floor")

    print("\nD7: is the count above the floor constant round small loops?")
    ks = np.arange(a.rad + 1, W.nth - a.rad - 1, a.grid)
    js = np.arange(a.rad + 1, W.nv - a.rad - 1, a.grid)
    stable = 0; total = 0; hist = {}
    for j in js[::3]:
        for k in ks[::3]:
            c = counts_on_loop(W, k, j, a.rad, floor, steps=16)
            total += 1
            if c is not None and len(c) == 1:
                stable += 1
                hist[c[0]] = hist.get(c[0], 0) + 1
    print(f"  constant on {stable} of {total} loops ({100*stable/max(total,1):.0f}%)")
    print(f"  counts where constant: {dict(sorted(hist.items()))}")

    print("\nScanning for monodromy...")
    hits = []; ok = 0; tot = 0
    for j in js:
        for k in ks:
            r = monodromy_abs(W.field, loop_points(k, j, a.rad, a.steps), floor)
            tot += 1
            if r["ok"]:
                ok += 1
                if r["order"] and r["order"] > 1:
                    hits.append({"k": float(k), "j": float(j),
                                 "theta_deg": float(k * DEG_PER_PROJ),
                                 "row": int(W.vs[min(int(j), len(W.vs) - 1)]),
                                 "order": int(r["order"]), "n": int(r["n"]),
                                 "changes": int(r["pairing_changes"]),
                                 "cycles": [list(c) for c in r["cycles"]],
                                 "min_gap": r.get("min_gap")})
        print(f"  row j={j:.0f}: {len(hits)} hits so far", flush=True)
    print(f"\n  well-defined on {ok} of {tot} loops")
    print(f"  NONTRIVIAL monodromy on {len(hits)}")
    for h in hits[:30]:
        print(f"    theta {h['theta_deg']:5.2f} deg, row {h['row']}: order {h['order']}, "
              f"{h['n']} pts, {h['changes']} changes, gap {h['min_gap']:.2e}")
    os.makedirs("out", exist_ok=True)
    json.dump({"sigma": a.sigma, "nsig": a.nsig, "floor": floor, "rad": a.rad,
               "noise": nz, "ok": ok, "total": tot, "hits": hits},
              open(f"out/ct_walnut_s{a.sigma}_n{a.nsig}_r{a.rad}_st{a.steps}.json", "w"), indent=1)


if __name__ == "__main__":
    main()
