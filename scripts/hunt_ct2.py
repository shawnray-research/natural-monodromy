"""
Second pass at the measured cone-beam scan, on a SIGNAL-DOMINATED parameter grid.

The first pass localized generators to boxes of 1e-8 that had the correct
A_1^2/A_1^2 structure, one birth-birth crossing and one death-death crossing,
and that survived replacing the bilinear interpolant with a bicubic one. They
were still not real, and the reason is measurable:

    change in the profile per gantry step (0.5 deg)   0.00529
    change in the profile per row step   (100 um)     0.00338
    air noise of the same smoothed profile           0.00390

Consecutive samples differ by LESS than the noise in v and by only 1.4x the
noise in theta. A two-parameter family sampled that finely is, at the scale of
one step, a family of noise realisations, so a codimension-2 point localized
inside one cell describes the noise rather than the walnut. That is why those
candidates survived a change of interpolant, which re-reads the same noisy
samples, and died under every coarsening of the sampling.

The fix is not more smoothing, which destroys the structure along with the
noise. It is to take parameter steps large enough that the field moves by more
than the noise between them:

    theta stride 4 projections = 2.0 deg     -> about 5x the noise per step
    row   stride 8 rows        = 800 um      -> about 7x the noise per step

The generator is then a property of the specimen, because between adjacent
samples the specimen has changed the profile by much more than the measurement
error. Loops still enclose a single cell, so they remain small in the sense that
matters: the feature inventory does not reorganise around them.
"""
import sys, os, json, argparse
import numpy as np
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.general import (field_diagram, monodromy_abs, loop_points,
                          localize_abs, rank_pairing_key, wall_certificate)
from mono.spectra import compactify
from hunt_ct import Walnut, DEG_PER_PROJ, UM_PER_ROW, VSTEP


class Coarse(Walnut):
    """Measured profiles on a stride grid, so each step beats the noise."""

    def __init__(self, sigma_pix=15.0, kstride=4, jstride=8, koff=0, joff=0):
        super().__init__(sigma_pix)
        self.S = self.S[koff::kstride, joff::jstride]
        self.nth, self.nv = self.S.shape[0], self.S.shape[1]
        self.kstride, self.jstride = kstride, jstride
        self.koff, self.joff = koff, joff

    def theta_deg(self, k):
        return (self.koff + k * self.kstride) * DEG_PER_PROJ

    def row(self, j):
        return self.vs[0] + (self.joff + j * self.jstride) * VSTEP


def step_snr(W):
    nz = W.noise()
    dk, dj = [], []
    for k in range(2, W.nth - 2, max(1, W.nth // 12)):
        for j in range(2, W.nv - 2, max(1, W.nv // 12)):
            f0 = W.profile(k, j)
            dk.append(np.abs(W.profile(k + 1, j) - f0).mean())
            dj.append(np.abs(W.profile(k, j + 1) - f0).mean())
    return nz, float(np.mean(dk)) / nz, float(np.mean(dj)) / nz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sigma", type=float, default=15.0)
    ap.add_argument("--nsig", type=float, default=3.0)
    ap.add_argument("--ks", type=int, default=4)
    ap.add_argument("--js", type=int, default=8)
    ap.add_argument("--rad", type=float, default=0.45)
    ap.add_argument("--steps", type=int, default=160)
    a = ap.parse_args()

    W = Coarse(a.sigma, a.ks, a.js)
    nz, sk, sj = step_snr(W)
    floor = a.nsig * nz
    print(f"grid {W.nth} x {W.nv}  (theta stride {a.ks*DEG_PER_PROJ:.1f} deg, "
          f"row stride {a.js*VSTEP*UM_PER_ROW:.0f} um)")
    print(f"  noise {nz:.5f}, floor {floor:.5f}")
    print(f"  per-step field change / noise:  theta {sk:.1f}x,  row {sj:.1f}x")
    if min(sk, sj) < 3:
        print("  WARNING: a step is not clearly above the noise")

    ks = np.arange(a.rad + 1, W.nth - a.rad - 1, 0.5)
    js = np.arange(a.rad + 1, W.nv - a.rad - 1, 0.5)
    print(f"\nscanning {len(ks)}x{len(js)} loop centers, radius {a.rad}")

    hits = []; ok = 0; tot = 0
    for j in js:
        for k in ks:
            r = monodromy_abs(W.field, loop_points(k, j, a.rad, a.steps), floor)
            tot += 1
            if r["ok"]:
                ok += 1
                if r["order"] and r["order"] > 1:
                    hits.append({"k": float(k), "j": float(j), "order": int(r["order"]),
                                 "n": int(r["n"]), "changes": int(r["pairing_changes"])})
    print(f"  well-defined on {ok} of {tot};  nontrivial on {len(hits)}")

    print("\nlocalizing and diagnosing wall types:")
    good = []
    for h in hits:
        L = localize_abs(W.field, h["k"] - 0.5, h["k"] + 0.5,
                         h["j"] - 0.5, h["j"] + 0.5, floor)
        if L is None or L["box"] > 1e-4 or L["changes"] != 2:
            continue
        u, v = L["u"], L["v"]
        # wall types, tracker independent
        kinds = []
        for rad in (0.05, 0.005):
            B = []; D = []; K = []
            for t in np.linspace(0, 2*np.pi, 480, endpoint=False):
                d = field_diagram(W.field(u + rad*np.cos(t), v + rad*np.sin(t)))
                p = d["points"]; p = p[(p[:,1]-p[:,0]) >= floor]
                B.append(np.diff(np.sort(p[:,0])).min())
                D.append(np.diff(np.sort(p[:,1])).min())
                K.append(rank_pairing_key(p))
            ch = [i for i in range(480) if K[i] != K[i-1]]
            kk = []
            for i in ch:
                kk.append("BB" if min(B[i-1],B[i]) < min(D[i-1],D[i]) else "DD")
            kinds.append(kk)
        if sorted(kinds[-1]) != ["BB", "DD"]:
            print(f"   theta {W.theta_deg(u):6.2f} row {W.row(v):7.1f}: "
                  f"wall types {kinds[-1]} -> not A_1^2/A_1^2")
            continue
        print(f"   theta {W.theta_deg(u):6.2f} deg, row {W.row(v):7.1f}: "
              f"order {L['order']}, box {L['box']:.1e}, walls {kinds[-1]}  <== candidate")
        good.append({"u": u, "v": v, "order": L["order"], "box": L["box"],
                     "theta_deg": W.theta_deg(u), "row": W.row(v),
                     "walls": kinds[-1]})

    json.dump({"sigma": a.sigma, "nsig": a.nsig, "floor": floor, "rad": a.rad,
               "ks": a.ks, "js": a.js, "snr_theta": sk, "snr_row": sj,
               "n_hits": len(hits), "good": good},
              open(f"out/ct2_k{a.ks}_j{a.js}.json", "w"), indent=1)
    print(f"\n{len(good)} candidates with the correct wall structure.")


if __name__ == "__main__":
    main()
