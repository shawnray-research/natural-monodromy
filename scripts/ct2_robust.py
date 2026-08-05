"""
The decisive test for the signal-dominated candidates.

A generator that is a property of the walnut must be recoverable from a
DIFFERENT SUBSET of the measurements. The stride grid makes that a clean
experiment: changing the offset selects entirely different projections and
entirely different detector rows, and changing the stride changes both the
subset and the spacing. A generator of the specimen stays at the same physical
(theta, row); a generator of the noise does not survive at all.

This is the test the first pass failed. Its candidates were localized to 1e-8
inside a single cell whose two corners differed by less than the noise, and they
vanished under every resampling.
"""
import sys, os, json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.general import (field_diagram, monodromy_abs, loop_points,
                          localize_abs, rank_pairing_key, wall_certificate)
from hunt_ct import DEG_PER_PROJ, UM_PER_ROW, VSTEP
from hunt_ct2 import Coarse, step_snr

TARGETS = [(26.33, 714.0), (7.71, 957.0)]


def to_index(W, theta_deg, row):
    k = (theta_deg / DEG_PER_PROJ - W.koff) / W.kstride
    j = ((row - W.vs[0]) / VSTEP - W.joff) / W.jstride
    return k, j


def walls_at(W, u, v, floor, rad=0.005, steps=480):
    B = []; D = []; K = []
    for t in np.linspace(0, 2*np.pi, steps, endpoint=False):
        d = field_diagram(W.field(u + rad*np.cos(t), v + rad*np.sin(t)))
        if d is None:
            return None
        p = d["points"]; p = p[(p[:, 1]-p[:, 0]) >= floor]
        if len(p) < 2:
            return None
        B.append(np.diff(np.sort(p[:, 0])).min())
        D.append(np.diff(np.sort(p[:, 1])).min())
        K.append(rank_pairing_key(p))
    ch = [i for i in range(steps) if K[i] != K[i-1]]
    return [("BB" if min(B[i-1], B[i]) < min(D[i-1], D[i]) else "DD") for i in ch]


def probe(sigma, ks, js, koff, joff, nsig=3.0, label=""):
    W = Coarse(sigma, ks, js, koff, joff)
    nz, sk, sj = step_snr(W)
    floor = nsig * nz
    rows = []
    for (th, rw) in TARGETS:
        k, j = to_index(W, th, rw)
        if not (1 < k < W.nth - 2 and 1 < j < W.nv - 2):
            rows.append(None); continue
        L = localize_abs(W.field, k - 0.7, k + 0.7, j - 0.7, j + 0.7, floor)
        if L is None or L["changes"] != 2:
            rows.append(None); continue
        wl = walls_at(W, L["u"], L["v"], floor)
        rows.append({"theta": W.theta_deg(L["u"]), "row": W.row(L["v"]),
                     "order": L["order"], "box": L["box"],
                     "walls": sorted(wl) if wl else None})
    print(f"  {label:34s} snr {sk:.1f}/{sj:.1f}x  ", end="")
    for (t, r) in zip(TARGETS, rows):
        if r is None:
            print(f"| {'lost':>26s} ", end="")
        else:
            print(f"| th {r['theta']:6.2f} row {r['row']:6.1f} o{r['order']} ", end="")
    print()
    return rows


def main():
    print("Baseline and then DIFFERENT SUBSETS of the same measured scan.")
    print("targets:", [f"theta {t:.2f} deg, row {r:.0f}" for t, r in TARGETS], "\n")
    allr = {}
    allr["base"] = probe(15.0, 4, 8, 0, 0, label="stride 4/8, offset 0/0 (base)")

    print("\n different GANTRY ANGLES (offset in theta):")
    for ko in (1, 2, 3):
        allr[f"koff{ko}"] = probe(15.0, 4, 8, ko, 0, label=f"stride 4/8, theta offset {ko}")

    print("\n different DETECTOR ROWS (offset in v):")
    for jo in (2, 4, 6):
        allr[f"joff{jo}"] = probe(15.0, 4, 8, 0, jo, label=f"stride 4/8, row offset {jo}")

    print("\n different STRIDES (different subset AND spacing):")
    for (ks, js) in ((3, 6), (5, 10), (6, 12), (3, 10), (5, 6)):
        allr[f"s{ks}_{js}"] = probe(15.0, ks, js, 0, 0, label=f"stride {ks}/{js}")

    print("\n different SMOOTHING:")
    for sg in (12.0, 18.0, 22.0):
        allr[f"sig{sg}"] = probe(sg, 4, 8, 0, 0, label=f"sigma {sg} px, stride 4/8")

    print("\n different NOISE FLOOR:")
    for ns in (2.0, 4.0, 5.0):
        allr[f"ns{ns}"] = probe(15.0, 4, 8, 0, 0, nsig=ns, label=f"floor {ns} sigma")

    json.dump(allr, open("out/ct2_robust.json", "w"), indent=1, default=str)

    print("\n\nSUMMARY: how often is each target recovered, and how far does it move?")
    for i, (th, rw) in enumerate(TARGETS):
        got = [v[i] for v in allr.values() if v[i] is not None]
        n = len(allr)
        if not got:
            print(f"  target theta {th}: recovered 0/{n}")
            continue
        dth = [abs(g["theta"] - th) for g in got]
        drw = [abs(g["row"] - rw) for g in got]
        ords = {g["order"] for g in got}
        wl = [tuple(g["walls"]) if g["walls"] else None for g in got]
        print(f"  target theta {th:5.2f}, row {rw:5.0f}: recovered {len(got)}/{n}")
        print(f"     orders {ords},  wall types {set(wl)}")
        print(f"     theta spread max {max(dth):.3f} deg,  row spread max {max(drw):.1f} "
              f"({max(drw)*UM_PER_ROW:.0f} um)")


if __name__ == "__main__":
    main()
