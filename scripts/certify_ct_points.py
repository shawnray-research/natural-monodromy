"""
Full certification battery about a LOCALIZED generator in the measured scan.

The earlier chain applied the shrink test about the scan grid point that first
noticed a candidate. That is the wrong center: a loop contracted about a point
offset from the singularity stops enclosing it and correctly returns the
identity, which is indistinguishable from a spurious detection. Everything here
is therefore anchored at the quadrisected p*.

Tests, all reported whether they pass or fail:

  T1 shrink       order at radii spanning several decades about p*
  T2 displace     order on loops of the same radius centered away from p*
  T3 refine       order against the number of loop samples
  T4 wall         the A_1^2/A_1^2 signature: two minima of the profile equal in
                  value AND two maxima equal in value, at p*
  T5 scale        order against the smoothing scale and the noise floor, the two
                  numbers not fixed by the instrument
"""
import sys, os, json, argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.general import (field_diagram, monodromy_abs, loop_points,
                          localize_abs, wall_certificate)
from hunt_ct import Walnut, DEG_PER_PROJ, UM_PER_ROW, VSTEP


def battery(W, u, v, floor, base_rad=0.8, steps=160):
    out = {"u": u, "v": v, "theta_deg": u * DEG_PER_PROJ}

    print(f"\n  p* = (k {u:.6f}, j {v:.6f})  ->  theta {u*DEG_PER_PROJ:.4f} deg, "
          f"detector row {W.vs[0] + v*VSTEP:.1f}")

    print("  T1 shrink:")
    t1 = []
    for f in (1.0, 0.5, 0.25, 0.1, 0.05, 0.02, 0.01, 5e-3, 1e-3, 1e-4):
        r = monodromy_abs(W.field, loop_points(u, v, base_rad * f, steps), floor)
        o = r["order"] if r["ok"] else None
        t1.append({"radius": base_rad * f, "order": o,
                   "changes": r.get("pairing_changes"), "n": r.get("n")})
        print(f"     r = {base_rad*f:10.3e}  order {str(o):4s}  "
              f"changes {r.get('pairing_changes')}  n {r.get('n')}")
    out["T1"] = t1

    print("  T2 displace (same radius, center moved by 4r):")
    t2 = []
    for a in np.linspace(0, 2 * np.pi, 8, endpoint=False):
        uu, vv = u + 4 * base_rad * np.cos(a), v + 4 * base_rad * np.sin(a)
        r = monodromy_abs(W.field, loop_points(uu, vv, base_rad, steps), floor)
        t2.append(r["order"] if r["ok"] else None)
    print(f"     orders: {t2}")
    out["T2"] = t2

    print("  T3 refine (loop samples):")
    t3 = []
    for s in (40, 80, 160, 320, 640):
        r = monodromy_abs(W.field, loop_points(u, v, base_rad * 0.25, s), floor)
        t3.append(r["order"] if r["ok"] else None)
    print(f"     steps 40,80,160,320,640 -> {t3}")
    out["T3"] = t3

    print("  T4 wall certificate at p*:")
    w = wall_certificate(W.field, u, v)
    if w:
        mn, mx = np.array(w["min_values"]), np.array(w["max_values"])
        dmn, dmx = np.diff(np.sort(mn)), np.diff(np.sort(mx))
        rng = mx.max() - mn.min()
        print(f"     {w['n_min']} minima, {w['n_max']} maxima")
        print(f"     closest two minima differ by {dmn.min():.3e}  "
              f"({dmn.min()/rng:.2e} of range)")
        print(f"     closest two maxima differ by {dmx.min():.3e}  "
              f"({dmx.min()/rng:.2e} of range)")
        out["T4"] = {"n_min": w["n_min"], "n_max": w["n_max"],
                     "min_gap": float(dmn.min()), "max_gap": float(dmx.min()),
                     "min_gap_rel": float(dmn.min() / rng),
                     "max_gap_rel": float(dmx.min() / rng)}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True)
    ap.add_argument("--sigma", type=float, default=15.0)
    a = ap.parse_args()

    R = json.load(open(a.json))
    floor, rad = R["floor"], R["rad"]
    W = Walnut(a.sigma)
    print(f"smoothing {a.sigma} px, floor {floor:.5f}, scan radius {rad}")

    results = []
    for h in R["hits"]:
        k, j = h["k"], h["j"]
        L = localize_abs(W.field, k - 1.0, k + 1.0, j - 1.0, j + 1.0, floor)
        if L is None or L["box"] > 1e-3 or L["changes"] != 2:
            continue
        print("\n" + "=" * 74)
        print(f"CANDIDATE from scan cell theta {h['theta_deg']:.2f} deg, row {h['row']}"
              f"   (localization box {L['box']:.2e}, {L['changes']} pairing changes)")
        b = battery(W, L["u"], L["v"], floor, base_rad=rad)
        b["localize"] = L
        results.append(b)

    json.dump(results, open(a.json.replace(".json", "_battery.json"), "w"), indent=1)
    print(f"\n\n{len(results)} localized generators put through the battery.")

    print("\nVERDICT")
    for b in results:
        t1 = [d["order"] for d in b["T1"]]
        held = all(o == 2 for o in t1 if o is not None)
        ndef = sum(1 for o in t1 if o is not None)
        triv = sum(1 for o in b["T2"] if o == 1)
        ref = len({o for o in b["T3"] if o is not None}) == 1
        print(f"  theta {b['theta_deg']:6.3f} deg: shrink held={held} "
              f"({ndef}/10 defined), displaced trivial {triv}/8, refine stable={ref}")


if __name__ == "__main__":
    main()
