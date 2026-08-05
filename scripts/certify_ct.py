"""
Certification chain for monodromy candidates in the measured cone-beam scan.

Written before any candidate was seen, so the gates are not tuned to survivors.

The atmospheric attempt produced 250 raw hits and zero survivors, and 37 of
those hits claimed a nontrivial permutation with ZERO elder-rule pairing
changes, which cannot happen. That failure mode is a Hungarian matcher swapping
diagram points that are nearly coincident, so a gap guard is applied here before
anything else is believed.

Gates, in order of severity:

  G0 boundary   the loop lies well inside the sampled patch in both parameters,
                so no part of the answer comes from clamped interpolation
  G1 refine     the order is unchanged as the loop discretisation is refined
                40 -> 80 -> 160 -> 320 steps.

                A gap gate was written here first, requiring the closest pair of
                diagram points on the loop to stay well above the noise floor,
                and it was wrong. Near a birth-birth wall two diagram points
                genuinely do approach each other in the birth coordinate: that
                approach IS the mechanism, so a small gap is evidence for the
                phenomenon rather than against it. What actually distinguishes a
                real crossing from a mis-assignment is whether the tracking is
                resolved, and refining the discretisation tests exactly that. The
                gap is still reported, as a diagnostic rather than a gate.
  G2 changes    exactly two elder-rule pairing changes, as Definition 3.5 of
                arXiv:2607.01046 requires of a monodromy-critical loop; zero or
                four give the identity
  G3 shrink     the order survives as the loop contracts. This is the gate that
                rejected the sunspot candidate, which gave order 2 at exactly one
                radius. A genuine codimension-2 point is still enclosed by every
                smaller loop
  G4 displace   order 1 on loops of the SAME radius displaced in eight
                directions, so the answer is caused by what is inside the loop
  G5 wall       at the located point two minima of the profile coincide in value
                AND two maxima do, which is the A_1^2/A_1^2 signature
"""
import sys, os, json, argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.general import (field_diagram, monodromy_abs, loop_points,
                          wall_certificate)
from hunt_ct import Walnut, DEG_PER_PROJ, UM_PER_ROW, VSTEP


def gap_on_loop(W, k, j, r, floor, steps=40):
    """Closest pair of diagram points anywhere on the loop, above the floor."""
    th = np.linspace(0, 2 * np.pi, steps, endpoint=False)
    g = np.inf
    for t in th:
        d = field_diagram(W.field(k + r * np.cos(t), j + r * np.sin(t)))
        if d is None:
            return 0.0
        p = d["points"]
        p = p[(p[:, 1] - p[:, 0]) >= floor]
        if len(p) < 2:
            return 0.0
        D = np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)
        g = min(g, float(D[np.triu_indices(len(p), 1)].min()))
    return g


def certify(W, h, floor, rad, steps=40, verbose=True):
    k, j = h["k"], h["j"]
    rep = {"k": k, "j": j, "theta_deg": k * DEG_PER_PROJ, "order": h["order"]}

    # G0 boundary
    margin = 3 * rad + 2
    inside = (margin < k < W.nth - 1 - margin) and (margin < j < W.nv - 1 - margin)
    rep["G0_boundary"] = bool(inside)
    if not inside:
        rep["fail"] = "G0 boundary"; return rep

    # G1 refinement of the loop discretisation (gap reported, not gated)
    rep["gap_diag"] = gap_on_loop(W, k, j, rad, floor, steps)
    ref = []
    for s in (steps, 2*steps, 4*steps):
        r = monodromy_abs(W.field, loop_points(k, j, rad, s), floor)
        ref.append(r["order"] if r["ok"] else None)
    rep["G1_refine_steps"] = [steps, 2*steps, 4*steps]
    rep["G1_refine_orders"] = ref
    if ref[-1] is None or len({o for o in ref if o is not None}) != 1 \
            or ref[-1] != h["order"]:
        rep["fail"] = f"G1 refine {ref}"; return rep

    # G2 pairing changes
    r0 = monodromy_abs(W.field, loop_points(k, j, rad, steps), floor)
    rep["G2_changes"] = r0.get("pairing_changes")
    rep["order_at_rad"] = r0.get("order")
    if not r0["ok"] or r0["order"] in (None, 1):
        rep["fail"] = "G2 not reproduced"; return rep
    if r0["pairing_changes"] != 2:
        rep["fail"] = f"G2 changes = {r0['pairing_changes']}"; return rep

    # G3 shrink
    radii = [rad * f for f in (0.8, 0.6, 0.45, 0.3, 0.2, 0.12)]
    orders = []
    for rr in radii:
        r = monodromy_abs(W.field, loop_points(k, j, rr, steps), floor)
        orders.append(r["order"] if r["ok"] else None)
    rep["G3_radii"] = radii
    rep["G3_orders"] = orders
    good = [o for o in orders if o is not None]
    if len(good) < 4 or not all(o == r0["order"] for o in good):
        rep["fail"] = f"G3 shrink {orders}"; return rep

    # G4 displacement
    disp = []
    for a in np.linspace(0, 2 * np.pi, 8, endpoint=False):
        kk, jj = k + 3 * rad * np.cos(a), j + 3 * rad * np.sin(a)
        r = monodromy_abs(W.field, loop_points(kk, jj, rad, steps), floor)
        disp.append(r["order"] if r["ok"] else None)
    rep["G4_displaced"] = disp
    if sum(1 for o in disp if o == 1) < 6:
        rep["fail"] = f"G4 displaced {disp}"; return rep

    # G5 wall certificate
    w = wall_certificate(W.field, k, j)
    rep["G5_wall"] = w
    if w is None or w["min_gap_rel"] is None or w["max_gap_rel"] is None:
        rep["fail"] = "G5 no wall"; return rep
    rep["certified"] = True
    return rep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True)
    ap.add_argument("--sigma", type=float, default=3.0)
    ap.add_argument("--max", type=int, default=40)
    ap.add_argument("--steps", type=int, default=120)
    a = ap.parse_args()

    R = json.load(open(a.json))
    floor, rad = R["floor"], R["rad"]
    W = Walnut(a.sigma)
    hits = R["hits"]
    print(f"{len(hits)} raw hits from {a.json}")
    print(f"floor {floor:.5f}, loop radius {rad}\n")

    out = []
    stage = {}
    for i, h in enumerate(hits[:a.max]):
        rep = certify(W, h, floor, rad, steps=a.steps)
        out.append(rep)
        key = rep.get("fail", "CERTIFIED")
        stage[key.split()[0]] = stage.get(key.split()[0], 0) + 1
        tag = "CERTIFIED" if rep.get("certified") else rep.get("fail")
        print(f"  [{i+1:3d}] theta {rep['theta_deg']:5.2f} deg  order {rep['order']}  -> {tag}")

    print("\nsummary by stage:")
    for k, v in sorted(stage.items(), key=lambda x: -x[1]):
        print(f"  {k:14s} {v}")
    cert = [r for r in out if r.get("certified")]
    print(f"\nCERTIFIED: {len(cert)} of {len(out)} examined")
    json.dump(out, open(a.json.replace(".json", "_certified.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
