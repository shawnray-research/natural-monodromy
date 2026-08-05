"""
Deep certification of monodromy in a single named real coastline.

Robustness battery (each is a way the finding could be an artifact; all must pass):
  R1 curve sampling density      -- singularity must converge, order must be stable
  R2 loop discretization         -- order must be stable as the loop is refined
  R3 smoothing scale             -- must survive over a band of scale-space scales
  R4 measurement noise           -- must survive perturbation of the raw lon/lat vertices
  R5 causal toggle               -- shrink onto it (persists) / displace off it (vanishes)
  R6 explicit bitangent circles  -- the A_1^2/A_1^2 geometry, reported by hand-checkable numbers
"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.shapes import (load_natural_earth_rings, project_local, condition,
                         is_simple)
from mono.core import circle_loop, radial_diagram, resample_closed, smooth_closed
from mono.scan import scan_plane, monodromy_on_loop
from mono.localize import localize, full_certificate
from mono.robust import find_nearest_monodromy

ISLANDS = ["Ireland", "Iceland", "Sri Lanka", "Cuba", "Japan", "New Zealand",
           "Jamaica", "Cyprus", "Taiwan", "Haiti", "Dominican Rep."]


def get_ring(name, rings):
    cands = [r for r in rings if r["name"] == name]
    if not cands:
        return None
    return max(cands, key=lambda r: len(r["lonlat"]))


def find_best(M, grid=72):
    lo = M.min(axis=0); hi = M.max(axis=0)
    pad = 0.12 * (hi - lo)
    res = scan_plane(M, (lo[0] - pad[0], hi[0] + pad[0]),
                        (lo[1] - pad[1], hi[1] + pad[1]),
                     nx=grid, ny=grid, loop_steps=40, loop_frac=0.8, tau=0.05)
    out, seen = [], []
    for h in res["hits"]:
        if any((h["x"] - s[0]) ** 2 + (h["y"] - s[1]) ** 2 < (2 * res["radius"]) ** 2
               for s in seen):
            continue
        r = res["radius"]
        loc = localize(M, h["x"] - r, h["x"] + r, h["y"] - r, h["y"] + r,
                       depth=30, steps=48)
        if loc is None:
            continue
        cx, cy, hw, order = loc
        cert = full_certificate(M, cx, cy, r / 30, steps=256)
        c3 = cert["C3"]
        if (cert["C1_order"] and cert["C1_order"] > 1
                and cert["C2_pairing_changes"] == 2
                and c3["min_gap_rel"] < 1e-6 and c3["max_gap_rel"] < 1e-6):
            seen.append((cx, cy))
            out.append({"x": cx, "y": cy, "cert": cert, "scan_radius": r})
    return out, res


def bitangent_report(M, p):
    """R6: the explicit A_1^2/A_1^2 geometry at p."""
    d = radial_diagram(M, np.asarray(p, float))
    f = d["f"]
    mins, maxs = d["mins"], d["maxs"]
    vmin = f[mins]; vmax = f[maxs]
    om = np.argsort(vmin); oM = np.argsort(vmax)
    # closest pair of minima, closest pair of maxima
    i1, i2 = om[0], om[1]
    gaps = np.abs(vmin[om][1:] - vmin[om][:-1])
    kmin = int(np.argmin(gaps))
    a1, a2 = om[kmin], om[kmin + 1]
    gapsM = np.abs(vmax[oM][1:] - vmax[oM][:-1])
    kmax = int(np.argmin(gapsM))
    b1, b2 = oM[kmax], oM[kmax + 1]
    return {
        "inner_circle_radius": float(np.sqrt(0.5 * (vmin[a1] + vmin[a2]))),
        "inner_tangency_points": [M[mins[a1]].tolist(), M[mins[a2]].tolist()],
        "inner_radius_mismatch": float(abs(np.sqrt(vmin[a1]) - np.sqrt(vmin[a2]))),
        "outer_circle_radius": float(np.sqrt(0.5 * (vmax[b1] + vmax[b2]))),
        "outer_tangency_points": [M[maxs[b1]].tolist(), M[maxs[b2]].tolist()],
        "outer_radius_mismatch": float(abs(np.sqrt(vmax[b1]) - np.sqrt(vmax[b2]))),
        "n_minima": len(mins), "n_maxima": len(maxs),
    }


def battery(name, lonlat, base_smooth=0.010, n_samp=1400):
    P = project_local(lonlat)
    M = condition(P, n=n_samp, smooth_frac=base_smooth)
    if not is_simple(M, step=2):
        return {"name": name, "status": "self-intersecting after smoothing"}
    hits, _ = find_best(M)
    if not hits:
        return {"name": name, "status": "no certified monodromy"}
    h = max(hits, key=lambda z: z["cert"]["C1_n_diagram_points"])
    cx, cy = h["x"], h["y"]
    r0 = h["scan_radius"]
    rep = {"name": name, "status": "ok", "n_certified": len(hits),
           "point": [cx, cy], "base_smooth": base_smooth, "n_samp": n_samp,
           "cert": h["cert"], "R6_bitangent": bitangent_report(M, [cx, cy])}

    # R1 sampling density
    r1 = []
    for N in (700, 1400, 2800, 5600):
        MM = condition(P, n=N, smooth_frac=base_smooth)
        b = find_nearest_monodromy(MM, (cx, cy), 2.0 * r0)
        r1.append({"N": N, "order": None if b is None else b["order"],
                   "drift": None if b is None else b["drift"]})
    rep["R1_sampling"] = r1

    # R2 loop discretization
    r2 = []
    for steps in (32, 64, 128, 256, 1024, 4096):
        res = monodromy_on_loop(M, circle_loop([cx, cy], r0 / 30, steps))
        r2.append({"steps": steps, "order": res["order"], "ok": res["ok"]})
    rep["R2_loop_steps"] = r2

    # R3 smoothing scale
    r3 = []
    for sf in (0.006, 0.008, 0.010, 0.013, 0.017, 0.022):
        MM = condition(P, n=n_samp, smooth_frac=sf)
        b = find_nearest_monodromy(MM, (cx, cy), 3.0 * r0)
        r3.append({"smooth_frac": sf, "order": None if b is None else b["order"],
                   "drift": None if b is None else b["drift"]})
    rep["R3_smoothing"] = r3

    # R4 measurement noise on the raw vertices (in km, before smoothing)
    rng = np.random.default_rng(11)
    r4 = []
    for eps_km in (0.0, 0.5, 1.0, 2.0, 5.0, 10.0):
        orders, drifts = [], []
        for trial in range(5):
            PN = P + rng.normal(0, eps_km, P.shape)
            MM = condition(PN, n=n_samp, smooth_frac=base_smooth)
            b = find_nearest_monodromy(MM, (cx, cy), 3.0 * r0)
            orders.append(None if b is None else b["order"])
            drifts.append(None if b is None else b["drift"])
        r4.append({"noise_km": eps_km, "orders": orders, "drifts": drifts,
                   "survived": sum(1 for o in orders if o and o > 1)})
    rep["R4_noise"] = r4

    # R5 toggle
    small = r0 / 30
    shrink = []
    for f in (1.0, 0.3, 0.1, 0.03, 0.01, 0.003):
        res = monodromy_on_loop(M, circle_loop([cx, cy], small * f, 384))
        shrink.append({"radius": small * f, "order": res["order"]})
    displace = []
    for dx, dy in [(3, 0), (-3, 0), (0, 3), (0, -3), (2, 2), (-2, -2)]:
        res = monodromy_on_loop(M, circle_loop(
            [cx + dx * small, cy + dy * small], small, 384))
        displace.append({"offset": [dx * small, dy * small], "order": res["order"]})
    rep["R5_shrink"] = shrink
    rep["R5_displace"] = displace
    return rep


def main():
    rings = load_natural_earth_rings("data/ne10_countries.geojson")
    reports = []
    for name in ISLANDS:
        r = get_ring(name, rings)
        if r is None:
            print(f"{name}: not found")
            continue
        rep = battery(name, r["lonlat"])
        reports.append(rep)
        print("=" * 78, flush=True)
        print(f"{name}   ({len(r['lonlat'])} measured vertices)  -> {rep['status']}",
              flush=True)
        if rep["status"] != "ok":
            continue
        c = rep["cert"]
        b = rep["R6_bitangent"]
        print(f"  certified singularities on this coastline : {rep['n_certified']}")
        print(f"  hero point (local km coords)              : "
              f"({rep['point'][0]:.6f}, {rep['point'][1]:.6f})")
        print(f"  C1 permutation order                      : {c['C1_order']}  "
              f"cycles {c['C1_cycles']} over {c['C1_n_diagram_points']} diagram points")
        print(f"  C2 elder-rule pairing changes             : {c['C2_pairing_changes']}"
              f"   (monodromy-critical iff exactly 2)")
        print(f"  R6 inner bitangent circle radius          : {b['inner_circle_radius']:.6f}"
              f"  (radius mismatch {b['inner_radius_mismatch']:.2e})")
        print(f"     outer bitangent circle radius          : {b['outer_circle_radius']:.6f}"
              f"  (radius mismatch {b['outer_radius_mismatch']:.2e})")
        print(f"     minima/maxima of the distance function : "
              f"{b['n_minima']}/{b['n_maxima']}")
        print(f"  R1 sampling  : " + ", ".join(
            f"N={x['N']}->{x['order']}(drift {x['drift']:.1e})"
            if x['drift'] is not None else f"N={x['N']}->None"
            for x in rep["R1_sampling"]))
        print(f"  R2 loop steps: " + ", ".join(
            f"{x['steps']}->{x['order']}" for x in rep["R2_loop_steps"]))
        print(f"  R3 smoothing : " + ", ".join(
            f"{x['smooth_frac']}->{x['order']}(drift {x['drift']:.1e})"
            if x['drift'] is not None else f"{x['smooth_frac']}->None"
            for x in rep["R3_smoothing"]))
        print(f"  R4 noise     : " + ", ".join(
            f"{x['noise_km']}km:{x['survived']}/5" for x in rep["R4_noise"]))
        print(f"  R5 shrink    : " + ", ".join(
            f"{x['radius']:.1e}->{x['order']}" for x in rep["R5_shrink"]))
        print(f"  R5 displace  : " + ", ".join(
            str(x["order"]) for x in rep["R5_displace"]))
    with open("out/hero_reports.json", "w") as f:
        json.dump(reports, f, indent=1, default=float)


if __name__ == "__main__":
    os.makedirs("out", exist_ok=True)
    main()
