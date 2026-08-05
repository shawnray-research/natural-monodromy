"""
Re-run the measured-data scans with TRANSPORT tracking.

Every negative in notes/RECORD.md, N4 to N9, was obtained with a tracker that has since
been shown to return the identity on systems whose true monodromy is an n-cycle.
Either outcome here is a result: if a measured field lights up, the obstruction
was the tracker; if it stays dark with a tracker that provably works, D10 becomes
a finding rather than a hypothesis.
"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.general import monodromy_transport, monodromy_abs, loop_points
from hunt_ct import Walnut, DEG_PER_PROJ, VSTEP
from hunt_ct2 import Coarse, step_snr

def scan(W, floor, rad, grid, steps=160, label=""):
    ks = np.arange(rad+1, W.nth-rad-1, grid)
    js = np.arange(rad+1, W.nv-rad-1, grid)
    tr_hits=[]; as_hits=[]; tr_ok=0; as_ok=0; tot=0
    for j in js:
        for k in ks:
            P = loop_points(k, j, rad, steps)
            rt = monodromy_transport(W.field, P, floor)
            ra = monodromy_abs(W.field, P, floor)
            tot += 1
            if rt["ok"]:
                tr_ok += 1
                if rt["order"] and rt["order"] > 1:
                    tr_hits.append({"k":float(k),"j":float(j),"order":int(rt["order"]),
                                    "cycles":[list(c) for c in rt["cycles"]],
                                    "n":rt["n"],"jump":rt["max_jump"]})
            if ra["ok"]:
                as_ok += 1
                if ra["order"] and ra["order"] > 1:
                    as_hits.append({"k":float(k),"j":float(j),"order":int(ra["order"])})
    print(f"\n{label}")
    print(f"  loops scanned {tot}")
    print(f"  TRANSPORT : well-defined {tr_ok}, nontrivial {len(tr_hits)}")
    print(f"  assignment: well-defined {as_ok}, nontrivial {len(as_hits)}")
    if tr_hits:
        from collections import Counter
        c = Counter(h["order"] for h in tr_hits)
        print(f"  transport order histogram: {dict(sorted(c.items()))}")
        for h in tr_hits[:12]:
            print(f"     theta {W.theta_deg(h['k']) if hasattr(W,'theta_deg') else h['k']*DEG_PER_PROJ:7.2f} "
                  f"order {h['order']} n {h['n']} cycles {h['cycles']} jump {h['jump']:.1f}")
    return tr_hits, as_hits

print("WALNUT CONE-BEAM SCAN, transport vs assignment")
print("="*70)
W1 = Coarse(15.0, 4, 8)
nz, sk, sj = step_snr(W1)
print(f"signal-dominated grid (D10 satisfied): steps {sk:.1f}x / {sj:.1f}x noise")
t1,a1 = scan(W1, 3.0*nz, 0.45, 0.5, 160, "grid stride 4/8, radius 0.45")

W2 = Walnut(15.0)
fl2 = 3.0*W2.noise()
t2,a2 = scan(W2, fl2, 0.8, 2.0, 120, "fine grid (D10 NOT satisfied), radius 0.8")

json.dump({"coarse_transport":t1,"coarse_assign":a1,
           "fine_transport":t2,"fine_assign":a2},
          open("out/rerun_transport.json","w"), indent=1)
