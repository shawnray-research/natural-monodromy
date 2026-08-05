"""
Vineyard monodromy in a rotating detonation combustor.

The theorem says a rotating wave carried through a stationary envelope has
monodromy of order n over T/n, non-degenerate when the envelope is not
2*pi/n-periodic. An RDC is that object physically: n detonation fronts circulate
in an annulus at a fixed speed, and the injector pattern and feed asymmetry are
fixed in the laboratory frame and modulate the fronts as they pass.

  ring       the annulus, a genuine S^1, no compactification
  field      measured luminosity, raw high-speed video, no threshold mask
  loop       T/n, the time for the pattern to advance one front spacing. Forced
             by the device
  envelope   the injector pattern, stationary in the laboratory frame

Data: Bohon et al., Zenodo 18886925, CC BY.

The tracker is transport, following crests by continuity in azimuth. Assignment
on (birth, death) is the wrong instrument here and returns the identity on
exactly this class.
"""
import sys, os, json, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.core import extended_persistence_circle, perm_order, cycle_type
from mono.general import monodromy_transport, monodromy_abs, rank_pairing_key


def load_rdc(path):
    """
    Raw RDC record -> (azimuthal field, times).

    Each frame carries `flat`, the annulus already unwrapped to
    (radius x azimuth), together with the ring center, the inner and outer radii
    and an absolute timestamp. The field used here is the mean over the radial
    band of that unwrapped image, which is luminosity as a function of azimuth.
    """
    import scipy.io as sio
    d = sio.loadmat(path)
    V = d["V"][0, 0]
    IM = np.squeeze(V["IM"])
    n = len(IM)
    F0 = np.squeeze(IM[0]["flat"])
    A = np.empty((n, F0.shape[1]), dtype=np.float32)
    t = np.empty(n)
    for i in range(n):
        A[i] = np.squeeze(IM[i]["flat"]).astype(np.float32).mean(axis=0)
        t[i] = float(np.ravel(IM[i]["t_rel"])[0])
    r = np.ravel(IM[0]["r"])
    return A, t, F0.shape, r


def wave_count_and_period(A):
    """Dominant azimuthal wavenumber n, and the frames for one front spacing."""
    S = np.abs(np.fft.rfft(A - A.mean(axis=1, keepdims=True), axis=1)).mean(axis=0)
    n = int(np.argmax(S[1:20]) + 1)
    ph = np.unwrap(np.angle(np.fft.rfft(A, axis=1)[:, n]))
    slope = np.polyfit(np.arange(len(ph)), ph, 1)[0]
    if abs(slope) < 1e-9:
        return n, None, S
    return n, abs(2 * np.pi / slope), S


def diag(f, floor):
    r = extended_persistence_circle(np.asarray(f, dtype=float))
    if r is None:
        return None
    p = np.array([[f[r["mins"][i]], f[r["maxs"][j]]] for (i, j, _) in r["pairs"]])
    keep = (p[:, 1] - p[:, 0]) >= floor
    return p[keep], np.array([int(r["maxs"][j]) for (i, j, _) in r["pairs"]])[keep]


def transport_order(A, i0, span, floor, nth, steps=None):
    """Track crests by continuity in azimuth around the loop of `span` frames."""
    idx = np.linspace(i0, i0 + span, steps or int(span) + 1)
    prof, pos = [], []
    for t in idx:
        a = int(np.floor(t)); b = min(a + 1, len(A) - 1); w = t - a
        f = (1 - w) * A[a] + w * A[b]
        d = diag(f, floor)
        if d is None or len(d[0]) < 2:
            return None
        o = np.argsort(d[1])
        prof.append(d[0][o]); pos.append(d[1][o])
    k = {len(p) for p in pos}
    if len(k) != 1:
        return {"ok": False, "counts": sorted(k)}
    k = k.pop()
    from scipy.optimize import linear_sum_assignment
    cur = list(range(k)); worst = 0
    for t in range(1, len(pos) + 1):
        a, b = pos[t - 1], pos[t % len(pos)]
        D = np.abs(a[:, None] - b[None, :]); D = np.minimum(D, nth - D)
        ri, ci = linear_sum_assignment(D)
        worst = max(worst, D[ri, ci].max())
        m = {int(i): int(j) for i, j in zip(ri, ci)}
        cur = [m[c] for c in cur]
    seps = [np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)[np.triu_indices(len(p), 1)].min()
            for p in prof if len(p) > 1]
    return {"ok": True, "order": perm_order(cur), "cycles": cycle_type(cur),
            "n": k, "max_jump": float(worst), "min_sep": float(min(seps)),
            "closure": float(np.abs(prof[0] - prof[-1]).max())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mat", required=True)
    ap.add_argument("--nth", type=int, default=720)
    ap.add_argument("--nsig", type=float, default=5.0)
    a = ap.parse_args()

    A, t, shp, radii = load_rdc(a.mat)
    dt = float(np.median(np.diff(t)))
    print(f"{A.shape[0]} frames at {1/dt/1e3:.1f} kHz, unwrapped annulus {shp[0]} radial "
          f"x {shp[1]} azimuthal, ring radii {radii[0]:.0f} and {radii[1]:.0f} px")
    a.nth = A.shape[1]
    print(f"azimuthal field {A.shape}  (frames x azimuth)")

    n, span, S = wave_count_and_period(A)
    print(f"\ndominant azimuthal wavenumber n = {n}")
    print(f"  harmonic amplitudes 1..8: "
          f"{np.array2string(S[1:9]/S[1:9].max(), precision=2)}")
    if span is None:
        print("  no steady rotation detected"); return
    print(f"  frames for the pattern to advance one front spacing (T/n): {span:.2f}")

    # noise floor from the frame-to-frame difference in a quiet azimuth band
    hf = A[1:] - A[:-1]
    nz = float(np.median(np.abs(hf)) * 1.4826)
    floor = a.nsig * nz
    print(f"  noise (robust sigma of frame differences) {nz:.4f}, "
          f"floor {a.nsig} sigma = {floor:.4f}")
    dA = float(np.median(np.abs(A[1:] - A[:-1])))
    print(f"  D10: median |change| per frame {dA:.4f} = {dA/nz:.1f}x noise")

    print("\nmonodromy over the forced loop T/n, transport tracking:")
    res = []
    for t0 in np.linspace(0, len(A) - 2 * span - 2, 6):
        r = transport_order(A, t0, span, floor, a.nth, steps=max(60, int(span * 4)))
        if r is None or not r.get("ok"):
            print(f"  start frame {t0:7.0f}: {r}"); continue
        print(f"  start frame {t0:7.0f}: order {r['order']}, cycles "
              f"{[c for c in r['cycles'] if len(c)>1]}, n={r['n']}, "
              f"min sep {r['min_sep']:.4f}, closure {r['closure']:.4f}, "
              f"jump {r['max_jump']:.0f}")
        res.append(r)

    print("\ncontrol, the FULL rotation n*(T/n), which must give the identity:")
    for t0 in np.linspace(0, len(A) - n * span - 2, 3):
        r = transport_order(A, t0, n * span, floor, a.nth, steps=max(120, int(n*span*3)))
        if r and r.get("ok"):
            print(f"  start frame {t0:7.0f}: order {r['order']}, n={r['n']}")
        else:
            print(f"  start frame {t0:7.0f}: {r}")


if __name__ == "__main__":
    main()
