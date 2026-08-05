"""
The tests that decide whether the two localized generators are real.

Everything so far has been computed on a field that is BILINEARLY interpolated
between the measured samples in (theta, v). Bilinear interpolation is only C^0:
its derivatives jump across cell boundaries, so it can in principle manufacture
codimension-2 coincidences that the underlying measurement does not have. If a
generator is an artifact of the interpolant it must move or vanish when the
interpolant is replaced.

  R1 interpolant   Catmull-Rom bicubic in (theta, v) instead of bilinear. This
                   is a genuinely different C^1 reconstruction of the same
                   samples. A real generator survives and stays put.
  R2 scale         the smoothing sigma and the noise floor are the only two
                   numbers in the pipeline not fixed by the instrument.
  R3 subsample     drop every second gantry angle and re-localize, so the
                   generator is recovered from half the measurements.
"""
import sys, os, json
import numpy as np
from scipy.ndimage import gaussian_filter1d

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.general import monodromy_abs, loop_points, localize_abs, wall_certificate
from mono.spectra import compactify
from hunt_ct import Walnut, DEG_PER_PROJ, VSTEP


class WalnutCubic(Walnut):
    """Same measured samples, Catmull-Rom bicubic in (theta, v)."""

    @staticmethod
    def _w(t):
        return np.array([-0.5*t**3 + t**2 - 0.5*t,
                         1.5*t**3 - 2.5*t**2 + 1.0,
                         -1.5*t**3 + 2.0*t**2 + 0.5*t,
                         0.5*t**3 - 0.5*t**2])

    def profile(self, k, j):
        k = float(np.clip(k, 1, self.nth - 3 - 1e-9))
        j = float(np.clip(j, 1, self.nv - 3 - 1e-9))
        k0, j0 = int(k), int(j)
        wk, wj = self._w(k - k0), self._w(j - j0)
        blk = self.S[k0-1:k0+3, j0-1:j0+3, :]
        return np.einsum("a,b,abu->u", wk, wj, blk)


class WalnutSub(Walnut):
    """Every second gantry angle dropped; k is still in original index units."""

    def __init__(self, sigma_pix, offset=0):
        super().__init__(sigma_pix)
        self.S = self.S[offset::2]
        self.nth = self.S.shape[0]
        self.offset = offset

    def profile(self, k, j):
        return super().profile((k - self.offset) / 2.0, j)


PTS = [(46.086646, 12.152450), (49.617361, 63.760125)]


def order_at(W, u, v, floor, rad, steps=160):
    r = monodromy_abs(W.field, loop_points(u, v, rad, steps), floor)
    return (r["order"] if r["ok"] else None), r.get("pairing_changes")


def main():
    out = {}
    print("=" * 74)
    print("R1  DIFFERENT INTERPOLANT: Catmull-Rom bicubic instead of bilinear")
    print("=" * 74)
    Wb, Wc = Walnut(15.0), WalnutCubic(15.0)
    floor = 3.0 * Wb.noise()
    r1 = []
    for (u, v) in PTS:
        print(f"\n  bilinear p* = ({u:.6f}, {v:.6f})   theta {u*DEG_PER_PROJ:.4f} deg")
        L = localize_abs(Wc.field, u - 0.5, u + 0.5, v - 0.5, v + 0.5, floor)
        if L is None:
            print("    bicubic: NO generator in the same neighborhood -> artifact")
            r1.append({"u": u, "v": v, "cubic": None}); continue
        d = np.hypot(L["u"] - u, L["v"] - v)
        print(f"    bicubic p* = ({L['u']:.6f}, {L['v']:.6f})  order {L['order']}  "
              f"changes {L['changes']}  box {L['box']:.1e}")
        print(f"    moved by {d:.2e} index units "
              f"= {d*DEG_PER_PROJ*3600:.2f} arcsec in theta-equivalent")
        os_, ch = [], []
        for f in (1.0, 0.25, 0.05, 0.01, 1e-3):
            o, c = order_at(Wc, L["u"], L["v"], floor, 0.8 * f)
            os_.append(o); ch.append(c)
        print(f"    bicubic shrink orders {os_}  changes {ch}")
        w = wall_certificate(Wc.field, L["u"], L["v"])
        mn, mx = np.sort(w["min_values"]), np.sort(w["max_values"])
        print(f"    bicubic wall: min gap {np.diff(mn).min():.2e}, "
              f"max gap {np.diff(mx).min():.2e}")
        r1.append({"u": u, "v": v, "cubic_u": L["u"], "cubic_v": L["v"],
                   "moved": float(d), "orders": os_, "changes": ch})
    out["R1"] = r1

    print("\n" + "=" * 74)
    print("R2  SCALE: smoothing sigma and noise floor")
    print("=" * 74)
    r2 = []
    for (u, v) in PTS:
        print(f"\n  p* near theta {u*DEG_PER_PROJ:.4f} deg")
        row = []
        for sig in (10.0, 12.0, 15.0, 18.0, 22.0):
            W = Walnut(sig)
            fl_base = W.noise()
            line = []
            for ns in (2.0, 3.0, 4.0, 5.0):
                fl = ns * fl_base
                L = localize_abs(W.field, u - 0.5, u + 0.5, v - 0.5, v + 0.5, fl)
                if L is None:
                    line.append("-")
                else:
                    o, _ = order_at(W, L["u"], L["v"], fl, 0.02)
                    line.append(f"{o}" if o else "-")
            print(f"    sigma {sig:5.1f} px:  floors 2,3,4,5 sigma -> {line}")
            row.append({"sigma": sig, "orders": line})
        r2.append(row)
    out["R2"] = r2

    print("\n" + "=" * 74)
    print("R3  SUBSAMPLE: half the gantry angles")
    print("=" * 74)
    r3 = []
    for off in (0, 1):
        Ws = WalnutSub(15.0, off)
        fl = 3.0 * Ws.noise()
        print(f"\n  keeping projections {off}, {off+2}, {off+4}, ... "
              f"(1.0 degree steps)")
        for (u, v) in PTS:
            L = localize_abs(Ws.field, u - 0.6, u + 0.6, v - 0.6, v + 0.6, fl)
            if L is None:
                print(f"    theta {u*DEG_PER_PROJ:7.4f} deg: not recovered")
                r3.append(None); continue
            d = np.hypot(L["u"] - u, L["v"] - v)
            o, c = order_at(Ws, L["u"], L["v"], fl, 0.05)
            print(f"    theta {u*DEG_PER_PROJ:7.4f} deg: recovered at "
                  f"({L['u']:.5f},{L['v']:.5f}), moved {d:.2e}, order {o}, changes {c}")
            r3.append({"off": off, "u": L["u"], "v": L["v"], "moved": float(d),
                       "order": o})
    out["R3"] = r3

    json.dump(out, open("out/ct_robust.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
