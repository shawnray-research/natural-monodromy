"""
Fetch a (gantry angle, detector row) patch of the measured cone-beam walnut scan.

The archive is 4.25 GB of full 2368 x 2240 detector panels. Only a band of rows
is needed, and the TIFF pixel data is contiguous from a known offset, so each
entry's deflate stream is read only as far as the last wanted row and then
abandoned. That turns the download into a few hundred megabytes.

Output: A[theta, v, u], the measured attenuation -log(I / I_air), with I_air
taken per row from the detector margins, which removes the row-dependent
detector gain. No flat field is distributed with this scan, and a per-row air
normalization is the honest substitute: it cannot invent structure along u,
which is the direction persistence is computed in.
"""
import sys, os, io, time
import numpy as np
sys.path.insert(0, "/private/tmp/claude-501/-Users-shawn-Documents-Research-ISTA-monodromy/ffd0f93b-6d80-4af2-a6fa-78520735c625/scratchpad")
from zget import opened

NU, NV = 2240, 2368
DATA_OFFSET = 10256          # first strip offset, strips are contiguous
V0, V1, VSTEP = 600, 1100, 2
U0, U1 = 300, 1900
NPROJ = 121                  # 0.0 to 60.0 degrees at 0.5 degree steps

def main():
    z = opened()
    need = DATA_OFFSET + V1 * NU * 2
    vs = np.arange(V0, V1, VSTEP)
    out = np.zeros((NPROJ, len(vs), U1 - U0), dtype=np.float32)
    t0 = time.time()
    for k in range(NPROJ):
        name = f"20201111_walnut_raw_data/20201111_walnut_{k+1:04d}.tif"
        with z.open(name) as f:
            buf = f.read(need)
        px = np.frombuffer(buf, dtype="<u2", count=V1 * NU,
                           offset=DATA_OFFSET).reshape(V1, NU)
        band = px[vs].astype(np.float32)
        air = np.median(np.concatenate([band[:, :60], band[:, -60:]], axis=1),
                        axis=1, keepdims=True)
        out[k] = -np.log(np.maximum(band[:, U0:U1], 1.0) / air)
        if k % 10 == 0:
            el = time.time() - t0
            print(f"  {k+1}/{NPROJ}  {el:6.1f}s", flush=True)
    os.makedirs("data/ct", exist_ok=True)
    np.save("data/ct/walnut_patch.npy", out)
    np.save("data/ct/walnut_patch_v.npy", vs)
    print("saved", out.shape, f"{out.nbytes/1e6:.0f} MB, {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
