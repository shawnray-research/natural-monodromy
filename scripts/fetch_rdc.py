"""
Fetch one raw high-speed video record from the rotating detonation combustor
dataset of Bohon et al., Zenodo 18886925, CC BY.

Raw rather than processed: the deposited processing pipeline applies a
luminosity threshold mask before its polar transform, and a masked field is not
the continuous scalar field persistence needs. The polar transform is done here
instead, from the raw frames.
"""
import io, os, sys, zipfile, urllib.request

URL = "https://zenodo.org/records/18886925/files/data.zip?download=1"
WANT = "CE2029_HSV_CB076_A0199_Phi1005_20180913_1_5000.mat"

class HttpFile(io.RawIOBase):
    def __init__(self, url):
        self.url = url; self.pos = 0
        with urllib.request.urlopen(urllib.request.Request(url, method="HEAD")) as f:
            self.size = int(f.headers["Content-Length"])
    def seekable(self): return True
    def readable(self): return True
    def seek(self, o, w=0):
        self.pos = o if w == 0 else (self.pos + o if w == 1 else self.size + o)
        return self.pos
    def tell(self): return self.pos
    def readinto(self, b):
        n = len(b)
        if self.pos >= self.size or n == 0:
            return 0
        hi = min(self.pos + n - 1, self.size - 1)
        r = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{hi}"})
        with urllib.request.urlopen(r) as f:
            d = f.read()
        b[:len(d)] = d; self.pos += len(d); return len(d)

def main():
    z = zipfile.ZipFile(io.BufferedReader(HttpFile(URL), 1 << 22))
    name = [n for n in z.namelist() if n.endswith(WANT)][0]
    out = os.path.join("data/rdc", WANT)
    print(f"fetching {WANT}  ({z.getinfo(name).file_size/1e6:.0f} MB)")
    with z.open(name) as fin, open(out, "wb") as fo:
        done = 0
        while True:
            chunk = fin.read(1 << 22)
            if not chunk:
                break
            fo.write(chunk); done += len(chunk)
            print(f"  {done/1e6:7.1f} MB", flush=True)
    print("saved", out)

if __name__ == "__main__":
    main()
