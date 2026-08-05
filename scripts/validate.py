"""Validation of the monodromy detector on known ground-truth cases."""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.core import (radial_diagram, circle_loop, monodromy_permutation,
                       perm_order, cycle_type, resample_closed, smooth_closed)

np.set_printoptions(precision=4, suppress=True)


def ellipse(n=2000, a=2.0, b=1.0):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.column_stack([a * np.cos(t), b * np.sin(t)])


print("=" * 72)
print("TEST 0 : timing + sanity of extended persistence on the ellipse")
print("=" * 72)
M = ellipse(2000)
t0 = time.time()
for _ in range(200):
    r = radial_diagram(M, np.array([0.13, 0.07]))
dt = (time.time() - t0) / 200
print(f"  one diagram (N=2000): {dt*1e3:.2f} ms")
print(f"  #minima={r['n_crit']}, diagram points={len(r['pairs'])}, kinds={r['kinds']}")
print(f"  points (b,d) = \n{r['points']}")

print()
print("=" * 72)
print("TEST 1 : ellipse center. Theory says the center IS an A_1^2/A_1^2 point")
print("         (inscribed circle r=b bitangent at 2 minima, circumscribed")
print("         circle r=a bitangent at 2 maxima) but with only 2 births and")
print("         2 deaths the pairing flips at ALL FOUR quadrant walls, so the")
print("         net monodromy must be TRIVIAL.")
print("=" * 72)
for rad in [0.02, 0.05, 0.1, 0.2]:
    loop = circle_loop([0.0, 0.0], rad, 256)
    perm, info = monodromy_permutation(M, loop)
    print(f"  radius={rad:5.3f}  perm={perm}  order={perm_order(perm)}  info={ {k:v for k,v in info.items() if k!='stack'} }")

print()
print("=" * 72)
print("TEST 2 : count the pairing flips around the ellipse center explicitly")
print("=" * 72)
loop = circle_loop([0.0, 0.0], 0.05, 720)
prev = None
flips = 0
for k, p in enumerate(loop):
    d = radial_diagram(M, p)
    pr = tuple(sorted((i, j) for (i, j, _) in d["pairs"]))
    if prev is not None and pr != prev:
        flips += 1
        print(f"    t={2*np.pi*k/len(loop):6.3f} rad :  pairing {prev} -> {pr}")
    prev = pr
print(f"  total pairing changes around the loop: {flips}  (theory: 4 -> identity)")
