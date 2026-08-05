"""
The two questions left open: what the knot means, and what survives at N > 3.

KNOT. The braid sigma_1 sigma_2^-1 is the classical pseudo-Anosov 3-braid. Under
the standard surjection from the 3-strand braid group modulo its center onto
SL(2,Z), sigma_1 -> [[1,1],[0,1]] and sigma_2 -> [[1,0],[-1,1]], so

    sigma_1 sigma_2^-1  ->  [[2,1],[1,1]],

the Anosov matrix whose mapping torus on the once-punctured torus is the
figure-eight knot complement. That is why the figure-eight ORBIT is celebrated in
braid terms: its trajectory braid is that map. The vineyard braid is its square,
so it has the same stable foliation and its dilatation is squared, which means
its topological entropy is exactly twice that of the trajectories. The knot is
then not decorative: it is the closure of the square of the orbit's own
pseudo-Anosov braid.

N > 3. The wall between bodies i and j is

    h_i - h_j = sum over k not i,j of [ K(r_ik) - K(r_jk) ] = 0,

kernel-free only when the sum has one term, which is N = 3. But as the bandwidth
falls, K decays so fast that the sum is dominated by each body's NEAREST
neighbor, so the wall tends to

    d_i = d_j,      d_i = distance from body i to its nearest neighbor,

which is kernel-free for every N. So the clean statement survives at N > 3 as a
small-bandwidth limit rather than an identity, and at N = 3 the limit is already
exact.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def sl2(word):
    S1 = np.array([[1, 1], [0, 1]], dtype=float)
    S2 = np.array([[1, 0], [-1, 1]], dtype=float)
    M = np.eye(2)
    for w in word:
        A = S1 if abs(w) == 1 else S2
        M = M @ (np.linalg.inv(A) if w < 0 else A)
    return M


def dil(M):
    return max(abs(np.linalg.eigvals(M)))


def main():
    print(__doc__)
    print("=" * 76)
    print("1. The knot, and what the braid actually is")
    print("=" * 76)
    traj = [1, -2]                      # sigma_1 sigma_2^-1
    vine = [-2, 1, -2, 1]               # (sigma_2^-1 sigma_1)^2
    Mt, Mv = sl2(traj), sl2(vine)
    print(f"  trajectory braid  s1 s2^-1        -> SL(2,Z) matrix\n{Mt.astype(int)}")
    print(f"  trace {np.trace(Mt):.0f}, dilatation {dil(Mt):.6f}, "
          f"golden ratio squared {( (1+5**0.5)/2 )**2:.6f}")
    print(f"\n  vineyard braid    (s2^-1 s1)^2   -> SL(2,Z) matrix\n{Mv.astype(int)}")
    print(f"  trace {np.trace(Mv):.0f}, dilatation {dil(Mv):.6f}")
    print(f"\n  is the vineyard dilatation the square of the trajectory one? "
          f"{abs(dil(Mv) - dil(Mt)**2):.2e}")
    print(f"  topological entropy: trajectories {np.log(dil(Mt)):.6f}, "
          f"vineyard {np.log(dil(Mv)):.6f}, ratio {np.log(dil(Mv))/np.log(dil(Mt)):.6f}")
    try:
        import snappy
        M = snappy.Manifold("m004")
        print(f"\n  m004, the once-punctured-torus bundle with monodromy [[2,1],[1,1]]:")
        print(f"    volume {float(M.volume()):.6f}, identified as {M.identify()[:2]}")
    except Exception as e:
        print(f"  snappy: {e}")
    print("\n  So the same pseudo-Anosov map appears twice over: as the monodromy of")
    print("  the figure-eight knot fibration, and as the trajectory braid of the")
    print("  figure-eight orbit. The vineyard braid is its square, with exactly")
    print("  twice the entropy, and closes to that same knot.")

    print("\n" + "=" * 76)
    print("2. N > 3: the small-bandwidth limit")
    print("=" * 76)
    rng = np.random.default_rng(5)
    ts = np.linspace(0, 1, 6000)

    def walls(P_of_s, N, K):
        out = []
        for a in range(N):
            for b in range(a + 1, N):
                g = []
                for s in ts:
                    P = P_of_s(s)
                    r = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)
                    g.append(sum(K(r[a, k]) - K(r[b, k])
                                 for k in range(N) if k not in (a, b)))
                g = np.array(g)
                out += [ts[i] for i in range(len(ts)-1)
                        if np.sign(g[i]) != np.sign(g[i+1])]
        return np.sort(np.array(out))

    def nn_walls(P_of_s, N):
        """Walls of the kernel-free limit: equal nearest-neighbor distance."""
        out = []
        for a in range(N):
            for b in range(a + 1, N):
                g = []
                for s in ts:
                    P = P_of_s(s)
                    r = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)
                    r += np.eye(N) * 1e9
                    g.append(r[a].min() - r[b].min())
                g = np.array(g)
                out += [ts[i] for i in range(len(ts)-1)
                        if np.sign(g[i]) != np.sign(g[i+1])]
        return np.sort(np.array(out))

    for N in (3, 4, 5):
        A = rng.normal(0, 1, (N, 2)); B = rng.normal(0, 1, (N, 2))
        P_of_s = lambda s: A*np.cos(2*np.pi*s) + B*np.sin(2*np.pi*s)
        ref = nn_walls(P_of_s, N)
        print(f"\n  N = {N}: kernel-free limit gives {len(ref)} walls")
        print(f"    {'sigma':>7} {'walls':>7} {'max distance to the limit walls':>34}")
        for sig in (0.60, 0.35, 0.20, 0.12, 0.07):
            w = walls(P_of_s, N, lambda r: np.exp(-r**2/(2*sig**2)))
            if len(w) == len(ref):
                print(f"    {sig:7.2f} {len(w):7d} {np.abs(w-ref).max():34.3e}")
            else:
                d = max(min(abs(x - ref)) for x in w) if len(w) else np.nan
                print(f"    {sig:7.2f} {len(w):7d} {d:34.3e}  (count differs)")


if __name__ == "__main__":
    main()
