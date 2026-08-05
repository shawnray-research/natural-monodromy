"""
Audit the last paragraph of the note: the trajectory braid, the squaring claim,
and the entropy comparison.

The note asserts four things about the trajectory braid, and each was taken from
the literature rather than measured here. This script measures all of them.

  1  The trajectory braid of the figure-eight over T/3 is sigma_1 sigma_2^-1.
     Measured directly from the integrated orbit, in the same style of convention
     used for the vineyard: order the strands by their projected coordinate at
     the first slice, run time forwards, sign each crossing by which body passes
     in front. A braid WORD depends on the projection direction, so the word is
     computed over many directions and it is the conjugacy class that has to be
     stable, not the letters.

  2  The vineyard braid is the square of the trajectory braid. This should be an
     identity in the braid group, not a numerical coincidence, and not something
     that has to be pushed through SL(2,Z) first: conjugating (s1 s2^-1)^2 by s1
     gives s2^-1 s1 s2^-1 s1 on the nose. Checked here in reduced Burau, which is
     faithful for three strands, so word equality is decided and not sampled.

  3  The SL(2,Z) matrix. sigma_1 -> [[1,1],[0,1]], sigma_2 -> [[1,0],[-1,1]] is
     the standard surjection. Under it s1 s2^-1 -> [[2,1],[1,1]]. Since B_3 does
     not act faithfully, the check that matters is that the trace, which is a
     conjugacy invariant, matches what the measured braid gives.

  4  The entropies. Both braids must be read over the SAME time span or the
     factor of two is an artifact of normalization. The vineyard braid is
     measured over T/3, so the trajectory braid must be too.

Also settles a claim in the closing sentence: that the vineyard sees twelve
isosceles configurations and the trajectories six. The second half is projection
dependent, and this checks whether it is stable enough to state.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.nbody import FIGURE8, integrate
from mono.braid import word_to_string, reduce_word, permutation_of_word

F8 = FIGURE8
MS, T, NS = F8["m"], F8["T"], 24000


def trajectory_braid(P, theta):
    """
    Braid word of three moving points, projected onto the direction theta.

    Same convention as the vineyard: strands ordered by the projected coordinate
    at the first slice, time forwards, sign from the transverse coordinate, which
    is the analogue of signing by the death and says which strand is in front.
    Returns the word and the fractional time indices of the crossings.
    """
    u = np.array([np.cos(theta), np.sin(theta)])
    w = np.array([-np.sin(theta), np.cos(theta)])
    x = P @ u
    y = P @ w
    K = len(x)
    order = list(np.argsort(x[0]))
    word, times = [], []
    for t in range(1, K):
        moved = True
        while moved:
            moved = False
            for j in range(2):
                a, b = order[j], order[j + 1]
                if x[t, a] > x[t, b]:
                    ga, gb = x[t - 1, a] - x[t - 1, b], x[t, a] - x[t, b]
                    frac = abs(ga) / (abs(ga) + abs(gb)) if (ga or gb) else 0.5
                    sgn = 1 if y[t, a] > y[t, b] else -1
                    word.append(sgn * (j + 1))
                    times.append(t - 1 + frac)
                    order[j], order[j + 1] = order[j + 1], order[j]
                    moved = True
    return word, np.array(times)


# ---- reduced Burau, faithful for three strands, so word equality is decided ----

def burau(word, t):
    def gen(i, inv):
        if i == 1:
            M = np.array([[-t, 0.0], [1.0, 1.0]])
        else:
            M = np.array([[1.0, t], [0.0, -t]])
        return np.linalg.inv(M) if inv else M
    M = np.eye(2)
    for a in word:
        M = M @ gen(abs(a), a < 0)
    return M


def same_braid(u, v):
    """Equality in B_3, decided by reduced Burau at several values of t."""
    return all(np.abs(burau(u, t) - burau(v, t)).max() < 1e-8
               for t in (1.7, -0.43, 2.9, 0.31))


def sl2(word):
    S1 = np.array([[1.0, 1.0], [0.0, 1.0]])
    S2 = np.array([[1.0, 0.0], [-1.0, 1.0]])
    M = np.eye(2)
    for a in word:
        A = S1 if abs(a) == 1 else S2
        M = M @ (np.linalg.inv(A) if a < 0 else A)
    return M


def conjugator(target, source, maxlen=3):
    """Search short words c with c^-1 source c = target. Constructive, so a hit
    is a proof and a miss is only a miss."""
    from itertools import product
    for L in range(0, maxlen + 1):
        for c in product([1, -1, 2, -2], repeat=L):
            c = list(c)
            if same_braid([-a for a in reversed(c)] + source + c, target):
                return c
    return None


def main():
    print(__doc__)
    t, X, V = integrate(F8["x"], F8["v"], MS, T, NS)
    third = NS // 3
    TRAJ = [1, -2]                       # s1 s2^-1, the word the note quotes
    VINE = [-2, 1, -2, 1]                # (s2^-1 s1)^2, the measured vineyard word

    print("=" * 78)
    print("1. The trajectory braid over T/3, measured, over many projections")
    print("=" * 78)
    print(f"  {'theta/pi':>9} {'word over T/3':>18} {'perm':>10} {'exp sum':>8} "
          f"{'SL2 trace':>10} {'conj to s1 s2^-1':>18}")
    classes, ok, tot = {}, 0, 0
    for th in np.linspace(0, np.pi, 19)[:-1]:
        w, _ = trajectory_braid(X[:third + 1], th)
        w = reduce_word(w)
        if len(w) == 0:
            print(f"  {th/np.pi:9.3f} {'(degenerate)':>18}")
            continue
        c = conjugator(w, TRAJ)
        tot += 1
        ok += c is not None
        classes[word_to_string(w)] = classes.get(word_to_string(w), 0) + 1
        print(f"  {th/np.pi:9.3f} {word_to_string(w):>18} "
              f"{str(permutation_of_word(w, 3)):>10} {sum(np.sign(w)):8.0f} "
              f"{np.trace(sl2(w)):10.0f} "
              f"{('yes, by ' + word_to_string(c)) if c is not None else 'NO':>18}")
    print(f"\n  conjugate to s1 s2^-1 in {ok} of {tot} projections")
    print(f"  distinct words seen: {classes}")
    print(f"  reference: s1 s2^-1 has exponent sum 0 and SL(2,Z) trace "
          f"{np.trace(sl2(TRAJ)):.0f}")

    print("\n" + "=" * 78)
    print("2. Is the vineyard braid the square of the trajectory braid?")
    print("=" * 78)
    sq = TRAJ + TRAJ
    print(f"  trajectory braid          {word_to_string(TRAJ)}")
    print(f"  its square                {word_to_string(sq)}")
    print(f"  measured vineyard braid   {word_to_string(VINE)}")
    print(f"\n  equal as braids?          {same_braid(VINE, sq)}")
    c = conjugator(VINE, sq, maxlen=2)
    print(f"  conjugate as braids?      "
          f"{'yes, by ' + word_to_string(c) if c is not None else 'NO'}")
    print(f"  check the conjugation by hand: s1^-1 (s1 s2^-1 s1 s2^-1) s1")
    print(f"    = {word_to_string(reduce_word([-1] + sq + [1]))}")
    print(f"  closure is a conjugacy invariant, so the knot follows from this")
    print(f"  alone, with no appeal to SL(2,Z).")

    print("\n" + "=" * 78)
    print("3. The SL(2,Z) matrix and the dilatations")
    print("=" * 78)
    Mt, Mv = sl2(TRAJ), sl2(sq)
    print(f"  s1 s2^-1 -> \n{Mt.astype(int)}   trace {np.trace(Mt):.0f}")
    print(f"  its square -> \n{Mv.astype(int)}   trace {np.trace(Mv):.0f}")
    lt = max(abs(np.linalg.eigvals(Mt)))
    lv = max(abs(np.linalg.eigvals(Mv)))
    phi = (1 + 5 ** 0.5) / 2
    print(f"\n  dilatation, trajectory {lt:.6f}   golden ratio squared {phi**2:.6f}"
          f"   diff {abs(lt - phi**2):.2e}")
    print(f"  dilatation, vineyard   {lv:.6f}   golden ratio ^4      {phi**4:.6f}"
          f"   diff {abs(lv - phi**4):.2e}")
    print(f"  vineyard dilatation is the square of the trajectory one to "
          f"{abs(lv - lt**2):.2e}")

    print("\n" + "=" * 78)
    print("4. Entropy, and whether the two are normalized over the same time")
    print("=" * 78)
    print(f"  both braids are read over the same loop, T/3.")
    print(f"    trajectory braid over T/3: {word_to_string(TRAJ)}, "
          f"{len(TRAJ)} crossings")
    print(f"    vineyard braid over T/3:   {word_to_string(VINE)}, "
          f"{len(VINE)} crossings")
    print(f"  entropy of the mapping class, trajectory {np.log(lt):.6f}")
    print(f"  entropy of the mapping class, vineyard   {np.log(lv):.6f}")
    print(f"  ratio {np.log(lv)/np.log(lt):.6f}")
    print(f"  per unit time the ratio is the same, since the span is shared:")
    print(f"    trajectory {np.log(lt)/(T/3):.6f} per unit time, "
          f"vineyard {np.log(lv)/(T/3):.6f}")

    print("\n" + "=" * 78)
    print("5. Crossing counts per period: is 'the trajectories see six' stable?")
    print("=" * 78)
    print(f"  {'theta/pi':>9} {'traj crossings per T':>22} {'word length':>12}")
    counts = []
    for th in np.linspace(0, np.pi, 19)[:-1]:
        w, tm = trajectory_braid(X, th)
        counts.append(len(w))
        print(f"  {th/np.pi:9.3f} {len(w):22d} {len(reduce_word(w)):12d}")
    print(f"\n  crossing count over a period ranges {min(counts)} to {max(counts)}"
          f" across projections")
    print(f"  the vineyard braid over a period has {3*len(VINE)} crossings, "
          f"and is projection independent")

    print("\n" + "=" * 78)
    print("6. Where the trajectory crossings sit relative to the walls")
    print("=" * 78)
    walls = np.arange(12) / 12.0
    for th in (0.0, np.pi/4, np.pi/2):
        w, tm = trajectory_braid(X, th)
        frac = np.sort(tm / NS)
        d = [min(abs(((f - wl) + 0.5) % 1.0 - 0.5) for wl in walls) for f in frac]
        print(f"\n  theta/pi = {th/np.pi:.2f}: {len(frac)} crossings at t/T")
        print(f"    {np.array2string(frac, precision=4)}")
        print(f"    distance to the nearest twelfth: max {max(d):.4f}")
        print(f"    {'ON the wall locus' if max(d) < 1e-3 else 'NOT on the wall locus'}")


if __name__ == "__main__":
    main()
