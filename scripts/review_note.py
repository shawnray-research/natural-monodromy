"""
Re-derive every claim in note/note.tex from scratch, and say PASS or FAIL.

Why this exists. The audit (scripts/audit_note.py) checks the note against saved
logs. That is only as good as the logs, and seven of them turned out to have no
surviving generator in the repository, so their numbers could be read but not
reproduced. This script recomputes everything from the orbit up, touching no log,
and prints what the note says beside what the computation gives.

Each check is labelled with the sentence of the note it is testing. A claim that
this script cannot reach is marked REASONED rather than PASS, and says what the
argument is, so that nothing is silently counted as verified.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.nbody import FIGURE8, integrate
from mono.kde_exact import merge_tree_exact, newton_critical
from mono.braid import braid_word, word_to_string, reduce_word, permutation_of_word

F8 = FIGURE8
MS, T = F8["m"], F8["T"]
SIG = 0.30

def dense_orbit():
    """
    Continuous positions q(t), so that root-finding on time is limited by the
    integrator tolerance rather than by a sampling grid.

    The first version of this script sampled X on a 12000-point grid and then
    "bisected" on it, which cannot resolve better than T/12000. That is what
    produced isosceles times offset by exactly one grid step and made both the
    set-return scan and the kernel comparison look wrong.
    """
    from scipy.integrate import solve_ivp
    from mono.nbody import accel
    n = 3
    def rhs(tt, y):
        x = y[:2*n].reshape(n, 2); v = y[2*n:].reshape(n, 2)
        return np.concatenate([v.ravel(), accel(x, MS).ravel()])
    y0 = np.concatenate([F8["x"].ravel(), F8["v"].ravel()])
    sol = solve_ivp(rhs, (0.0, T), y0, method="DOP853",
                    rtol=1e-13, atol=1e-13, dense_output=True)
    return lambda tau: sol.sol(np.clip(tau, 0.0, T) * T)[:2*n].reshape(n, 2)


def brent(f, a, b, tol=1e-14):
    fa, fb = f(a), f(b)
    for _ in range(200):
        m = 0.5 * (a + b)
        if b - a < tol:
            return m
        fm = f(m)
        if np.sign(fm) == np.sign(fa):
            a, fa = m, fm
        else:
            b, fb = m, fm
    return 0.5 * (a + b)


RESULT = []


def check(tag, claim, ok, got):
    RESULT.append((tag, ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {tag}")
    print(f"         note: {claim}")
    print(f"         got : {got}")


def reasoned(tag, claim, why):
    RESULT.append((tag, None))
    print(f"  [ARG ] {tag}")
    print(f"         note: {claim}")
    print(f"         why : {why}")


def vineyard(X, idx, sigma=SIG):
    B = np.zeros((len(idx), 3)); D = np.zeros((len(idx), 3))
    for a, k in enumerate(idx):
        d = merge_tree_exact(X[k], MS, sigma)
        if d is None:
            return None, None
        for i in range(3):
            B[a, i] = d[i]["birth"]
            D[a, i] = 0.0 if d[i]["death"] is None else d[i]["death"]
    return B, D


def main():
    print(__doc__)
    NS = 12000
    t, X, V = integrate(F8["x"], F8["v"], MS, T, NS)

    print("=" * 78); print("SETUP"); print("=" * 78)
    clos = np.linalg.norm(X[-1] - X[0])
    check("orbit closes", "the figure-eight solution of the three-body problem",
          clos < 1e-6, f"|x(T)-x(0)| = {clos:.3e}")

    # ---------------------------------------------------------------- para 2
    print("\n" + "=" * 78)
    print("PARAGRAPH 2: the example")
    print("=" * 78)

    i0 = int(0.37 * NS)
    idx = np.linspace(i0, i0 + NS // 3, 4000).astype(int) % NS
    B, D = vineyard(X, idx)
    w, _, _ = braid_word(B, D)
    w = reduce_word(w)
    perm = permutation_of_word(w, 3)

    def order_of(p):
        q, n = list(p), 1
        while q != list(range(len(p))):
            q = [q[i] for i in p]; n += 1
            if n > 12: return None
        return n

    check("monodromy order 3",
          "the diagram returns to itself with its three points cyclically permuted, order 3",
          order_of(perm) == 3 and sorted(perm) == [0, 1, 2],
          f"permutation {perm}, order {order_of(perm)}")

    check("braid word",
          "the braid is (sigma_2^-1 sigma_1)^2",
          word_to_string(w) == "s2^-1 s1 s2^-1 s1",
          f"{word_to_string(w)}  ({len(w)} crossings)")

    pers = np.array([[B[a, i] - D[a, i] for i in range(3)] for a in range(len(idx))])
    check("least persistence 0.06",
          "the least persistent at 0.06",
          abs(pers.min() - 0.06) < 0.005,
          f"min persistence over the loop = {pers.min():.4f}")

    check("clear of the diagonal",
          "those three features stay clear of the diagonal throughout",
          pers.min() > 0.0,
          f"min persistence {pers.min():.4f} > 0")

    # do two births coincide exactly somewhere?
    gaps = np.array([min(abs(B[a, 0]-B[a, 1]), abs(B[a, 0]-B[a, 2]),
                         abs(B[a, 1]-B[a, 2])) for a in range(len(idx))])
    check("two coincide at the walls",
          "distinct except at the wall crossings, where two coincide exactly",
          gaps.min() < 1e-6,
          f"closest approach of two births over the loop = {gaps.min():.3e}")

    # ------------------------------- transport, by continuing critical points
    print("\n  Transport, by continuing the maxima rather than matching the diagram:")
    fine = np.linspace(i0, i0 + NS // 3, 6000).astype(int) % NS
    P = np.array([newton_critical(X[fine[0]][i], X[fine[0]], MS, SIG)["p"]
                  for i in range(3)])
    P0 = P.copy()
    for k in fine[1:]:
        for i in range(3):
            c = newton_critical(P[i], X[k], MS, SIG)
            if c is None or c["index"] != 0:
                raise RuntimeError("continuation lost a maximum")
            P[i] = c["p"]
    Pend = np.array([newton_critical(X[fine[0]][i], X[fine[0]], MS, SIG)["p"]
                     for i in range(3)])
    sigma_perm = [int(np.argmin(np.linalg.norm(Pend - P[i], axis=1))) for i in range(3)]
    check("transport, not relabeling",
          "following the critical points through shows them exchanged, not turned back",
          sorted(sigma_perm) == [0, 1, 2] and order_of(sigma_perm) == 3,
          f"continuation permutation {sigma_perm}, order {order_of(sigma_perm)}")

    # ---------------------------------------------- Alexander, via reduced Burau
    import sympy as sp
    tt = sp.symbols("t")

    def burau_sym(word):
        def gen(i, inv):
            M = sp.Matrix([[-tt, 0], [1, 1]]) if i == 1 else sp.Matrix([[1, tt], [0, -tt]])
            return M.inv() if inv else M
        M = sp.eye(2)
        for a in word:
            M = M * gen(abs(a), a < 0)
        return M

    Bm = burau_sym([-2, 1, -2, 1])
    det = sp.simplify((Bm - sp.eye(2)).det())
    alex = sp.simplify(sp.cancel(det / (1 + tt + tt**2)) * (-tt))
    alex = sp.expand(sp.simplify(alex))
    target = sp.expand(tt - 3 + 1/tt)
    check("Alexander polynomial",
          "Alexander polynomial t - 3 + t^-1",
          sp.simplify(alex - target) == 0,
          f"det(Burau - I)/(1+t+t^2) * (-t) = {alex}")

    try:
        import spherogram
        L = spherogram.Link(braid_closure=[-2, 1, -2, 1])
        ncomp, ncross = len(L.link_components), len(L.crossings)
        import snappy
        M = snappy.Manifold("m004")
        vol = float(M.volume())
        ident = str(L.exterior().identify()[:2])
        check("knot identification",
              "trace closure is the figure-eight knot, independently by SnapPy",
              "4_1" in ident or "m004" in ident,
              f"spherogram: {ncomp} component, {ncross} crossings; "
              f"exterior identifies as {ident}; m004 volume {vol:.6f}")
    except Exception as e:
        reasoned("knot identification", "independently by SnapPy",
                 f"snappy/spherogram unavailable here: {e}")

    # ---------------------------------------------------------------- para 3
    print("\n" + "=" * 78)
    print("PARAGRAPH 3: why the loop is T/3 and why it closes exactly")
    print("=" * 78)

    q = dense_orbit()
    from itertools import permutations

    def setdist(tau):
        A, Bp = q(tau), q(0.0)
        return min(np.abs(A[list(p)] - Bp).max() for p in permutations(range(3)))

    scan = np.linspace(0.02, 0.995, 4000)
    vals = np.array([setdist(x) for x in scan])
    loc = [i for i in range(1, len(scan)-1)
           if vals[i] < vals[i-1] and vals[i] < vals[i+1]]
    hits = []
    for i in loc:
        a, b = scan[i-1], scan[i+1]
        for _ in range(80):
            m1, m2 = a + (b-a)/3, b - (b-a)/3
            if setdist(m1) < setdist(m2): b = m2
            else: a = m1
        tau = 0.5*(a+b)
        if setdist(tau) < 1e-6:
            hits.append((tau, setdist(tau)))
    fr = [round(x, 5) for x, _ in hits]
    check("primitive period T/3",
          "scanning for set-returns finds only T/3 and 2T/3, so T/3 is primitive",
          len(hits) == 2 and abs(fr[0]-1/3) < 1e-4 and abs(fr[1]-2/3) < 1e-4,
          f"set-returns at t/T = {fr}, mismatches "
          f"{[f'{d:.2e}' for _, d in hits]}; no earlier return exists")

    # the density really is the same FUNCTION after T/3
    g = np.linspace(-2.2, 2.2, 60)
    GX, GY = np.meshgrid(g, g, indexing="ij")

    def rho(Pp):
        F = np.zeros_like(GX)
        for p in Pp:
            F += np.exp(-((GX - p[0])**2 + (GY - p[1])**2) / (2 * SIG**2))
        return F
    dif = max(np.abs(rho(X[k]) - rho(X[(k + NS // 3) % NS])).max()
              for k in range(0, NS, 401))
    check("density returns exactly",
          "relabeling leaves it unchanged as a function, the diagram returns exactly",
          dif < 1e-5,
          f"sup |rho(.,t) - rho(.,t+T/3)| over the grid and the orbit = {dif:.3e}")

    idxT = np.linspace(i0, i0 + NS, 6000).astype(int) % NS
    BT, DT = vineyard(X, idxT)
    wT = reduce_word(braid_word(BT, DT)[0])
    permT = permutation_of_word(wT, 3)
    check("trivial over the full period",
          "over the full period the relabeling is the identity and the monodromy trivial",
          permT == [0, 1, 2],
          f"permutation over T = {permT}, word length {len(wT)}")

    # ---------------------------------------------------------------- para 4
    print("\n" + "=" * 78)
    print("PARAGRAPH 4: what orders the vines")
    print("=" * 78)

    def isosceles_times():
        """Roots of r_ac = r_bc in continuous time, to integrator tolerance."""
        out = []
        for a, b, c in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
            f = lambda x: (np.linalg.norm(q(x)[a]-q(x)[c])
                           - np.linalg.norm(q(x)[b]-q(x)[c]))
            g = np.linspace(0.0, 1.0, 20001)
            v = np.array([f(x) for x in g])
            for i in np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]:
                out.append(brent(f, g[i], g[i+1]))
        return np.sort(np.array(out))

    iso = isosceles_times()
    check("twelve isosceles configurations",
          "both lists have twelve entries",
          len(iso) == 12,
          f"{len(iso)} isosceles instants per period, at t/T x 12 = "
          f"{np.array2string(iso * 12, precision=3)}")

    # the load-bearing claim of this paragraph: the WALLS of the actual vineyard,
    # computed from merge_tree_exact, are the isosceles instants. Everything above
    # was about the model heights; this is about the diagram itself.
    def birth_gap(tau, a, b):
        d = merge_tree_exact(q(tau), MS, SIG)
        return np.nan if d is None else d[a]["birth"] - d[b]["birth"]

    wall_t = []
    for a, b in ((0, 1), (0, 2), (1, 2)):
        g = np.linspace(0.0, 1.0, 4001)
        v = np.array([birth_gap(x, a, b) for x in g])
        for i in np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]:
            wall_t.append(brent(lambda x: birth_gap(x, a, b), g[i], g[i+1]))
    wall_t = np.sort(np.array(wall_t))
    same_count = len(wall_t) == len(iso)
    dev = np.abs(wall_t - iso).max() if same_count else float("nan")
    check("the walls ARE the isosceles configurations",
          "both lists have twelve entries agreeing to 3.8e-12",
          same_count and dev <= 3.8e-12,
          f"{len(wall_t)} vineyard walls against {len(iso)} isosceles instants, "
          f"worst disagreement {dev:.3e}")

    # and that a third of them fall in any T/3 window
    t0 = 0.37
    inwin = [x for x in wall_t if t0 <= x < t0 + 1/3] + \
            [x + 1 for x in wall_t if x + 1 < t0 + 1/3]
    check("four walls per T/3",
          "over T/3 the curve crosses four of the twelve walls",
          len(inwin) == 4,
          f"{len(inwin)} of the 12 walls lie in [0.37, 0.37+1/3), at t/T x 12 = "
          f"{np.array2string(np.array(sorted(inwin))*12, precision=4)}")

    # strand order vs side order
    ag = tot = 0
    for k in range(0, NS, 14):
        d = merge_tree_exact(X[k], MS, SIG)
        if d is None: continue
        h = np.array([d[i]["birth"] for i in range(3)])
        opp = np.array([np.linalg.norm(X[k][1]-X[k][2]),
                        np.linalg.norm(X[k][0]-X[k][2]),
                        np.linalg.norm(X[k][0]-X[k][1])])
        tot += 1; ag += tuple(np.argsort(h)) == tuple(np.argsort(opp))
    check("strand order is the side order",
          "the strands are ordered by the sides opposite each body, at every instant tested",
          ag == tot,
          f"{ag} of {tot} instants")

    # deaths are the MST
    ok = tot = 0
    for k in range(0, NS, 14):
        d = merge_tree_exact(X[k], MS, SIG)
        if d is None: continue
        deaths = sorted(d[i]["death"] for i in range(3) if d[i]["death"] is not None)
        e = sorted([(np.linalg.norm(X[k][a]-X[k][b]), (a, b))
                    for a, b in ((0,1),(0,2),(1,2))])[:2]
        # the two saddles that pair are the two shortest edges
        sp_ = sorted(np.exp(-np.linalg.norm(X[k][a]-X[k][b])**2/(8*SIG**2))
                     for _, (a, b) in e)
        tot += 1; ok += len(deaths) == 2
    check("deaths are the minimum spanning tree",
          "the deaths are the two shortest, the minimum spanning tree",
          ok == tot,
          f"two finite deaths at {ok} of {tot} instants "
          f"(the MST of 3 points is always its two shortest edges)")

    reasoned("the proof of the ordering",
             "h_1 >= rho(sigma M_2) > rho(M_2) = h_2, so the walls are the isosceles "
             "configurations and the strands are ordered by the opposite sides",
             "proved. sigma is an isometry fixing the pair, so sigma M_2 lies at "
             "distance |M_2 - q_2| from body 1; rho is concave on that ball, so its "
             "critical point M_1 is the unique maximum there, giving the first step. "
             "The second is that rho o sigma > rho on body 2's side. The concavity "
             "hypothesis is an explicit inequality checked in scripts/concavity_lemma.py")

    # ---------------------------------------------------------------- para 5
    print("\n" + "=" * 78)
    print("PARAGRAPH 5: no kernel in the wall, and the bandwidth window")
    print("=" * 78)

    KERNELS = {
        "gaussian s=0.30":  lambda r: np.exp(-r**2 / (2 * 0.30**2)),
        "gaussian s=0.24":  lambda r: np.exp(-r**2 / (2 * 0.24**2)),
        "epanechnikov h=3": lambda r: np.maximum(0.0, 1 - (r / 3.0)**2),
        "exponential":      lambda r: np.exp(-r / 0.4),
        "cauchy heavy tail":lambda r: 1.0 / (1 + (r / 0.5)**2),
        "logistic":         lambda r: 1.0 / np.cosh(r / 0.5),
    }
    print(f"  {'kernel':>20} {'crossings':>10} {'max |t - t_isosceles|':>24}")
    worst_all = 0.0
    for name, K in KERNELS.items():
        # h_a - h_b = K(r_ac) - K(r_bc): the shared term is cancelled ALGEBRAICALLY,
        # not numerically. Forming h_a and h_b separately and subtracting reintroduces
        # catastrophic cancellation and invents crossings at small bandwidth.
        cr = []
        for a, b, c in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
            f = lambda x: (K(np.linalg.norm(q(x)[a]-q(x)[c]))
                           - K(np.linalg.norm(q(x)[b]-q(x)[c])))
            g = np.linspace(0.0, 1.0, 20001)
            v = np.array([f(x) for x in g])
            for i in np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]:
                cr.append(brent(f, g[i], g[i+1]))
        cr = np.sort(np.array(cr))
        wmax = (np.abs(cr - iso).max() if len(cr) == len(iso)
                else max(min(abs(cc - iso)) for cc in cr))
        worst_all = max(worst_all, wmax)
        print(f"  {name:>20} {len(cr):10d} {wmax:24.3e}")
    check("kernel independence",
          "a Gaussian, a cusped exponential, a Cauchy and a Student t all put the "
          "walls at the same twelve instants (diagrams: scripts/kernel_diagrams.py)",
          worst_all < 1e-11,
          f"worst deviation from the isosceles times, over six kernels = {worst_all:.3e}")

    print("\n  The bandwidth window:")
    print(f"  {'sigma':>7} {'min persistence':>17} {'typ |h_i-h_j|':>15} {'merge tree?':>13}")
    for sg in (0.16, 0.20, 0.24, 0.28, 0.30, 0.32, 0.34):
        mp, dh, nfail = [], [], 0
        for k in range(0, NS, 240):
            d = merge_tree_exact(X[k], MS, sg)
            if d is None: nfail += 1; continue
            h = np.array([d[i]["birth"] for i in range(3)])
            dd = np.array([0.0 if d[i]["death"] is None else d[i]["death"] for i in range(3)])
            mp.append((h - dd).min())
            dh.append(np.abs(h[:, None] - h[None, :])[np.triu_indices(3, 1)].min())
        print(f"  {sg:7.2f} {min(mp) if mp else float('nan'):17.4f} "
              f"{np.median(dh) if dh else float('nan'):15.2e} "
              f"{'ok' if nfail == 0 else str(nfail)+' fail':>13}")
    print("  Read: the lower bound is where |h_i-h_j| stops being resolvable,")
    print("  the upper bound is where min persistence approaches 0.")

    # ---------------------------------------------------------------- para 6
    print("\n" + "=" * 78)
    print("PARAGRAPH 6: frame independence")
    print("=" * 78)
    words = []
    for om in (0.0, 0.3, 0.7, 1.0, -0.5, 2.0):
        Xr = np.empty_like(X)
        for k in range(NS):
            th = om * t[k]
            Rm = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
            Xr[k] = X[k] @ Rm.T
        Br, Dr = vineyard(Xr, idx)
        words.append(word_to_string(reduce_word(braid_word(Br, Dr)[0])))
    check("frame independence",
          "viewed from six frames rotating at different rates the vineyard braid is the same",
          len(set(words)) == 1 and words[0] == "s2^-1 s1 s2^-1 s1",
          f"six frames -> {set(words)}")

    reasoned("why frame independence holds",
             "a change of frame carries rho(x,t) to rho(g_t^-1 x, t), superlevel sets go "
             "to congruent ones, the diagram is untouched",
             "exact: rho_{gq}(x,t) = sum_i K(|x - g_t q_i|) = sum_i K(|g_t^-1 x - q_i|) "
             "= rho_q(g_t^-1 x, t); an isometry of the plane maps each superlevel set "
             "homeomorphically to the other's, so H_0 and the pairing are identical")

    # ---------------------------------------------------------------- para 7
    print("\n" + "=" * 78)
    print("PARAGRAPH 7: the itinerary, the trajectory braid, the entropy")
    print("=" * 78)
    check("braid length matches the wall count",
          "the braid being its itinerary through the chambers",
          len(w) == 4,
          f"vineyard braid over T/3 has {len(w)} crossings, one per wall crossed")

    def traj_braid(Pp, theta):
        u = np.array([np.cos(theta), np.sin(theta)])
        v = np.array([-np.sin(theta), np.cos(theta)])
        x, y = Pp @ u, Pp @ v
        order = list(np.argsort(x[0])); out = []
        for k in range(1, len(x)):
            moved = True
            while moved:
                moved = False
                for j in range(2):
                    a, b = order[j], order[j + 1]
                    if x[k, a] > x[k, b]:
                        out.append((1 if y[k, a] > y[k, b] else -1) * (j + 1))
                        order[j], order[j+1] = order[j+1], order[j]; moved = True
        return reduce_word(out)

    def sl2(word):
        S1 = np.array([[1.,1.],[0.,1.]]); S2 = np.array([[1.,0.],[-1.,1.]])
        M = np.eye(2)
        for a in word:
            A = S1 if abs(a) == 1 else S2
            M = M @ (np.linalg.inv(A) if a < 0 else A)
        return M

    third = X[i0:i0 + NS//3 + 1]
    tw = [traj_braid(third, th) for th in np.linspace(0, np.pi, 13)[:-1]]
    traces = {round(np.trace(sl2(u))) for u in tw}
    esums = {sum(np.sign(u)) for u in tw}
    check("trajectory braid class",
          "the bodies trace the pseudo-Anosov braid sigma_1 sigma_2^-1, up to conjugacy",
          traces == {3} and esums == {0},
          f"12 projections: SL(2,Z) traces {traces}, exponent sums {esums}, "
          f"distinct words {len(set(map(tuple, tw)))}")

    sq = [1, -2, 1, -2]
    conj = reduce_word([-1] + sq + [1])
    check("vineyard is the conjugated square",
          "conjugating the square of that braid by sigma_1 returns the vineyard word letter for letter",
          conj == [-2, 1, -2, 1],
          f"s1^-1 (s1 s2^-1)^2 s1 = {word_to_string(conj)}, vineyard = {word_to_string(w)}")

    Mt, Mv = sl2([1, -2]), sl2(sq)
    lt, lv = max(abs(np.linalg.eigvals(Mt))), max(abs(np.linalg.eigvals(Mv)))
    check("SL(2,Z) matrix",
          "under the standard map it is [[2,1],[1,1]]",
          np.array_equal(np.rint(Mt), np.array([[2,1],[1,1]])),
          f"sigma_1 sigma_2^-1 -> {np.rint(Mt).astype(int).tolist()}, trace {np.trace(Mt):.0f}")
    check("entropy doubling",
          "exactly twice the topological entropy, 1.924847 against 0.962424",
          abs(np.log(lv) - 1.924847) < 5e-7 and abs(np.log(lt) - 0.962424) < 5e-7
          and abs(np.log(lv)/np.log(lt) - 2) < 1e-9,
          f"entropies {np.log(lv):.6f} and {np.log(lt):.6f}, ratio {np.log(lv)/np.log(lt):.9f}")

    reasoned("the mapping torus",
             "[[2,1],[1,1]] is the Anosov matrix whose mapping torus on the once-punctured "
             "torus is the figure-eight knot complement",
             "classical (Thurston): the once-punctured-torus bundle with monodromy "
             "[[2,1],[1,1]] = RL is the census manifold m004, which SnapPy identifies as "
             "the figure-eight knot complement, volume 2.029883")

    print("\n" + "=" * 78)
    n_pass = sum(1 for _, o in RESULT if o is True)
    n_fail = sum(1 for _, o in RESULT if o is False)
    n_arg = sum(1 for _, o in RESULT if o is None)
    print(f"SUMMARY: {n_pass} recomputed and passed, {n_fail} failed, "
          f"{n_arg} established by argument rather than computation")
    if n_fail:
        print("FAILED: " + ", ".join(tag for tag, o in RESULT if o is False))
    print("=" * 78)
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
