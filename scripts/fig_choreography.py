"""
Figure: a natural example of vineyard monodromy, from the three-body problem.

Every panel is computed, not illustrated. The braid diagram is drawn from the
measured crossing times, and the knot in the last panel is the closure that was
identified by conjugacy and confirmed by the Alexander polynomial.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mono.style import use_paper_style, tidy, SHAPE, BB, DD, LOOP, MUTED
use_paper_style(usetex=True)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from mono.nbody import FIGURE8, integrate
from mono.kde import kde_grid
from mono.kde_exact import merge_tree_exact, exact_maxima, exact_saddle
from mono.braid import braid_word, reduce_word, initial_order

COL = [BB, DD, "#159947"]
SIG = 0.30


def compute():
    F8 = FIGURE8
    m, T, N, NS = F8["m"], F8["T"], 3, 12000
    t, X, V = integrate(F8["x"], F8["v"], m, T, NS)
    i0 = int(0.37 * NS); i1 = i0 + NS // N
    K = 3000
    idx = np.linspace(i0, i1, K).astype(int)
    B = np.zeros((K, N)); D = np.zeros((K, N))
    for a, k in enumerate(idx):
        d = merge_tree_exact(X[k], m, SIG)
        for i in range(N):
            B[a, i] = d[i]["birth"]
            D[a, i] = 0.0 if d[i]["death"] is None else d[i]["death"]
    w, _, cr = braid_word(B, D)
    # trajectory braid for comparison
    wt, _, crt = braid_word(X[idx][:, :, 0], X[idx][:, :, 1])
    return dict(t=t, X=X, m=m, T=T, i0=i0, i1=i1, idx=idx, B=B, D=D,
                word=w, cross=cr, tword=wt)


def _rounded(ax, pts, r, **kw):
    """Polyline through pts with quadratic fillets of radius r at the corners."""
    from matplotlib.path import Path
    import matplotlib.patches as mpatches
    pts = [np.asarray(p, dtype=float) for p in pts]
    verts = [pts[0]]
    codes = [Path.MOVETO]
    for i in range(1, len(pts) - 1):
        a, b, c = pts[i - 1], pts[i], pts[i + 1]
        va, vc = a - b, c - b
        na, nc = np.linalg.norm(va), np.linalg.norm(vc)
        if na < 1e-12 or nc < 1e-12:
            continue
        rr = min(r, 0.5 * na, 0.5 * nc)
        verts += [b + va / na * rr, b, b + vc / nc * rr]
        codes += [Path.LINETO, Path.CURVE3, Path.CURVE3]
    verts.append(pts[-1]); codes.append(Path.LINETO)
    ax.add_patch(mpatches.PathPatch(Path(verts, codes), fill=False, **kw))


def _strand_tracks(crossings, K, init_order, n, x0, x1, half):
    """Strand polylines and crossing events, starting from the CORRECT order."""
    order = list(init_order)
    tracks = {i: [(x0, float(order.index(i)))] for i in range(n)}
    events = []
    span = x1 - x0
    for (tt, a, b, sgn) in crossings:
        events.append([x0 + 0.10 * span + 0.80 * span * (tt / K),
                       int(a), int(b), sgn, None, None])
    for ev in events:
        x, a, b = ev[0], ev[1], ev[2]
        pa, pb = order.index(a), order.index(b)
        assert abs(pa - pb) == 1, f"non-adjacent crossing at positions {pa},{pb}"
        ev[4], ev[5] = pa, pb
        for i in range(n):
            tracks[i].append((x - half, float(order.index(i))))
        order[pa], order[pb] = order[pb], order[pa]
        for i in range(n):
            tracks[i].append((x + half, float(order.index(i))))
    for i in range(n):
        tracks[i].append((x1, float(order.index(i))))
    return tracks, events, order


def _paint(ax, tracks, events, colors, half, n, lw=3.0, ms=13):
    for i in range(n):
        P = np.array(tracks[i])
        ax.plot(P[:, 0], P[:, 1], color=colors[i], lw=lw,
                solid_capstyle="round", zorder=3)
    for (x, a, b, sgn, pa, pb) in events:
        over = a if sgn > 0 else b
        ax.plot([x], [0.5 * (pa + pb)], "o", color="white", ms=ms,
                zorder=4, mec="none")
        y0 = pa if over == a else pb
        y1 = pb if over == a else pa
        ax.plot([x - half, x + half], [y0, y1], color=colors[over], lw=lw,
                solid_capstyle="round", zorder=5)


def draw_braid(ax, crossings, K, title, colors, init_order):
    """
    Braid diagram. init_order MUST be mono.braid.initial_order(births): the word
    is recorded relative to that order, and starting from list(range(n)) makes
    strands appear to cross through one another.
    """
    n = 3
    tracks, events, _ = _strand_tracks(crossings, K, init_order, n, 0.0, 1.0, 0.05)
    _paint(ax, tracks, events, colors, 0.05, n)
    ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.45, n - 0.55)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(title, loc="left", pad=6)


def draw_closure(ax, crossings, K, colors, init_order):
    """
    The measured braid, closed up. Each strand end at height p is carried around
    the OUTSIDE of the diagram back to the left end at the same height, with the
    arcs nested so they neither cross each other nor re-enter the diagram. The
    first version ran a straight segment across the top of each arc, which cut
    straight through the braid, and the arcs did not meet the strand ends.
    """
    n = 3
    XL, XR = 0.12, 0.80
    tracks, events, _ = _strand_tracks(crossings, K, init_order, n, XL, XR, 0.042)
    for p in range(n):
        off = 0.045 + 0.042 * p
        ybot = -0.42 - 0.30 * p
        _rounded(ax, [(XR, p), (XR + off, p), (XR + off, ybot),
                      (XL - off, ybot), (XL - off, p), (XL, p)],
                 r=0.055, color="0.45", lw=1.7, zorder=1,
                 capstyle="round", joinstyle="round")
    _paint(ax, tracks, events, colors, 0.042, n)
    m = 0.045 + 0.042 * (n - 1)
    ax.set_xlim(XL - m - 0.035, XR + m + 0.035)
    ax.set_ylim(-0.42 - 0.30 * (n - 1) - 0.20, n - 0.50)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def main():
    R = compute()
    X, idx, B, D = R["X"], R["idx"], R["B"], R["D"]
    K = len(idx)
    fig = plt.figure(figsize=(13.6, 7.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.92], hspace=0.42, wspace=0.30,
                          left=0.05, right=0.975, top=0.87, bottom=0.06)

    # (a) the orbit
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(X[:, 0, 0], X[:, 0, 1], color=MUTED, lw=1.0, zorder=1)
    P0 = X[R["i0"]]
    for i in range(3):
        seg = X[R["i0"]:R["i1"] + 1, i]
        ax.plot(seg[:, 0], seg[:, 1], color=COL[i], lw=2.0, zorder=3)
        ax.plot(*P0[i], "o", color=COL[i], ms=11, mec="white", mew=1.3, zorder=5)
    ax.set_aspect("equal"); tidy(ax, 4, 4)
    ax.set_title(r"\textbf{(a)}\quad the figure-eight choreography:" "\n"
                 r"\phantom{\textbf{(a)}\quad}three bodies, one curve, over $T/3$",
                 loc="left", pad=8)

    # (b) the density and its critical points
    axb = fig.add_subplot(gs[0, 1])
    lim = 1.7 * np.abs(X).max()
    F, Xg, Yg = kde_grid(P0, R["m"], SIG, lim, 400)
    axb.contourf(Xg, Yg, F, levels=18, cmap="Blues")
    axb.contour(Xg, Yg, F, levels=18, colors="white", linewidths=0.4, alpha=0.7)
    mx = exact_maxima(P0, R["m"], SIG)
    for i in range(3):
        axb.plot(*mx[i]["p"], "o", color=COL[i], ms=10, mec="white", mew=1.3, zorder=5)
    sad = []
    for i in range(3):
        for j in range(i + 1, 3):
            s = exact_saddle(0.5 * (P0[i] + P0[j]), P0, R["m"], SIG)
            if s is None:
                continue
            if not any(np.linalg.norm(s["p"] - q) < 1e-5 for q in sad):
                sad.append(s["p"])
    # The number of saddles of the field is NOT constant round the loop: it is 2
    # when the bodies are near collinear and 3 when they are spread, the third
    # appearing together with a central local minimum. Only two of them are ever
    # merge saddles, because three maxima merge into one component in exactly two
    # steps, so the H_0 vineyard is unaffected. This instant has two, and the
    # assertion keeps the panel honest if the instant is ever changed.
    assert len(sad) == 2, f"instant has {len(sad)} saddles, caption says two"
    for q in sad:
        axb.plot(*q, "x", color="black", ms=9, mew=2.0, zorder=5)
    axb.set_aspect("equal"); tidy(axb, 4, 4)
    # Caption states what is actually there. The maxima are NEAR the bodies, not
    # at them: neighboring bumps pull them off by up to 2.9e-02 at this sigma.
    # And there are two saddles, not one per pair: bodies 1 and 2 have none
    # between them. Both agree with the Morse count #max - #saddle = 1 for a
    # density decaying to zero, and with the merge tree's two finite deaths.
    axb.set_title(r"\textbf{(b)}\quad the mass density $\rho(x;t)$: its three"
                  "\n" r"\phantom{\textbf{(b)}\quad}maxima ($\bullet$) and the two merge "
                  r"saddles ($\times$)", loc="left", pad=8)

    # (c) the vineyard
    ax3 = fig.add_subplot(gs[0, 2], projection="3d")
    ts = np.linspace(0, 1, K)
    # Vertical guides at the three STARTING diagram points. The vineyard is
    # closed, so the diagram at t = T/3 is the same set of three points; each
    # vine therefore has to land on one of these guides, and the panel is only
    # doing its job if the reader can see that it lands on a DIFFERENT one.
    for i in range(3):
        ax3.plot([B[0, i]] * 2, [D[0, i]] * 2, [0, 1], color="0.62", lw=0.8,
                 ls=(0, (3, 3)), zorder=1)
    for i in range(3):
        ax3.plot(B[:, i], D[:, i], ts, color=COL[i], lw=2.2, zorder=3)
        ax3.scatter(B[0, i], D[0, i], 0, color=COL[i], s=30, edgecolor="white",
                    linewidth=0.6, depthshade=False, zorder=4)
        ax3.scatter(B[-1, i], D[-1, i], 1, facecolor="white", edgecolor=COL[i],
                    s=34, linewidth=1.6, depthshade=False, zorder=4)
    ax3.set_xlabel(r"birth", fontsize=9, labelpad=8)
    ax3.set_ylabel(r"death", fontsize=9, labelpad=8)
    ax3.set_zlabel(r"$t$", fontsize=9, labelpad=-2, rotation=0)
    ax3.set_zticks([0, 1]); ax3.set_zticklabels([r"$0$", r"$T/3$"])
    ax3.xaxis.set_major_locator(MaxNLocator(3)); ax3.yaxis.set_major_locator(MaxNLocator(3))
    ax3.tick_params(labelsize=6.5, pad=2)
    ax3.set_box_aspect((1, 1, 0.95), zoom=0.86)
    ax3.view_init(elev=18, azim=-62)
    ax3.set_title(r"\textbf{(c)}\quad the closed vineyard: each vine ($\bullet$) ends"
                  "\n" r"\phantom{\textbf{(c)}\quad}on a \emph{different} start point ($\circ$)",
                  loc="left", pad=0)

    # (d) trajectory braid, (e) vineyard braid, (f) the closure
    axd = fig.add_subplot(gs[1, 0])
    wt, _, crt = braid_word(X[idx][:, :, 0], X[idx][:, :, 1])
    init_t = initial_order(X[idx][:, :, 0])
    draw_braid(axd, crt, K,
               r"\textbf{(d)}\quad braid of the \emph{trajectories}: "
               r"$\sigma_1\sigma_2^{-1}$", COL, init_t)
    axe = fig.add_subplot(gs[1, 1])
    init_v = initial_order(R["B"])
    draw_braid(axe, R["cross"], K,
               r"\textbf{(e)}\quad braid of the \emph{vineyard}: "
               r"$(\sigma_2^{-1}\sigma_1)^2$, conjugate to the square of (d)", COL, init_v)
    axf = fig.add_subplot(gs[1, 2])
    draw_closure(axf, R["cross"], K, COL, init_v)
    axf.set_title(r"\textbf{(f)}\quad the same braid, closed up: this is the"
                  "\n" r"\phantom{\textbf{(f)}\quad}figure-eight knot $4_1$, "
                  r"$\Delta(t)=t-3+t^{-1}$", loc="left", pad=6)

    fig.suptitle(r"Vineyard monodromy of order 3 in the three-body problem, and the knot "
                 r"it closes to", fontsize=14, y=0.965)
    out = "figs/fig3_choreography.pdf"
    fig.savefig(out); fig.savefig(out.replace(".pdf", ".png"))
    print("wrote", out)
    print(f"  vineyard braid crossings: {len(R['word'])}, trajectory braid: {len(wt)}")
    P0 = np.column_stack([B[0], D[0]]); P1 = np.column_stack([B[-1], D[-1]])
    Cm = np.linalg.norm(P0[:, None, :] - P1[None, :, :], axis=2)
    from scipy.optimize import linear_sum_assignment as _lsa
    ri, ci = _lsa(Cm)
    print(f"  panel (c) closure: vine i ends on start point {list(map(int, ci))}, "
          f"max mismatch {Cm[ri, ci].max():.2e}")


if __name__ == "__main__":
    os.makedirs("figs", exist_ok=True)
    main()
