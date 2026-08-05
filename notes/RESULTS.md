# Vineyard monodromy occurs in Newton's equations

Shawn Ray, 29 July 2026. Every number is reproduced by a script in this
repository and by a dated log in `out/verify/`. Nothing is quoted from memory.

---

## The claim

Existing examples of vineyard monodromy are engineered: *Braiding Vineyards*
realizes a prescribed braid by constructing a manifold to order. The question is
whether the phenomenon occurs in a system that was **not constructed to exhibit
it**.

It does. Two independent statements, in order of how much I think they matter.

## 1. Standard vineyard tracking is blind to a whole class of monodromy

This is checkable in an afternoon and depends on no naturalness question.

There are two distinct mechanisms, and they need different trackers.

- **Pairing monodromy.** The critical points stay put; the elder-rule pairing
  reassigns around a codimension-2 point, so a minimum that was paired with one
  maximum returns paired with another. This is the `A_1^2/A_1^2` mechanism.
- **Transport monodromy.** The critical points themselves move and exchange
  places over the loop, while staying distinguishable throughout.

Frame-to-frame Hungarian matching of `(birth, death)` coordinates, which is what
essentially every vineyard implementation does including all of mine, **fails on
the second**. On rotating waves whose true monodromy is an `n`-cycle it returned
order 1 in **11 of 11** cases, with the closest pair of diagram points never
below `9.57e-03`, so this is not a near-collision artifact. The vines exchange in
the death coordinate while remaining separated in birth, and distance-based
matching follows the wrong branch through the exchange.

Tracking instead by transport, that is by continuity of the critical points in
the domain, recovers order `n` in all 11.

Choosing the tracker to match the mechanism is not optional. On the measured
X-ray scan below, the two disagree in both directions on the *same* grid:
assignment reports 16 nontrivial loops, the mechanism-appropriate pair tracker
reports 71. Every candidate tested from either was then rejected by
certification, so neither disagreement changes the verdict there, but a bare
count from the wrong tracker means nothing.

An earlier version of this note claimed the transport tracker had confirmed the
X-ray negative. That was a category error and is withdrawn: transport follows
maxima by position, the X-ray scan is pairing monodromy where the maxima do not
move, and a tracker blind to the mechanism confirms nothing.

## 2. A mechanism that forces transport monodromy, and Newton realizes it

**Symmetry destroys it; transport creates it.** If two extrema are exchanged by
an isometry of the *instantaneous* field, they carry equal critical values and
their paired saddles carry equal values, so their diagram points **coincide**.
Permuting coincident points is not a nontrivial permutation. Measured with an
exact analytic merge tree, for `k` identical bumps exchanged by a rotation:

```
k = 2, 3, 4, 5:   spread of birth values      = 0.000e+00
                  closest pair of diagram pts = 0.000e+00
```

This rejects the obvious first guesses, rigid rotating waves, k-armed spirals,
the rotating regular N-gon. Symmetry, which looks like the natural source of a
forced permutation, is exactly the wrong place to look.

**Statement.** Let `A` be a fixed envelope on the circle and `h` have period
`2*pi/n`. Put `f(theta, t) = A(theta) h(theta - omega t)`. Then `f` is exactly
`T/n`-periodic, no critical point is created or destroyed, the crest set rotates
rigidly by `2*pi/n`, and the monodromy is a single `n`-cycle, non-degenerate
exactly when `A` is not `2*pi/n`-periodic. Verified for `n = 3, 4, 5`; a vine
traced continuously travels `119.43` of `120` degrees and lands on its
neighbor's starting position **and** starting value.

The figure-eight three-body choreography is a worked instance of the same
principle, with the exchange driven by a time shift rather than a spatial
isometry: `x_i(t + T/N) = x_{i+1}(t)` acts on time, so at every instant the
bodies occupy different points of the curve and stay distinguishable.

```
integrator DOP853, rtol = atol = 1e-13, 12000 steps
orbit closure |x(T) - x(0)|                     7.389e-08
configuration identical after T/3, mismatch     2.981e-08,  relabelling [1, 2, 0]
configuration identical after T,   mismatch     0.00e+00,   relabelling [0, 1, 2]

                      sigma = 0.20    0.25    0.30
vineyard braid        s2^-1 s1 s2^-1 s1   (4 crossings)   order 3
trajectory braid      s1 s2^-1            (2 crossings)   order 3

Alexander polynomial of the closure  t - 3 + 1/t,  coefficients [1, -3, 1]
figure-eight knot 4_1 has            t - 3 + 1/t,  coefficients [1, -3, 1]   MATCH
```

**The vineyard braid is not the square of the trajectory braid.** I claimed it
was, then withdrew the claim because the supporting comparison used two different
strand labellings. It is now settled, and against the original claim.

*Proposition.* For any choreography, both braids have the same underlying
permutation over `T/N`: a single `N`-cycle. For the trajectories this is the
definition, `x_i(t + T/N) = x_{i+1}(t)`. For the vineyard, diagram point `i`
belongs to body `i`, and at `t = T/N` the configuration is identical as a set, so
body `i` sits at body `i+1`'s former position carrying its former birth and
death.

*Corollary.* If the vineyard braid were the square of the trajectory braid up to
conjugacy then, under `B_N -> S_N`, an `N`-cycle would equal the square of an
`N`-cycle. That square is an `N`-cycle for odd `N` but splits into two cycles of
length `N/2` for even `N`. **So the relation is false for every even `N`**, with
no computation.

*And the `N = 3` evidence was vacuous.* At `N = 3` the square of a 3-cycle is a
3-cycle, so the permutation invariant cannot distinguish anything, and both
exponent sums are `0`, so `e(V) = 2 e(T)` reads `0 = 0`. Neither conjugation
invariant available at `N = 3` can detect a difference. The agreement that
prompted the claim was never evidence for it.

Nothing is asserted for odd `N >= 5`, and nothing needs to be: a relation that
fails for every even `N` is not a relation.

**The loop is `T/3`, not `T`.** `T/3` is the primitive period in the space of
*unlabelled* configurations, and the full period `T` is its triple cover. Over `T`
the relabelling is the identity, the braid is `(s2^-1 s1)^6` with 12 crossings,
and the monodromy is order 1, converged at 4000 and 12000 loop samples with
closure mismatch `0.00e+00`. The `T/3` loop is still forced by the orbit rather
than chosen, but it must be named correctly.

**One further caveat, found while checking the above.** The figure-eight passes
through six configurations per period, at `t/T = 0, 1/6, 1/3, 1/2, 2/3, 5/6`, at
which two bodies are mirror images and therefore carry *exactly* equal density
values, so two diagram points coincide (`0.000e+00`). A `T/3` loop **based** at
one of those instants begins at a degenerate diagram and its braid is not well
defined. Eleven basepoints were tested. Exactly one of them, `t/T = 0`, sits on a
symmetric instant, and it is exactly the one that fails, returning order 2 with 3
crossings. The other ten all avoid the symmetric instants and all return order 3,
with word `(s2^-1 s1)^2` or its conjugate `(s1 s2^-1)^2`. So the count is ten out
of ten among admissible basepoints, and the single exception is predicted by the
degeneracy rather than left over as an anomaly. Passing *through* those
instants mid-loop is fine and expected, since a coincidence of diagram points is
what a vine crossing is. Only the basepoint matters.

## What is chosen, and what is not

Not chosen: the orbit, the loop `T/3`, the permutation, the braid, the knot.

**Chosen: the kernel bandwidth.** Three point masses are three delta functions;
the field is a Gaussian kernel density estimate at a bandwidth `sigma` that I
picked. The result is stable over `sigma = 0.20` to `0.30`, and it fails outside
that: `sigma = 0.35` degenerates at `N = 3`, and `sigma = 0.12` mis-assigns at
`N = 5` where two diagram points sit `2e-4` apart. Any statement about this
example has to carry the bandwidth with it.

## What I am not claiming

- **Not a priority claim.** I have not surveyed the literature well enough to say
  "first". The claim is that these examples exist and reproduce.
- **Not "natural" in the sense of measured.** The figure-eight is an exact
  solution of a physical law, not an experimental record.
- **Not a rare specimen.** The theorem above is a recipe: it manufactures
  monodromy of any order from a two-line formula. That makes the right claim
  *ubiquity*, not rarity. Monodromy is what transport-recurrent systems do.
- **Not the `A_1^2/A_1^2` singularity of arXiv:2607.01046.** That theorem is
  about the planar radial transform. This monodromy is forced by transport, and I
  have not shown the mechanisms coincide.

## Measured data: nine negatives, and why

Documented in `notes/RECORD.md`. The reason is now precise and is not the tracker,
which was the obvious suspect and has been ruled out by re-running with the
mechanism-appropriate one.

Making the two parameters physical makes them **sampled**, and between adjacent
samples the field moves by less than the measurement noise, so a codimension-2
point localized inside one cell describes the noise. In a cone-beam X-ray scan of
a walnut, consecutive detector rows differ by `0.0034` against a noise level of
`0.0039`. On a grid coarsened until every step beat the noise by `3.3x`, the
exact pair tracker found 71 loops at order 2, every one a single transposition,
and all of the ones tested collapsed to order 1 under shrinking and showed wall
types `DD,DD` or `BB,BB` rather than the required one birth-birth and one
death-death. Since the same diagnostic returns `BB,DD` at the coastline
generators, it is discriminating rather than under-powered. That is recorded as
D10.

The coastline result, by contrast, survives. **Both trackers were run at both
certified generators of the Japan outline, and they disagree, which is the point:**

```
exact PAIR tracker    order 2, single transposition, radii 1e-3 .. 1e-6   both points
TRANSPORT tracker     order 1 at every radius 1e-3 .. 1e-7               both points
                      measured transport jump 0.0
```

The coastline is *pairing* monodromy: the critical points do not move at all, and
the transport jump of 0.0 says so directly. A tracker that follows maxima by
position must return the identity there, and does. Reporting only the favourable
number would have hidden the best evidence in the note, which is that the two
trackers separate the two mechanisms exactly as the framework predicts.

The wall-type test agrees. At both Japan generators the diagram has 8 points and
shows **2 pairing changes, one birth-birth and one death-death**, stable at radii
`1e-3`, `1e-4`, `1e-5`. That is the `A_1^2/A_1^2` signature. An earlier run of this
diagnostic returned 0 and 0 and was wrong: it used a different curve construction
which resolves only 4 diagram points instead of 8, so the relevant walls were not
in the feature set at all.

## Reproducing

```
python3 scripts/verify_headline.py       # orbit, braid, Alexander polynomial
python3 scripts/rotating_wave_test2.py   # symmetry degenerates
python3 scripts/modulated_wave.py        # transport gives order n
python3 scripts/rerun_transport.py       # assignment vs transport on measured data
```


---

# Appendix: the rotating-wave theorem in full

Result of 29 July 2026. This replaces the "symmetry implies monodromy" statement,
which is false in a specific and fixable way, with one that holds and is verified.

## What is false

> "f_{t+T} = f_t o g for an isometry g, so the diagram returns identically but the
> tracked points are permuted by the action of g. If g moves the extrema, the
> permutation is provably nontrivial."

The last step fails. If an isometry of the INSTANTANEOUS field exchanges two
extrema, then those extrema carry equal critical values, their paired saddles
carry equal values too, and their diagram points **coincide**. Vineyard monodromy
is a permutation of diagram points; permuting coincident points is not a
nontrivial permutation of anything.

Measured, with the exact analytic merge tree, for k identical bumps exchanged by
a rotation:

```
k = 2, 3, 4, 5:   spread of birth values among the k bumps = 0.000e+00
                  closest pair of diagram points           = 0.000e+00
```

Exactly zero, to machine precision, as it must be. **Symmetry is what kills it.**
This rejects rigid rotating waves and k-armed spirals directly, and it is the
same degeneracy that makes the rotating N-gon useless (RECORD, N-series).

## What is true

The exchange must be realised by **transport**, not by symmetry: the features
must swap places by moving along the loop while staying at different phases of
one process, so that they remain distinguishable at every instant. That is what
the figure-eight choreography does, and it is why it works:

```
body values at t = 0    [1.00395, 1.00395, 1.00773]
body values at t = T/3  [1.00773, 1.00395, 1.00773 -> shifted]
spread                   3.779e-03   (nonzero, so the permutation is visible)
```

## The theorem

Let A be a fixed positive envelope on the circle and h have period 2*pi/n. Put

    f(theta, t) = A(theta) * h(theta - omega t),     omega = 2*pi/T.

Then:

1. **The loop closes exactly and is the system's own.** Since h has period
   2*pi/n, f(theta, t + T/n) = f(theta, t) identically. The loop length T/n is
   the physical period, not an analyst's choice. Verified closure ~1e-15.
2. **D7 holds by construction.** No crest is created or destroyed; there are n
   maxima at every instant.
3. **The monodromy is the n-cycle.** Over T/n the crest set rotates rigidly by
   2*pi/n, so each maximum is carried to the position its neighbor occupied,
   and the induced permutation of features is a single n-cycle.
4. **It is non-degenerate exactly when A is not 2*pi/n-periodic.** The n maxima
   then sit at different heights of a fixed envelope, so their diagram points are
   distinct. A constant recovers the degenerate case of the false theorem.

## Verification

Tracking maxima by continuity of position, which is unambiguous because positions
move smoothly, against the naive matcher:

```
  n    mod   seed | by transport         cycles | by matching   min sep
  3   0.15      1 |            3    [(0, 1, 2)] |           1   1.23e-02
  3   0.15      2 |            3    [(0, 1, 2)] |           1   8.19e-02
  3   0.35      1 |            3    [(0, 1, 2)] |           1   1.73e-02
  4   0.15      1 |            4 [(0, 1, 2, 3)] |           1   3.62e-02
  4   0.35      2 |            4 [(0, 1, 2, 3)] |           1   1.86e-01
  5   0.15      1 |            5 [(0,1,2,3,4)]  |           1   9.57e-03
  5   0.35      2 |            5 [(0,1,2,3,4)]  |           1   4.18e-02
```

A single vine traced continuously over one period:

```
  t=0.00  position 118.89 deg  value +0.84621
  t=0.40  position 169.65 deg  value +1.00512
  t=0.80  position 213.84 deg  value +0.96394
  t=1.00  position 238.32 deg  value +0.83160
  travel 119.43 deg of a 120 deg slot; index map [1, 2, 0]
```

It lands on maximum 1's starting position AND its starting value.

## The methodological finding, which matters independently

**The standard matcher misses this monodromy entirely.** Frame-to-frame Hungarian
matching of diagram points reported order 1 in 11 of 11 cases where the true
answer is order n, with the closest pair of diagram points never below 9.6e-03,
so this is not a near-collision artifact. The vines swap by passing each other in
the death coordinate while remaining separated in birth, and distance-based
matching follows the wrong branch through the exchange.

Consequence: for this class the features must be tracked by **transport**, that
is by continuity of the critical points in the domain, not by matching diagram
points. Any pipeline that tracks vineyards by assignment alone will report the
identity on a system that has order-n monodromy. This is the same lesson as
mono/exact.py, now with a concrete class of false negatives rather than false
positives.

## Where this is measured

The requirement is a rotating wave carried through a stationary spatial
inhomogeneity, so that the rotational symmetry is broken in the lab frame while
the temporal period is exact. Candidates, in order of how well the loop is
forced and how well the field is resolved:

- spinning azimuthal modes in annular combustors, where the burner pattern is the
  fixed envelope, imaged by OH* chemiluminescence. Thermoacoustic instability is
  a first-order engineering problem, so the loop is not only natural but watched
- rotating detonation engines
- wavy vortex flow in Taylor-Couette between eccentric or machined cylinders
- spiral waves anchored to a tissue heterogeneity, optical mapping
- rotating convection in an annulus with a fixed thermal boundary pattern

The gate to apply before touching any of these, from D10: the field must change
by more than the measurement noise between consecutive samples of the loop.

---

# What the monodromy is measuring

The permutation on its own is a restatement of the choreography property. The
content is in what the vineyard's walls turn out to be.

## The walls are isosceles configurations

For N equal masses and any radial kernel K, the density peak belonging to body i
is `h_i = sum_j K(r_ij)`. For three bodies the shared term cancels:

```
h_1 - h_2 = K(r_13) - K(r_23)
```

so if K is monotone decreasing then `h_1 = h_2` exactly when `r_13 = r_23`, that
is when body 3 is equidistant from bodies 1 and 2. The vineyard orders its
strands by h, so:

> **A crossing of strands i and j is an isosceles configuration with apex at the
> third body.**

Measured on the figure-eight: 12 vineyard crossings per period, 12 isosceles
configurations per period, agreeing to **3.8e-12 of a period**, and the crossing
pair is the pair not at the apex in all 12. They sit at `t/T = k/12` exactly,
evenly spaced, which the orbit's symmetry forces.

## The bandwidth cannot move a wall

The condition `r_ik = r_jk` contains no kernel. Six kernels of completely
different shape give the same crossing times:

```
Gaussian sigma 0.30      12 crossings   max |t - t_isosceles|  9.4e-13
Gaussian sigma 0.60      12                                    5.4e-14
exponential exp(-2r)     12                                    2.2e-14
inverse 1/(1+r^2)        12                                    1.1e-14
Epanechnikov-like        12                                    5.8e-15
heavy tail 1/(1+r)       12                                    6.7e-15
```

This answers the standing objection that the bandwidth is an analyst's choice.
It is, and it is invisible to the answer. What the bandwidth controls is only
whether the crowding difference is numerically resolvable: the typical
`|h_i - h_j|` runs from 3.7e-04 at sigma 0.22 to 4.0e-02 at sigma 0.34, which is
exactly why the stability window has the lower bound it has. The upper bound is
where the density peaks merge. Neither bound is topological.

## The braid is an invariant of the shape, not of the motion

Birth and death are built from mutual distances alone, so they are unchanged by
any rigid motion of the configuration, including a time-dependent one. Viewing
the same orbit from a frame rotating q times per period:

```
q =  0.0   vineyard  s2^-1 s1 s2^-1 s1     trajectory  s1 s2^-1
q =  0.5   vineyard  s2^-1 s1 s2^-1 s1     trajectory  s1 s2^-1
q =  1.0   vineyard  s2^-1 s1 s2^-1 s1     trajectory  s1 s2^-1 s1^-1 s2^-1 s1^-1 ...
q =  2.0   vineyard  s2^-1 s1 s2^-1 s1     trajectory  s2^-1 s1^-1 s1^-1
q = -1.5   vineyard  s2^-1 s1 s2^-1 s1     trajectory  s1 s1 s2
q =  3.0   vineyard  s2^-1 s1 s2^-1 s1     trajectory  s2^-1 s1^-1 s1^-1 s2^-1
```

with the mutual distances identical to 6.7e-16 in every frame. The vineyard braid
is invariant; the trajectory braid picks up a twist per turn of the frame.

So the two braids are invariants of different objects. The trajectory braid is an
invariant of the motion in the plane. **The vineyard braid is an invariant of the
curve in shape space**, the quotient of configuration space by translations and
rotations. Persistent homology of the mass density performs that reduction for
free, and the isosceles loci it lands on are the standard walls of that quotient:
three codimension-one strata that meet where all three distances are equal, at
the equilateral configurations.

This also supplies, in correct form, the claim withdrawn in T4. The vineyard
braid does carry something the trajectory braid does not, but the difference is
not a crossing count at N=5. It is that one is frame-independent and the other is
not, and that is a proof rather than a numerical comparison.

## What the monodromy adds

The order is 3 because the choreography's time shift permutes the three bodies
cyclically, and over `T/3` the shape curve advances by 4 of the 12 isosceles
walls, one third of the circuit. Over the full period it crosses all 12 and
returns to the identity. So the monodromy is measuring how the orbit's own
time-shift symmetry acts on the chambers of the isosceles arrangement in shape
space, and the braid records the itinerary through them.

---

# Rigor pass: the objections, answered

Written against the two claims that have to hold: that this is monodromy, and
that the isosceles reading is right. Log at `out/verify/rigor.log`.

**Setup, stated once.** Filtration: superlevel sets of the planar mass density
`rho(x,t) = sum_i K(||x - q_i(t)||)`, degree 0. Birth is the value at a local
maximum, death the value at the saddle where its component merges into an older
one, and the one essential class dies at the infimum of `rho`, which is 0. The
domain is the plane, no boundary treatment is needed since `rho` decays to 0.
Braid convention: strands ordered left to right by birth value at the first
slice, time increasing, crossing sign taken from the death coordinate, closure
the ordinary trace closure.

**The closure is exact, not numerical.** The choreography property
`x_i(t + T/3) = x_{i+1}(t)` is exact for the true solution, and the density is a
SUM over bodies, so relabelling leaves it unchanged as a function:
`rho(., t + T/3) = rho(., t)` identically. The diagram is a function of `rho`, so
it returns exactly. The measured residual is 1.2e-09 and saturates there as the
integrator tolerance is tightened from 1e-9 to 1e-13, which places it in the
precision of the published initial conditions rather than in the integration.

**T/3 is primitive.** Scanning the whole period for times at which the
configuration returns as a SET finds exactly two, `s/T = 0.33333` with
relabelling [1,2,0] and `s/T = 0.66667` with [2,0,1]. There is no shorter return,
so the order is 3 and not a multiple of something smaller.

**The permutation is not bookkeeping.** The diagram carries no labels, so a
permutation of it is nontrivial only if its points are distinct, which is exactly
what fails for the rotating N-gon. Birth values at `t/T = 0.37` are
`[1.032995, 1.000124, 1.032904]`, spread 3.3e-02. Vines are followed by
continuity of critical points in the domain, which is the original Vines and
Vineyards definition; the birth-death plane is the image, not the bookkeeping.

**Three features, never near the diagonal.** Over a T/3 loop sampled at 600
instants the count is 3 throughout, and the minimum persistence is 0.42 at
sigma 0.22, 0.21 at 0.26 and 0.060 at 0.30. Pairwise separation touches zero only
at the isosceles instants, which is what a crossing is.

**The isosceles identification is exact in one direction, empirical in the
other.** At an isosceles configuration with apex k, the configuration is
symmetric under reflection in the perpendicular bisector, and that reflection is
an isometry of the whole density exchanging the other two bodies. Their maxima
are therefore exchanged by an isometry and their birth values are EXACTLY equal,
whatever the kernel and however far the maxima sit from the bodies. Constructed
on random isosceles triangles: `|birth_i - birth_j|` is 0.000e+00 at maximum
displacements up to 0.30.

That is the direction `isosceles => wall`. The converse, that every wall is
isosceles, is not proved here; what is measured is that on the figure-eight the
two lists have the same length, 12 and 12, and agree to 3.8e-12. The correct
statement is therefore: **every isosceles configuration is a wall, for any
monotone radial kernel, and on this orbit there are no others.**

**Invariance is exact.** 200 random rotations and translations change the diagram
by 4.4e-16, since `rho` composed with an isometry has identical critical values.
Scaling is different: scaling the configuration by 0.8 changes the diagram by
1.2e-01. So the honest statement is that the braid is an invariant of the curve
in configuration space modulo TRANSLATIONS AND ROTATIONS. The wall locus
`r_ik = r_jk` is separately scale invariant, so the walls, but not the diagram,
descend to the shape sphere.

**The bandwidth window, and why it has the bounds it has.**

```
 sigma   min persistence   typ |h_i-h_j|    verdict
  0.16          0.8052        7.5e-06       walls not resolvable
  0.20          0.5519        5.3e-04       walls not resolvable
  0.24          0.3065        5.4e-03       ok
  0.28          0.1222        2.3e-02       ok
  0.30          0.0598        4.0e-02       ok
  0.32          0.0190        6.4e-02       feature near the diagonal
  0.34          0.0008        9.5e-02       feature near the diagonal
  0.36       degenerate
```

Two competing constraints, neither topological: below the window the crowding
differences fall towards numerical resolution, above it a feature reaches the
diagonal. The honest range is 0.24 to 0.30, narrower at the top than earlier
statements in this file. The wall LOCATIONS are unaffected throughout, since
`r_ik = r_jk` contains no sigma.

---

# The mode count, which is not the component count

A real hazard, and it had to be checked rather than assumed. The number of modes
of an isotropic Gaussian mixture is not the number of components:
Carreira-Perpinan and Williams conjectured that a homoscedastic mixture of M
components in more than one dimension has at most M modes, and report a
counterexample due to Duistermaat, communicated to them privately: three equal
Gaussians at the vertices of an equilateral triangle in the plane, with a fourth
mode at the center for a range of variances.

That is this setup exactly. Worse, `exact_maxima` in `mono/kde_exact.py` seeds
Newton from the body positions and returns one maximum per body, so it cannot
find a mode at the centroid. A missing maximum changes sublevel-set
connectivity, so the elder-rule pairing among the three that were found could be
wrong, and the braid would inherit it.

**The finder.** Newton from every cell of a dense grid covering the configuration
plus a margin, deduplicated, classified by the sign of the Hessian eigenvalues.
No seeding from the bodies.

**Validated against the counterexample.** On the equilateral triangle of
circumradius 1 the finder reports

```
sigma 0.700   3 max, 3 saddles, 1 minimum   center is a minimum
sigma 0.710   4 max, 3 saddles, 0 minima    center is a MAXIMUM
sigma 0.735   4 max, 3 saddles, 0 minima
sigma 0.740   1 max, 0 saddles, 0 minima
```

so it does see the fourth mode. The window is only about 0.025 wide in sigma,
and a first scan in steps of 0.04 stepped straight over it, which is the reason
this needed a fine scan rather than a coarse one.

**Result on the figure-eight.** Exactly three modes everywhere tested: 40 frames
by 16 bandwidths from 0.20 to 0.35, 640 combinations, mode count 3 in 634 of
them. The six exceptions are all at sigma 0.35 and are 2 modes, not 4: two peaks
merging, which is what ends the usable window at the top. There is no
Duistermaat mode anywhere near it.

The structural reason is that the figure-eight is nowhere near equilateral. The
relative spread of its three mutual distances runs from 0.7244 to 0.8513 over the
period, and the fourth mode lives at spread 0.

**A correction this forces.** The entry T3c reported 2 or 3 saddles depending on
the instant, and read the 3-saddle case as three saddles plus a central minimum.
That was wrong. The exhaustive search gives 3 maxima, 2 saddles and 0 minima at
every instant; the midpoint-seeded count double-counted, because seeding at the
midpoint of a far-apart pair converges to a saddle belonging to a different pair
and a 1e-5 dedup tolerance then keeps both. The merge tree itself was never
affected, since it uses two finite deaths throughout, but the accounting of why
was wrong and is corrected here.

## Two things the mode-count check left asserted, now measured

**How far from equilateral does the fourth mode survive?** The note said the orbit
"never comes near" the equilateral shape, which was true but unmeasured, since
the window had only been mapped at spread 0. Deforming an equilateral triangle by
pulling one vertex radially and stepping sigma at 0.005, a fifth of the width of
the window at spread 0:

```
relative spread   4-mode window in sigma   width
   0.0000              0.710 to 0.735      0.030
   0.0025              0.720 to 0.735      0.020
   0.0050              0.725 to 0.740      0.020
   0.0100              0.735 to 0.740      0.010
   0.0173              0.745 only          0.005
   0.0247              none
   0.0392 and beyond   none
```

The window shrinks monotonically and closes by a relative spread of about 0.02.
The orbit's minimum is 0.7244, more than thirty-five times further out. The
sigma step is fine enough that a narrower window cannot slip between samples, and
by 0.0173 the window is already down to a single sample.

**The Morse count as a self-check.** For a smooth function on the plane decaying
at infinity, `#max - #saddle + #min` must be 1, so anything else means a critical
point was missed or misclassified. Across the 640 in-window combinations the
count is 1 every time, at grid resolutions 24, 32 and 44. That is a stronger
statement than "three modes", because it would catch a single missed critical
point; it would not catch a maximum and a saddle missed together, so it is
necessary rather than sufficient.

**A reporting bug, corrected.** An earlier version of `mode_count.py` printed the
Morse count as the Cartesian product of the count SETS collected across all
frames rather than per frame. With maxima {2,3} and saddles {1,2} at sigma 0.38
it therefore printed {0,1,2}, which reads as the finder failing at large
bandwidth. It was not: the per-frame count is 1 at every frame and every
resolution tried. The display is fixed and the log regenerated.

---

# What the diagram retains, and why the braid uses only part of it

Another crack at the interpretation, prompted by asking how much of the
configuration the diagram actually keeps. I had been calling the birth values a
measure of "crowding", which treats them as a lossy summary. The answer is more
interesting than either "lossy" or "complete".

**The strand order is the ordering of the triangle's sides.** From
`h_1 - h_2 = K(r_13) - K(r_23)` and its two companions, the ordering of the peak
heights is the ordering of the sides opposite each body. Checked at 858 instants
along the orbit, in agreement at **858 of 858**. The vineyard braid is literally a
braid of the three sides of the triangle, ordered by length, and its walls are the
isosceles configurations because that is where two sides are equal.

**In exact arithmetic the heights determine the shape.** Writing
`u = K(r_23), v = K(r_13), w = K(r_12)` and subtracting the self term
`g_i = h_i - K(0)`, the system

```
g_1 = w + v,   g_2 = w + u,   g_3 = v + u
```

is linear, so `u = S - g_1` with `S = (g_1+g_2+g_3)/2`, and `K` inverts. Using the
density at the bodies this recovers the three mutual distances to **4.9e-09**, and
a triangle is determined by its sides. So the map from shape to diagram is a
bijection in principle.

*I got this wrong twice before getting it right, and both are worth recording.*
The first version forgot the self term `K(0)`, which is what makes the system
linear, and the recovery was off by 1.5. The second version measured the
Jacobian, printed determinants of 1e-17 and condition numbers of infinity, and
concluded underneath them that the map was "nonsingular and well conditioned",
which the numbers directly contradicted.

**In practice the inversion is useless, for a structural reason.** With a Gaussian
of width sigma the influence of a separation r on a peak height is
`exp(-r^2 / 2 sigma^2)`, so at the working bandwidths:

```
sigma   shortest side   longest side   K(short)   K(long)
 0.22       0.7976         1.9324      1.4e-03    1.8e-17
 0.26       0.7976         1.9324      9.1e-03    1.0e-12
 0.30       0.7976         1.9324      2.9e-02    9.8e-10
 0.34       0.7976         1.9324      6.4e-02    9.7e-08
```

The longest side moves the peak heights in the twelfth to thirtieth decimal
place, against heights of order 1. Recovering it needs more precision than double
arithmetic carries, which is why the Jacobian determinant runs from 1e-8 down to
1e-17 and the condition number to infinity. The metric content is there and is
not retrievable.

### The two measured claims, now proved

Both were the same statement, and it has a three-line proof. Write L for the
perpendicular bisector of q_i q_j, sigma for the reflection in L, H_i and H_j for
the two half planes, M_i for the maximum belonging to body i and h_i = rho(M_i).

> **Proposition.** Let K be strictly decreasing. Suppose sigma M_j lies in the
> basin of body i, so h_i >= rho(sigma M_j), and M_j is interior to H_j. If
> r_ik < r_jk then h_i > h_j.
>
> *Proof.* sigma fixes the pair {q_i, q_j}, so the only term that moves belongs
> to body k: `rho o sigma (x) - rho(x) = K(|x - sigma q_k|) - K(|x - q_k|)`. For
> x interior to H_j, `|x - sigma q_k| < |x - q_k|`, since sigma q_k lies in H_j
> and L is the perpendicular bisector of q_k and sigma q_k. K strictly
> decreasing gives rho o sigma > rho on the interior of H_j, so
> `h_i >= rho(sigma M_j) = rho o sigma (M_j) > rho(M_j) = h_j`. []

The two corollaries are exactly the two claims. The strand order is the order of
the opposite sides, since h_i > h_j iff r_ik < r_jk and r_jk is the side opposite
body i. And h_i = h_j forces r_ik = r_jk, so the walls are the isosceles
configurations and nothing else. The reflection argument the note already had is
the equality case of this one.

No derivatives, no bandwidth, no Gaussian, which is why the four-kernel
measurement now reads as a confirmation rather than as the evidence.

**The first attempt was wrong and the failure is worth keeping.** It took the
hypothesis to be that h_i is the largest value of rho on the whole half plane
H_i, which is tidier and false: when r_ik < r_jk, body k lies in H_i too and its
peak is frequently the taller one. Measured margins:

```
sigma    0.22     0.24     0.26     0.28     0.30     0.32
margin  -2.7e-5  -1.5e-4  -5.3e-4  -1.7e-3  -3.8e-3  -7.3e-3     FAILS everywhere
```

Only the single value rho(sigma M_j) is ever needed, and sigma M_j sits near q_i
by construction. With that hypothesis the margin is 0 or -2.2e-16 at every
bandwidth tried, and it vanishes exactly at the walls, which is the equality
case: there sigma carries M_j onto M_i.

Log: `out/verify/wall_proof.log`. The remaining gap is that the hypothesis is
checked along the orbit rather than proved, so the note says "as it does here".

### Can the concavity hypothesis be made kernel-independent? Yes.

(C) was written for a Gaussian, which made it look Gaussian-specific. It is not.
For a radial kernel the Hessian of `K(|x-q|)` has closed-form eigenvalues, `K''(r)`
along the radius and `K'(r)/r` across it:

```
Hess K(|x-q|) = K''(r) P_radial + (K'(r)/r) P_tangential
lambda_max(Hess rho)(x) <= sum_j max( K''(r_j), K'(r_j)/r_j )              (C')
```

so rho is strictly concave at x as soon as (C') is negative. For the Gaussian,
`K'(r)/r = -K/s^2` and `K'' = K(r^2/s^4 - 1/s^2)`, and K'' is the larger for every
r > 0, so (C') collapses to `sum_j K_j r_j^2 < s^2 sum_j K_j`, which is (C). The
criterion was always general; only my statement of it was not.

```
kernel              worst (C')   worst true lambda   concave?
gaussian s=0.30        -5.5767            -5.5768       yes
cauchy s=0.35          -9.7697            -9.8843       yes
student t3 s=0.35      -6.5216            -6.6179       yes
```

The bound is nearly tight, so little is lost by avoiding eigenvalues.

**Cusped kernels are a separate and easier case.** If K has a corner at r = 0, as
`exp(-r/l)` does, the body's own term pulls with magnitude `|K'(0+)| = 1/l` from
every direction while the others pull with `K(r_ij)/l < 1/l`. The maximum
therefore sits exactly ON the body: measured displacement 0.00e+00 at three
bandwidths. So `sigma M_j = q_i = M_i`, the hypothesis holds with equality, and no
concavity is needed at all.

So the note's "wherever that concavity holds" is a hypothesis about the
configuration and the bandwidth, not about the kernel, and it is checkable for
any radial kernel from the formula above.

### Can the minimum-spanning-tree claim be proved? Not yet, but the evidence is much stronger.

The merge tree of rho is the single-linkage tree of the ultrametric
`m(i,j) = max over paths from q_i to q_j of min rho`, so the two merges are always
the two largest m. The claim is the further statement that this ordering is the
reverse ordering of the Euclidean distances. That is not automatic: a ridge is
raised both by the pair being close and by the third body sitting near their
midpoint, and those can compete.

Searched for a crossover. **Zero disagreements in 37127 clean random triangles**
across three bandwidths, with 39767 dropped as degenerate (fewer than two saddles,
or a saddle value under 1e-12 where the ridge is beneath the arithmetic).

**The first version of this search was wrong and reported 259 disagreements.** It
attributed each saddle to the nearest pair midpoint, which in an elongated
triangle collapses two saddles onto one pair. The clearest "counterexample" had
sides 3.99, 1.79, 2.20 and listed a single pair for two saddles, which is what
gave it away. Asking directly whether Newton from each pair midpoint reaches an
index-1 critical point removes the artifact entirely.

So the claim stands as measured, now over 37000 random triangles rather than 858
instants of one orbit. The note still says "measured rather than proved", which
remains the right thing for it to say.

Log: `out/verify/two_questions.log`.

### Finishing it: concavity replaces the basin hypothesis

The gap was `h_i >= rho(sigma M_j)`, argued through basins. Concavity closes it
in one line, and the hypothesis becomes an explicit inequality.

Two exact facts. First, sigma is an isometry fixing the pair {q_i, q_j}, so

```
|sigma M_j - q_i| = |sigma M_j - sigma q_j| = |M_j - q_j| = delta_j
```

and sigma M_j lies in the closed ball B(q_i, delta_j). Second, a critical point of
a strictly concave function on a convex set is its unique maximum there. M_i is a
critical point and lies in that ball, so

```
h_i = rho(M_i) = max over B(q_i, delta_j) of rho  >=  rho(sigma M_j)
```

with no mention of basins, gradient ascent, or which peak belongs to whom.

An explicit criterion, no eigenvalues needed. For a Gaussian,
`Hess rho (x) = (1/s^2) sum_j K_j [ (x-q_j)(x-q_j)^T/s^2 - I ]`, so for any unit u

```
u' Hess rho u = (1/s^2) sum_j K_j [ ((x-q_j).u)^2/s^2 - 1 ]
             <= (1/s^2) [ sum_j K_j |x-q_j|^2/s^2 - sum_j K_j ]
```

and rho is strictly concave at x as soon as

```
sum_j K_j |x - q_j|^2  <  s^2 sum_j K_j                                    (C)
```

the K-weighted mean square distance to the bodies being below s^2. Near a body
almost all the weight sits on that body, where the distance is nearly zero, so (C)
is comfortable.

> **Theorem.** Let K be the Gaussian of width s and suppose (C) holds throughout
> each ball B(q_i, delta), delta = max_j |M_j - q_j|. Then h_i > h_j precisely
> when r_ik < r_jk. Hence the walls of the vineyard are exactly the isosceles
> configurations, and the strand order is the order of the opposite sides.

Measured margins. `slack` is the margin in (C); `lam_max` is the true top Hessian
eigenvalue, showing how much (C) gives away by avoiding eigenvalues; `concl` is
the theorem's conclusion checked end to end.

```
case                        n        slack      lam_max    delta       concl
figure-eight, sigma=0.22  300    4.506e-02   -1.924e+01   0.0054   0.000e+00
figure-eight, sigma=0.24  300    4.993e-02   -1.505e+01   0.0126   0.000e+00
figure-eight, sigma=0.28  300    5.137e-02   -8.357e+00   0.0463   0.000e+00
figure-eight, sigma=0.30  300    4.517e-02   -5.577e+00   0.0795   0.000e+00
figure-eight, sigma=0.32  300    3.176e-02   -3.030e+00   0.1337   0.000e+00
figure-eight, sigma=0.36  210   -2.114e-02   -3.466e-01   0.2775   0.000e+00
random, min sep > 0.75   1500    4.285e-02   -6.509e+00   0.0557  -2.2e-16
random, min sep > 1.00   1500    8.380e-02   -1.044e+01   0.0054  -0.0e+00
random, min sep > 1.50   1500    8.999e-02   -1.111e+01   0.0000  -2.2e-16
```

(C) holds across the whole working window with roughly a factor of two in hand,
against `s^2 = 0.09` at the working bandwidth. At s = 0.36 the sufficient
criterion goes negative while the true eigenvalue is still -3.5e-01, so rho is
concave and only the bound has given out; that bandwidth is past the note's upper
limit anyway. And (C) gives out exactly where the diagram does: sweeping an
isosceles triangle tighter at s = 0.30, the slack falls to 1.6e-02 at base 0.75
and at base 0.70 there are no longer three peaks at all.

Log: `out/verify/concavity_lemma.log`. The note now gives the concavity reason in
place of the basin assertion.

### How far the last hypothesis can be closed

The proposition needs `h_i >= rho(sigma M_j)`, which the note reports as checked.
`scripts/basin_lemma.py` splits that into three steps and proves the first two.

**Proved, for every configuration.** For a Gaussian, `K'(r) = -(r/sigma^2)K(r)`,
so `grad rho(x) = -(1/sigma^2) sum_j (x - q_j) K(|x - q_j|)` and a critical point
satisfies `sum_j (M - q_j) w_j = 0` with `w_j = K(|M - q_j|)`. Every maximum is
therefore a **weighted centroid of the bodies**,

```
M_i = sum_j w_j q_j / sum_j w_j ,     |M_i - q_i| <= sum_{j!=i} w_j r_ij / sum_j w_j
```

exactly, not to first order. Verified to 1.1e-15 over 4300 configurations. Since
sigma is an isometry fixing the pair {q_i, q_j},

```
|sigma M_j - q_i| = |sigma M_j - sigma q_j| = |M_j - q_j| = delta_j
```

so sigma M_j sits exactly as far from body i as M_j sits from body j. That much
is unconditional, and it is what the note now says.

**Still measured.** That this puts sigma M_j inside body i's basin. The whole
remaining gap is one length against another, and both are measurable:

```
case                       n     identity      bound        hyp    headroom
figure-eight, sigma=0.24  300    4.74e-16   -7.9e-17   -2.2e-16      0.0364
figure-eight, sigma=0.28  300    8.93e-16   -5.2e-17   -2.2e-16      0.1340
figure-eight, sigma=0.30  300    8.93e-16   -4.8e-17   -2.2e-16      0.2304
figure-eight, sigma=0.32  300    9.42e-16   -4.9e-17   -2.2e-16      0.3872
random triangles, 0.30   4000    1.12e-15   -2.2e-16   -2.2e-16      0.1453
```

`hyp` is `h_i - rho(sigma M_j)`, zero to roundoff because the walls are the
equality case. `headroom` is the largest peak displacement over the distance from
the body to the nearest saddle. It is a factor of a few, not an order of
magnitude: 4 per cent at sigma = 0.24, 23 per cent at the working 0.30, 39 per
cent at 0.32. Comfortable, and worth not overstating. Holding on 4000 random
triangles as well as on the orbit says it is not a property of this orbit.

**Not done.** A uniform bound on the peak displacement against the basin radius,
over all triangles. That is what would close the gap completely, and it is a real
lemma rather than a loose end.

A discarded intermediate is worth recording. A first attempt certified (*) by
checking that rho(M_i) exceeds the max of rho on a sphere of radius 1.05 max
delta about q_i. That radius puts the sphere barely outside M_i whenever body i
carries the largest displacement, so the margin is zero by construction, and the
check reported roundoff-level failures on random triangles that were nothing of
the kind.

Log: `out/verify/basin_lemma.log`.

### Kernel independence, redone on actual diagrams

The original six-kernel check ran through the model heights `h_i = sum_j K(r_ij)`,
the density read AT the bodies. That is the right quantity for the algebra but it
is not the birth value, and `mono/kde_exact` hard-codes Gaussian derivatives so it
cannot produce a diagram for any other kernel. Two things came out of redoing it
with GUDHI, which needs no derivatives (`scripts/kernel_diagrams.py`).

The claim holds, on real diagrams. Thresholding the birth gap does not work off a
cubical complex, since births are quantised to the grid and a fixed threshold
finds hundreds of spurious minima, including for the Gaussian control. Asking the
question directly does work: at an isosceles instant two births must coincide, and
anywhere else they must not.

```
kernel                  max gap ON twelfths   min gap OFF   separation
gaussian s=0.30               8.500e-09        6.193e-07         73x
cauchy s=0.35                 1.932e-08        1.413e-02     731256x
exponential cusp 0.35         2.482e-07        1.475e-03       5942x
student t3 s=0.35             2.451e-08        7.946e-03     324127x
quartic compact h=1.1     no three-peak diagram
```

And one of the original six does not survive. A compactly supported kernel wide
enough to reach all three bodies makes rho unimodal: the Epanechnikov with
bandwidth 3 turns the density into a single paraboloid with ONE maximum, not
three. As a monotone ordering function it is fine, which is why it passed the
model-height test; as a density kernel it has no three-peak diagram to be
independent about. The note said "from compactly supported to heavy tailed" and
now names the four kernels that actually have diagrams.

**The conclusion, which is the useful part.** The diagram's ordering information
is decided by differences of order 1e-2 and is robust; its metric information is
below machine precision. The braid is built from the ordering alone. So the braid
is not a coarsening of the diagram out of convenience: it is precisely the part
of the diagram that survives, and the part that is discarded is the part that
could not have been trusted. That is a reason to prefer the braid to the diagram
here, rather than a concession.

**And this is special to three bodies.** The heights give N numbers and the shape
needs N(N-1)/2, equal only at N = 3. For four bodies the heights cannot determine
the shape even in principle.

**Measured afterwards, and it sharpens the conclusion rather than softening it.**
The account above rests on an inversion formula, which invites the objection that
the formula is the problem and not the diagram. `scripts/completeness_test.py`
tests the diagram directly instead, with no formula anywhere in it: take every
pair of instants around the orbit, and compare the distance between their birth
triples with the distance between their triangles taken up to congruence.

```
sigma   instants   pairs     min |dB|/|dS|   median    smallest |dB| among
                                                       well separated shapes
0.24      1000     499500      1.294e-04    3.528e-02        3.220e-04
0.30      1000     499500      5.515e-04    1.766e-01        3.667e-03
```

No collisions: no two clearly different shapes on this orbit have nearly equal
diagrams, and the lower Lipschitz ratio stays positive. So the map from shape to
diagram **is injective here**, and the earlier verdict was too pessimistic in one
direction and too optimistic in the other. It is injective, and it is useless to
invert, and those are consistent: injectivity only needs the ratio to stay above
zero, while inversion needs it bounded away from zero by more than the arithmetic
can resolve. 1.3e-04 clears the first bar and not the second.

Two claims retired for the record. The idealized-height inversion is exact, but
it is exact about the density read AT the bodies, not about the birth values,
which sit at maxima displaced towards the other two. Feeding true births through
it leaves 33 to 35 per cent relative error at every bandwidth in the working
window, and no convergence: below sigma about 0.18 the merge tree stops returning
anything at all, so the apparent zeros at small bandwidth in an earlier run were
empty tables, not agreement. The docstring of `complete_invariant.py` has been
corrected to say so.

What none of this touches is the braid, which reads the ordering and never the
values. That is why the ordering result, 858 of 858, is the load-bearing one.

Logs: `out/verify/completeness_test.log`, `out/verify/complete_invariant.log`.

---

# Closing the four remaining gaps

## Cross-check against an independent implementation

Every number ran through `mono/`, so a bug in the two core routines would have
propagated everywhere. Checked against GUDHI 3.12, which uses a cubical complex
on a sampled grid, a different algorithm and a different implementation.

*The planar merge tree.* Superlevel `H_0` of rho is sublevel `H_0` of `-rho`.
GUDHI's values converge to the exact merge tree as the grid refines, on both
births and deaths:

```
grid    max |birth diff|   max |death diff|   spacing
 400        2.5e-04            5.3e-05        0.0150
 800        9.2e-05            2.1e-05        0.0075
1600        2.8e-05            4.3e-06        0.0037
```

*The circle machinery,* used for the coastline, the X-ray scan and the rotating
waves. Against GUDHI's periodic cubical complex, on random trigonometric
polynomials: **every finite pair agrees to 0.00e+00.** The one apparent
discrepancy was mine: `mono` reports the essential class as a pair, global
maximum to global minimum, while GUDHI reports it separately with infinite death,
so the two lists are offset by one until that is accounted for.

## What the deaths encode

The births order the strands; the crossing signs come from the deaths, and those
were unexplained. For superlevel `H_0` of a kernel density the merge structure is
single linkage, so the prediction is that the deaths are the minimum spanning
tree. Measured at 858 instants:

```
essential class is the body opposite the LONGEST side        858 of 858
the two finite deaths are the saddles on the two SHORTEST    858 of 858
```

So the diagram splits cleanly. The births are the ordering of all three sides;
the deaths are the two shortest, the minimum spanning tree of the three bodies.
The braid's strand order comes from the first and its crossing signs from the
second, which completes the account of what the braid is made of.

## The twelve walls, and what is collinear among them

```
isosceles instants   12   at t/T = k/12
collinear instants    6   at t/T = k/6
of the 6 collinear instants, how many are also isosceles:  6 of 6
```

So six of the twelve walls are Euler configurations, collinear with one body
equidistant from the other two, and the other six are proper isosceles triangles.
The vineyard braid counts all twelve. The trajectory braid's six crossings are
equal-x events at the ODD twelfths, which is projection dependent and does not
coincide with the syzygies at all.

## The symmetry group, and what it does and does not force

Searching over spatial isometries, time shifts and time reversal:

```
identity   t -> +t + {0, 4, 8}/12
reflect y  t -> +t + {2, 6, 10}/12
reflect x  t -> -t + {2, 6, 10}/12
rotate pi  t -> -t + {0, 4, 8}/12          group order 12
```

An earlier version searched shifts only, found six, and asserted twelve
underneath the number six. With time reversal included it is genuinely twelve.

But it does not give the result I wanted from it. Every shift appearing is an
even multiple of `T/12`, so each group element maps Euler configurations to Euler
configurations and isosceles-only to isosceles-only. The group forces each set of
six to be evenly spaced at `T/6`; it does not mix the two classes, so the
interleaving that puts a wall at every `T/12` is a separate fact and not a
consequence of these symmetries. Stated as measured rather than derived.

---

# The mechanism of the monodromy, corrected and sharpened

Chasing the braid verification turned up something more important than the
verification. The order-3 claim is right, but for a reason I had not established,
and I nearly convinced myself it was wrong along the way.

## The vines genuinely collide, six times per period

At each of the six collinear (Euler) instants the two outer bodies are mirror
images, so they have equal peak heights AND equal merge saddles: their diagram
points **coincide exactly**. Classifying all twelve walls:

```
 t/T    tied pair   |P_i - P_j| at the wall   collinear?   what happens
 0/12      0,1            0.00e+00             0.0e+00     vines MEET
 1/12      0,2            1.05e-01             1.1e+00     labels swap only
 2/12      0,2            1.39e-08             2.1e-09     vines MEET
 3/12      1,2            1.05e-01             1.1e+00     labels swap only
 ...                                                       (alternating)
                     six collisions, six label swaps
```

At the six proper isosceles instants the tie involves the tallest peak, which
carries the essential class at death 0. The two points are then `(h, 0)` and
`(h, s)`, distinct; the points stay put and only the body labels swap. No vine
crossing there.

At the six Euler instants the two points are identical, the vines meet, and how
they continue is not decided by the diagram. Following the bodies through, which
move smoothly, the two points are **exactly exchanged**:

```
t/T = 0,  offset -400 steps:  body0 (1.005040, 0.529297)   body1 (1.003065, 0.468503)
          offset +400 steps:  body0 (1.003065, 0.468503)   body1 (1.005040, 0.529297)
   cross costs 0.00000, bounce costs 0.12165  ->  they CROSS
```

Over `T/3` there are two such collisions, giving transpositions that compose to a
3-cycle. **Order 3 confirmed, and now explained.**

## Which corrects the account of why matching fails

I had written that diagram-plane matching is "blind" to this monodromy. The
mechanism is sharper and worse. At an exact collision both resolutions have zero
matching cost, so the assignment is degenerate and the solver returns whichever
it happens to pick, which is the bounce. Tracking the diagram points that way
returns the IDENTITY, at every basepoint and every sampling, with residual 1e-8.
That is not a near-miss being mishandled: it is a tie the method cannot break in
principle. The critical points break it, because they pass through smoothly.

*I nearly reported the identity as the answer.* Two intermediate computations
gave order 1, and both were indexed by body, so they were reporting the
relabelling rather than the vine permutation; a third was the degenerate match
above. The claim survived, but it needed this to be established rather than
assumed.

## The braid, and a caveat that comes with it

Because the vines meet, the vineyard is not generic: three curves in
`(birth, death, t)` that intersect do not form a braid without a resolution
convention. Two natural conventions, taking the death comparison just before or
just after each event, give mirror words, `(s2^-1 s1)^2` and `(s2 s1^-1)^2`.

An independent extraction confirms this is the only difference: an event-list
algorithm, finding every birth crossing first and then reading generator indices
off the positions held at that instant, returns the same generator sequence
`2,1,2,1` with every sign flipped, over both `T/3` and the full period.

The knot does not care. SnapPy identifies the closure of BOTH conventions as

```
m004, 4_1, K2_1, K4a1        hyperbolic volume 2.029883
```

which is the figure-eight knot, matching the known volume, because 4_1 is
amphichiral. The full period closes to the 3-component link L12a1882, volume
21.631599.

## The chamber itinerary

The six chambers are the six orderings of the side lengths, and the walls between
them make the Cayley graph of `S_3`, a hexagon. The orbit's itinerary over one
period is

```
A B E F D C A B E F D C
```

every chamber exactly twice, in the same cyclic order both times: the hexagon
traversed twice.

## What happens at N > 3

`h_i - h_j = sum over k not i,j of [K(r_ik) - K(r_jk)]`. At `N = 3` that sum has
one term and the wall is `r_ik = r_jk`, with no kernel in it. At `N > 3` it is a
weighted sum and the wall moves with the kernel. Comparing a Gaussian against a
heavy tail on the same path:

```
N = 3   12 walls each, max shift 0.000e+00
N = 4   24 walls each, max shift 4.4e-02
N = 5   36 walls each, max shift 6.6e-02
```

So the kernel-free geometry is exactly a three-body phenomenon. For four bodies
and up the braid depends on the smoothing and is no longer an invariant of the
shape curve alone.

---

# The two research questions, answered

## What the knot means

The braid `sigma_1 sigma_2^-1` is the classical pseudo-Anosov 3-braid. Under the
standard map from the 3-strand braid group modulo its center onto SL(2,Z),

```
sigma_1 -> [[1,1],[0,1]]      sigma_2 -> [[1,0],[-1,1]]
sigma_1 sigma_2^-1 -> [[2,1],[1,1]]      trace 3, dilatation 2.618034
```

which is the golden ratio squared, and `[[2,1],[1,1]]` is the Anosov map whose
mapping torus on the once-punctured torus is the figure-eight knot complement.
SnapPy: `m004`, volume `2.029883`, identified as `4_1`. That is the reason the
figure-eight ORBIT is celebrated in braid terms: its trajectory braid is exactly
that map.

The vineyard braid is its square:

```
(sigma_2^-1 sigma_1)^2 -> [[2,3],[3,5]]   trace 7, dilatation 6.854102
dilatation is the square, to 0.00e+00
topological entropy   trajectories 0.962424   vineyard 1.924847   ratio 2.000000
```

So the knot is not decorative. The vineyard braid is the square of the orbit's own
pseudo-Anosov braid, it has the same stable foliation, its **topological entropy
is exactly twice** that of the trajectories, and it closes to the very knot whose
fibration that braid generates.

### Audited afterwards, and one claim withdrawn

All of the above was taken from the literature rather than measured, so
`scripts/trajectory_braid.py` measures it. Three of the four parts hold, and the
fourth does not.

Holds. The trajectory braid over T/3, computed directly from the integrated orbit
in the same convention as the vineyard, is **conjugate to `sigma_1 sigma_2^-1` in
18 of 18 projection directions**, with exponent sum 0 and SL(2,Z) trace 3 in every
one. The literal word is not stable (`s1 s2^-1`, `s2 s1^-1` and
`s1^-1 s2^-1 s1 s2` all appear); the conjugacy class is.

Holds, and more cleanly than stated. The squaring needs no representation at all:

```
s1^-1 (s1 s2^-1 s1 s2^-1) s1  =  s2^-1 s1 s2^-1 s1
```

which is the measured vineyard word letter for letter. Confirmed in reduced
Burau, faithful on three strands, so this is decided rather than sampled. Closure
is a conjugacy invariant, so the knot follows from this line alone. The two braids
are conjugate, not equal, and the note now says so.

Holds. Both braids are read over the same T/3, so the entropy ratio of exactly 2
is not an artifact of normalization. Per unit time as well: 0.456420 against
0.912839.

**Withdrawn.** "The doubling is the same factor of two as twelve walls against
six, the vineyard seeing every isosceles configuration and the trajectories half
of them." The trajectory crossing count is projection dependent: 6 per period for
most directions but **12** for a band near theta = pi/2, and the crossings sit on
the wall locus only in the single projection along the long axis of the eight
(there they are the odd twelfths, the proper isosceles configurations, to 1e-4).
The vineyard count of 12 is projection independent; the trajectory count of 6 is
not, so the comparison was not an invariant one. The sentence is gone from the
note, replaced by the point the check actually established: writing the
trajectory braid down requires a frame and a projection, and the vineyard braid
requires neither.

Log: `out/verify/trajectory_braid.log`.

## What replaces the clean story at N > 3

The wall between bodies i and j is `sum over k not i,j of [K(r_ik) - K(r_jk)] = 0`.
The shared term `K(r_ij)` has already cancelled, so as the bandwidth falls the
surviving sum is dominated by the single body nearest to the pair, and the wall
tends to

```
r_{i k*} = r_{j k*},     k* = argmin over k not i,j of min(r_ik, r_jk)
```

that is, **i and j equidistant from the body nearest to them**: the isosceles
condition again, with the apex now chosen as the nearest neighbor of the pair.
This is kernel-free for every N. Measured against Gaussian kernels of falling
width:

```
N   predicted walls   sigma 0.80    0.40      0.20      0.10      0.05
3        12           0.0e+00    0.0e+00   0.0e+00   0.0e+00   (numerics fail)
4        24           3.3e-02    7.2e-03   6.7e-04   1.7e-04   (numerics fail)
5        44           2.5e-02    7.8e-03   2.2e-03   1.7e-04   0.0e+00
```

At `N = 3` there is only one candidate for `k*`, so the limit is already an
identity and holds at every bandwidth, which is why that case is kernel-free
outright. For `N > 3` it is a genuine limit, and the braid is an invariant of the
shape only asymptotically.

*A wrong guess, recorded.* My first proposal was that the walls tend to equal
NEAREST-NEIGHBOUR distance, `d_i = d_j`. That is wrong, and it fails already at
`N = 3`: it ignores that the shared term cancels, so when i and j are each other's
nearest neighbor the dominant terms disappear and the subleading ones decide. It
predicted 18 walls at `N = 3` where there are 12.

At `N = 6` the predicted count is 120 and between 64 and 112 are realised, with
the distance to the predicted walls still falling like the bandwidth. Some
predicted walls are simply not crossed by the path, so the limit locus is an upper
bound on what a given orbit sees. Stated as measured.

---

# Literature check

Every attribution the note makes or relies on, verified against the sources
rather than recalled.

**The orbit.** Found numerically by Moore, *Braids in classical gravity*, Phys.
Rev. Lett. 70 (1993) 3675, who classified periodic three-body orbits by braid
type; existence proved by Chenciner and Montgomery, *A remarkable periodic
solution of the three-body problem in the case of equal masses*, Annals of
Mathematics 152 (2000) 881. Their statement is that the bodies "follow, with time
shift equal to 1/3 of the period, the same path on the plane", which is exactly
the relabelling this note uses.

**The braid and the knot.** The figure-eight knot complement is a once-punctured
torus bundle with monodromy `RL = [[2,1],[1,1]]`, standard. `B_3` modulo its
center is `PSL(2,Z)` with `sigma_1 -> [[1,1],[0,1]]` and
`sigma_2 -> [[1,0],[-1,1]]`, so `sigma_1 sigma_2^-1 -> [[2,1],[1,1]]`. That the
orbit's trajectory braid is this map is why Moore's classification puts it where
it does.

**Duistermaat.** Reported by Carreira-Perpinan and Williams as a private
communication, not a paper of his own: three equal Gaussians at the vertices of
an equilateral triangle give four modes for a window of variances. Their own
statement of the conjecture it refutes is that a homoscedastic mixture of M
components in more than one dimension has at most M modes.

**Braiding Vineyards**, arXiv:2504.11203, SODA 2026. Their definition of
monodromy is that following the family for a period "permutes the set of points
in a non-trivial way", which is the definition used here, and what is verified is
that the diagram POINTS exchange, not merely the labels. Their theorem is that
any knot or link CAN appear, proved by constructing a space from the given link.
So the relation to this note is realizability against occurrence: they show a
knot can be made to appear, and here one appears without being made to. The note
therefore claims no novelty for a knot appearing in a vineyard, only for this one
turning up in a solution of Newton's equations.

**The A_1^2/A_1^2 theorem**, *The Singular Source of Vineyard Monodromy*,
arXiv:2607.01046. Not relied on here. The mechanism in this note is transport
along a forced loop, and no claim is made that the two coincide.
