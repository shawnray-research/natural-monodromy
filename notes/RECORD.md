# Falsification log: hunting natural vineyard monodromy

Every negative is recorded with *which necessary ingredient failed*. The three
ingredients (from arXiv:2607.01046 and Braiding Vineyards):

1. **L**: a genuine closed loop in the parameter space, not imposed by us.
2. **S**: a degeneracy the loop encircles, meaning a birth-birth wall *and* a
   death-death wall crossed, i.e. an A₁²/A₁²-type codim-2 point.
3. **P**: the elder-rule pairing must change **exactly twice** (0 or 4 give the identity).

---

## N1. Aharonov-Bohm ring, electron density over one flux quantum

*Setup.* 1D ring conductor with static disorder potential `V(θ)`, threaded by flux
`α = Φ/Φ₀`. `H(α) = ½(-i∂_θ - α)² + V(θ)`. Field is `ρ(θ;α) = Σ_occ |ψ_n|²`.
Loop is `α ∈ [0,1)`, forced by gauge invariance.

*Result.* **0 of 48** (potential, filling) combinations showed monodromy.
Periodicity verified to about `1e-11`, so the loop is exactly closed and the setup is
sound.

*Which ingredient failed: **S**.* The loop is impeccable (**L** holds) but the field never
reorganizes. The density is *pinned to the static potential wells*: peaks of `ρ` sit in
minima of `V` and neither move much nor reorder as flux threads the ring. The
Aharonov-Bohm modulation of the density is small compared to the peak separations, so
the vineyard never crosses a birth-birth wall, let alone a death-death wall too.

*Lesson (drives the next step).* A genuine closed loop is necessary but nowhere near
sufficient. **The field must undergo large-amplitude reorganization around the loop.**
A one-parameter knob that only perturbs a field pinned by a fixed potential cannot do
it. What *does* force large reorganization is proximity to a **degeneracy**: around a
conical intersection the eigenstate rotates within the degenerate subspace, so the
density morphs from `|φ₁|²` to `|φ₂|²` and back, an O(1) reorganization, exactly the
regime where both walls can be crossed. This is Edelsbrunner's hunch sharpened: not
"quantum density on a loop" but "quantum density on a loop **encircling a degeneracy**".

---

## N2. Quantum density encircling a diabolical point (conical intersection)

*Setup.* Two-parameter quantum ring `H(u,v) = -½∂²_θ + V₀ + uV₁ + vV₂`, real symmetric,
so eigenvalue degeneracies are codimension 2 and isolated. Located true diabolical points
(gap driven to between `0.0e+00` and `1e-6`). Field is `|ψ_n(θ)|²`, a measurable
probability density. Loop is a circuit in the `(u,v)` control plane, encircling versus not
encircling.

*Result.* The predicted mechanism is **confirmed**: encircling gives
`reorganization = 2.00`, with the density becoming essentially disjoint from its starting
state as it morphs from `|φ_a|²` to `|φ_b|²` and back, the Berry-phase eigenvector
rotation. Not encircling gives `0.02` to `0.57`. But **monodromy order = 1** in every
case.

*Which ingredient failed: **P**, by a mechanism not previously articulated.*
Diagnostics: `survivors = 1, diagonal = 4`. Almost every persistence point is
**annihilated into the diagonal** during the circuit and a fresh set is born. A vine that
does not survive one period cannot permute with anything; with a single survivor there is
nothing left to permute.

*Lesson, corrected after the successful scan below.* The first reading was
"reorganization must be intermediate". The scan that followed refutes that: certified
monodromy points in the *same* quantum system have `reorg` as small as 0.003. The correct
lesson is sharper and more useful:

**The codimension-2 point you must encircle is a degeneracy of the FIELD'S CRITICAL
VALUES, not a degeneracy of the Hamiltonian's spectrum.**

- N1 failed because the parameter space was **one**-dimensional: there is no codim-2
  point to encircle at all, so nothing forces a permutation.
- N2 failed because a conical intersection is the wrong codim-2 point. It is a
  degeneracy of *eigenvalues*, which rotates the eigenvector and therefore *replaces*
  the density wholesale (`reorg = 2.00`, `survivors = 1`, `diagonal = 4`), so features are
  annihilated into the diagonal rather than permuted.
- What is needed is a point where a **birth-birth wall** (two minima of the field exchange
  value) crosses a **death-death wall** (two maxima exchange value): the exact analogue of
  A₁²/A₁², in the control plane rather than the observation plane.

This corrects the natural intuition, Edelsbrunner's included, that conical intersections
are where quantum monodromy should live. They are the canonical monodromy generator for
*eigenvectors*; they are the wrong object for the persistence diagram of a *density*.

---

## N3. PHT over the direction circle, as a source of permutation monodromy

*Setup.* Degree-0 PHT of a filled planar region over the whole direction circle. The
attraction was that the loop is **forced**: it is the transform's own parameter space, so
there is no center, no radius, and no point to search for first. Persistence computed
exactly from boundary tangencies, validated on hand-computable controls (a convex region
gives exactly one diagram point at every direction; a tilted arch gives two births and one
merge).

*Result.* A spiral starts with 5 diagram points and exactly **1** survives the full
circuit, the essential class. Order 1. Same for every number of turns from 1.5 to 5.

*Which ingredient failed: none of L, S or P. The target was misdefined.*
Arya et al.'s "non-trivial geometric monodromy" for the spiral means *no global sections
exist*, and their text says so explicitly: the non-essential classes "cannot be followed
over the whole circle". Nothing permutes. Chambers et al.'s "monodromy of order k" is a
permutation. These are different phenomena, and the brief's litmus test, *reproduce the
order-2 spiral result*, asks for something that does not exist.

*What my run actually did.* It reproduced Arya's obstruction exactly: 4 of 5 classes
undefinable over the circle. That is the correct answer, not a pipeline failure.

*Lesson, and it is structural.* For a connected region there is one global component at
every direction, and every other component is an arm the sweeping line catches and then
absorbs. Only the essential vine is defined over the whole circle, so the direction circle
**cannot** produce permutation monodromy for connected planar regions. Permutation
monodromy needs the opposite regime: a loop small enough that the vines survive it,
encircling a codimension-2 point where the pairing turns over. See `notes/RECORD.md`.

---

## N4. Integral-field spectroscopy: measured spectra over a two-parameter sky plane

*Setup.* MaNGA 127-fibre IFU (`manga-7443-12703`), 74 x 74 spaxels, a full spectrum at
every spaxel, so the parameter space is physical position on the sky and is dense for free.
Field = the measured spectrum, compactified to the circle. Loops of radius 0.7 spaxels,
small enough that features should survive, which N3 showed is the required regime.

Two things made this look like the right dataset, and both are genuine advantages over the
coastlines:

* **The smoothing scale is not a free parameter.** Structure narrower than the instrumental
  line spread function cannot be real, so smoothing to the LSF (sigma = 1.02 px, shipped in
  SPECRES) is the physically correct operation, not a choice.
* **The data carries its own error bars.** IVAR gives the variance of every flux value, so
  "is this feature real" is measurable rather than arguable.

*Result.* **0 of 484 loops returned a well-defined answer**, at a 3-sigma cut and at a
5-sigma cut alike. Cause in every case: the number of significant features is not constant
around the loop.

*First diagnosis, wrong.* The initial scan filtered by persistence relative to the largest
persistence in the slice, which has no units. That gave 22 to 43 critical points per
spectrum and 15 to 24 pairing changes per loop, and I read it as noise domination.

*Correct diagnosis.* Using the error bars instead, features surviving a 5-sigma cut are not
noise: 6-sigma excursions in a 110-element window are impossible by chance. They are real
spectral structure. The count still varies, and it varies at **every** brightness and
**every** threshold:

```
S/N > 10, 8-sigma cut : counts 3..8
S/N > 20, 8-sigma cut : counts 4..8
S/N > 30, 8-sigma cut : counts 4..8
S/N > 40, 8-sigma cut : counts 4..8
```

*Which ingredient failed: **the survival requirement**, and for a physical reason rather
than an instrumental one.* Different parts of a galaxy genuinely have different numbers of
detectable lines, because ionization state, velocity dispersion blending neighboring lines,
and continuum shape all change with position. The spectra are not different amplitudes of
one structure; they are different structures. Raising S/N cannot fix that, because it is not
noise.

*Lesson, and it is a new requirement.* **D7: the field must present the SAME set of
features everywhere on the loop, with only their amplitudes varying.** Permutation
monodromy is a statement about features exchanging places, so the features have to exist
throughout in order to exchange anything. A system whose feature *inventory* changes across
the parameter space can never supply it, no matter how good the data is.

*Where that points.* A measurement in which one fixed set of spectral lines is present at
every pixel, always far above the noise, while their shapes and relative amplitudes vary
strongly. Solar spectropolarimetry is the clean case: Hinode SOT/SP records the same Fe I
630.15 and 630.25 nm pair at every pixel of a raster at S/N of order 10^3, and Zeeman
splitting changes the profile shape drastically between umbra, penumbra and quiet Sun. Same
inventory everywhere, large amplitude variation, dense two-dimensional raster. It also
satisfies the "somebody would be wrong" clause, since inversion codes fit these profiles
pixel by pixel and profile-component identity across a raster is a live difficulty.

---

## N5. Hinode SOT/SP spectropolarimetry of a solar active region

*Setup.* Chosen because N4 named exactly what was missing. The same Fe I 630.15 / 630.25 nm
pair is recorded at every pixel of the raster, always far above the noise, while the profile
shape varies strongly with the magnetic field. Level 1 Stokes profiles, 81 slit positions x
512 along-slit x 112 wavelengths, from the strong-field patch of the 2007-02-10 15:00:08 map
(field up to 4054 G, continuum 0.49 to 1.95 Ic). Noise measured directly from the flat
continuum edge of each profile: **S/N ~ 534 per wavelength point**.

*Result: D7 is satisfied for the first time.* With a 12-sigma significance cut, 96 per cent
of pixels carry exactly the same number of features, and **84 per cent of loops return a
well-defined answer**, against 0 of 484 for the IFU. The machinery finally has something it
can track.

*But zero nontrivial monodromy*, across six configurations:

```
Stokes I, 12 sigma, r=1.5 : 84% well defined, 0 monodromy
Stokes I,  8 sigma, r=1.5 : 60%              0
Stokes I,  8 sigma, r=3.0 : 45%              0
Stokes I,  5 sigma, r=2.0 : 25%              0
Stokes V, 12 sigma, r=1.5 : 19%              0
Stokes V,  5 sigma, r=2.0 :  1%              0
```

*Which ingredient failed: **D3**, richness, and it failed because of D7.* The stable
inventory is exactly **3** diagram points, one of which is the essential class produced by
compactification, leaving **2** real vines. The mechanism needs at least 3 births and 3
deaths with an elder birth outside the swapping pair: the worked example in
arXiv:2607.01046 needs b0, b1, b2 and D_I, D_II, D_III. Two real vines cannot supply it, for
the same reason the ellipse center cannot: too few features and the pairing simply flips
back.

*Lesson, and it is the sharpest constraint yet.* **Stability and richness pull in opposite
directions.** The significance cut that stabilizes the feature inventory is the same cut
that strips the inventory below the minimum complexity monodromy requires:

```
12 sigma : 96% of pixels stable,  3 features
 8 sigma : 86%                    3
 5 sigma : 62%                    3 to 5
 3 sigma : 27%                    3 to 11
```

There is no threshold at which both hold. **D9: the field must carry at least four features
that are simultaneously far above the noise and present everywhere on the loop.** For a
measured one-dimensional spectrum that demands an intrinsically dense, strong-lined
spectrum, not merely a high signal to noise ratio on a two-line profile. Two Fe I lines at
S/N 534 are not enough; what is needed is many strong lines, for instance a molecular band
head, a rich absorption forest, or a spectral range wide enough to hold several strong
transitions at once.

---

## N6. Sunspot umbra: one candidate found, and rejected by its own certification

*Setup.* Same instrument as N5, better region. Hinode SOT/SP map 2007-05-03 13:15:07, a
proper sunspot with an umbra down to 0.13 Ic. 105 slit positions x 512 along-slit
downloaded around the umbra; 3202 pixels below 0.5 Ic. Inside the dark core the kG field
splits each Fe I line into resolved Zeeman components, so the feature inventory rises:

```
pixels below 0.25 Ic, 8-sigma cut : counts {3:24, 4:555, 5:144, 6:53, 7:3}, mode 4 (71%)
```

**This is the first region satisfying D7 and D9 at once**, which is why it was worth
scanning.

*Result.* Five configurations gave zero monodromy. One configuration
(5-sigma, dark < 0.25 Ic, radius 1.5) produced a single hit at slit 66, along-slit 212:
order 2 over 6 points. But in that same configuration only **5 of 274 loops were even well
defined**, against 26 per cent at 8-sigma and 49 per cent at 12-sigma. A lone hit inside a
2 per cent island is the profile of an artifact, so it was certified rather than reported.

*Certification, and it fails.*

```
[1] threshold      order 2 ONLY at ksig=5; order 1 at ksig=10; every other cut fails
[2] shrink loop    radius 0.6, 0.9, 1.2 -> order 1;  only radius 1.5 -> order 2
[3] displace       every displaced loop leaves the core: toggle cannot even be run
[4] steps          order 2 at 24, 48, 64, 96, 144 steps, but order 1 at 32
```

Test [2] is decisive. Genuine monodromy persists as the loop shrinks onto the singularity:
the Japan generator held order 2 across three and a half orders of magnitude of radius.
This one appears at exactly one radius and vanishes on either side. **Rejected.**

*Which ingredient failed: none of the data requirements; the candidate was not real.*
D7 and D9 were met for the first time and the machinery ran. What the umbra does not
supply is a configuration where the pairing genuinely turns over twice, and the one
apparent case was a tracking artifact at a marginal threshold.

*Value of this negative.* It is the certification chain catching a false positive that a
less careful pass would have reported, in a project that has already had to retract an
overclaim once. The protocol works in the direction that matters.

---

## N7. Measured climate field on a latitude circle, over the annual cycle

*Setup.* The strongest structural setup yet, on two counts. The domain is a **genuine S^1**,
the longitude circle, needing no compactification: NCEP daily sea-level-pressure climatology,
365 days x 73 latitudes x 144 longitudes, longitude wrapping exactly to 360. And the loop is
**astronomically forced**, the annual cycle, closing exactly because a day-of-year
climatology is periodic by construction.

*Result.* The critical-point count on a latitude circle runs from **2 to 36** over the year at
every latitude tested. Filtering by persistence never stabilises it:

```
lat    thr=0.10      thr=0.20      thr=0.30      thr=0.40
70.0   2-7  (53%)    2-5  (60%)    1-4  (79%)    1-3  (79%)
60.0   2-9  (41%)    1-6  (73%)    1-5  (84%)    1-4  (86%)
40.0   3-12 (26%)    2-8  (58%)    1-6  (63%)    1-4  (39%)
```

The best case leaves 1 to 4 features at 86 per cent stability, and by then the inventory is
too thin for the mechanism, which needs at least three minima and three maxima.

*Which ingredient failed: **D7**, and this time for a PHYSICAL reason.* In the galaxy data
(N4) the inventory changed because sensitivity varied across the field. Here it changes
because the atmosphere itself changes: winter has vigorous planetary-wave activity and summer
has little, so the number of pressure highs and lows on a latitude circle genuinely differs
between January and July. No threshold can repair that, because it is not noise.

*Lesson.* The annual cycle is a perfect loop attached to a field whose feature inventory is
seasonal. **A forced loop is worth nothing if the thing going round it is not the same kind of
object at both ends.** Across N4, N5 and N7 the pattern is now unambiguous: for measured
*fields*, stability and richness are in direct conflict, and the conflict is sometimes
instrumental and sometimes physical but always present.

*What this leaves.* The one setting where D7 holds automatically is the one the choreography
work stumbled into: when the features are tied to **persistent physical objects** rather than
to level sets of a field. A fixed set of tracked particles cannot change its inventory. So the
remaining experimental target is measured particle tracking in a periodically driven system,
where the objects persist by construction and the kernel bandwidth can be fixed by the imaging
point-spread function rather than chosen, which is the same argument that made the Hinode
smoothing scale legitimate.

---

## N8. Measured atmospheric field, SMALL loops in a two-parameter physical control space

*Setup.* This corrects the structural error in N7. Taking the whole annual cycle forces the
field to reorganise completely, which is what destroys the inventory. The coastline result
worked because its loops were *small*: the field barely changes, D7 holds for free, and the
monodromy comes from encircling a codimension-2 point rather than from reorganisation. So:

  - domain: the longitude circle at fixed latitude, a genuine S^1 needing no compactification
  - field: geopotential height, NCEP daily climatology, 365 days x 17 pressure levels
  - parameters: **(day of year, pressure level)**, both physical, day periodic, level
    interpolated in log-pressure
  - loops: small, 3 days by 0.08 in log-p
  - the field is zonally truncated to planetary wavenumbers 1..K, which is the standard
    definition of the planetary-wave field rather than an arbitrary smoothing

*Small loops did fix D7*, partially: raw 2.5-degree data gave constant counts on only 4 per
cent of loops, and truncation to K = 3 or 4 raised that to 50 to 60 per cent.

*Result.* 250 raw hits across eight latitudes and three truncations. **Zero survive
certification.**

```
raw hits                                        250
  of which claim order > 1 with ZERO pairing
  changes, which is impossible                   37   <- pure tracking artifact
survive C2 (exactly two pairing changes)          36
survive boundary guard + shrink + toggle + C3      0
```

*Which ingredient failed: the candidates were not real.* The 37 hits with zero pairing
changes are the diagnostic: a permutation with no elder-rule change cannot happen, so the
Hungarian matcher was swapping near-coincident diagram points. Many others sat at 11 hPa, the
top of the interpolation range, and were boundary artifacts. Nothing survived shrinking the
loop, which is the test that matters.

*Lesson.* Small loops are necessary but not sufficient. They buy D7, but the codimension-2
crossing still has to be there, and in a measured field at this resolution the genuine
crossings are outnumbered by matching artifacts among near-coincident points. The 36 hits
that passed C2 and then failed the shrink test are the honest measure of how badly a
C2-only filter would have misled.

---

## N9. Measured cone-beam X-ray projections, and the sampling criterion

*Setup, chosen from the synthesis rather than by browsing for data.* The
recurring obstruction was that measured FIELDS trade inventory stability against
richness. The stated escape was a field whose features are tied to persistent
OBJECTS. An X-ray projection is exactly that: each dense structure contributes a
peak at `s_i(theta) = x_i cos theta + y_i sin theta`, the structures are
permanent, and the peaks cross as the gantry turns.

  data        University of Helsinki cone-beam scan of a walnut, 721 projections
              at 0.5 deg over 360 deg, 50 um detector pixels, Zenodo 4279549.
              Raw detector counts, not a reconstruction. A 121 x 250 x 1600 patch
              was pulled out of the 4.25 GB archive by reading each deflate
              stream only as far as the wanted rows
  domain      detector coordinate, one-point compactified. Exact here rather than
              a device: the attenuation really is zero in air at both ends
  parameters  (gantry angle, detector row), both physical, both swept by the
              machine
  field       -log(I / I_air), air taken per row from the panel margins

*This setting finally cleared D9.* At 15 detector pixels of smoothing, 285 um at
the object, the persistence spectrum has a 7.5x gap: six features from 1.09 down
to 0.027, then nothing above 0.004. No previous measured field managed that.
Hinode gave three features. D7 reached 47 per cent of small loops.

*And it produced candidates with the correct structure.* Two localized to boxes
of 1e-8, held order 2 over four decades of loop radius with exactly two pairing
changes at every radius, were stable from 40 to 640 loop samples, and survived
replacing the bilinear interpolant with a Catmull-Rom bicubic one, moving 0.04
index units.

*They were still not real, and the reason is measurable.*

```
change in the profile per gantry step (0.5 deg)   0.00529
change in the profile per row step   (100 um)     0.00338
air noise of the same smoothed profile            0.00390
```

Consecutive samples differ by LESS than the noise across rows and by 1.4x the
noise in angle. A codimension-2 point localized to 1e-8 inside a cell whose
corners differ by less than the measurement error is a description of the noise.
That is exactly why those candidates survived a change of interpolant, which
re-reads the same noisy samples, and why they died under every coarsening:
subsampling rows lost them at both offsets, and blurring the panel isotropically
erased them.

*Second pass on a signal-dominated grid.* Strides of 2.0 deg and 800 um put both
parameters at 3.3 to 3.7x the noise per step. That gave three candidates that
localized to 7.5e-9 with exactly two pairing changes AND the correct wall
structure. They were then tested against eighteen different subsets of the same
scan, changing which projections and which rows are used:

```
best candidate (theta 26.33 deg, row 714)  recovered in  4 of 18 subsets
                                           order 2 and BB/DD walls every time
                                           spread 0.76 deg, 157 um
second candidate (theta 7.71 deg, row 957) recovered in  0 of 18 subsets
```

**Not certified.**

*Two things worth keeping.*

**D10, the sampling criterion.** For a measured two-parameter family, a
codimension-2 point localized inside one grid cell means something only if the
field changes by more than the measurement noise across that cell. Below that,
localization converges onto noise, and it does so with every appearance of
success: tiny boxes, stable orders over decades of radius, coincidences to 1e-11.
This retroactively explains the sunspot candidate and the atmospheric hits.

**The wall-type test.** Counting elder-rule pairing changes is NOT sufficient to
identify an A_1^2/A_1^2 point in noisy data. A loop that crosses a single wall
out and back also gives exactly two changes, and composes to the identity. The
crossings must be checked to be one birth-birth and one death-death, which is
done by asking, at each change, whether it is two births or two deaths that
coincide. This is tracker-independent and it is now the sharpest gate available.

*The deeper lesson, which is the point of the whole attempt.* Making the two
parameters physical is what makes an example natural, since the machine sweeps
them rather than the analyst drawing a loop. But physical parameters are
SAMPLED, and sampling noise is a new obstruction that the coastline never faced:
there the two parameters are a point in the plane and the field is evaluated
exactly at any point, with no interpolation and no noise between samples. The
naturalness that was wanted and the certifiability that was achieved pull against
each other for a reason, and D10 is the statement of it.

---

## N10. Exceptional points. Refuted at the mechanism, before any data

Proposed as the "golden" target: track the measured intensity |psi(x,y)|^2 of a
microwave billiard around the experimental loop in (slit width, stub position)
that encircles an exceptional point. The attraction was that encircling an EP
provably interchanges the two coalescing modes, and that the enclosed
degeneracy is attested by physics rather than by our own machinery, so the
singularity certificate would be "inherited, not earned".

This is testable analytically at zero cost, so it was tested before obtaining
any data. EP normal form H = [[0,1],[eps,0]], eps = (s-s0) + i(delta-delta0),
mode psi = phi_1 + lambda phi_2 with lambda = sqrt(eps), observable |psi|^2.

**Result 1: over ONE encirclement the vineyard is not closed.** lambda runs from
+sqrt(r) to -sqrt(r), so the field returns as the OTHER mode's intensity. Measured
across four random complex mode pairs, the field moves by 78 to 95 per cent of
its own maximum. Vineyard monodromy is not defined on a vineyard that does not
close. The mode interchange is holonomy of the eigenvectors, not of a diagram.

**Result 2: over TWO encirclements the vineyard closes and D7 fails.** The
diagram cardinality varies around the loop in 4 of 4 runs, [3,4,5], [4,5,6,7,9],
[4,5], [4,5,6,7,8]. Features are created and destroyed, not permuted. This is
exactly N2: at a spectral degeneracy the features annihilate.

**Result 3: with real mode profiles the diagram cannot even see the swap.** The
two coalescing modes are then mirror images, the field changes by 91 per cent
and the diagram returns to itself to 6.7e-16. Persistence is not injective, and
here that erases the entire effect.

*The structural error, which generalizes.* An exceptional point is a degeneracy
of the SPECTRUM: of the eigenvalues and eigenvectors of an operator. Vineyard
monodromy needs a degeneracy of the CRITICAL VALUES of one scalar field. These
are different loci in the same control plane, and certifying that a loop encloses
the first says nothing about whether it encloses the second. So the step offered
as the strongest, inheriting the certificate from physics, is the one that fails.

*The filter this gives, which is the useful part.* **The objects that the forced
symmetry relabels must be the critical points of a single scalar field, not the
eigenvectors of an operator.** This is why the choreography works: the bodies are
relabelled by x_i(t + T/N) = x_{i+1}(t), and the bodies ARE the maxima of the
mass density. It is why conical intersections (N2) and exceptional points (N10)
both fail. Applied to the proposed fallback portfolio it rejects the optoacoustic
EP3 platform and the coupled-pendulum mode profiles for the same reason, without
computation.

*A third thing called monodromy.* Molecular quantum monodromy, also on that
fallback list, is real and measured, but it is monodromy of the lattice of
quantum numbers in a joint spectrum (Cushman and Duistermaat), not a permutation
of persistence points. Treating it as vineyard monodromy would repeat exactly the
conflation recorded in notes/RECORD.md.

---

## P1. A positive: rotating waves through a fixed envelope

Not a failure. Recorded here because it came out of testing the proposed
"symmetry implies monodromy" theorem, which is false, and the corrected version
is verified. Details in notes/RESULTS.md.

Symmetry-exchanged extrema carry EQUAL values (spread 0.000e+00 for k=2,3,4,5),
so their diagram points coincide and the permutation is invisible. What works is
transport: a rotating wave carried through a stationary envelope has an exactly
closed loop of length T/n, D7 by construction, and monodromy of order n, provided
the envelope is not 2*pi/n-periodic. Verified order n for n=3,4,5.

Separately, and important on its own: frame-to-frame Hungarian matching of
diagram points MISSES this monodromy, reporting order 1 in 11 of 11 cases where
the truth is order n, with the closest pair of diagram points never below
9.6e-03. Vineyards in this class must be tracked by transport, not by assignment.

---

## T1. The tracker audit, and two mechanisms

Prompted by a correct objection: every measured-data negative N4 to N9 was
obtained with `linear_sum_assignment` on (birth, death), the same tracker later
shown to return order 1 on rotating waves whose true monodromy is an n-cycle. If
that tracker is unreliable, the negatives are unreliable.

Resolving it required separating two mechanisms that had been conflated:

- **pairing monodromy**: critical points stationary, the elder-rule pairing
  reassigns around a codimension-2 point. The A_1^2/A_1^2 case. Coastlines,
  quantum ring, and every measured-field scan in this project.
- **transport monodromy**: critical points move and exchange places. The
  choreography and the rotating wave.

A transport tracker following maxima by position is structurally BLIND to
pairing monodromy, because there the critical points do not move at all
(measured transport jump 0.0 at the Japan generator). Applying it to the
measured-field scans was a category error, and the first re-run, which reported
0 nontrivial loops on the walnut and was briefly taken as confirmation, proved
nothing.

The correct instrument for pairing monodromy is the combinatorial pair tracker
in `mono/exact.py`, which follows vines by critical-point identity and resolves
pairing changes by which side preserved its labels, with no assignment on
(birth, death) anywhere.

**Result of the audit.**

```
Japan coastline, both certified generators, exact pair tracker
    order 2, single transposition, stable at loop radius 1e-3 .. 1e-6
    -> the coastline positive SURVIVES an assignment-free tracker

walnut, signal-dominated grid (D10 satisfied, steps 3.7x / 3.3x noise)
    assignment tracker      16 nontrivial of 957 well-defined
    exact pair tracker      71 nontrivial of 3363, ALL order 2, all transpositions
    of the 8 tested: shrink collapses to order 1 in all 8
                     wall types DD,DD or BB,BB or none; never BB and DD
    -> N9 stands, and now means something: the obstruction is not the tracker
```

So assignment gives false negatives on transport monodromy AND false positives
on measured pairing monodromy, and neither failure changes the verdict on the
measured data once the right tracker is used. D10 is a finding rather than a
hypothesis.

Also corrected here: the earlier claim that the vineyard braid is not the
trajectory braid, which compared permutations in two different strand labellings.
At N=3 both braids have permutation [1,2,0] and exponent sum 0, so neither
invariant separates them. The comparison at N=5 has not been redone with a
conjugacy invariant and the claim is withdrawn until it is.

---

## T2. Three corrections closing out the audit

**Japan wall-type, resolved.** An earlier run reported 0 birth-order and 0
death-order changes at the Japan generators, which would have meant the strongest
measured candidate was not an A_1^2/A_1^2 point. It was an artifact of the curve
construction: that run built the outline with resample_closed + smooth_closed and
resolved only 4 diagram points, where certify_hero uses
`condition(P, n=1400, smooth_frac=0.010)` and resolves 8. On the correct curve
both generators show **2 pairing changes, walls BB and DD**, at radii 1e-3, 1e-4
and 1e-5. The diagnostic is therefore discriminating, not under-powered: it
passes the coastline and rejects the walnut.

**Full period, resampled.** The previous `full_period_and_writhe.log` reported 7
crossings and order 2 under a heading asserting triviality. Two faults: the index
range ran off the integration array so it covered 0.63 T, and the loop was
sampled at the same K as a loop three times shorter. Corrected, from a
non-degenerate basepoint: **12 crossings, (s2^-1 s1)^6, perm [0,1,2], order 1**,
identical at K = 4000 and 12000, closure mismatch 0.00e+00.

**Basepoint dependence, new.** The figure-eight has six instants per period,
t/T = 0, 1/6, 1/3, 1/2, 2/3, 5/6, at which two bodies are mirror images and carry
exactly equal density values, so two diagram points coincide at 0.000e+00. A T/3
loop based at one of them starts from a degenerate diagram and returns order 2
with 3 crossings. Of the eleven basepoints tested, exactly one, t/T = 0, lies on a
symmetric instant, and it is precisely the one that fails. The other ten avoid the
symmetric instants and all return order 3, word (s2^-1 s1)^2 or its conjugate:
ten of ten among admissible basepoints, with the lone exception predicted rather
than anomalous. Passing through a
coincidence mid-loop is fine, since that is what a vine crossing is; only the
basepoint matters. Any statement of the result must name a basepoint as well as
a bandwidth.

---

## T3. Figure 3 panels (d), (e), (f) were drawn wrong

Found by inspection, 29 July. The braid WORDS were always correct and are
verified independently; only the pictures were wrong.

**Bug 1, wrong initial strand order.** `braid_word` records crossings relative to
the order `argsort(births[0])`, which is `[1,2,0]` for the vineyard and `[2,0,1]`
for the trajectories. `draw_braid` assumed `list(range(3))`. With the wrong
order, two of the four vineyard crossings came out NON-ADJACENT, drawn as two
strands swapping straight through a third. That is not an Artin generator and not
a braid diagram. Fixed by adding `mono.braid.initial_order` and passing it to the
drawing routines, with an assertion that every drawn crossing is adjacent.

**Bug 2, the closure arcs.** Each arc was drawn as two semicircles joined by a
horizontal segment across the top, which ran straight through the braid, and the
semicircles were centered on the strand height rather than attached to the strand
ends. Replaced by arcs that leave the strand end, route around the OUTSIDE of the
diagram, nested so they cross neither each other nor the braid, and return to the
same height on the left.

**Panel (c) did not show what its caption claimed.** It said "closed vineyard"
while showing only the three starting points, so neither the closure nor the
permutation was visible. Now the three starting diagram points carry dashed
vertical guides, starts are filled and ends are open circles, and the reader can
see each vine land on a different start point. Verified numerically: vine i ends
on start point [1, 2, 0] with maximum mismatch 2.88e-08.

The large excursions of the vines to death = 0 in panel (c) are real, not an
artifact: for a density decaying to zero the essential class genuinely dies at 0,
and a vine moves onto that level exactly when its body becomes the global
maximum. Those excursions are the birth-birth wall crossings.

## T3b. Panels (a) and (b), checked properly

(a) verified: the drawn segment is exactly T/3 (4000 of 12000 steps); every
position of bodies 1 and 2 lies on body 0's orbit to 5.06e-09 and 7.73e-11, so
drawing body 0's full orbit as the single grey curve is correct; the three T/3
segments together cover the whole curve to 4.30e-08.

(b) caption was wrong twice, though the picture was right. The three maxima are
NEAR the bodies, not at them, displaced by 2.89e-02, 1.58e-04 and 2.90e-02
because neighboring bumps pull them. And there are TWO saddles, not one per
pair: bodies 1 and 2 have no saddle between them. Both are confirmed by
#max - #saddle = 3 - 2 = 1, the Morse count for a density decaying to zero, and
by the merge tree using exactly two finite deaths plus one essential class.
Caption now reads "its three maxima and two saddles".

Colour convention across panels, checked and intentional: (a) and (b) index by
body; (c) plots per-body vines; (d), (e), (f) color each strand by its vine
label, and vine label i IS body i because braid_word is fed arrays indexed by
body. Panels (d) and (e) start with different colors on top because they order
strands by different quantities, x-coordinate and birth value respectively. That
is correct, not an inconsistency. Panels (e) and (f) share one initial order.

## T3c. The saddle count is NOT constant, and why it does not matter

Checked after asking whether panel (b)'s "two saddles" holds all the way round.
It does not. Sampled at 500 instants, at sigma = 0.20, 0.25 and 0.30, over both
the T/3 loop and the full period:

```
maxima            always 3
saddles           2 OR 3   (3 at roughly a third of instants)
finite deaths     always 2
```

**CORRECTED, see the mode-count section of notes/RESULTS.md.** The saddle count
is always 2, not 2 or 3. The midpoint-seeded search used here double-counts: a
seed at the midpoint of a far-apart pair converges to a saddle belonging to a
different pair, and the dedup tolerance keeps both. An exhaustive search from a
grid gives 3 maxima, 2 saddles and 0 minima at every instant. The conclusion of
this entry, that the H_0 vineyard is unaffected, still holds.

The field genuinely gains a third saddle when the three bodies spread out, and it
arrives together with a central local minimum, so Euler characteristic is
preserved: 3 - 3 + 1 = 1, matching 3 - 2 = 1 in the collinear case. The
transition is an ordinary fold.

**The H_0 vineyard is untouched.** Three maxima merge into one component in
exactly two steps whatever the geometry, so only two saddles are ever merge
saddles and the third never enters the merge tree. That is why the finite-death
count is 2 at every instant and every bandwidth tested, so D7 holds for the
vineyard and the braid is unaffected.

Panel (b) now says "the two merge saddles" rather than "two saddles", which is
true at every instant rather than only at the one drawn, and the script asserts
that the plotted instant really has two so the caption cannot drift.

---

## T4. The squaring claim, closed

"The vineyard braid is not the trajectory braid squared" was asserted, then
withdrawn when its only support turned out to compare permutations written in two
different strand labellings. Now settled, and the answer is that the ORIGINAL
claim was right but for none of the reasons given, and the evidence that
originally suggested the opposite was vacuous.

**Proposition.** For any choreography both braids have underlying permutation a
single N-cycle over T/N. Trajectories: that is the definition. Vineyard: diagram
point i belongs to body i, and at T/N the configuration is identical as a set, so
body i sits at body i+1's former position carrying its former birth and death.

**Corollary.** Squaring would force N-cycle = (N-cycle)^2 under B_N -> S_N. That
square is an N-cycle for odd N and two cycles of length N/2 for even N, so the
relation is FALSE for every even N. Verified for N = 4, 6, 8, 10.

**The N=3 evidence was vacuous.** The square of a 3-cycle is a 3-cycle, so the
permutation invariant is blind at N=3; and both exponent sums are 0, so
e(V) = 2 e(T) reads 0 = 0. Neither available conjugation invariant can detect a
difference at N=3. Confirmed on the real orbit: cycle lengths [3] and [3],
exponent sums 0 and 0.

Nothing is claimed for odd N >= 5 and nothing needs to be.

A methodological note worth keeping: a first attempt tested this on kinematic
choreographies built from random closed curves, N points at equal time offsets
along one curve, which is legitimate since the relabelling is kinematic rather
than dynamical. It was abandoned because the merge trees were unreliable for
widely scattered points: several runs returned a vineyard permutation of cycle
type [1,2] at N=3, which the proposition above says is impossible. The numerics
were measuring themselves. The structural argument needs none of it.


---

# Appendix A: what a measured example would have to satisfy (D1 to D10)

## What the data must look like

The theory needs a scalar field over a parameter space in which a loop can encircle a
codimension-2 point. So the dataset must be, at minimum, three-dimensional:

```
f(x ; u, v)        x       the domain the field lives on   (1D profile, spectrum, or S^1)
                   (u,v)   a genuinely two-dimensional PHYSICAL parameter space
```

## Hard requirements, each learned from a previous failure

| # | requirement | why, and which failure taught it |
|---|---|---|
| D1 | the parameter space is genuinely 2D | N1: a one-parameter knob (AB flux) has no codim-2 point to encircle at all |
| D2 | dense sampling in **both** parameters, roughly 30 x 30 or better | a loop must be resolvable; most experiments scan one knob finely and the other coarsely, so this is the binding constraint |
| D3 | the field has at least 3 local minima and 3 local maxima | the ellipse-center validation: 2 and 2 flips the pairing at all four walls and gives the identity |
| D4 | features survive one circuit | N2: at a conical intersection the density is replaced wholesale, `survivors = 1`, so nothing is left to permute |
| D5 | the measurement is genuinely experimental | this is the whole point; reanalysis and simulation do not count for the headline claim |
| D6 | the loop must be SMALL enough that the features survive it | N3: the direction circle is forced but too big; for a connected region only the essential vine is defined over the whole circle |
| D7 | the field must present the SAME feature inventory everywhere on the loop, varying only in amplitude | N4: in an IFU the number of detectable lines genuinely changes across the galaxy, so 0 of 484 loops were well defined at any brightness or threshold |
| D8 | somebody must actually track these features and be wrong | otherwise the example is true and pointless, which is what sank the coastlines |
| D9 | at least four features simultaneously far above the noise AND present everywhere on the loop | N5: Hinode satisfies D7 at last (84% of loops well defined) but the cut that stabilizes the inventory leaves only 3 points, one of them the essential class |

## Ranked candidate classes

Ranked by (impact) x (chance the data actually satisfies D2).

1. **Spectral imaging**: a full spectrum measured at every pixel, so the parameter space
   *is* physical space and is automatically dense. Sub-cases: astronomical integral-field
   spectroscopy (MaNGA, MUSE), solar spectropolarimetry (Hinode SOT/SP, IRIS), mass
   spectrometry imaging, Raman/IR microscopy, STEM-EELS spectrum imaging.
2. **Solar spectra around a sunspot**: as above, but the loop encircles a genuine physical
   singularity, a magnetic flux tube. The loop is still chosen, but what it encircles is
   not.
3. **Atmospheric field on a latitude circle**: domain is the longitude circle S^1, which is
   a *forced* closed domain; parameters (latitude, pressure level). Planetary-wave
   structure gives many extrema.
4. **ARPES over the Brillouin zone**: highest prestige, because the BZ is a torus and its
   non-contractible loops are forced by the crystal, not chosen. Risk: an energy
   distribution curve may have too few peaks for D3.

## Status

Being worked through in order. Findings and rejections recorded here.

## D10. The parameter steps must beat the measurement noise

Added after N9, and the only requirement here that is about the MEASUREMENT
rather than about the field.

If the two parameters are physical they are sampled by an instrument, and
between adjacent samples the field moves by some amount that must be compared
with the measurement noise. When the step change is below the noise, a
codimension-2 point can be localized inside one cell to arbitrary precision and
will exhibit every sign of being genuine: a box of 1e-8, a stable order over four
decades of loop radius, value coincidences to 1e-11, survival under a change of
interpolant. All of it describes the noise realisation stored in those samples.

The operational form: measure the mean change of the field per unit step in each
parameter, divide by the noise of the same smoothed field measured where there is
no signal, and require at least 3. Then verify by re-deriving the generator from a
DIFFERENT subset of the samples, changing which ones are used, and requiring it to
reappear at the same physical coordinates.

This is what separates a natural example from a certifiable one. The coastline
evaluates its field exactly at any point of the parameter plane, so it has no D10
problem and no natural loop either. A scanner sweeps its parameters physically,
which is what makes the loop natural, and samples them, which is what makes D10
bite.

---

# Appendix B: two different things are called monodromy

This should have been pinned down before any data was touched. It was not, and it is the
reason the PHT route behaved strangely.

## The two notions

**Chambers, Fillmore, Stephenson, Wintraecken (Braiding Vineyards; The Singular Source).**
Monodromy of order `k`: transport the diagram points once around the base loop and read off
the permutation; `k` is its order. Nontrivial means `k >= 2`. This requires vines that
**survive the whole circuit and exchange places**.

**Arya, Giunti, Hickok, Kanari, McGuire, Turner (PHT of star-shaped objects).** *Trivial
geometric monodromy* means a set of global **sections** of the PHT bundle exists. Their
own words, Example 2.24, on the spiral:

> "There is one section encoding the essential connected component of M. However, there
> cannot be a single section encoding the other connected component(s) of M over the whole
> circle because it would not satisfy the conditions to be a function. Furthermore, a union
> of sections cannot encode the other connected component(s) because it would violate
> continuity."

So their spiral is nontrivial because its non-essential classes **cannot be followed around
the circle at all**. Nothing permutes. It is an obstruction to *existence of sections*, not
a permutation.

## What my computation actually shows, and why it is right

Exact filled-region `PHT_0`, validated on hand-computable controls (convex region gives
exactly one diagram point at every direction; a tilted arch gives two births and one merge):

```
spiral, 2.5 turns   diagram sizes [3,4,5,6]   5 points at theta=0,  1 survives the circuit
spiral, 5.0 turns   diagram sizes [5,6,7]     5 points at theta=0,  1 survives the circuit
ellipse (convex)    diagram sizes [1]         1 point,              1 survives
```

The one survivor is the essential class. Four of the five spiral classes cannot be followed
around the circle. **That is precisely Arya's nontrivial geometric monodromy, reproduced.**
The permutation is trivial because there is nothing left to permute, and that is not a
failure of the pipeline.

## Consequence for the search

The brief says "a spiral has non-trivial geometric monodromy in its degree-0 PHT" (Arya's
sense) and then sets the target as "a non-trivial permutation of the off-diagonal points
(order-k monodromy, k >= 2)" (Chambers' sense). The spiral satisfies the first and not the
second, so the proposed litmus test, *reproduce the known order-2 spiral result*, is
testing for something that does not exist. There is no order-2 spiral result to recover.

This kills the PHT-over-the-direction-circle route **as a source of permutation
monodromy**, and it does so for a structural reason rather than a numerical one: for a
connected region, at every direction there is one global component and every other
component is an arm that the sweeping line catches and then absorbs, so only the essential
vine is defined over the whole circle.

**What permutation monodromy actually needs** is a loop small enough that the vines
involved persist all the way around it, encircling a codimension-2 point where the pairing
turns over. That is the Chambers and Fillmore setting, and it is what the radial-transform
work did exhibit genuinely. The defect there was never the mathematics; it was that the
loop was chosen and the example meant nothing.

So the physical target is now sharp:

> **a measured two-parameter physical control space, containing a codimension-2 degeneracy
> of the field's critical values, where the features involved survive a circuit around it,
> and where somebody actually tracks those features and would be wrong.**

`notes/RECORD.md` D1 to D5 already encode the data requirements. This adds D6: the loop
must be small enough for the vines to survive it, which is exactly the opposite of what the
direction circle offers.

## Also fixed along the way

- Rasterizing a region and sweeping pixels is wrong: a diagonal front adds pixels touching
  nothing already added, so a convex region reported 47 diagram points instead of 1.
  Eight-connectivity got it to 7. Computing `H_0` exactly from boundary tangencies gives 1.
- The "region locally above" test used `p + eps*omega` with `eps = 0.0206` against a ribbon
  half-width of `0.025`. Replaced by the sign of the inward normal, which needs no epsilon.
- Labels must be carried on **every** tangency, not only those currently in a pair.
  Otherwise the global extremum switching arms looks like destruction and the essential
  vine is wrongly killed, which is how the bug was caught: the essential class exists at
  every direction, so any computation that lets it die is wrong.
- Tracking a critical point by position still fails where the global extremum genuinely
  jumps between arms of a star (263 samples in one step, at 5 of 359 steps). The vine is
  continuous in the diagram while its location is not, so position-based tracking has a
  hard limit and the elder-rule swap logic has to carry those steps.

---

## N11. Rotating detonation combustor. The right structure, the wrong wavenumber

First test of the rotating-wave prediction against measured data. An RDC is the
theorem's object physically: detonation fronts circulate in an annulus at fixed
speed, the annulus is a genuine S^1 needing no compactification, the loop T/n is
set by the device rather than by an analyst, and the injector pattern is
stationary in the laboratory frame and modulates each front as it passes, which
is exactly the symmetry-breaking envelope the theorem requires.

Data: Bohon et al., Zenodo 18886925, CC BY. Raw high-speed video rather than the
deposited processed version, because that pipeline applies a luminosity
threshold mask and a masked field is not the continuous field persistence needs.
Each raw frame carries the annulus already unwrapped to 50 radial by 360
azimuthal samples, with timestamps at 87.5 kHz and the ring radii.

```
BD0038   5001 frames   dominant azimuthal wavenumber n = 1
CE2029   5000 frames   dominant azimuthal wavenumber n = 1
BD0041                 two COUNTER-rotating waves, a standing pattern
```

**Both co-rotating conditions run a single front.** With n = 1 there is one crest
on the circle, hence one diagram point, hence nothing to permute: the monodromy
is order 1 by construction, not by failure. The harmonic amplitudes, 1.00, 0.62,
0.33, 0.23 and falling, are the spectrum of one sharp pulse, not of several
fronts. The third record is counter-rotating, which is a standing wave with no
net transport, so it can only ever serve as the negative control.

D10 fails independently: the median frame-to-frame change of the azimuthal field
is 0.7 times the noise, on 8-bit data whose maximum value is 12.

**The requirement this sharpens.** The theorem needs n >= 2 crests, since the
permutation is of n diagram points. Most physical rotating waves prefer the
lowest azimuthal wavenumber, because it is the least damped, so n = 1 is the
common case and is useless here. A measured test needs a system that runs
robustly at n >= 2. Wavy vortex flow in Taylor-Couette is the obvious candidate,
where the azimuthal wavenumber is typically 4 to 6 and set by the aspect ratio
rather than by chance.

Not a refutation of anything. The prediction was not tested, because neither
available record has more than one crest.
