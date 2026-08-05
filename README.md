# Vineyard monodromy in Newton's equations

The figure-eight solution of the three-body problem carries a closed vineyard
with monodromy of order 3. The write-up is one page, `note/note.pdf`.

Sum a radial kernel at each body and follow the degree-zero persistence of its
superlevel sets over a third of the period. The diagram returns to itself with
its three points cyclically permuted. The braid is `(s2^-1 s1)^2` and its closure
is the figure-eight knot. The walls are the isosceles configurations, twelve per
period, and that has a proof rather than only a measurement.

## Layout

```
note/             the write-up, source and PDF
figs/             the figure
mono/             the library
scripts/          everything that produces a number
out/verify/       a log behind every number quoted anywhere
notes/RESULTS.md  the same findings at length
notes/RECORD.md   the working record: what was tried and what failed
```

`data/` and `papers/` are not tracked. The first is too large for a repository
and the second is not mine to redistribute. `scripts/fetch_walnut.py` and
`scripts/fetch_rdc.py` pull the measured inputs back from Zenodo.

## Reproducing

Two scripts check the write-up as a whole:

```
python3 scripts/review_note.py    re-derives every claim from the initial
                                  conditions, touching no saved log
python3 scripts/audit_note.py     every number in the note against its log
```

The rest, by what they establish:

```
validate.py             the machinery against ground truth
verify_headline.py      orbit, braid, Alexander polynomial
rigor.py                the adversarial checks, including primitivity of T/3
crosscheck_gudhi.py     both core routines against GUDHI
crosscheck_braid.py     braid extraction and closure, by other implementations

what_it_means.py        the walls are the isosceles configurations
why_twelfths.py         why there is a wall at every T/12
wall_proof.py           the proof that the walls are those and nothing else
concavity_lemma.py      the hypothesis that proof needs, as an explicit inequality
basin_lemma.py          the weighted-centroid identity behind it
two_questions.py        that hypothesis for kernels other than the Gaussian, and
                        a search for a counterexample to the spanning tree

kernel_diagrams.py      kernel independence, on actual diagrams
frame_test.py           the braid is an invariant of the shape, not the motion
deaths_and_symmetry.py  what the deaths encode
mode_count.py           the mode count, against Duistermaat's counterexample
complete_invariant.py   what the diagram retains of the shape
completeness_test.py    and whether it determines it

trajectory_braid.py     the trajectory braid, its conjugacy class and entropy
knot_meaning_and_N.py   the pseudo-Anosov reading, and what happens at N > 3
braid_square_closed.py  the squaring claim
rerun_transport.py      assignment tracking against transport tracking
rotating_wave_test2.py  why symmetry degenerates monodromy
modulated_wave.py       and why transport produces it

fig_choreography.py     rebuild the figure
```

The rest is the measured-data work, which produced negatives rather than a
second example: `hunt_ct*`, `certify_ct*`, `ct*_robust`, `rdc_analyze`,
`fetch_*`, `ep_mechanism2`, `certify_hero`. They are here because
`notes/RECORD.md` cites their results, and a negative that cannot be re-run is
not a negative. The one that matters is D10, the sampling criterion: making the
two parameters physical makes them sampled, and between adjacent samples the
field can move by less than the measurement noise, which manufactures monodromy
that is not there.

## The audit

`audit_note.py` is a whitelist with arithmetic. Every quantity the note asserts
is declared against the logged value it came from, and the rounding or bound is
checked numerically. It also scans the note for any number with no log behind it.
It passes 29 of 29 with none undeclared. Run it after any change to
`note/note.tex`.

An earlier version pasted whole logs into one blob and asked whether each number
appeared somewhere in it. That is too weak: short strings match incidentally, and
it let the note claim agreement to `1e-14` when the worst measured value was
`9.37e-13`.
