# Research directions after PR 160

**Assessment, 2026-11-04.** The upstream merge supplies finite good-chart,
collar and formal-comparison ingredients. It does not close uniform bad
regions, the crossover, or a nontrivial continuum theory. This document
records a proposed research program, not a claim that any of these
problems has been solved.

## Bad regions and the shared environment

**Proved finite interface.** `ENVIRONMENT_COMPARISON.md` derives environment
covariance decay from uniform conditional label influences and shows how
the existing hard-island cross-ratio lemma would feed those influences.
The criterion requires all other labels frozen, including arbitrarily large
bad clusters; the present collar estimates cover neither this hypothesis
nor locality of physical-observable conditional means.

**Sketch / next.** Start with a pair of labelled blocks, freeze the
complement, and explicitly disintegrate every block constraint and all
path-link branches. Identify which constraints intersect a proposed collar,
which interface variables are shared, and whether the four sector weights
have one common positive kernel. Work first with a fixed calm annulus and
islands separated by it, tracking normalization and cases of empty sectors.
A successful conditional factorization would be useful even if the uniform
all-label Dobrushin criterion fails.

A second route is a graded expansion: Bertini, Cirillo and Olivieri's
[graded cluster expansion](https://arxiv.org/abs/1110.6028) combines mixing
on a good subset with sufficiently sparse, well-separated bad regions.
**Literature assessment:** its hypotheses and conclusions concern the
models and settings of that paper; applying it to constrained gauge fibres
requires a new verification, including boundary-conditioned probabilities,
random geometry, gauge constraints and uniform estimates on bad regions.
The hidden-phase example in `ENVIRONMENT_OBSTRUCTION.md` explains why
unconditional joint rarity is insufficient by itself.

## Multiscale and crossover

**Sketch / next.** Build a scale-indexed ledger for good-field localization,
bad-region weights, the effective action, source derivatives and errors
under actual blocking. At each step specify the metric, physical units,
background and constraint domains. A single-scale Gaussian Schur bound is
available in `WEIGHTED_SCHUR.md`; the missing result is a nonlinear,
gauge-covariant map preserving suitable locality norms and controlling its
large-field terms. Sum errors only after their propagation constants have
been bounded. State separately the ultraviolet small-coupling induction,
the finite crossover window and the infrared mixing input; do not assume
the effective action stays Wilson or that a strong-coupling expansion
is available at crossover.

**Literature assessment:** Dimock's [small-field](https://arxiv.org/abs/1108.1335)
and [large-field](https://arxiv.org/abs/1212.5562) expositions develop
Balaban-style methods for a three-dimensional scalar model. They are
valuable guides to localization, extraction and large-field bookkeeping,
not four-dimensional Yang--Mills results. Extract an explicit list of
additional gauge-theory hypotheses before importing any theorem.

## Continuum, reconstruction and nontriviality

**Sketch / next.** Choose a countable separating family of gauge-invariant
observables and a physical scale, then establish uniform distributional
bounds and subsequential convergence of all required joint correlations.
Check Euclidean covariance, symmetry, reflection positivity, clustering,
and the growth/regularity assumptions of a specifically stated
reconstruction theorem, including the strengthened hypotheses discussed in
[Osterwalder and Schrader II](https://projecteuclid.org/journals/communications-in-mathematical-physics/volume-42/issue-3/Axioms-for-Euclidean-Greens-functions-II-with-an-Appendix-by/cmp/1103899050.pdf).
Fixed positive gradient-flow time is a useful intermediate observable,
but flow is nonlocal and its removal/local reconstruction needs its own
bounds. Prove a nonzero renormalized observable and a physical-unit gap on
a dense reconstructed subspace, then extend by the actual reconstruction
argument. The positive-time ratio criterion is conditional and does not
provide these inputs (`POSITIVE_TIME_TIGHTNESS.md`). A lattice gap with a
rate that vanishes in physical units is not enough.

## Experiments and formalization

**Assessment / proposed experiments.** First implement exact small binary
label models to check matrix orientation, possible influence row sums and
the environment comparison bound, including the hidden-phase example.
These are algebraic diagnostics, not evidence that the Yang--Mills
hypotheses hold. A later branch-complete constrained SU(2) pilot should
compare explicit interface integrals and conditional odds across frozen
geometries; only launch after defining an independently checkable test,
resource estimate and stop criterion. No GPU is warranted for the present
finite binary identity. Larger B200 work requires a CPU reference and a
credible computational advantage, within the user's resource limits.

**Sketch / formalization.** The finite environment comparison gives a
precise candidate for Lean after its mathematical statement is independently
checked. Independently, finish the measurable conditional-coupling and
cost-integral identifications in the existing G1 bridge. Neither task
substitutes for the gauge-specific uniform estimates.
