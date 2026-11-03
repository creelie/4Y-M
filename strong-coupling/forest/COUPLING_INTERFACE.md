# The coupling bridge: minimality is enough for the recurrence

**Latest status (2026-11-03, proved / kernel checked).** The trimmed
`WassersteinBridge.lean` supplied by James is now a repository source and
compiles against Lean/Mathlib v4.33.0. Added `local_nonimprovement_of_cost`
and `cost_minimal_comparison` connect the existing abstract cost argument
to actual finite discrepancy integrals. `CouplingMinimizer.lean` proves
that fixed-marginal probability couplings on a compact metric space form
a nonempty compact set and that every continuous real cost attains its
minimum there. Both modules are imported by `YangMills.lean`.
Nothing in this interface proves an SU(2) coupling construction, (I1),
a lattice gap, or a continuum theorem.

## 1. The stationary-coupling route and the minimiser route

**Proved (finite algebra, already used mathematically by the paper).**
Let Gamma be a set of couplings with fixed marginals, and for pi in Gamma
write d(pi,i) for the expected coordinate discrepancy. For positive weights,

    cost(pi) = sum_i w_i d(pi,i).

Suppose pi minimises this cost and a local update at i gives another
coupling pi_i in Gamma. If discrepancies away from i are unchanged, then

    cost(pi_i)-cost(pi) = w_i[d(pi_i,i)-d(pi,i)] >= 0,
    d(pi,i) <= d(pi_i,i).

An expected update estimate

    d(pi_i,i) <= sum_j C_ij d(pi,j) + b_i

therefore supplies the recurrence d(pi)<=C d(pi)+b. Invariance of pi under
each coupled update, or under a random scan, is unnecessary for this
particular proof. The update must preserve the two marginals, which is
a separate requirement from preserving pi itself.

With approximate minimality cost(pi)<=cost(pi_i)+epsilon_i, the precise
additional source is epsilon_i/w_i. This weight-dependent loss matters
for exponentially decaying weights. It cannot be silently replaced by
a spatially uniform error independent of w_i.

## 2. What the formal source proves

**Proved (kernel-checked finite interface, Lean/Mathlib v4.33.0):**

- `cost_update`: the exact change of cost under an update fixing the other
  coordinate discrepancies;
- `recurrence_with_error`: the recurrence with epsilon_i/w_i;
- `recurrence_of_minimality`: the exact recurrence without a stationarity
  assumption;
- `comparison_of_minimality`: the boundary comparison from the existing
  `DobrushinComparison.lean`, with the abstract update and minimality
  hypotheses retained explicitly.

The type of couplings in this source is abstract. In particular, the file
neither defines a probability measure nor assumes hidden instances giving
optimal transport. Its purpose is to expose the hypotheses needed to link
actual couplings to the already verified arithmetic comparison.
The code contains no `sorry` or additional axioms. The source passed
`lake env lean -j 2 YangMills/StrongCoupling/CouplingMinimality.lean`
against the pinned toolchain. The new module is included in `YangMills.lean`. The entire repository default target was not
rebuilt during this pass.

## 3. What remains to instantiate it

**Sketch / open formalisation.** For finite compact metric spin spaces,
the ordinary mathematical argument in the paper selects a cost-minimising
coupling using compactness, then uses measurable conditional couplings to
construct each pi_i with the required marginals and expected discrepancy.
The global continuous-cost minimiser now has a Lean implementation in
`CouplingMinimizer.lean`. It uses weak compactness, continuous marginal
maps and continuity of integration on compact spaces. The measurable
conditional-coupling construction remains open. The older random-scan
bridge is an alternative conditional interface; the implemented trimmed
bridge instead uses cost minimality.

The new minimality interface avoids asking for a coupling simultaneously
invariant under every individually selected update kernel. It does not
avoid the measurable-coupling construction. Approximate optimal updates
can be used, but their error must be tracked or sent to zero in a justified
limit. The environment covariance of `ANNEALED.md` still needs a separate
recursion even after a generic coupling theorem is fully formalised.

## 4. Assessment of the supplied bridge (mathematical review)

The displayed weighted resolvent and source-floor implications have the
stated hypotheses in the right places: positive weights are needed for
entrywise convergence from the column criterion; nonnegative weights
suffice for the scalar comparison; separate subsolutions do not give a
stability estimate. The measure recurrence is conditional on the supplied
updates and coupling, and the random-scan argument only needs finite
discrepancy expectations, rather than invariance under each update.
This paragraph originally assessed the larger supplied text mathematically.
The trimmed minimality bridge has since been kernel checked; the omitted
larger resolvent/random-scan source has not been compiled in this pass.

## 5. Newly checked measure-level chain and remaining inputs

**Proved (kernel checked, Lean/Mathlib v4.33.0).** The new bridge verifies
coupling marginals under updates, integrated update bounds, real weighted
comparison, observable bounds and the additive decay floor. The two added
cost theorems require finite coordinate integrals before and after updates,
positive weights and preservation of other coordinate expectations.
`CouplingMinimizer.exists_minimizer` constructs a minimum of the integral
of any continuous real cost among actual probability couplings; it does
not merely assume a compact abstract coupling type.

**Sketch / open.** Instantiate the continuous cost with the specification's
weighted coordinate discrepancies, identify its ordinary integral with the
finite extended-integral cost used by the bridge, and construct measurable
coupled conditional updates with the desired one-step estimates. Existence
of a global minimiser is distinct from measurable selection of local
optimal couplings as the boundary varies. No such selection, SU(2)
instantiation or infinite-volume theorem is asserted here.

**Validation.** All four modules in this dependency chain were compiled
from source against the pinned toolchain in an isolated build. The compact
minimiser file has only unused-section-variable linter warnings. No `sorry`
or additional axioms were introduced. The entire default target was not
rebuilt. Source snapshots and compiler output are retained privately.
