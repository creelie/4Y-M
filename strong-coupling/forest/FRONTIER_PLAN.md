# Research frontier and next efforts (2026-11-03)

**Assessment.** The frontier is the passage from controlled good regions
to the full constrained measure, then across scales and into the
continuum. James's ideas sharpen the formal interface but do not supply
the missing dynamical estimates. This note records priorities, not new
results or an experiment execution log.

**Implementation follow-up (2026-11-03).** The latest `ROADMAP.md` section
records the first pass through this plan: the measure bridge and compact
cost-minimiser existence are kernel checked; finite collar/fold, flow and
spectral interfaces are recorded with explicit limits. Current formal
status is in `COUPLING_INTERFACE.md` section 5. The priorities below were
the starting plan; no open Yang--Mills dynamical wall is declared closed.

## 1. Current status

- **A -- proved, stated fixed-lattice regime:** existing strong-coupling
  results remain intact. G1's arithmetic and abstract cost-minimality
  interface are kernel checked. The actual measure-theoretic construction
  remains unfinished (`COUPLING_INTERFACE.md`, `A4_LEAN.md`).
- **B -- proved finite reduction / computation:** the complete finite-box
  interacting correction is nonzero and independently checked, including
  exact small-box arithmetic. Spatial flow and uniform large-volume and
  continuum control remain open (`B_WICK.md`).
- **C/D/I -- proved conditional results / open mixture problem:** rarity
  holds in the stated torus setting. Covariance decay now applies to an
  actual good-window conditional law under explicit calmness and collar
  assumptions (`INDEPENDENT_LABELS.md`). Bad branches, comparison between
  environments and connected island weights remain open. Rarity alone
  does not control the environment covariance
  (`ENVIRONMENT_OBSTRUCTION.md`).
- **D8/E -- open:** multiscale composition and nontrivial continuum
  construction remain the largest walls. Neither a finite-lattice gap
  nor a conditional continuum theorem closes them (`D_STATUS.md`).

**Assessment of James's contribution.** The complete trimmed bridge
integrates an update estimate, derives a discrepancy recurrence and bounds
observables. It still assumes local non-improvement and the one-step
estimate. Connecting it to our checked minimality module is concrete
formalisation work; proving the requisite objects exist remains a separate
obligation. The supplied bridge has not been compiled in this review.
No new E argument was included. See `JAMES_IDEAS_REVIEW.md`.

## 2. Next CPU experiments

All entries below are **sketches / proposed diagnostic computations**, not
claims of uniform bounds. Execution requires a concrete, reviewed plan
specifying the model, observables, validation, resource budget and stopping
criterion.

1. **Bad-branch transport -- highest priority.** Follow every branch of
   the explicit example in `BAD_BRANCHES.md`, then integrate over
   neighbouring variables and perturb the coarse constraint. Measure
   normalisation and expectations of bounded Lipschitz observables.
   **Expected information:** whether integration regularises the apparent
   path-Jacobian singularity, and which variables must move together for
   a viable transport estimate. A singular coordinate density alone does
   not settle this. Validate branch completeness and integration near
   turning points before interpreting response estimates; growth under
   refinement would be a diagnostic, not a proof of divergence.
2. **Two-island connected weights.** In a tractable constrained model,
   compute log Z_11 + log Z_00 - log Z_10 - log Z_01 versus separation,
   retaining hard domains and every branch. State the exterior, coarse
   constraints and label events explicitly and keep them consistent across
   the four partition functions. **Expected information:** a candidate
   mechanism or obstruction for extending the continuous-source estimate
   of `CONNECTED_NORMALIZERS.md` to actual label changes. Stable integration
   and cancellation checks must precede any decay interpretation; a finite
   toy model cannot certify the required Yang--Mills bound.
3. **B's spatial-flow correction.** Extend the complete Wick calculation,
   including nonlinear flow contributions at the required perturbative
   order. Check small cases through independent formulas and alternative
   gauge trees before increasing volume. **Expected information:** whether
   the complete flowed observable has the cancellations and connected
   structure needed for a uniform analytical estimate. A finite coefficient
   is not a continuum extrapolation or a mass estimate.

## 3. Conditional GPU candidate

**Sketch / proposed computation.** The strongest candidate is B's larger
tensor contractions after the spatial-flow calculation passes CPU checks
and profiling establishes a useful GPU advantage. Existing finite boxes
run efficiently on CPUs; repeating them on B200s adds little. Bad-island
integration is a possible later candidate only after a validated
formulation and cost comparison exist.

Before submission, write down the mathematical question, expected output,
validation, budget and stopping criterion; reread and defend the proposal.
Use only the private audited submission wrapper. B200 GPUs only, at most
four GPUs and 32 CPUs in aggregate across parallel jobs; GPU jobs at most
one hour and CPU jobs at most 24 hours. Prefer short pilots and stop if
validation fails or expected information cannot justify further cost.
Operational tooling, job logs and raw data remain outside version control.
No jobs were submitted for this planning note.

## 4. Formalisation and analytical priorities

All entries are **sketch / future work** with undischarged hypotheses.

- **G1 formalisation:** compile and reconcile James's bridge with existing
  modules; connect global cost minimality to its local `hmin`; prove the
  marginal-preservation and finite-integral interfaces, minimiser existence
  and measurable conditional coupling selection. Instantiate the
  specification-dependent one-step estimate. A passing conditional bridge
  does not discharge those inputs.
- **D/I analysis:** use the full private-link submersion or multiple charts
  to handle bad branches; derive connected hard-label normaliser bounds
  and the environment recursion, including collar costs and changed
  conditioning laws.
- **D8 analysis:** control spatial locality and accumulated errors through
  renormalisation. Finite form comparison alone does not supply locality;
  see `FRONTIER_LEMMAS.md` for the weighted-operator and summable-error
  interfaces.
- **E analysis:** identify renormalised, suitably smeared observables with
  nonvanishing spectral weight, physical-energy tightness and convergent
  correlations, together with the reconstruction hypotheses. The
  conditional spectral criterion in `FRONTIER_LEMMAS.md` identifies useful
  inputs, not an existence theorem.

**Priority assessment.** Bad-branch transport and connected island weights
are the immediate analytical targets. The measure-level G1 bridge is the
parallel formalisation target. G1 need not block independent work on D/I,
B or E. (I1), D7, the required D5, (D), D8 and (E) remain open.
