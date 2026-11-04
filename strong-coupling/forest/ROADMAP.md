# Roadmap after the A/B/C directions (November 2026)

This file records:
- where every direction of `DIRECTIONS.md` ended up;
- how the directions relate to the two remaining walls, (D) the cube estimate and (E) existence;
- which path we bet on, and why.

The detail behind each line is in the linked files. The split of (D) into pieces D1–D8 is in `D_STATUS.md`.

## 1. Track A: strong coupling, fixed β

| item | outcome | where |
|---|---|---|
| A1 SU(3) blocks | **done.** β_W ≤ 0.77 (was 0.76); α₃(0.77) ≤ 0.99978 and m ≥ 1.1·10⁻⁴ there. Certified cover at r = 77/50 with 2825 boxes; constants on 826 boxes | `../gauge-fixed/blocks/su3/README.md` |
| A2 SU(2) 0.351 → 0.355 | **open.** Earlier Monte Carlo suggested (1+β)/4; this is only first order. The exact zero-field attainable test is [1+I2(β_W)/I1(β_W)]/4; neither the full response supremum nor an improved endpoint is proved | `FRONTIER_LEMMAS.md`, `../gauge-fixed/blocks/feedback/README.md` |
| A3 larger blocks | **done, negative.** Built-up blocks reach a ceiling of about 0.369, against 0.366: not worth proving | `../gauge-fixed/blocks/shapes/README.md` |
| A4 Lean port onto LGT | **scoped.** The critical path is a W₁ (Kantorovich) Dobrushin theorem, which the Lean ecosystem lacks. LGT itself proves only β_W < 0.0104 (total variation) | `A4_LEAN.md` |

Current proved ranges:
- SU(2): β_W ≤ 0.351.
- SU(3): β_W ≤ 0.77.

**Relation to (D): none.** Track A is about fixed coupling, while (D) is posed along β_W → ∞. Pushing the range by
0.01 leaves the gap to the scaling region near β_W ≈ 2.2 essentially unchanged. A4 raises trust in proved results; it
does not add content.

## 2. Track B: the flowed-slice route

| item | outcome | where |
|---|---|---|
| B1 3D next order | **done** | `slice/THREE_D.md` |
| B2 lattice S′ at O(g²) | **done in the free vacuum**, uniformly in a. The vertex bound is proved (Lemma 3). Numerically, the old `leak_norm` drift sits in the normalisation (even) sector, while the hard leak is flat with empirical c ≈ 1.7. The O(g²) vacuum correction is open | `slice/sprime/LATTICE_SPRIME.md` |
| B3 stochastic quantization (CCHS) | **done** (literature). Spacetime-flowed slices are smoothed observables; the 3D measure the route needs does not exist | `slice/CCHS.md` |
| B4 3D unflowed Monte Carlo | not pursued (low priority) | |

**Relation to (D): a parallel route, not a prerequisite.** If it worked, it would go around the cube estimate.
Everything it has produced so far is perturbative. Turning it into a statement about the measure needs the same
nonperturbative control of large fields and of all scales that D7 and D8 ask for. It also needs a 3D measure that
nobody has constructed (B3).

## 3. Track C: the cube estimate

| item | outcome | where |
|---|---|---|
| C1 sharper cube constants | not pursued. Sharper local constants change nothing without the composition step | |
| C2 rare bad fibres | **half done.** "Rare" is proved on the torus by chessboard estimates, with explicit thresholds. "Harmless" cannot come from disagreement percolation at weak coupling | `rare/RARE.md` |
| C3 tree bounds, continuum | only reachable through the 3D route; not started | |

**Relation to (D): part of it.**
- C2's "harmless" half is D7.
- C3 is D8 and beyond.

## 4. (D), the cube estimate

The estimate: D_c(n, α) for block variables, with a supremum over boundary conditions, uniformly along
β_W ∼ b log(1/a). Pieces, from `D_STATUS.md`:

- **Proved:**
  - D1: the Gaussian version at every scale (n ≥ 10⁸²⁷/α);
  - D2: one composition step;
  - D3: decay of the constrained law on a calm fibre, at rate 2.5·10⁻⁴;
  - D6: large fields are rare.
- **Open:**
  - D4: the mass off the chart on fibres with no calm completion;
  - D5: fibre decay ⇒ cube coefficient, with the supremum over boundary conditions;
  - D7: large fields are harmless;
  - D8: the multiscale composition, the coarsest averages and the crossover.

New since v4.29.0 (`D_STATUS.md` §3) concerns D7. The single-link W₁ coefficients are bounded uniformly in the field,
by β_W per entry and 4.5 β_W per row (Lean). So a bad island of n sites amplifies influence by at most
(4.5 β_W)^{4n}, and its rarity beats this once β_W ≥ 1.1·10⁷. If three caveats are met, D7 reduces to D5:
1. an annealed comparison theorem;
2. rarity under an arbitrary exterior;
3. the transport of the island bound to chart coordinates.

## 5. (E), existence

(E) asks for a continuum limit that satisfies the Osterwalder–Schrader axioms and has a gap in physical units
(`OPEN.md` §3–4). The paper proves (D) ⇒ a gap ≥ log(1/α)/(2ℓ) for every continuum limit. What is left after (D) is
the existence of a non-trivial limit with the axioms, which is the Clay problem itself. It is not attacked; it waits
on (D).

## 6. The bet

**Do (D) itself, starting with D5. Keep B as a cheap diagnostic, not a route.**

- **Why not B.** B postpones the wall rather than avoiding it. Its upgrade from perturbation theory to the measure is
  D7 + D8 under another name, and it also needs a 3D measure that does not exist. B is still useful as a fast source
  of yes/no answers on which estimates should be true. It should run when a numerical check would settle a question
  raised by D.
- **Why D5.** It is the next step from a proved theorem (`chartdecay`). It is a finite argument with known constants,
  and the only new work is controlling the chart's boundary terms. Each proved piece of D narrows the problem
  permanently. Through §4, D5 may bring D7 with it.
- **The real risk is D8.** Composition across all scales is Bałaban's territory: that program ran for years without
  reaching a gap, and the proved Gaussian constant (10⁸²⁷) shows the distance. The realistic goal is to make D5 and D7
  theorems and isolate D8 cleanly. Closing (D) is not the expected outcome.

## 7. Order of work

1. **D5 on calm exteriors:** the cube coefficient from Theorem `chartdecay`(ii), with the boundary terms of the chart.
   **Done for links in no average (`D5.md`):** C_G ≤ 3·10¹³ plus an off-chart error. Global calmness is
   volume-limited, so D5 must be local, and then it merges with D7 through the missing local off-chart bound (L).
2. **Annealed Dobrushin comparison:** random coefficients dominated in finite intersections by a Bernoulli field
   (caveat 1). A self-contained probability lemma; with §4 it turns D7 into D5.
   **Done, negative (`ANNEALED.md`):** one-point rarity gives only an O(q) floor (approximate uniqueness). A bootstrap
   lemma reduces the rest to quenched bounds (I1) for laws conditioned on the large-field set. Those are Bałaban
   territory, not a probability lemma.
3. **Local rarity under an arbitrary exterior** (caveat 2), in the constrained-law form (L) of `D5.md` §4.
   **Partly done (`ISLANDS.md` §1):** Proposition L (sketch) gives the off-chart bound in window form with
   volume-independent constants, by adding the private links of the straddling pairs to the free set. Rarity of calm
   windows under an arbitrary exterior in finite-intersection form is still open (only the density bound `localcalm`).
4. **Island bound in chart coordinates** (caveat 3). The section S has row sums ≤ 14 and range one.
   **Done, crude (`ISLANDS.md` §2):** with a one-block collar, an island of n blocks costs at most
   (4.5β_W C_G)^{324n}. With item 3 the thresholds balance at β_W ≈ 10¹⁶⁴, independently of the volume, but only
   given (I1), and only for a gap in lattice units.
5. **D4 and D8:** multiscale analysis, no short route visible.
   **Scoped (`D4_D8.md`):**
   - D4 in window form is the large-field part of (I1).
   - The base case of D8 lands at an effective coupling of about 10–700 (heuristic, one loop). That is far below every
     weak-coupling threshold here (≥ 10³⁸), so the crossover is unavoidable.
6. **Side items, if cheap:** the A2 second-order proof (0.351 → 0.355), and the Lean W₁ Dobrushin theorem (A4).
   Neither bears on (D).
   **Outcome:**
   - A2 is parked. It needs a new Arb certificate (`../gauge-fixed/blocks/feedback/README.md`).
   - The arithmetic half of the Lean W₁ Dobrushin theorem is proved
     (`YangMills/StrongCoupling/DobrushinComparison.lean`, no `sorry`). The coupling half is still open.

### After §7: the calm-data link (`CALM_LINK.md`)

- The work here is merged with the OpenGluon calm-data note. The constant for the first kind of link is C_* ≤ 7.4·10¹². The second kind
  has C_b ≤ 9.2·10¹³ plus a shell term.
- The shell term, averaged over v, reduces to the localised off-chart covariance (sketch, §2).
- A window bootstrap gives that covariance with decay rate α_G/2, down to a floor, at β_W ≳ 10⁶⁴ (§3).
- What is left of the calm-data cube coefficient is one statement, (R1): one-point calm rarity under an arbitrary
  exterior. All three routes to it fail (§4), and it is the large-field half of a renormalisation-group expansion.
  **D is blocked there**, and at D8. Work moves to other paths.

### After the calm-data link: B2 closed, A2 corrected, lead (b) assessed (November 2026)

- **B2.** The `leak_norm` drift is located (numerical, `slice/sprime/LATTICE_SPRIME.md`, "The split measurement"):
  - the cubic leak is entirely hard and flat for τ/a² ≥ 4, with c ≈ 1.7 against N(2);
  - the drift is the two-particle normalisation, which runs like the coupling and cancels in S′.

  What is left is the O(g²) vacuum.
- **A2.** The 4.1% slack in Lemma `blockK` is the direct y₃ response d₃ = C d₁, not the feedback loop (that is 0.15%).
  The supremum over test functions lies between 0.2719 and 0.2849. So the "3.4% margin" for β_W = 0.355 is not
  established (`../gauge-fixed/blocks/feedback/README.md`, Correction).
- **Lead (b), Dobrushin–Pechersky.** Its Lyapunov condition needs field strength to be attenuated from a block's
  boundary inward by a factor Δ^χ. A constant abelian-embedded field shows no attenuation at the linearised level. The
  same example shows that (R1) cannot hold under *every* exterior without a nonabelian mechanism, so (R1) should be
  posed averaged over exteriors (`CALM_LINK.md` §4, heuristic).

### The large-field wall (November 2026, `LARGE_FIELD.md`)

- **New idea (sketch):** condition on a fixed product of link balls in fibre gauge, not on plaquette-level calm sets.
  The bad side has no thin lens (proved); the good side is convex, so Bakry–Émery, Herbst and Helffer–Sjöstrand apply
  given a deterministic Hessian bound (input (a)); at ω = ∅ that is Theorem `chartdecay`. Rarity is then needed only on the torus, annealed, at a fixed
  threshold. Island threshold 4.5·10³⁷ at the chart radius of Theorem `chartdecay` (was 8.3·10¹⁶³).
- **Proved:** a box-local Peierls lemma for excess action under any exterior. It is marginal for (R1): forced curvature
  has attenuation 1 (heuristic).
- **Unchanged:** D8. Chessboard rarity covers O(log β) scales only (≤ 19 at β = 10⁶⁴); RG-improved rarity is Bałaban's.
- **Next:** the relative chart for ω ≠ ∅ (`LARGE_FIELD.md` §3a), then input (b).

### Relative-chart audit and local witness (November 2026, `RELATIVE_CHART.md`)

- **Proved (deductions from the paper):** the canonical reconstruction R(c) is in the full fibre of c. A link
  farther than r from R(c) forces a plaquette defect > r/50 in a fixed block collar, for
  r ≤ 50·10^(-8)ρ_G. This supplies a specified version of input (b), without an exterior-uniform rarity claim.
- **Proved (limitations):** the chart product uses tangent coordinates, not physical links. Sandwiching domains
  does not transfer covariance decay; h is not an inhomogeneous fibre projection. Full-fibre disintegration also
  retains incoming averages at the window boundary.
- **Sketch / open:** construct a chart solving all touching constraints and check its domain and Hessian;
  prove constrained bad-side coefficients and the environment bootstrap hypothesis. The single-scale reduction
  to two deterministic checks is not yet justified. The earlier β ≈ 4.5·10³⁷ is conditional ledger arithmetic,
  with an unchecked witness threshold and finite-intersection losses. D8 is unchanged.

### Touching-constraint graph and witness losses (November 2026, `TOUCHING_CHART.md`)

- **Proved (deduction):** around an assumed calm completion, all touching averages
  admit a smooth local graph after freeing incoming private links. The graph
  has finite range, path-link Lipschitz constant <44, and an explicit
  product of inverse path-link determinants as its Haar disintegration density.
- **Proved (deduction):** input (b) gives finite-intersection torus rarity
  [6*2401*z(r/50)]^(|S|/2401) when the bracket is below one.
- **Sketch / open:** background selection, the actual label-domain geometry,
  the graph-density Hessian, constrained bad-side response and the environment
  bootstrap. The graph is a local parametrisation, not relative chart decay.
  D7, (I1), (D) and D8 remain open; no working island threshold is certified.

### Rough coarse backgrounds: an exact obstruction (`ROUGH_BACKGROUND.md`)

- **Proved:** on a fine torus of side four there is a canonical reconstruction
  R(c) with every relative physical link label good and an exact four-link
  fibre curve with Q''(0)=-4. Its constrained action Hessian is negative.
- **Computation:** an eigensolver found the obstruction; the final four-link
  certificate uses integer arithmetic only (`d5/rough_certificate.py`).
- **Sketch / next:** add a coarse-curvature hypothesis or label before seeking
  relative chart decay. Background selection, the label domain and the
  environment term remain separate checks. D8 is unchanged.

### A specified canonical graph law (`CANONICAL_CHART.md`, `COARSE_LABELS.md`)

- **Proved (deduction using the paper's `smallfibre` and covariance argument):**
  private-link fixed points around R(c), with the frozen exterior as input,
  select a measurable background and solve all touching constraints.
  The graph density has a local product of inverse determinants.
- **Proved (deduction):** on the stated product graph domain, with
  kappa_R+176t<=10^(-10) and lambda>=2*10^11, the effective Hessian is at
  least lambda/104. With kappa_R+352t<=10^(-10), the same law has weighted
  covariance decay at a conservative coarse-distance rate 10^(-14).
  This is a theorem about a specified restricted chart law, not (I1).
- **Proved:** coarse curvature has a fine witness at eta/105. A coverage
  count improves the relative-link finite-intersection exponent denominator
  from 2401 to 256; including coarse-curvature labels gives denominator 1296.
- **Computation:** rational bounds with explicit Taylor remainders, finite
  overlap checks, and a nonabelian determinant diagnostic accompany the notes.
- **Sketch / open:** identify or replace the physical label-conditioned
  domain, control bad-side constrained response, and prove the environment
  comparison. D7 and (I1) are open; no island threshold is certified.
  D8 and the continuum existence question are unchanged.

### Match the law by labelling independent links (`INDEPENDENT_LABELS.md`)

- **Proved:** the original physical all-link label domain is nonconvex in
  non-path coordinates even when c=I (`LABEL_DOMAIN.md`, analytic BCH proof).
- **Proved:** if a pair's non-path links are within r of R(c), every path-link
  solution is within 14r of c_k, including all cut-off branches. Hence the
  path solution is unique in this regime; no component is discarded.
- **Proved (deduction):** label only non-path physical links, plus coarse
  curvature. On W with W^(+5) good, use A_W plus incoming path links as
  the free set and freeze its complement. The full conditional law
  equals the product graph law. Under eta+352r<=10^(-10) and
  lambda>=2*10^11, the previous weighted covariance estimate therefore
  holds for this actual enlarged good-window law.
- **Proved:** bad non-path data can have three path-link solutions with the
  same average (`BAD_BRANCHES.md`). Bad-side estimates must include branches.
- **Computation:** quaternion endpoint/midpoint probes support the domain
  obstruction; rational inequalities and numerical roots check the bad-branch
  example. They are labelled separately from the analytic arguments.
- **Sketch / next:** constrained bad-side transport (including branches),
  island/collar comparison, and the environment covariance. The five-block
  collar has not been entered into a working island ledger. (I1), D7, D5
  and (D) remain open; D8 is unchanged.

### Coupling minimality interface, checked in Lean (`COUPLING_INTERFACE.md`)

- **Proved (kernel checked):** `CouplingMinimality.lean` derives the finite
  discrepancy recurrence and weighted boundary comparison from minimality
  of coupling cost. Individual-update or random-scan stationarity is not
  needed for this route. Approximate minimality costs epsilon_i/w_i.
- **Computation / validation:** checked against the pinned Lean/Mathlib
  v4.33.0 successfully. The supplied larger bridge text
  was assessed mathematically but was not compiled as a repository file.
- **Sketch / open:** instantiate the abstract coupling/update hypotheses;
  this formal interface does not establish (I1) or any new gauge-theory gap.
  Next research target remains bad-side transport including path branches,
  then island/environment comparison. D8 is unchanged.

### Analytical interfaces and finite-box corrections (November 2026)

- **A, proved deduction + exact computation:** the zero-exterior attainable
  test response is [1+I2(beta_W)/I1(beta_W)]/4. The older (1+beta)/4 is its
  first-order approximation, not a certified optimum. No endpoint changes.
- **B, proved finite reduction + computation:** `B_WICK.md` includes the full
  finite-box interacting-measure, observable and Haar corrections. One B200
  pilot was independently reproduced on CPUs; an exact rational small-box
  check gives a nonzero normalised order-g^2 correction and exact two-tree
  invariance. Individual sectors depend on gauge coordinates. This is not
  the infinite-volume flowed-vacuum correction of B2.
- **C, proved scope limitation:** overlap-H witness information alone cannot
  avoid the 1/H exponent loss (`FRONTIER_LEMMAS.md`). Better constants require
  extra gauge geometry, not just another union bound.
- **D/I, proved deductions:** `CONNECTED_NORMALIZERS.md` gives exponential
  locality of connected log-normalisers for continuous source changes
  inside the good product chart. `ENVIRONMENT_OBSTRUCTION.md` shows that
  even all finite-intersection rarity bounds and zero quenched covariances
  do not alone control the environment covariance. Hard label flips and
  noncalm islands remain outside the new estimate.
- **D8, proved finite interface:** quadratic-form perturbations survive
  constrained minimisation, but do not imply locality; a weighted operator
  estimate and scale-summable errors are sufficient additional inputs in
  the stated finite setting (`FRONTIER_LEMMAS.md`). No RG construction.
- **E, proved conditional lemma:** positive spectral weight and a uniform
  physical-energy first moment imply a nonzero limiting correlation if
  correlations converge. A physical gap alone does not imply nontriviality.
  Neither tightness nor continuum existence is established.

The remaining priorities are bad-island constrained transport, connected
label/island normalisers and the environment bootstrap, then multiscale
composition. (I1), D7, the required D5, (D), D8 and (E) remain open.

### Review of James's supplied ideas (2026-11-03)

- **Assessment:** `JAMES_IDEAS_REVIEW.md` records both supplied texts. The
  first is useful as a continuum/nontriviality checklist and conditional
  coupling infrastructure, rather than a new dynamical estimate.
- **Proved (existing kernel-checked interface):** the second text's useful
  minimality refinement is already in `CouplingMinimality.lean`. Its
  proposed measure bridge could instantiate that interface; the pasted
  code is truncated and explicitly uncompiled.
- **Sketch / next:** actual measurable coupled updates and cost-minimiser
  existence for G1; bad-branch transport and connected hard-label
  normalisers for D/I; physical-energy tightness and nonvanishing
  renormalised correlations for E. These tracks can proceed independently.
  The E ideas mentioned in the second text's final sentence are absent
  from the supplied file. No new item is closed; D8 is unchanged.

**Third-text follow-up (assessment, 2026-11-03):** the complete trimmed
bridge is now supplied and reviewed in `JAMES_IDEAS_REVIEW.md` section 3.
It assumes local non-improvement (`hmin`); deriving it from a global
cost-minimising measure still needs the existing finite interface to be
instantiated with actual updates and finite coordinate integrals. The
source remains uncompiled in this review. No new E argument is included;
the workflow prompt remains truncated. G1 and the dynamical walls remain
open.

### Frontier and next experiment/proof priorities (2026-11-03)

**Assessment / future directions:** `FRONTIER_PLAN.md` consolidates the
A--E/D8/I status, James's conditional bridge and the next efforts. CPU
priorities are branch-complete bad-side transport diagnostics, connected
two-island partition weights, and the complete spatial-flow correction.
Larger B tensor contractions are a conditional B200 candidate only after
CPU validation and a cost comparison. G1's measure-level construction is
the next formalisation target; bad-island/environment bounds, all-scale
locality and continuum nontriviality remain analytical obligations. The
note records expectations and resource limits; no experiments were run
and no new theorem is claimed.

### Frontier implementation pass (2026-11-03)

- **A/G1, proved / kernel checked:** James's trimmed measure bridge is now
  `WassersteinBridge.lean`, with added actual-measure cost-minimality
  connection theorems. `CouplingMinimizer.lean` proves nonempty compact
  fixed-marginal probability couplings and existence of a minimum for
  every continuous real cost on a compact metric space. All four modules
  in the dependency chain compile against v4.33.0; the complete default
  target was not rebuilt. Measurable conditional selection, discrepancy
  instantiation, specification-dependent update bounds and limit arguments
  remain open. No new strong-coupling endpoint is claimed.
- **B, proved finite reduction + computation:** `NONLINEAR_FLOW_CHECK.md`
  supplies the nonlinear full-link flow-jet formulas. The exact compact
  one-variable calibration has a normalised correction three times the
  linear-flow-only answer. CPU checks also validate local spatial-action
  vertices on a 3x2x2x2 box and show that freezing tree links changes the
  flow. The complete flowed Wick covariance and B2 remain open.
- **C/D/I, proved finite-measure lemma + open application:**
  `ISLAND_KERNEL.md` proves that a common positive collar kernel's uniform
  cross-ratio bound survives arbitrary separate hard island weights and
  branch mixtures, with a binary-label influence consequence. The actual
  constrained gauge factorisation and uniform calm-interface estimates
  are not yet supplied; the environment recursion remains open.
- **D/I, proved scalar local result + computation:** `TRANSVERSE_FOLD.md`
  solves a transverse variable through the explicit scalar path fold.
  Independent quadrature representations agree, but hard-boundary fluxes
  and uniform SU(2) bad-side transport remain uncontrolled.
- **D8, proved finite interface:** `WEIGHTED_SCHUR.md` bounds locality and
  perturbation amplification under one Gaussian elimination. No all-scale
  Yang--Mills estimates or large-field RG construction are provided.
- **E, proved conditional spectral lemma:**
  `POSITIVE_TIME_TIGHTNESS.md` converts a fixed-physical-time half-time
  correlation-ratio lower bound into tight normalised spectral weight and
  nonzero normalised correlations. The Yang--Mills ratio bound, local-field
  construction, convergence and reconstruction remain open.

**Computation / resource assessment:** bounded reviewed CPU jobs supplied
kernel checks and finite diagnostics. No GPU was allocated: the numerical
pilots take seconds on CPUs. Submission tooling, exact source bundles,
operational logs and raw results remain outside Git. The next targets are
measurable conditional couplings, actual island/collar factorisation,
hard-boundary transport and complete flowed Wick coefficients. (I1), D7,
the required D5, (D), D8 and (E) remain open.

### Upstream contribution comparison (2026-11-03)

**Assessment / prepared contribution:** `UPSTREAM_STATUS.md` compares
research head `6e0608b` with `creelie/4Y-M:main` at `fdd2dc5`.
Upstream already has the endpoint improvements and newer paper/release
work; a full merge has seven conflicting paths. A focused six-file Lean
contribution is pushed as `romanov360:upstream/g1-coupling` (`e7a52e4`),
based directly on current upstream and cleanly mergeable. Its proofs and
dependency pins match the checked sources. No PR was opened because
GitHub CLI lacks API authentication; the note links the prepared compare
form. Analytical notes need a separate review against upstream's newer
forest-tilt and coarse-defect results. Private tooling/data are excluded.


### Complete upstream contribution (2026-11-03)

**Integration status / validation:** the user requested the full public
research continuation. `upstream/all-research` merges this research history
with upstream `main` at `fdd2dc5`, superseding the six-file candidate as
the intended contribution. `RESEARCH_INDEX.md` maps all fronts, scientific
programs, validation scope and remaining external dependencies. Both paper
sources and release artifacts retain upstream's current versions; all
unique research notes, scripts, small certificates and formal modules are
included. Alternate certificate run logs are retained separately.

Eight standalone mathematical programs behind the newer Wick, nonlinear
flow, branch, scalar collar and spectral diagnostics are now supplied in
`wick/` and `d5/`; raw data and submission tooling remain excluded. Earlier
A2/A3 exact-response and ceiling language is explicitly superseded by the
Bessel correction. The formal and analytical scopes remain those of the
preceding implementation pass: no new (I1), D8 or continuum theorem is
claimed by this integration. The complete contribution is open as
[draft PR #160](https://github.com/creelie/4Y-M/pull/160), with 129 changed
files. `UPSTREAM_STATUS.md` records the reconciliation and earlier scope.

### After upstream merge (2026-11-04)

**Proved finite criterion / open application:** `ENVIRONMENT_COMPARISON.md`
connects uniform conditional label odds to environment covariance through
a Dobrushin resolvent. The necessary all-background influence bounds and
locality of conditional physical observables are not proved for the gauge
fibres. `RESEARCH_DIRECTIONS_2026_11_04.md` records literature and next
steps for bad regions, multiscale crossover and continuum reconstruction.
No new Yang--Mills mass gap or continuum construction is claimed.

### Interfaces and boundary obstruction (2026-11-04)

- **Proved abstract interfaces:** `CROSSOVER_LEDGER.md` gives an exact
  scale-dependent error ledger and a small-gain split criterion;
  `CONTINUUM_GAP_INTERFACE.md` gives sufficient semigroup convergence and
  density hypotheses for a physical gap to survive reconstruction. Neither
  supplies the needed Yang--Mills inputs.
- **Proved scalar obstruction:** `TRANSVERSE_FOLD.md` now gives a submersion
  with a hard boundary where the normalizer and bounded-observable response
  are not uniformly Lipschitz. A full-rank constraint alone cannot close bad
  transport; boundary transversality or another mechanism must be proved.
- **Open:** uniform gauge-island/environment estimates, the nonlinear
  crossover and a nontrivial reconstructed continuum theory.

### Boundary-preserving transport criterion (2026-11-04)

**Proved conditional criterion:** `BOUNDARY_TRANSPORT.md` derives a fibre
response identity and a Wasserstein bound from a boundary-preserving lift
with uniformly controlled speed and score oscillation. Its active-stratum
rank test shows why the scalar hard-boundary counterexample cannot satisfy
those hypotheses. **Open:** whether the actual gauge labels admit compatible
uniform lifts and score bounds; the bulk private-link submersion alone does
not establish them. No gauge-island or continuum result is claimed.

### Three-front decomposition (2026-11-04)

- **Proved finite obstruction:** `RARE_LABEL_INFLUENCE.md` shows that even
  arbitrarily rare labels can have conditional influence tending to one.
  Graded control or conditional collar estimates remain to be proved.
- **Proved scalar recursion / open gauge input:** `CROSSOVER_COUNT.md`
  quantifies the number of steps and accumulated one-loop-size errors and
  lists five separate D8 obligations, including matching physical units.
  Its weighted refinement proves how marginal growth changes the error
  summability threshold, and gives a scalar example where that threshold
  is sharp. Actual gauge RG constants and cancellations remain unproved.
- **Proved symmetry obstruction:** `CONTINUUM_SYMMETRY_TEST.md` shows that
  tightness and exact hypercubic symmetry do not imply continuum rotation
  invariance. The required rotation-defect and reconstruction bounds remain
  open. `OS_LIMIT_LEDGER.md` proves a countable diagonal extraction with
  positivity and distinguishes an ordinary variance witness from a
  nonzero reflection-norm witness. Its uniform bounds and the analytic
  reconstruction hypotheses remain open. None of these results closes D,
  D8 or the Clay problem.

### Flat-holonomy correction to input (b) (2026-11-04)

**Proved gauge-field counterexample:** `FLAT_HOLONOMY_OBSTRUCTION.md`
exhibits a flat SU(2) configuration in block tree gauge with links far from
identity. Thus the unrestricted absolute version of bad-link-implies-large-
plaquette is false, not merely unchecked. A same-endpoints path comparison
is proved there. A relative witness and its annealed rarity were already proved in
`RELATIVE_CHART.md` and `COARSE_LABELS.md`; it was incorrect to list that
witness itself as open. Matching its labels to the conditional law and
controlling bad-side response and the environment remain open.

**Proved conditional counting:** `WITNESS_COUNTING.md` turns any local
finite-multiplicity curvature witness plus joint plaquette rarity into a
joint bad-label rarity bound, with explicit entropy and overlap factors.
The sharper existing gauge result is in `COARSE_LABELS.md`; both are
annealed rather than uniform in the exterior.

### SU(2) shell hazard (2026-11-04)

**Proved:** `SU2_RADIAL_HAZARD.md` gives a uniform O(1+kappa) radial
hazard bound on SU(2), sufficient in order for the fixed-threshold shell
ledger in `CALM_LINK.md` if its conditioning and flux assumptions are
verified. The sharper Gaussian-scale estimate, the actual maximum-shell
flux, and large-field displacement remain open.
