# PR #160: scientific programs and progress toward the open inputs

This summarises the full public contribution in
[creelie/4Y-M#160](https://github.com/creelie/4Y-M/pull/160).
The PR adds rigorous finite results, formal verification and better-defined
partial inputs. The main continuum problems remain open. Detailed status
and hypotheses are in [ROADMAP.md](strong-coupling/forest/ROADMAP.md),
[D_STATUS.md](strong-coupling/forest/D_STATUS.md) and [RESEARCH_INDEX.md](strong-coupling/forest/RESEARCH_INDEX.md).

## Programs 1–8: the eight standalone scientific programs

All eight support finite calculations or conditional arguments. Generated
outputs and private operational tooling are excluded from the contribution.

1. **[wick.py](strong-coupling/forest/wick/wick.py) — computation.** Computes the complete
   finite-box perturbative correction to Wilson-observable correlations,
   including action, observable and Haar-measure contributions.
2. **[jet.py](strong-coupling/forest/wick/jet.py) — independent computation.** Reproduces those
   coefficients using explicit covariance derivative formulas, with lower
   memory requirements than the automatic-differentiation implementation.
3. **[exact.py](strong-coupling/forest/wick/exact.py) — exact rational computation.** Computes
   the smallest-box coefficients and verifies their equality in two
   maximal gauge trees. This does not kernel-verify the analytic Wick
   reduction or establish an infinite-volume result.
4. **[flow_calibration.py](strong-coupling/forest/wick/flow_calibration.py) — finite model
   calibration.** Checks an exactly solvable single-variable SU(2) flow
   against Haar quadrature. Omitting nonlinear flow changes its normalised
   order-g-squared correction by a factor of three.
5. **[lattice_flow_jets.py](strong-coupling/forest/wick/lattice_flow_jets.py) — computation.**
   Tests nonlinear spatial-flow formulas on a small lattice, including the
   exponential-coordinate metric correction and movement of initially
   zero tree links. It checks the local vector field, not the full flowed
   interacting covariance.
6. **[branch_controls.py](strong-coupling/forest/d5/branch_controls.py) — computation.** Examines
   SU(2) constraint Jacobians near path folds and hard-label boundaries,
   checking derivatives against independent finite differences. A singular
   path-only Jacobian does not establish divergent integrated transport.
7. **[frontier_diagnostics.py](strong-coupling/forest/d5/frontier_diagnostics.py) — computation.**
   Tests transverse smoothing of a scalar fold, connected collar bounds
   after separate hard island restrictions, and finite spectral-tightness
   examples. These scalar diagnostics do not prove SU(2) transport bounds.
8. **[frontier_calibration.py](strong-coupling/forest/d5/frontier_calibration.py) — exact and
   finite calibrations.** Checks supporting examples concerning Bessel
   responses, energy shifts, rarity overlap, locality and nontriviality.

Usage and dependencies are documented in [wick/README.md](strong-coupling/forest/wick/README.md)
and [d5/README.md](strong-coupling/forest/d5/README.md). All eight ran successfully from their
public locations. Both floating Wick methods agreed with the exact
rational small-box coefficients; the separate Lean core block certificate
also compiled. These checks do not certify the required uniform
Yang--Mills estimates.

## Progress by research item

| Item | Progress and remaining limit |
|---|---|
| **A1 — SU(3) blocks** | **Already upstream:** the certified improvement to beta_W <= 0.77. This PR preserves it and adds supporting research material; it claims no new endpoint improvement. |
| **A2 — SU(2) response** | **Proved deduction / exact computation:** corrected the purported exact response using the zero-field Bessel test. Its supremum over all fields and Lipschitz tests remains open; no improved coupling threshold is proved. |
| **A3 — larger blocks** | **Computation / heuristic:** added shape calculations and response diagnostics. Proposed finite-coupling ceilings are not rigorous impossibility results. |
| **A4 — formalisation** | **Kernel checked:** finite comparison arithmetic, cost-minimality reasoning, measure-level update lemmas and compact-space continuous-cost minimiser existence. Measurable SU(2) conditional updates and the specification-dependent estimates remain missing. |
| **B1 — three dimensions** | **Computation / power-counting assessment:** evidence supporting finite next-order corrections without the four-dimensional ultraviolet logarithm. No nonperturbative theorem. |
| **B2 — slice-ratio route** | **Proved finite/perturbative reductions plus computation:** free-vacuum lattice analysis, exact finite-box interacting corrections and validated nonlinear-flow formulas. The full flowed interacting estimate and its uniform limits remain open. |
| **B3 — stochastic quantisation** | **Literature assessment:** clarified why existing state-space results do not supply the needed slice measure or correlation estimate. |
| **B4 — further unflowed 3D Monte Carlo** | **Not developed further:** retained existing diagnostics; no new uniform slice-ratio estimate follows. |
| **C1 — sharper cube constants** | **Not pursued as a separate advance:** local constants alone do not resolve the composition problem. No new certified full cube estimate. |
| **C2 — rarity** | **Proved deductions / computation:** torus rarity analysis, curvature witnesses and sharper overlap counting. Rarity still does not establish harmlessness. |
| **C3 — continuum tree bounds** | **Open / no substantive closure:** the contribution does not construct the required continuum tree bounds. |
| **D1 — Gaussian estimate at every scale** | **Already proved in the paper:** preserved the existing Gaussian result. This PR adds no nonlinear all-scale extension. |
| **D2 — one composition step** | **Already proved in the paper:** preserved the existing constrained-energy and locality results. Finite weighted refinements are recorded under D8. |
| **D3 — calm-fibre decay** | **Proved conditional deduction:** constructed a canonical relative chart and matched an actual good-window conditional law with covariance decay under explicit calmness and collar assumptions. This improves the usable good-region input. |
| **D4 — off-chart mass without a calm completion** | **Proved obstructions / open estimate:** identified rough-background and multiple-branch failures of simpler reductions. Local scalar regularisation is available, but the required uniform bad-side bound remains open. |
| **D5 — fibre decay to uniform cube contraction** | **Proved conditional reductions / open target:** supplied calm-case estimates and corrected chart/domain matching. Transport across exteriors, collar costs and the required uniform cube coefficient remain open. |
| **D6 — large fields are rare** | **Proved torus setting / deductions:** preserved the paper's rarity result and added local witnesses, overlap bounds and finite-intersection analysis. No uniform rarity under every exterior is inferred. |
| **D7 — large fields are harmless** | **Proved conditional lemmas / open application:** a common collar-kernel bound survives separate hard island weights and branch mixtures. Actual gauge-theory factorisation, hard-boundary transport and environment control remain open. |
| **D8 — multiscale composition** | **Proved finite interfaces:** weighted locality and Gaussian elimination estimates. No all-scale renormalisation construction. |
| **E — continuum existence** | **Proved conditional spectral criteria:** tightness and nontriviality inputs from fixed-physical-time correlation ratios. No continuum field construction or reconstruction theorem. |
| **I1 — overarching cube-condition input** | **Partial ingredients improved; open:** the chart, rarity, transport and comparison work above does not yet prove continuum-uniform cube contraction. |
| **I2 — overarching slice-ratio input** | **Partial ingredients improved; open:** free-vacuum reductions, exact finite-box corrections and positive-time spectral criteria do not yet prove the full interacting cutoff-uniform ratio. |
| **I3 — overarching tree-bound input** | **No substantive closure:** continuum tree bounds remain open. |

Here **I1/I2/I3** mean the cube condition, slice ratio and tree bounds, as
in [IDEAS.md](strong-coupling/forest/IDEAS.md). The local mixture bootstrap in
[ANNEALED.md](strong-coupling/forest/ANNEALED.md) also uses I1/I2 labels for its own intermediate
hypotheses; those remain open as scoped in the latest status notes.

## Relation to the Clay problem

**Assessment.** The PR strengthens the finite foundations and identifies
more precisely the estimates still needed. It proves neither continuum
Yang--Mills existence nor a continuum mass gap, and closes no central
continuum wall. The main remaining work is uniform bad-region and
environment control, multiscale composition through the crossover, and a
nontrivial continuum construction satisfying the reconstruction axioms.
