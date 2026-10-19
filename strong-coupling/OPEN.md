# The four open problems of the strong-coupling paper

The table `tab:gaps` of `prd/main.tex` lists four problems that separate the
results of the paper from the Clay mass gap. Normalisation as in the paper:
weight exp(Nβ Σ_p Re Tr Q_p), β_W = N²β. A problem counts as solved only when
every case below it is proved.

## 1. Intermediate fixed coupling

Uniqueness of the Gibbs measure, decay of correlations and a transfer-matrix
gap at a fixed lattice spacing, for every coupling up to the region where
Monte Carlo puts the crossover (β_W of order 2.2 for SU(2)).

| Range | Status | Where |
|---|---|---|
| SU(2), β_W < 2/9 | proved; the ceiling of every single-link Wasserstein–Dobrushin argument on the ungauged measure | `prd/`, Theorems E and F |
| every SU(N), \|β\| < 1/(8d) | proved by curvature; no curvature bound uniform over configurations goes further on the ungauged measure | `prd/`, Theorems B and C |
| SU(2), 2/9 ≤ β_W < 1/3 | proved by forest gauge fixing, not yet refereed; 1/3 is the ceiling of single-link Dobrushin after any forest gauge fixing | `gauge-fixed/forest.md` |
| every SU(N), \|β\| < 1/29 (SU(3): β_W < 9/29) | proved by the gauge-fixed Hessian bound 29/2, not yet refereed | `gauge-fixed/curvature.md` |
| SU(3), β_W ≤ 0.57 | proved by single-link Dobrushin in the forest gauge with an analytic spectral-gap bound for the one-link measure, exact tilted Haar moments (Schur–Weyl, Weyl integration, Bessel sums, ball arithmetic) and a blended corrector, not yet refereed; at β_W = 0.57 the remaining losses are the variance bound (0.050 in the proof against 0.018) and the gap bound (1.64 against a Galerkin value of 2.29, an upper bound for the true gap); Galerkin values at matrices of the allowed set rK put the limit of single-link Dobrushin estimated this way at or below β_W ≈ 0.77, beyond which blocks of several links are needed | `gauge-fixed/su3.md` |
| the best gauge-fixed Hessian constant h_F | 9 + √17 ≤ h_F ≤ 29/2; numerically 13.9740139454, attained by an abelian field; certificates on one plaquette and one vertex stop at 14.3839 | open |
| SU(2) from β_W = 1/3 to about 2.2 | no method; block (Dobrushin–Shlosman) criteria are the natural tool, and the paper's Theorem G gives one that is stable in every 't Hooft sector, but no box has been certified at these couplings | open |

## 2. The crossover to Bałaban's regime

Bałaban's renormalisation group controls the effective actions only at weak
coupling and in finite volume. Joining it to a strong-coupling statement
requires a criterion that one finite box can verify at a coupling where the
renormalisation group has already brought the effective action. Theorem G of
the paper is such a criterion at strong coupling; no criterion has been
verified in the window, and the hierarchical (Migdal–Kadanoff) models, where
confinement at every coupling is known, have no proved transfer to Z⁴.
Open.

## 3. The continuum limit with the Osterwalder–Schrader axioms

In two dimensions the continuum theory is known (Driver, Sengupta, Lévy). In
three and four dimensions no non-trivial limit with the axioms has been
constructed; lattice holonomies without smoothing converge to a trivial
(Haar white-noise) limit, so a smoothing that keeps reflection positivity
and a proof that the limit is not Gaussian are both needed. Open.

## 4. The gap in physical units

A positive gap at fixed spacing, as in problem 1, is a gap in lattice units.
The Clay problem asks for Δ_lat(a)/a ≥ m > 0 as a → 0 along the
renormalisation-group trajectory, together with the limit of problem 3
(Hypothesis R of `tex/main.tex`). The no-go theorem of `tex/main.tex` shows
that a gap bounded below in lattice units uniformly in the spacing forces
every renormalised two-point function to vanish, so the bound has to scale
with a. Open.
