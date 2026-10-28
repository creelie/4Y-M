# Yang–Mills: lattice obstructions, the lattice mass gap for every compact simple group, and Lean certificate

Deep Bhattacharjee, Ushashi Bhattacharya and James Montgomery Lockwood.

Archived on Zenodo; the concept DOI [10.5281/zenodo.23444338](https://doi.org/10.5281/zenodo.23444338)
resolves to the latest release.

| Release | Version DOI |
|---|---|
| v1.0.0 | [10.5281/zenodo.23444339](https://doi.org/10.5281/zenodo.23444339) |
| v2.0.0 | [10.5281/zenodo.23449141](https://doi.org/10.5281/zenodo.23449141) |
| v2.0.1 | [10.5281/zenodo.23449152](https://doi.org/10.5281/zenodo.23449152) |
| v2.0.2 | [10.5281/zenodo.23449466](https://doi.org/10.5281/zenodo.23449466) |
| v3.0.0 | [10.5281/zenodo.23484043](https://doi.org/10.5281/zenodo.23484043) |
| v3.0.1 | [10.5281/zenodo.23484057](https://doi.org/10.5281/zenodo.23484057) |
| v3.1.0 | [10.5281/zenodo.23486278](https://doi.org/10.5281/zenodo.23486278) |
| v3.1.1 | [10.5281/zenodo.23486348](https://doi.org/10.5281/zenodo.23486348) |
| v3.2.0 | [10.5281/zenodo.23487516](https://doi.org/10.5281/zenodo.23487516) |
| v3.2.1 | [10.5281/zenodo.23487558](https://doi.org/10.5281/zenodo.23487558) |
| v3.3.0 | [10.5281/zenodo.23490004](https://doi.org/10.5281/zenodo.23490004) |
| v3.3.1 | [10.5281/zenodo.23490028](https://doi.org/10.5281/zenodo.23490028) |
| v3.4.0 | [10.5281/zenodo.23541099](https://doi.org/10.5281/zenodo.23541099) |
| v3.5.0 | [10.5281/zenodo.23711288](https://doi.org/10.5281/zenodo.23711288) |
| v3.6.0 | [10.5281/zenodo.23717210](https://doi.org/10.5281/zenodo.23717210) |
| v3.7.0 | [10.5281/zenodo.23717915](https://doi.org/10.5281/zenodo.23717915) |
| v3.8.0 | [10.5281/zenodo.23718452](https://doi.org/10.5281/zenodo.23718452) |
| v3.9.0 | [10.5281/zenodo.23719806](https://doi.org/10.5281/zenodo.23719806) |
| v3.10.0 | [10.5281/zenodo.23720481](https://doi.org/10.5281/zenodo.23720481) |
| v3.11.0 | [10.5281/zenodo.23721334](https://doi.org/10.5281/zenodo.23721334) |
| v3.12.0 | [10.5281/zenodo.23721638](https://doi.org/10.5281/zenodo.23721638) |
| v3.13.0 | [10.5281/zenodo.23721860](https://doi.org/10.5281/zenodo.23721860) |
| v3.14.0 | [10.5281/zenodo.23724715](https://doi.org/10.5281/zenodo.23724715) |

`YangMills_masterzipfinal2.zip` holds the repository at v3.14.0 and the
Hexagon package of the paper *Mass Gap of Four-Dimensional Lattice
Yang–Mills Theory for Every Compact Simple Gauge Group*
([download](https://github.com/creelie/4Y-M/raw/main/YangMills_masterzipfinal2.zip)).

This repository contains three papers and a Lean 4 / Mathlib project
(`YangMills/`) checking their finite statements:

- *Spacing-Uniform Lattice Inequalities Cannot Produce the Yang–Mills Mass Gap*
  (`tex/main.tex`, compiled to `main.pdf`; RevTeX 4-2, PRD style);
- *Mass Gap of Lattice Yang–Mills Theory up to the Dobrushin Ceiling*
  (`strong-coupling/prd/`, PDF also at
  `preprintYangMillsDobrushinCeiling/YangMillsDobrushinCeiling.pdf`; see
  `strong-coupling/README.md`);
- *Mass Gap of Four-Dimensional Lattice Yang–Mills Theory for Every Compact
  Simple Gauge Group* (`strong-coupling/forest/`, REVTeX, 38 pages; PDF also
  at `preprintYangMillsForestGauge/YangMillsForestGauge.pdf`; the one-column
  version in the AMS article class, 58 pages, is `strong-coupling/forest/cmp/`;
  the working notes and checks are in `strong-coupling/gauge-fixed/`). It
  was first posted as Hexagon preprint 2610.00341 under the title *Mass Gap
  of Four-Dimensional Lattice Yang–Mills Theory in a Forest Gauge*.

## Status

**Proved (and checked in Lean, no `sorry`, standard axioms only):**

- No-go: if the lattice gap is bounded below in lattice units uniformly in the
  spacing, every renormalised two-point function vanishes at positive time in
  the limit (`SpectralNoGo.lean`).
- Critical scaling: if the lattice gap is at least `m_a · a` with `m_a → m`, the
  limiting two-point functions obey the ratio bound of a theory with gap `m`
  (`LimitRatio.lean`); this is the one-sided bound the no-go leaves open.
- Wilson action bounds `0 ≤ S ≤ 2N|P|`, oscillation `≤ 2Nβ|P|`, and the
  Holley–Stroock factor `e^{-2Nβ|P|} → 0` with the volume (`WilsonAction.lean`).
- Holley–Stroock perturbation for Poincaré and entropy inequalities on finite
  spaces (`HolleyStroock.lean`).
- Bianchi identity on a cube (`WilsonAction.lean`).
- One-loop running coupling attains every value at a spacing-independent scale
  (`RunningCoupling.lean`).
- A correlated local action whose variance per volume is unbounded, refuting
  the independence step of the cumulant argument (`CumulantCounterexample.lean`).

**Lattice SU(2) up to the Dobrushin ceiling (checked in Lean, no `sorry`,
standard axioms only):**

- `StrongCoupling/QxCoefficients.lean`: the closed form of every Taylor
  coefficient of 𝒟(κ) = (8+2κ²)I₁² − κ²I₀² − κI₁(2κ) by Vandermonde's
  identity, the bound 16^m ≤ 4m·C(2m,m)², and the nonnegativity of all
  coefficients (the transverse half of the sharp quarter).
- `StrongCoupling/ZonalBound.lean`: the identities α = 16/(2s+3)² and
  4α² + T = (2s−5)p(s)/(4(2s+3)⁴), the sign of the cubic p on (5/2, 4.004],
  the derivative and residual of the Riccati barrier κ/(3/2+√(κ²+25/4)), the
  monotonicity of T, and the final zonal inequality.
- `StrongCoupling/Thresholds.lean`: the Dobrushin row sum 6(d−1)|β|, the
  threshold β_W < 2/9 at d = 4, and the gap m = −log(9β_W/2) > 0.
- `StrongCoupling/GaugeFix.lean`: the line field of the forest gauge
  (`strong-coupling/gauge-fixed/`): one fixed link in every plaquette,
  straight lines, twelve partners per active link (so the Dobrushin row sum
  is 12|β| = 3β_W), and the per-link budget 29/2 of the gauge-fixed Hessian
  bound (seven theorems, `decide`, axioms `propext` at most).

**Forest gauge fixing (proved by hand, not yet refereed; third paper,
`strong-coupling/forest/`, notes in `strong-coupling/gauge-fixed/`):** four-dimensional SU(2) has a unique
free-boundary state with a transfer-matrix gap m ≥ −log(3β_W) for β_W < 1/3,
and every SU(N) has a gauge-fixed Hessian bound 29/2 giving Bakry–Émery,
decay and a gap for |β| < 1/29 (SU(3): β_W < 9/29). For SU(3) a single-link
Dobrushin bound in the same gauge, with a spectral-gap estimate for the
one-link measure and a corrector of degree two in the matrix entries,
evaluated with the exact moments of the one-link measure at every staple sum
(character expansion, rigorous error bounds, a cover of the staple sums by
2341 boxes), gives uniqueness, decay and a gap m ≥ −log(1.33319 β_W) for
β_W ≤ 3/4. The single-link route stops at or below β_W ≈ 0.77 with any
bound; past that, blocks of several links are needed. For SU(2), blocks made
of the three free links of one plaquette pass the single-link ceiling 1/3:
the free links split into such triples with period two, the middle link of
each block is integrated out, and the block constant is bounded analytically
with constants of the von Mises–Fisher partition function certified in ball
arithmetic (`strong-coupling/gauge-fixed/blocks/`). This gives uniqueness,
decay and a gap m ≥ −½ log α(β_W) for β_W ≤ 0.345, with α(1/3) ≤ 0.9561 and
α(0.345) ≤ 0.9973.
The same blocks work for SU(3) (`strong-coupling/gauge-fixed/blocks/su3/`):
the law of a block is controlled by the modulus Λ(r) of the cover, the
Hessian of log Z on r conv SU(3) and its third derivative, bounded on 841
boxes in ball arithmetic at r = 3/2 (eigenvalues in [0.056015, 0.247316],
third cumulant ≤ 0.182126). This gives α₃(3/4) ≤ 0.9653 and a gap
m ≥ −½ log α₃(β_W), so m ≥ 0.0177 at β_W = 3/4 instead of 1.1·10⁻⁴ from
single links. The blocks need only Λ(2β_W), not 2rΛ(r) < 1: a second cover
at r = 38/25 (1063 boxes, 2rΛ(r) ≤ 1.029967) with the constants on 784 boxes
(eigenvalues in [0.054978, 0.248280], third cumulant ≤ 0.183546) gives
α₃(0.76) ≤ 0.99950, so the SU(3) range extends to β_W ≤ 0.76, past the end of
the single-link route, with m ≥ 2.5·10⁻⁴ at β_W = 0.76.

**Every compact simple group (third paper, proved by hand, not refereed):**
in the forest gauge the least Ricci curvature of G is κ_G = h^∨/(2ℓ_ρ)
(Ric = −B/4), and for N|β| < 2κ_G/29 the theory has a unique free-boundary
state, exponential clustering and a transfer-matrix gap
m_G(β) = arccosh(1 + (κ_G/(N|β|) − 29/2)/12). The ranges in β_W = N²β run
from N²/29 for SU(N) to 124/29 for E₈ (exact table checked by
`strong-coupling/gauge-fixed/checks/groups/ricci.py`).

**The continuum step at one physical scale (conditional):** the paper
replaces the single-link Dobrushin condition by a condition D(n, α) on one
cube of n⁴ sites, averaged over its translates by 2Z⁴, and proves that it
gives a lattice gap ≥ log(1/α)/(2n). If D(n_k, α) holds along the continuum
tuning with n_k a_k ≤ ℓ for a fixed physical length ℓ, then every continuum
limit has a mass gap ≥ log(1/α)/(2ℓ). Two statements remain unproved: (D)
the cube estimate uniformly in the spacing, and (E) the existence of a
nontrivial Osterwalder–Schrader limit.

**What the cube estimate requires (proved):** D(n, α) holds at strong
coupling for every even n above an explicit size (for SU(2),
α = β_W/(1−3β_W)·(5/n − 2/n²)); if D(n, α) holds, then in the cube with the
identity outside the gauge-fixed links stay on average at distance at least
(1−α)D/3 from the identity (for SU(2), by four colours and tetrahedral
quaternions, at mean distance at least 1−α; for every compact simple group
with the Hilbert–Schmidt metric of an irreducible representation, at least
the fraction 1−α of the Haar disorder and mean distance (1−α)D/2), so for
every fixed n the condition fails once β is large; and for SU(2) there is an exterior, legal at every coupling, under
which the conditional law of one line of links has a ten-step gap at most
e^20 24^10 b^15 e^(−b/1280), so (D) cannot come from bounds on single lines
that are uniform in the exterior. None of this decides (D) along the
continuum tuning, where n grows with β. For SU(2) the condition needs
n⁴ ≥ 2(1−α) β_W/(5 L(β_W)), L(b) = 1 + log 24 + (3/2) log 9b: every
gauge-fixed link is the holonomy of a disc of at most (n/2)(n/2+1) plaquettes,
and a lattice Stokes inequality with the energy bound gives the size. The same
holds for every SU(N) with L_N(b) = 1 + ((N²−1)/2) log(81π²b) + N(N−1) log(12/11),
from the Jacobian of the exponential map and the volume of small balls. A
stationary flux configuration, every plaquette within 4π/L² of the identity
and a loop of side L/2 with holonomy −1, shows that the n² lost there is of
the right order. The continuum step also holds for observables smoothed by the
spatial gradient flow of Sakai and Sasaki, which is reflection positive and
leaves the transfer matrix unchanged up to unitary equivalence. Renormalised
by one condition at a single physical separation, these flowed fields have
continuum limits along subsequences for every tuning (Theorem "renormalised
flowed observables"), and a nontrivial limit exists if and only if the ratio
C(σ)/C(0) of the two-point function of the smeared packet at a fixed physical
separation to its reflected norm does not tend to zero (a point-field ratio at
coincidence would not do: it tends to zero for the free massless field); that
ratio condition is not proved. If the ratio is bounded below, renormalising at
separation σ gives tight spectral measures with explicit tail and moment
bounds, and a one-slab bound ⟨Ψ, (1−T)Ψ⟩ ≤ Ca on the transfer matrix would give
the ratio bound; neither is proved for Yang–Mills. For product test functions
the ratio lies between ratios of correlations of one spatially smeared slice
field at two fixed physical separations, so what is open is a lower bound
liminf D(r′)/D(r) > 0 for one pair r < r′; it holds for the free field. Tree
bounds on the truncated correlations of the flowed fields give all n-point
Schwinger functions as reflection-positive distributions along a
subsequence (Theorem "higher correlations"); they are not proved along the
continuum tuning. In a
periodic spatial box the normalising reflected norm C(0) is strictly positive
for product test functions, by Lüscher's strictly positive transfer matrix. At every n
the mean distance of the gauge-fixed links from the identity is at least
(1/9)·max(6β_W, 1)^(−1/2). For a Gaussian field the averaged cube coefficient
is exactly 1 at every size when the field is massless, and below
1 − (5/6)⁴(mℓ)²/((mℓ)²+160) at physical side ℓ with mass m. Finite counts:
`strong-coupling/gauge-fixed/checks/cube/cube_requirements.py`; Monte Carlo
probes that the proofs do not use: `influence_mc.py`, `eps_mc.py` in the same
directory (see `strong-coupling/OPEN.md`).

**Open:** the Clay problem itself (the four open problems of the
strong-coupling paper are listed in `strong-coupling/OPEN.md`). The paper isolates the residual statement
(Hypothesis R: `liminf Δ_lat(a)/a ≥ m > 0` with Osterwalder–Schrader convergence)
and records that no argument in the eight audited manuscripts, and no paper in
the literature, proves it.

## Build

```
lake exe cache get
lake build
lake env lean scripts/Axioms.lean   # prints the axioms of the main theorems
cd tex && latexmk -pdf main.tex
```
