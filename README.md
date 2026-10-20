# Yang–Mills: lattice obstructions, audit and Lean certificate

Deep Bhattacharjee and Ushashi Bhattacharya.

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

`YangMills_masterzipfinal2.zip` holds the repository at v3.2.1, the Hexagon
package of the forest-gauge paper and the journal packages of the papers
([download](https://github.com/creelie/4Y-M/raw/main/YangMills_masterzipfinal2.zip)).

This repository contains three papers and a Lean 4 / Mathlib project
(`YangMills/`) checking their finite statements:

- *Spacing-Uniform Lattice Inequalities Cannot Produce the Yang–Mills Mass Gap*
  (`tex/main.tex`, compiled to `main.pdf`; RevTeX 4-2, PRD style);
- *Mass Gap of Lattice Yang–Mills Theory up to the Dobrushin Ceiling*
  (`strong-coupling/prd/`, PDF also at
  `preprintYangMillsDobrushinCeiling/YangMillsDobrushinCeiling.pdf`; see
  `strong-coupling/README.md`);
- *Mass Gap of Four-Dimensional Lattice Yang–Mills Theory in a Forest Gauge*
  (`strong-coupling/forest/`, PDF also at
  `preprintYangMillsForestGauge/YangMillsForestGauge.pdf`; the journal
  version for Communications in Mathematical Physics, one column in the AMS
  article class, is `strong-coupling/forest/cmp/`; the working notes and
  checks are in `strong-coupling/gauge-fixed/`).

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
2341 boxes), gives uniqueness, decay and a gap m ≥ −log(1.3332 β_W) for
β_W ≤ 3/4. The single-link route stops at or below β_W ≈ 0.77 with any
bound; past that, blocks of several links are needed. For SU(2), blocks made
of the three free links of one plaquette pass the single-link ceiling 1/3:
the free links split into such triples with period two, the middle link of
each block is integrated out, and the block constant is bounded analytically
with constants of the von Mises–Fisher partition function certified in ball
arithmetic (`strong-coupling/gauge-fixed/blocks/`). This gives uniqueness,
decay and a gap m ≥ −½ log α(β_W) for β_W ≤ 0.345, with α(1/3) ≤ 0.9561 and
α(0.345) ≤ 0.9973.

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
