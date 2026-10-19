# Yang–Mills: lattice obstructions, audit and Lean certificate

Deep Bhattacharjee and Ushashi Bhattacharya.

Archived on Zenodo: release v1.0.0, [doi:10.5281/zenodo.23444338](https://doi.org/10.5281/zenodo.23444338).

This repository contains two papers and a Lean 4 / Mathlib project
(`YangMills/`) checking their finite statements:

- *Spacing-Uniform Lattice Inequalities Cannot Produce the Yang–Mills Mass Gap*
  (`tex/main.tex`, compiled to `main.pdf`; RevTeX 4-2, PRD style);
- *Mass Gap of Lattice Yang–Mills Theory up to the Dobrushin Ceiling*
  (`strong-coupling/prd/`, PDF also at
  `preprintYangMillsDobrushinCeiling/YangMillsDobrushinCeiling.pdf`; see
  `strong-coupling/README.md`).

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

**Forest gauge fixing (proved by hand, not yet refereed;
`strong-coupling/gauge-fixed/`):** four-dimensional SU(2) has a unique
free-boundary state with a transfer-matrix gap m ≥ −log(3β_W) for β_W < 1/3,
and every SU(N) has a gauge-fixed Hessian bound 29/2 giving Bakry–Émery,
decay and a gap for |β| < 1/29 (SU(3): β_W < 9/29). For SU(3) a single-link
Dobrushin bound in the same gauge, with a spectral-gap estimate for the
one-link measure and exact tilted Haar moments, gives uniqueness, decay and a
gap for β_W ≤ 0.49.

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
