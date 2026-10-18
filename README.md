# Yang–Mills: lattice obstructions, audit and Lean certificate

Deep Bhattacharjee and Ushashi Bhattacharya.

This repository contains the paper
*Spacing-Uniform Lattice Inequalities Cannot Produce the Yang–Mills Mass Gap*
(`tex/main.tex`, compiled to `main.pdf`; RevTeX 4-2, PRD style) and a Lean 4 /
Mathlib project (`YangMills/`) checking its finite statements.

## Status

**Proved (and checked in Lean, no `sorry`, standard axioms only):**

- No-go: if the lattice gap is bounded below in lattice units uniformly in the
  spacing, every renormalised two-point function vanishes at positive time in
  the limit (`SpectralNoGo.lean`).
- Wilson action bounds `0 ≤ S ≤ 2N|P|`, oscillation `≤ 2Nβ|P|`, and the
  Holley–Stroock factor `e^{-2Nβ|P|} → 0` with the volume (`WilsonAction.lean`).
- Holley–Stroock perturbation for Poincaré and entropy inequalities on finite
  spaces (`HolleyStroock.lean`).
- Bianchi identity on a cube (`WilsonAction.lean`).
- One-loop running coupling attains every value at a spacing-independent scale
  (`RunningCoupling.lean`).
- A correlated local action whose variance per volume is unbounded, refuting
  the independence step of the cumulant argument (`CumulantCounterexample.lean`).

**Open:** the Clay problem itself. The paper isolates the residual statement
(Hypothesis R: `Δ_lat(a) = m a + o(a)` with Osterwalder–Schrader convergence)
and records that no argument in the eight audited manuscripts, and no paper in
the literature, proves it.

## Build

```
lake exe cache get
lake build
lake env lean scripts/Axioms.lean   # prints the axioms of the main theorems
cd tex && latexmk -pdf main.tex
```
