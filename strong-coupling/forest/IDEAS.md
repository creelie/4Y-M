# Ideas tried for the three open inputs, and outcomes

This is a log of attempts, kept so that negative results are not repeated. The open inputs are:

- (I1) the cube condition;
- (I2) the slice ratio liminf D_k(r′)/D_k(r) > 0;
- (I3) the tree bounds.

## Done, positive

- **SU(2) plaquette blocks, β_W ≤ 0.345 → 0.351** (I1 at strong coupling, `gauge-fixed/blocks/sharp/`).
  - Two changes: interpolation of exponents in Lemma "one site" (ii), and the vMF Poincaré constant.
  - Certified in 128-bit interval arithmetic and rechecked in Lean.
  - It fails at 0.3511.
- **Composite free slice** (I2, `slice/`). The paper's magnetic slice field has the r > 0 ratio in the free theory, with or without the flow.
- **First-step energy vs. r > 0** (I2). The unflowed magnetic field has a divergent first-step energy but convergent r > 0 ratios. So the target is the r > 0 ratio, or the soft-part bound S′, not a first-step bound.
- **NLO leakage** (I2, `big-bang-monotile/yang-mills-continuum/LEAKAGE.md`).
  - The first-step energy has a one-loop logarithm, with coefficient 8/11 of 2b₀.
  - The r > 0 ratio sees only a UV-finite norm. No significant trend at t = 0.5√τ.

## Done, negative or limiting

- **Ceiling of the block method.**
  - The first-order coefficient K ≥ (1+β)/4 + O(β²) (Monte Carlo: Cov = 0.2500, 0.26255, 0.27492 at β = 0, 0.05, 0.1) gives α ≥ 10β(1+β).
  - That reaches 1 at β ≈ 0.366, so plaquette blocks with this Dobrushin matrix cannot go past about 0.366, whatever the constants.
- **3D Monte Carlo test of the unflowed ratio** (I2).
  - The unflowed slice correlator falls to 0.03–0.055 at one step and is noise (≲ 0.02) from two steps on, at every spacing tried.
  - The r > 0 unflowed ratio is not testable at this statistics.
  - Flowed ratios are measurable: D(.5)/D(.25) ≈ 0.15–0.30.
- **Electric slice fields** (I2). A field containing F₀ᵢ has a contact term, so its D(0) diverges. Only fields of one slice's spatial links (as in the paper) are suitable.

## Open, not yet tried

- **R2.** Sharpen the feedback constant C. The current C ≈ 0.038 is about twice the true feedback, and closing that gap would move 0.351 toward the 0.357–0.366 window. *Update: the measured response is (1+β)/4 in every field; the slack is the 1/(1−C) factor, not C (`../gauge-fixed/blocks/feedback/`). A second-order expansion around the product measure is the route.*
- **R6.** Larger blocks (two plaquettes, or a cube) to beat the 0.366 single-plaquette ceiling. *Update: ceilings measured, about 0.369 at best (`../gauge-fixed/blocks/shapes/`); not worth a proof.*
- **I11.** Rerun the 3D unflowed test with about 100× the statistics (multilevel or a reweighted estimator).
- **S′ at order g² on the lattice.** Show that the soft-part first-step bound holds with no logarithm in a lattice perturbative computation. The continuum computation says it does. *Update: done in the free vacuum, uniformly in a (`slice/sprime/LATTICE_SPRIME.md`).*
