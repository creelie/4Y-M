# PR #289, part 2: what input the framework needs

This continues `PROGRESS_SINCE_345304e.md`. Details are in `forest/A_GAUSSIAN_FIBRE_4D.md` §8.3, §8.4, §9.1, §10 and in `forest/H1_PATHCOV.md` (run 3).

## What this part accomplishes

- **The sufficient input is stated precisely (§8.4).** The framework needs exactly **(U)**: a uniform BP(ℓ/a, ϖ) for resampling cubes of fixed physical size ℓ. Thm `avgtransfer` turns (U) into a gap of at least 16ϖ/(625e²ℓ).
- **The tower route looks blocked in the UV (§8.3).** Under (LSC), the per-scale losses multiply over about k UV scales, which gives exponential growth. This is conditional: its step 1 (the optimal Gaussian C_B exceeds 1) is argued only by consistency.
- **A gap gives (U) in the free massive field (§10).**
  - Computed exactly to n = 24, ϖ_n converges to ϖ_∞ = (mℓ)²⟨1,(−Δ_D+(mℓ)²)⁻¹1⟩ on the unit 3-cube: 0.0196 at mℓ = 1 and 0.0719 at mℓ = 2.
  - The decay seen earlier at n ≤ 8 was an O(1/n) boundary correction.
  - A trial-field argument gives a lower bound that is uniform in n.
  - Higher chaos reduces to the linear case via Amit's Gaussian Gibbs-sampler theorem.
- **C_B^lin, the linear lower bound for §8.3 step 1:** 0.658, 0.724 and 0.769 for n = 2, 3, 4, volume-stable. It is below 1 and rising, so it is inconclusive.
- **Path covariance, run 3:** Prop 168 is volume-stable through β_W = 4. With Haar V the covariance is ≈ 0.
- **R_∠ volume test (§9.1):** the far-overlap floor scales like 1/V, so (D) is not refuted.

## Latest result and next steps

The L = 9 run finished: C_B^lin = 0.7685 at n = 4, against 0.7700 at L = 8, so the n = 4 value holds at about 0.769. The two volumes differ by 0.2%, more than at n = 2 and 3 because the cube is now half the box. The conclusion is unchanged: the bound rises with n but still extrapolates to below 1, so this test can't decide step 1 of the §8.3 no-go argument.

That finishes every numerical item in the plan. Before it can count as a proof, these remain:

1. **Free field, minimum over k:** prove it occurs as k → 0. So far that is only seen numerically on the scan grid.
2. **Free field, infinite volume:** write out the limiting argument that extends Amit's reduction (all observables to linear ones) to infinite volume.
3. **Yang–Mills, the main open problem:** prove the input (U) directly. The candidate route in §8.4 is a comparison between two scales, with a constant that doesn't depend on the lattice spacing a.
