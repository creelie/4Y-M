# PR #290: progress summary

Branch `research/next-steps-3`, 2026-11-09. Question: does a **uniform block-Poincaré constant (U)** hold? (U) says that resampling blocks of fixed physical side ℓ mixes the field at a rate ϖ that does not degrade as the lattice spacing a → 0, so that a gap ∝ ϖ/ℓ carries over to the continuum.

## Results

### Free field (theory, settled)
- Step 2 closed: chaos reduction by operator positivity, in any volume.
- The k-infimum sits at k → 0. This was checked numerically (35 groups, monotone in |k|) and **proved uniformly up to a factor 5**: ϖ_k ≥ ϖ_0(1+λ)/(1+Θλ) with Θ ≤ 5. The best constant found numerically is Θ* ≈ 1.06.
- The two-scale comparison is uniform in the cutoff (C → ≈ 0.97).
- Write-up: `forest/A_GAUSSIAN_FIBRE_4D.md` §10.1–10.3.

### Interacting, direct Monte Carlo of ϖ (3D)
Estimator: nested conditional heat-bath resampling of non-adjacent cubes, with a jackknife error estimate. Code and raw data: `forest/experiments/{u1bp,su2bp}/`.

| Theory | Result | Section |
|---|---|---|
| Compact U(1)₃ | ϖ ∝ (ma)^0.75 → 0 at fixed mℓ: **(U) fails**. Caveat: this theory has two scales, so it is not a single-scale limit. This shows the test can detect a failure. | §10.4 |
| SU(2)₃, β 5–15 | ϖ = O(1) with a slow coherent drift | §10.5–10.6 |
| SU(2)₃, fixed box Lg² = 12.8, β 5–25 | ϖ **levels off at ≈ 1.0** | §10.8 |

The SU(2)₃ fixed-box table (observable mean, ± ≈ 0.02):

| ℓg² | β = 5 | 10 | 15 | 20 | 25 | c₀ (c₀ + c₁a) |
|---|---|---|---|---|---|---|
| 1.6 | 1.44 | 1.27 | 1.22 | 1.19 | 1.13 | 1.080(17) |
| 2.4 | 1.42 | 1.23 | 1.15 | 1.14 | 1.09 | 1.022(16) |
| 3.2 | 1.40 | 1.18 | 1.10 | 1.09 | 1.06 | 0.975(15) |
| 4.8 | 1.29 | 1.12 | 1.08 | 1.08 | 1.04 | 0.983(15) |

- **The box size does not explain the drift.** At the same β, L = 16 and L = 32 agree within errors, so the drift is a lattice-spacing effect.
- **A plateau fits better than decay to zero.** c₀ + c₁a gives χ² = 9.4 on 15 dof, while c·aᵖ gives χ² = 31.7 on 15 dof (p-value ≈ 0.007).
- **The pre-registered ratio rule is inconclusive.** It gives β25/β10 = 0.89–0.92, between the 0.85 and 0.95 thresholds. The plateau fit predicts this value.
- **Status:** evidence for (U) in a single-scale gauge theory, not a proof.

### Erratum (§10.7)
The jackknife in `u1bp/analyze.py` and `su2bp/analyze.py` overstated errors by √20. This is now fixed, and the stored errors are corrected. Central values are unchanged. The fix strengthens both the U(1) failure and the SU(2) plateau preference.

## Not in this PR (next, separate PR)
- 4D SU(2) analogue of the fixed-box test, at modest β.
- Sharpen the free-field constant Θ from 5 to 1. The current trial-field method cannot reach this, since it is limited to Θ* ≈ 1.06.

## Where things are
- `strong-coupling/forest/A_GAUSSIAN_FIBRE_4D.md` §10 has the full write-up.
- `strong-coupling/CLAY_NEXT_DIRECTIONS.md` holds the dated log and next steps.
- `strong-coupling/forest/experiments/` holds the code, launchers, and raw `results*.jsonl`.
