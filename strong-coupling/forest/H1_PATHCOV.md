# C1. Prop 168 path-score covariance sums: first volume scan

Status: **diagnostic run 1 of direction C** (`../CLAY_NEXT_DIRECTIONS.md`). Result: **negative on method.**
The estimator used here measures its own noise floor, not the covariance sum. Nothing about (H1) follows.

## Setup

- Code: `research/experiments/pathcov/` (local, untracked). Heat bath + overrelaxation for free links, Metropolis near paths. Self-test passed: logZ to 1e-12, analytic vs finite-difference scores to 1e-10.
- Run: b = 4, conditioning V = id, coarse L ∈ {4, 6, 8}, β_W ∈ {0.5, 1.5, 2.3}. 200 thermalization sweeps, then T = 400 measurements, one every 2 sweeps. 10 CPUs, 76 min.
- Measured quantity: S = Σ_{e≠j} ‖Cov(a_e, a_j)‖_op (3×3 path scores), split into shells by |Δ|_∞.
- Two checks are built in:
  - a null, computed from a time-permuted pairing, which has zero true covariance;
  - the same sum on each half of the run.
- At b = 4, η_m ≥ 48, so the Prop 168 threshold is undefined (NaN in the output). This run was meant to test only the shape: volume stability and the decay exponent.

## Result

| L | β_W | S | null | S − null | S_half / S | S / Var |
|---|---|---|---|---|---|---|
| 4 | 0.5 | 0.119 | 0.119 | 0.000 | 1.410 | 8.0 |
| 4 | 1.5 | 4.39 | 4.35 | 0.04 | 1.412 | 8.0 |
| 4 | 2.3 | 15.35 | 13.91 | 1.44 | 1.410 | 8.8 |
| 6 | 0.5 | 0.268 | 0.267 | 0.001 | 1.412 | 18.1 |
| 6 | 1.5 | 9.85 | 9.84 | 0.01 | 1.413 | 18.1 |
| 6 | 2.3 | 34.55 | 31.47 | 3.08 | 1.414 | 19.9 |
| 8 | 0.5 | 0.476 | 0.475 | 0.001 | 1.413 | 32.2 |
| 8 | 1.5 | 17.51 | 17.45 | 0.06 | 1.414 | 32.2 |
| 8 | 2.3 | 61.39 | 56.05 | 5.35 | 1.414 | 35.3 |

## Reading

1. **S is the noise floor.**
   - At β_W ≤ 1.5, S matches the null to within 0.4%.
   - On half the data S grows by exactly √2, which is the 1/√T scaling of pure sampling noise.
   - S/Var grows like the number of terms summed (8 → 18 → 32 for L = 4, 6, 8, i.e. ∝ L²·const at this shell cut-off).
   - Summing |Cov| over ~L⁴ pairs, each estimated with error ~Var/√T, gives a bias of ~ L⁴ Var/√T. That is what we see.
   - So the "growth with volume" here is an artefact. It says nothing about Prop 168.
2. **β_W = 2.3 has a ~10% excess over the null** (1.4, 3.1, 5.3, roughly ∝ L²). The time-permuted null also destroys autocorrelation, so this excess may be an autocorrelation effect (slower decorrelation at weak coupling) and not a true covariance. It is not evidence either way until it is checked against a blocked-null estimate.
3. The same-base shell (the 3 links sharing j's base point) is at or near the null everywhere. This is the only shell small enough for the noise to be negligible, and it is tiny: ≤ 0.04 in absolute terms, ≤ 3% of Var.

## Fix for run 2 (before b = 8, random V, or bphidden)

- **Debias per entry.** Use split-sample products Cov_1(e,j)·Cov_2(e,j), which are unbiased for Cov², or subtract a per-shell null estimated with blocked permutations. Then sum per shell and report the shell profile with errors, not Σ|Ĉov|.
- **Sum by shell, then extrapolate.** Fit the decay exponent p only on shells where the signal exceeds 3σ.
- **Use more statistics where it matters.** Use T ≥ 4000 with integrated-autocorrelation estimates, and fix one or a few j (translation average) instead of all pairs.
- **The bphidden control is still worth running.** A sum that truly grows with volume must survive the debiasing. That is the check that the new estimator can see growth at all.
