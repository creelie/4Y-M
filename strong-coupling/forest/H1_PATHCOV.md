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

## Run 2 (debiased shell profile)

- The estimator is per-shell split-sample covariance with errors (shellQ ± err, Z = Q/err).
- It was run at b ∈ {4, 6}, L ∈ {4, 6, 8} (b = 6 at L = 4 only), β_W ∈ {1.5, 2.3}, T = 3000, with 10 CPUs for 7.9 h.

**Control first (passes).**
- The planted-covariance control at L = 4, 6, 8 behaves as it should.
  - With amp = 0 every shell has |Z| ≤ 1.1.
  - With amp = 0.1 every shell recovers the planted value to within 1–3%, at Z = 19–130.
  - This includes the growing shells: 6.29 vs 6.39 truth at L = 8.
- So the estimator is unbiased and can see volume growth. This also does the job of the bphidden check: a real growing sum would have shown up.

| b | L | β_W | shell Z (|Δ|_∞ = 0, 1, 2, …) | Σ shells, significant part |
|---|---|---|---|---|
| 4 | 4 | 1.5 | −0.5, −0.2, −0.3 | 0 |
| 4 | 6 | 1.5 | −0.2, 1.0, 0.1, 0.6 | 0 |
| 4 | 8 | 1.5 | −0.1, 1.4, 0.4, 1.5, −0.9 | 0 |
| 4 | 4 | 2.3 | 4.3, 2.3, 0.9 | 0.025 |
| 4 | 6 | 2.3 | 11.0, 9.5, −0.8, −1.0 | 0.56 |
| 4 | 8 | 2.3 | 21.9, 24.3, 0.1, 1.0, −0.6 | 0.58 |
| 6 | 4 | 1.5 | 0.5, −0.1, −0.1 | 0 |
| 6 | 4 | 2.3 | −0.5, −0.6, −0.4 | 0 |

**Reading.**
1. **The run-1 growth was pure estimator bias.** After debiasing, β_W = 1.5 shows no covariance in any shell at any L.
2. **β_W = 2.3, b = 4 has a real but short-range covariance.**
   - Shell 0 is 2.0e-4 and shell 1 is 9.5e-4, and both values agree between L = 6 and L = 8.
   - Shells ≥ 2 are zero within 1σ.
   - The significant sum, 0.56 → 0.58, is volume-stable. (L = 4 is too small, since its shells wrap.)
   - This is consistent with Prop 168: the path-score covariance sum is bounded uniformly in volume with finite range ≈ 1 block.
3. **At b = 6 the covariance vanishes at both β_W.** Coarser blocking kills even the shell-1 term at β_W = 2.3.

**Next for C.**
- Depth scan in β_W ∈ {2.3, 3.0, 4.0} at b = 4, L = 6, 8, to see where shell 2 turns on.
- Random-V stress test at b = 4, β_W = 2.3.
- b = 8 only if the depth scan shows growth.

## Run 3 (depth in β_W, b = 4)

| L | β_W | V | shellQ[0..2] | CSbound | CSbound_sig |
|---|---|---|---|---|---|
| 6 | 2.3 | Haar | ≈0 (all |Z|<1.4) | 0 | 0 |
| 8 | 2.3 | Haar | ≈0 (all |Z|<1) | 0.22 | 0 |
| 6 | 3.0 | id | .095, .364, .0034 | 14.05 | 14.05 |
| 8 | 3.0 | id | .094, .363, .0032 | 15.05 | 13.96 |
| 6 | 4.0 | id | .292, 1.129, .0101 | 25.79 | 24.62 |
| 8 | 4.0 | id | .292, 1.134, .0114 | 28.92 | 24.96 |

The significant shell sums agree between L = 6 and L = 8 to within errors, so under the stated decision rule **Prop 168 is volume-stable through β_W = 4**. The raw CSbound creeps up with L only through insignificant far shells; CSbound_sig is flat. With Haar-random V the covariance is consistent with zero at every shell. The bound grows with β_W (14 → 25), which matches weaker decay toward weak coupling. It does not show a volume divergence.
