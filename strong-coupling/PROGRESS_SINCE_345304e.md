# Progress since 345304e

Branch `research/next-steps-2`. These are the next steps from `CLAY_NEXT_DIRECTIONS.md`. Details are in `forest/A_GAUSSIAN_FIBRE_4D.md` §7–9 and `forest/H1_PATHCOV.md`.

## Direction A: abelian Gaussian fibre model

| Item | Result | Where |
|---|---|---|
| Lemma A4 | (D) on all of L² reduces to the first-chaos angle sum R_∠ = sup_x Σ_y ‖P_{V_x}P_{V_y}‖, using the Fock exponential law and polar decomposition. | §7 |
| ρ₁ (hypothesis X) | Measured 0.79–0.86, flat in n'/b (0.789, 0.828, 0.835 at n'/b = 1, 2, 3), volume-stable to 1e-3, offset spread ±0.02. | §8 |
| B under a flat ρ | Π C_B ≈ n_k^α with α = log₂(1/(1−ρ²)) ≈ 1.7. | §8.1 |
| (B′) | Cor bpuptuning needs Γ bounded, so polynomial growth empties the gap bound. This is forced: the model is free Maxwell (gapless). It tests (F), (X) and (D) per scale only, and the decay of ρ_j must come from non-abelian or confining input. | §8.2 |
| R_∠, L = 4 | 4.08 at (b, L, n') = (3, 4, 3) and 3.84 at (4, 4, 4). Face overlap 0.155, edge 0.07, and a 2–5e-3 floor at d_∞ = 2 (power-law signature). | §9 |
| R_∠, L = 6, n'/b = 2 | 11.2 at (3, 6, 6) and 10.1 at (4, 6, 8), insensitive to b. With m = 3 every cube is a neighbour, and dim C' = 63. | this file |
| Numerical fix | Removing the torus-harmonic part of the path currents makes the Pythagoras check exact (8e-14, previously 9e-4). ρ₁ moves by 3e-4, so §8 stands. | §9 |

**Open:** the volume comparison of R_∠ at fixed n'/b = 1, (3, 4, 3) vs (3, 6, 3), is still running. Bounded means (D) holds in the Gaussian model. Growth with L, which happens if the overlap decay exponent is ≤ 4, means (D) fails as framed and Lemma A4 needs a summable replacement.

## Direction C: path covariance (Prop 168)

- **Run 2 (debiased shell profile):**
  - The control passes: 1–3% recovery, null |Z| ≤ 1.1.
  - At β_W = 2.3, b = 4, only shells 0–1 are nonzero, and their sum (0.56 → 0.58 from L = 6 to 8) is volume-stable.
  - β_W = 1.5 and b = 6 show no covariance.
- **Run 3 (running):** β_W ∈ {3.0, 4.0} with V = id, plus Haar-random V at β_W = 2.3, at b = 4 and L = 6, 8. Equal shell sums at L = 6 and 8 would mean Prop 168 is volume-stable through β_W = 4.

## Net status

- A: (F) and (X) are measured, B is settled negatively for the Gaussian model, and (D) awaits the R_∠ volume test.
- C: short-range and volume-stable at β_W = 2.3, with the depth scan pending.
- The two pending results will be added in a follow-up.
