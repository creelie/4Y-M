# The true feedback of the SU(2) plaquette block (direction A2)

**Question.** Lemma 13 of `forest/main.tex` bounds the response of a plaquette block by

    K = (1+β) s / (1 − C),   C = 4βsA + (320/3) β³ M e^w,

and the Dobrushin column sum is α = 40βK = 10β_W K. At β_W = 0.351 this gives s = 0.25168, C = 0.03904, K = 0.28489 and
α = 0.99995. `DIRECTIONS.md` (A2) assumed C was about twice the true feedback. How large is the true response?

**Answer.**
- The true response is K_true ≈ (1+β)/4 = 0.2719 at every outside field tested.
- The proof's K is 4.8% higher. Almost all of the excess is the factor 1/(1 − C), which costs 4.1%. The single-link
  constant s = 0.25168 against the true 0.2500 costs 0.7%.
- The outside fields hardly change the response: the feedback the proof pays for with 1/(1 − C) is not visible in
  the data.
- If K could be brought to (1+β)/4, the column sum would reach 1 at β_W = 0.366, the first-order ceiling of A3. A
  bound within 3.4% of the truth reaches 0.355.
- Correction to `DIRECTIONS.md`:
  - The conditional feedback itself is not overestimated by a factor 2. βA(20β) = 0.0344 against the bound 0.0390.
  - The slack is that the true response does not grow by 1/(1 − C) at all: the y₁ → y₃ → y₁ loop that the Neumann
    series charges for cancels in the covariance.

## Method (`response.py`)

The block law on (S³)³ is μ_b(dy) ∝ exp(Σᵢ⟨yᵢ, bᵢ⟩ + 4β Re(y₁y₂y₃)). The derivative of μ_{b₁+th}(f) in t is
Cov_μ(f, ⟨y₁, h⟩). The response K_true(b) is the supremum of this covariance over functions f that are 1-Lipschitz for
Σᵢ|yᵢ − yᵢ′|. Each test family gives a lower bound on it:

| family | f | 1-Lipschitz because |
|---|---|---|
| A | ⟨y₁,k₁⟩ + ⟨y₂,k₂⟩ + ⟨y₃,k₃⟩ | sum of the three marginals |
| B | ⟨y₁,k₁⟩ + ⟨y₂y₃, k⟩ | \|d(y₂y₃)\| ≤ \|dy₂\| + \|dy₃\| |
| C | ⟨y₃,k₃⟩ + ⟨y₁y₂, k⟩ | the same |
| D | ⟨y₁y₂y₃, k⟩ | the same |

The best unit vectors are the normalised covariances, so each family's bound is the sum of its covariance norms.
Sampling uses exact heat-bath sweeps over the three links with von Mises–Fisher rejection on S³: 10⁵ parallel chains,
100 sweeps, 30 sweeps of burn-in. The outside fields are 0 or of the maximal size |b| = R = 20β allowed in the lemma.

## Results at β_W = 0.351 (β = 0.08775, `response_0.351.out`)

| outside fields | A | **B** | C | D | ⟨y₁,h⟩ alone |
|---|---|---|---|---|---|
| all zero | 0.2504–0.2505 | **0.2719–0.2720** | 0.0003 | 0.0002 | 0.2501 |
| b₂ maximal | 0.2588–0.2589 | **0.2717–0.2719** | 0.1065 | 0.0002 | 0.2500 |
| b₂, b₃ maximal | 0.2619–0.2654 | **0.2703–0.2717** | 0.104–0.106 | 0.038 | 0.2500 |
| b₂, b₃ maximal, aligned | 0.2619–0.2651 | **0.2703–0.2715** | 0.104–0.105 | 0.038 | 0.2500 |

Each range runs over the three directions h = e₀, e₁, (e₀+e₁)/√2. The last configuration has two of its three
directions; the run was cut off by a container restart. With 10⁵ chains the noise is about 0.003.

- Family B, the route y₁ → y₂y₃ that produces the (1+β) of the lemma, gives (1+β)/4 = 0.27194 to within noise in every
  configuration. The largest outside fields change it by at most about 0.6%.
- The single-link response is 1/4, as it must be: the proof's s = 0.25168 is the certified bound for the single-link
  constant over |b| ≤ R.

## What a proof would need

The present proof bounds the block response by a single-link response, times (1+β) for the y₂y₃ route, times a
Neumann series for the y₁ → y₃ → y₁ feedback. The data say the true response is (1+β)/4 + O(β²) uniformly in the
outside fields. A proof of that form needs the expansion of Cov_μ(f, ⟨y₁,h⟩) to second order in the plaquette coupling
4β, with a remainder bound uniform over |bᵢ| ≤ R:
- Expand μ_b around the product measure ⊗ᵢ ν_{bᵢ}. The first order gives the (1+β) term exactly.
- The second order, of size (4β)² ≈ 0.12 times a moment of the vMF laws, is where the bound must beat the present
  4.1%.
- The remainder can be bounded by the Poincaré inequality of the block law, which is already certified in
  `blocks/sharp/`.

Target: K ≤ 0.2817 at β_W = 0.355, i.e. 10β_W K < 1. The truth is about 0.2722 there, so there is a 3.4% margin for
the second-order and remainder terms. This is a finite computation of the same kind as `blocks/sharp/`, but it is not
done here.

## Files

| File | Contents |
|---|---|
| `response.py` | Monte Carlo of the block response (`python3 response.py 0.351 100000 100`) |
| `response_0.351.out` | Output of that run (11 of the 12 configurations) |
