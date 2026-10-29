import Mathlib

/-!
# Perturbations of the cube condition: the finite arithmetic

The proposition `prop:cubeopen` of the forest-gauge paper shows that the cube
condition `D(n, α)` survives a finite-range perturbation of the action with
Lipschitz constant `ε`, with `α' = α e^Ω (1 + 2Ω) + 96 s Ω` and
`Ω = 48 M D N_n ε`; the proposition `prop:smallfield` shows that, for
gauge-invariant perturbations in the small-field regime, the cube must satisfy
`m(m+1) η + κ D ≥ (1 - α) D / 3`. The analytic steps (Kantorovich duality,
the loop bound) are in the paper; this file checks the arithmetic that joins
them:

* `cube_constants` : `48 · 2 = 96`, `96 · 48 = 4608`, `4608 = 4 · 1152`;
* `lipschitz_times_diameter` : `Δ_V M ε = Ω` with `Δ_V = 48 N_n D`;
* `coefficient_total` : `2 Δ_V s M |V| ε = 96 s Ω N_n` with `|V| = 48 N_n`;
* `common_tilt` : the two terms of the common tilt add up to
  `e^Ω (1 + 2Ω) W`;
* `two_tilts` : `1 - e^{-2η} ≤ 2η`;
* `alpha_prime_ge` : `α' ≥ α + 96 s Ω`;
* `eps_lt_of_alpha_prime_lt_one` : `α' < 1` forces
  `ε < (1 - α) / (4608 s M D N_n)`;
* `cube_side_sq_le` : `m(m+1) ≤ 2 m²` for `m ≥ 1`;
* `cube_volume_ge` : `m(m+1) ≥ x ≥ 0` gives `m⁴ ≥ x² / 4`;
* `small_field_eps` : with `x = c D / η` the admissible `ε` is below
  `(1 - α) η² / (1152 s M D³ c²)`;
* `small_field_side` : the two bounds on `ε_V` combine to
  `m(m+1) η + κ D ≥ (1 - α) D / 3`.
-/

open Real

namespace YangMills.StrongCoupling

theorem cube_constants : 48 * 2 = 96 ∧ 96 * 48 = 4608 ∧ 4608 = 4 * 1152 := by
  norm_num

/-- With `Δ_V = 48 N_n D` and `Ω = 48 M D N_n ε`, `Δ_V M ε = Ω`. -/
theorem lipschitz_times_diameter (N D M ε : ℝ) :
    (48 * N * D) * M * ε = 48 * M * D * N * ε := by ring

/-- The total of the second terms: `2 Δ_V s M |V| ε = 96 s Ω N_n`. -/
theorem coefficient_total (N D M s ε : ℝ) :
    2 * (48 * N * D) * s * M * (48 * N) * ε = 96 * s * (48 * M * D * N * ε) * N := by
  ring

/-- The common tilt: the first term is at most `e^Ω (1 + Ω) W`, the second at
most `Ω e^Ω W`, and together at most `e^Ω (1 + 2Ω) W`. -/
theorem common_tilt (x y W Ω : ℝ)
    (hx : x ≤ exp Ω * (1 + Ω) * W) (hy : y ≤ Ω * exp Ω * W) :
    x + y ≤ exp Ω * (1 + 2 * Ω) * W := by
  nlinarith

/-- Two tilts of one measure: `1 - e^{-2η} ≤ 2η`. -/
theorem two_tilts (η : ℝ) : 1 - exp (-2 * η) ≤ 2 * η := by
  have := add_one_le_exp (-2 * η)
  linarith

/-- `α' = α e^Ω (1 + 2Ω) + 96 s Ω ≥ α + 96 s Ω`. -/
theorem alpha_prime_ge (α Ω s : ℝ) (hα : 0 ≤ α) (hΩ : 0 ≤ Ω) :
    α + 96 * s * Ω ≤ α * exp Ω * (1 + 2 * Ω) + 96 * s * Ω := by
  have h1 : 1 ≤ exp Ω := one_le_exp hΩ
  have h2 : 1 ≤ exp Ω * (1 + 2 * Ω) := by nlinarith
  nlinarith

/-- If `α' < 1` then `ε < (1 - α) / (4608 s M D N_n)`. -/
theorem eps_lt_of_alpha_prime_lt_one (α s M D N ε : ℝ) (hα : 0 ≤ α) (hε : 0 ≤ ε)
    (hs : 0 < s) (hM : 0 < M) (hD : 0 < D) (hN : 0 < N)
    (h : α * exp (48 * M * D * N * ε) * (1 + 2 * (48 * M * D * N * ε))
        + 96 * s * (48 * M * D * N * ε) < 1) :
    ε < (1 - α) / (4608 * s * M * D * N) := by
  have hΩ : 0 ≤ 48 * M * D * N * ε := by positivity
  have hge := alpha_prime_ge α (48 * M * D * N * ε) s hα hΩ
  have hpos : 0 < 4608 * s * M * D * N := by positivity
  rw [lt_div_iff₀ hpos]
  nlinarith

/-- `m(m+1) ≤ 2 m²` for `m ≥ 1`. -/
theorem cube_side_sq_le (m : ℕ) (hm : 1 ≤ m) : m * (m + 1) ≤ 2 * m ^ 2 := by
  nlinarith

/-- If `m(m+1) ≥ x ≥ 0` then `N_n = m⁴ ≥ x² / 4`. -/
theorem cube_volume_ge (m : ℕ) (hm : 1 ≤ m) (x : ℝ) (hx : 0 ≤ x)
    (h : x ≤ (m : ℝ) * (m + 1)) : x ^ 2 / 4 ≤ (m : ℝ) ^ 4 := by
  have h2 : (m : ℝ) * (m + 1) ≤ 2 * (m : ℝ) ^ 2 := by
    exact_mod_cast cube_side_sq_le m hm
  have h3 : x ≤ 2 * (m : ℝ) ^ 2 := le_trans h h2
  have h4 : x ^ 2 ≤ (2 * (m : ℝ) ^ 2) ^ 2 := pow_le_pow_left₀ hx h3 2
  nlinarith

/-- The admissible perturbation on a small-field cube: if
`ε < (1 - α) / (4608 s M D N_n)` and `N_n ≥ x² / 4` with `x = c D / η`, then
`ε < (1 - α) η² / (1152 s M D³ c²)`. -/
theorem small_field_eps (α s M D N ε c η : ℝ) (hα : α < 1)
    (hs : 0 < s) (hM : 0 < M) (hD : 0 < D) (hc : 0 < c) (hη : 0 < η)
    (hN : (c * D / η) ^ 2 / 4 ≤ N)
    (hε : ε < (1 - α) / (4608 * s * M * D * N)) :
    ε < (1 - α) * η ^ 2 / (1152 * s * M * D ^ 3 * c ^ 2) := by
  have hx : 0 < (c * D / η) ^ 2 / 4 := by positivity
  have hNpos : 0 < N := lt_of_lt_of_le hx hN
  have h1 : (1 - α) / (4608 * s * M * D * N)
      ≤ (1 - α) / (4608 * s * M * D * ((c * D / η) ^ 2 / 4)) := by
    apply div_le_div_of_nonneg_left (by linarith) (by positivity)
    have : 0 < 4608 * s * M * D := by positivity
    exact mul_le_mul_of_nonneg_left hN this.le
  have h2 : (1 - α) / (4608 * s * M * D * ((c * D / η) ^ 2 / 4))
      = (1 - α) * η ^ 2 / (1152 * s * M * D ^ 3 * c ^ 2) := by
    field_simp
    ring
  linarith [h1, h2.le, h2.ge]

/-- The small-field bound: `(1 - α) D / 3 ≤ ε_V ≤ m(m+1) η + κ D`. -/
theorem small_field_side (α D εV m η κ : ℝ)
    (hlow : (1 - α) * D / 3 ≤ εV) (hup : εV ≤ m * (m + 1) * η + κ * D) :
    (1 - α) * D / 3 ≤ m * (m + 1) * η + κ * D :=
  le_trans hlow hup

end YangMills.StrongCoupling
