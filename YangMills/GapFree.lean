import Mathlib

/-!
# The ascent without a gap hypothesis

Finite and elementary parts of Section `ssec:onehyp` of the paper (Theorem
`thm:ascentfree`, Theorem F, Theorem `thm:gapcloses`, Corollary `cor:basescale`):
the chord bound for the exponential on `[0, 1/20]`, the recursion
`y ↦ 16 (e^y − 1)²` and its doubly exponential decay, the constants `16.83`,
`0.8415`, `0.172`, `1.9·10⁻¹¹`, the Haar measure of a small ball of `SU(2)`, and
the arithmetic of the base scale.
-/

open Real

namespace YangMills.GapFree

/-- `e^{1/20} ≤ 1.051275`, from the Taylor bound of order four. -/
theorem exp_twentieth_le : Real.exp (1 / 20) ≤ 1.051275 := by
  have h := Real.exp_bound (x := (1 / 20 : ℝ)) (by rw [abs_of_nonneg] <;> norm_num)
    (n := 4) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at h
  rw [abs_of_nonneg (by norm_num : (0 : ℝ) ≤ 1 / 20)] at h
  have h2 := (abs_le.mp h).2
  norm_num at h2 ⊢
  linarith

/-- The chord of the exponential on `[0, 1/20]`: `e^y − 1 ≤ 1.0255 y`. -/
theorem exp_sub_one_le (y : ℝ) (h0 : 0 ≤ y) (h1 : y ≤ 1 / 20) :
    Real.exp y - 1 ≤ 1.0255 * y := by
  have hc := convexOn_exp.2 (Set.mem_univ (0 : ℝ)) (Set.mem_univ (1 / 20 : ℝ))
    (by linarith : 0 ≤ 1 - 20 * y) (by linarith : 0 ≤ 20 * y)
    (by ring : (1 - 20 * y) + 20 * y = 1)
  simp only [smul_eq_mul, mul_zero, zero_add, Real.exp_zero, mul_one] at hc
  have h3 : 20 * y * (1 / 20) = y := by ring
  rw [h3] at hc
  have h4 := exp_twentieth_le
  nlinarith

/-- One doubling of the time: `16 (e^y − 1)² ≤ 16.83 y²` and `≤ 0.8415 y` on `[0, 1/20]`. -/
theorem doubling_step (y : ℝ) (h0 : 0 ≤ y) (h1 : y ≤ 1 / 20) :
    16 * (Real.exp y - 1) ^ 2 ≤ 16.83 * y ^ 2 ∧ 16 * (Real.exp y - 1) ^ 2 ≤ 0.8415 * y := by
  have h := exp_sub_one_le y h0 h1
  have hp : 0 ≤ Real.exp y - 1 := by have := Real.add_one_le_exp y; linarith
  have hsq : (Real.exp y - 1) ^ 2 ≤ (1.0255 * y) ^ 2 := pow_le_pow_left₀ hp h 2
  constructor
  · nlinarith
  · nlinarith

/-- The recursion of Theorem `thm:ascentfree`(a): if `0 ≤ y_{i+1} ≤ 16 (e^{y_i} − 1)²` and
`0 ≤ y_0 ≤ 1/20`, then `y_i ≤ 1/20` and `16.83 y_i ≤ (16.83 y_0)^(2^i)` for every `i`. -/
theorem ascent_decay (y : ℕ → ℝ) (h1 : y 0 ≤ 1 / 20)
    (hpos : ∀ i, 0 ≤ y i)
    (hrec : ∀ i, y (i + 1) ≤ 16 * (Real.exp (y i) - 1) ^ 2) :
    ∀ i, y i ≤ 1 / 20 ∧ 16.83 * y i ≤ (16.83 * y 0) ^ (2 ^ i) := by
  intro i
  induction i with
  | zero => exact ⟨h1, by simp⟩
  | succ i ih =>
    obtain ⟨hi, hu⟩ := ih
    have hs := doubling_step (y i) (hpos i) hi
    have hr := hrec i
    refine ⟨by nlinarith [hpos i], ?_⟩
    have h2 : 16.83 * y (i + 1) ≤ (16.83 * y i) ^ 2 := by nlinarith
    have h3 : (16.83 * y i) ^ 2 ≤ ((16.83 * y 0) ^ (2 ^ i)) ^ 2 :=
      pow_le_pow_left₀ (by nlinarith [hpos i]) hu 2
    calc 16.83 * y (i + 1) ≤ (16.83 * y i) ^ 2 := h2
      _ ≤ ((16.83 * y 0) ^ (2 ^ i)) ^ 2 := h3
      _ = (16.83 * y 0) ^ (2 ^ (i + 1)) := by rw [← pow_mul, pow_succ]

/-- The constants of the ascent: `16·1.0255² ≤ 16.83`, `16.83/20 = 0.8415`,
`1.0255/20 ≤ 0.061·0.8415…` in the form used, and `(16.83·13823)⁻² ≤ 1.9·10⁻¹¹`. -/
theorem ascent_constants :
    (16 : ℝ) * 1.0255 ^ 2 ≤ 16.83 ∧ (16.83 : ℝ) / 20 = 0.8415 ∧
    (1.0255 : ℝ) ≤ 0.061 * 16.83 ∧ (1 : ℝ) / (16.83 * 13823) ^ 2 ≤ 1.9e-11 := by
  refine ⟨by norm_num, by norm_num, by norm_num, by norm_num⟩

/-- `log (1/0.8415) > 0.172`: the gap `δ_* ≥ 0.172/T_0` and `Δ ≥ 0.172/t_0`. -/
theorem log_gap_constant : (0.172 : ℝ) < Real.log (1 / 0.8415) := by
  rw [Real.lt_log_iff_exp_lt (by norm_num)]
  have h := Real.exp_bound (x := (0.172 : ℝ)) (by rw [abs_of_nonneg] <;> norm_num)
    (n := 3) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at h
  rw [abs_of_nonneg (by norm_num : (0 : ℝ) ≤ 0.172)] at h
  have h2 := (abs_le.mp h).2
  norm_num at h2 ⊢
  linarith

/-- The fixed point of the recursion lies above `1/20`: `16 (e^y − 1)² < y` at `y = 1/20`,
and the map exceeds the identity at `y = 0.059`. -/
theorem fixed_point_bracket :
    16 * (Real.exp (1 / 20) - 1) ^ 2 < (1 / 20 : ℝ) ∧
    (0.059 : ℝ) < 16 * (Real.exp 0.059 - 1) ^ 2 := by
  constructor
  · have h := exp_twentieth_le
    have hp : 0 ≤ Real.exp (1 / 20) - 1 := by have := Real.add_one_le_exp (1 / 20 : ℝ); linarith
    nlinarith
  · have hq := Real.quadratic_le_exp_of_nonneg (by norm_num : (0 : ℝ) ≤ 0.059)
    have h3 : (0.059 + 0.059 ^ 2 / 2 : ℝ) ≤ Real.exp 0.059 - 1 := by linarith
    have h4 : (0.059 + 0.059 ^ 2 / 2 : ℝ) ^ 2 ≤ (Real.exp 0.059 - 1) ^ 2 :=
      pow_le_pow_left₀ (by norm_num) h3 2
    norm_num at h4 ⊢
    linarith

/-- Theorem `thm:gapcloses`: the Haar measure of a ball of `SU(2)`,
`2√2/(3π³) ≥ 1/33`. -/
theorem haar_ball_constant : (1 : ℝ) / 33 ≤ 2 * Real.sqrt 2 / (3 * Real.pi ^ 3) := by
  have hs : (1.4142 : ℝ) ≤ Real.sqrt 2 := Real.le_sqrt_of_sq_le (by norm_num)
  have hpi : Real.pi < 3.1416 := Real.pi_lt_d4
  have hpi0 : 0 < Real.pi := Real.pi_pos
  have h3 : Real.pi ^ 3 < 3.1416 ^ 3 := pow_lt_pow_left₀ hpi hpi0.le (by norm_num)
  rw [div_le_div_iff₀ (by norm_num) (by positivity)]
  norm_num at h3 ⊢
  nlinarith

/-- Theorem `thm:gapcloses`: `6 + 4 log 33 + 6 log 8 ≤ 6 · 5.42`. -/
theorem action_constant :
    6 + 4 * Real.log 33 + 6 * Real.log 8 ≤ 6 * (5.42 : ℝ) := by
  have h33 : Real.log 33 < 3.5 := by
    rw [Real.log_lt_iff_lt_exp (by norm_num)]
    have h1 : Real.exp 3.5 = Real.exp 1 ^ 3 * Real.exp 0.125 ^ 4 := by
      rw [← Real.exp_nat_mul, ← Real.exp_nat_mul, ← Real.exp_add]; norm_num
    have hq := Real.quadratic_le_exp_of_nonneg (by norm_num : (0 : ℝ) ≤ 0.125)
    have h4 : (1 + 0.125 + 0.125 ^ 2 / 2 : ℝ) ^ 4 ≤ Real.exp 0.125 ^ 4 :=
      pow_le_pow_left₀ (by norm_num) hq 4
    have he : (2.7182818283 : ℝ) ^ 3 ≤ Real.exp 1 ^ 3 :=
      pow_le_pow_left₀ (by norm_num) Real.exp_one_gt_d9.le 3
    have hm := mul_le_mul he h4 (by positivity) (by positivity)
    rw [h1]; refine lt_of_lt_of_le ?_ hm; norm_num
  have h8 : Real.log 8 = 3 * Real.log 2 := by
    rw [show (8 : ℝ) = 2 ^ 3 by norm_num, Real.log_pow]; norm_num
  have h2 := Real.log_two_lt_d9
  rw [h8]; linarith

/-- Corollary `cor:basescale`: the arithmetic of the two cases. With `ε ≤ 10⁻³` and
`T⁴ ε ≤ 3/32`: for `T ≤ 15`, `(64 T² ε)² ≤ 0.384 < 0.632²`; for `T ≥ 16`,
`T³ ε ≤ 3/(32·16) < 0.006 < 0.15`. -/
theorem base_scale_arithmetic :
    (64 : ℝ) ^ 2 * (3 / 32) * 1e-3 ≤ 0.62 ^ 2 ∧ (0.62 : ℝ) < 0.632 ∧
    (3 : ℝ) / (32 * 16) < 0.006 ∧ (9.6 : ℝ) / 64 = 0.15 ∧
    (4 : ℝ) * (1 / 20) / 96 ^ 3 < 2.27e-7 := by
  refine ⟨by norm_num, by norm_num, by norm_num, by norm_num, by norm_num⟩

/-- Corollary `cor:basescale`: in the first case the bound `64 T² ε ≤ 0.62` follows from
`T⁴ ε ≤ 3/32` and `ε ≤ 10⁻³`. -/
theorem base_scale_small (T ε : ℝ) (hε : 0 ≤ ε) (h1 : ε ≤ 1e-3)
    (h2 : T ^ 4 * ε ≤ 3 / 32) : 64 * T ^ 2 * ε ≤ 0.62 := by
  have hsq : (64 * T ^ 2 * ε) ^ 2 ≤ 0.62 ^ 2 := by
    have : (64 * T ^ 2 * ε) ^ 2 = 64 ^ 2 * (T ^ 4 * ε) * ε := by ring
    rw [this]
    have h3 : 0 ≤ T ^ 4 * ε := by positivity
    nlinarith
  exact (pow_le_pow_iff_left₀ (by positivity) (by norm_num) two_ne_zero).mp hsq

end YangMills.GapFree
