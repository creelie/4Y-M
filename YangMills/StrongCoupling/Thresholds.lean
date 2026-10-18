import Mathlib

/-!
# The Dobrushin threshold for SU(2) at the sharp quarter

For `SU(2)` the conditional law of one link is a von Mises–Fisher law on `S³`
whose staple moves by `4|β|` times the displacement of a neighbouring link.
With the sharp chordal-Lipschitz bound `1/4` (`YangMills.StrongCoupling.zonal_final`
and `Dcoef_nonneg` give its two halves), each Dobrushin coefficient is at most
`4|β| · 1/4 = |β|`. A link lies in `2(d-1)` plaquettes, each with three other
links, so a row of the Dobrushin matrix has at most `6(d-1)` nonzero entries:

* `neighbours_count` : `3 · 2(d-1) = 6(d-1)`;
* `row_sum_le` : a row sum is at most `6(d-1)|β|`;
* `alpha_lt_one_iff` : `6(d-1)|β| < 1 ↔ |β| < 1/(6(d-1))`;
* `wilson_threshold_d4` : at `d = 4`, with `β_W = 4β`, the range is `β_W < 2/9`;
* `gap_pos`, `gap_d4` : the gap `m = -log(6(d-1)|β|)` is positive and equals
  `-log(9β_W/2)` at `d = 4`.
-/

open Finset Real

namespace YangMills.StrongCoupling

/-- The coefficient bound `4|β| · (1/4) = |β|`. -/
theorem coefficient_bound (β L : ℝ) (hL : L ≤ 1 / 4) :
    4 * |β| * L ≤ |β| := by
  have := abs_nonneg β
  nlinarith

theorem neighbours_count (d : ℕ) : 3 * (2 * (d - 1)) = 6 * (d - 1) := by ring

/-- A row of the Dobrushin matrix: at most `6(d-1)` entries, each at most `|β|`. -/
theorem row_sum_le {ι : Type*} (s : Finset ι) (C : ι → ℝ) (β : ℝ) (d : ℕ)
    (hcard : s.card ≤ 6 * (d - 1)) (hC : ∀ i ∈ s, C i ≤ |β|) :
    ∑ i ∈ s, C i ≤ 6 * (d - 1 : ℕ) * |β| := by
  calc ∑ i ∈ s, C i ≤ ∑ _i ∈ s, |β| := sum_le_sum hC
    _ = s.card * |β| := by rw [sum_const, nsmul_eq_mul]
    _ ≤ 6 * (d - 1 : ℕ) * |β| := by
      apply mul_le_mul_of_nonneg_right _ (abs_nonneg β)
      exact_mod_cast hcard

theorem alpha_lt_one_iff (β : ℝ) (d : ℕ) (hd : 2 ≤ d) :
    6 * (d - 1 : ℕ) * |β| < 1 ↔ |β| < 1 / (6 * (d - 1 : ℕ)) := by
  have hpos : (0 : ℝ) < 6 * (d - 1 : ℕ) := by
    have : (1 : ℝ) ≤ (d - 1 : ℕ) := by exact_mod_cast (show 1 ≤ d - 1 by omega)
    linarith
  rw [lt_div_iff₀ hpos, mul_comm]

/-- At `d = 4` with `β_W = 4β` the Dobrushin range is `β_W < 2/9`. -/
theorem wilson_threshold_d4 (βW : ℝ) (h0 : 0 ≤ βW) :
    6 * ((4 - 1 : ℕ) : ℝ) * |βW / 4| < 1 ↔ βW < 2 / 9 := by
  rw [abs_of_nonneg (by positivity)]
  norm_num
  constructor <;> intro h <;> linarith

theorem gap_pos (β : ℝ) (d : ℕ) (hβ : β ≠ 0) (h : 6 * (d - 1 : ℕ) * |β| < 1)
    (hd : 2 ≤ d) : 0 < -Real.log (6 * (d - 1 : ℕ) * |β|) := by
  have hpos : 0 < 6 * (d - 1 : ℕ) * |β| := by
    have : (1 : ℝ) ≤ (d - 1 : ℕ) := by exact_mod_cast (show 1 ≤ d - 1 by omega)
    have := abs_pos.mpr hβ
    positivity
  have := Real.log_neg hpos h
  linarith

theorem gap_d4 (βW : ℝ) (h0 : 0 ≤ βW) :
    -Real.log (6 * ((4 - 1 : ℕ) : ℝ) * |βW / 4|) = -Real.log (9 * βW / 2) := by
  rw [abs_of_nonneg (by positivity)]
  norm_num
  ring_nf

end YangMills.StrongCoupling
