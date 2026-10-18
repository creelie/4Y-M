import Mathlib

/-!
# The zonal half of the sharp quarter: the algebraic steps

Let `A = I₂(k)/I₁(k)`, which solves the Riccati equation `A' = 1 - A² - 3A/k`.
The zonal part of the chordal-Lipschitz bound is at most
`1/4 + 4A'² + A²/4 - A/(2k) - 2A²/k²`.

* `riccati_residual`, `riccati_residual_pos`, `hasDerivAt_Alo` : the function
  `A_lo = k/(3/2 + s)`, `s = √(k² + 25/4)`, is a strict subsolution of the
  Riccati equation, with residual `5(s - 5/2)/(2s(3/2+s)²)`;
* `alpha_identity`, `phi_identity` : at `A = A_lo` the bound equals
  `1/4 + (2s-5)(8s³+28s²-98s-695)/(4(2s+3)⁴)`;
* `cubic_neg`, `phi_neg` : that correction is negative for `5/2 < s ≤ √577/6`;
* `T_antitone`, `zonal_final` : for `A_lo ≤ A ≤ k/4`, `A' ≥ 0` and `9k² ≤ 88`
  the bound is strictly below `1/4`.
-/

namespace YangMills.StrongCoupling

/-- `α(s) = 1 - 3y - k²y²` at `y = 2/(2s+3)`, `k² = s² - 25/4`. -/
theorem alpha_identity (s : ℝ) (hs : 0 < s) :
    1 - 3 * (2 / (2 * s + 3)) - (s ^ 2 - 25 / 4) * (2 / (2 * s + 3)) ^ 2
      = 16 / (2 * s + 3) ^ 2 := by
  have : (2 * s + 3) ≠ 0 := by positivity
  field_simp; ring

theorem phi_identity (s : ℝ) (hs : 0 < s) :
    4 * (16 / (2 * s + 3) ^ 2) ^ 2 + (s ^ 2 - 25 / 4) * (2 / (2 * s + 3)) ^ 2 / 4
        - (2 / (2 * s + 3)) / 2 - 2 * (2 / (2 * s + 3)) ^ 2
      = (2 * s - 5) * (8 * s ^ 3 + 28 * s ^ 2 - 98 * s - 695) / (4 * (2 * s + 3) ^ 4) := by
  have : (2 * s + 3) ≠ 0 := by positivity
  field_simp; ring

theorem cubic_neg (s : ℝ) (h1 : 5 / 2 < s) (h2 : 36 * s ^ 2 ≤ 577) :
    8 * s ^ 3 + 28 * s ^ 2 - 98 * s - 695 < 0 := by
  have hs4 : s < 4.01 := by nlinarith
  nlinarith

theorem phi_neg (s : ℝ) (h1 : 5 / 2 < s) (h2 : 36 * s ^ 2 ≤ 577) :
    (2 * s - 5) * (8 * s ^ 3 + 28 * s ^ 2 - 98 * s - 695) / (4 * (2 * s + 3) ^ 4) < 0 := by
  apply div_neg_of_neg_of_pos
  · exact mul_neg_of_pos_of_neg (by linarith) (cubic_neg s h1 h2)
  · positivity

/-- Residual of the Riccati equation `A' = 1 - A² - 3A/k` along `A_lo = k/(3/2+s)`,
`s = √(k²+25/4)`, whose derivative is `(3s/2+25/4)/(s(3/2+s)²)`. -/
theorem riccati_residual (k s : ℝ) (hk : 0 < k) (hs : 0 < s) (hks : s ^ 2 = k ^ 2 + 25 / 4) :
    1 - (k / (3 / 2 + s)) ^ 2 - 3 * (k / (3 / 2 + s)) / k
        - (3 * s / 2 + 25 / 4) / (s * (3 / 2 + s) ^ 2)
      = 5 / 2 * (s - 5 / 2) / (s * (3 / 2 + s) ^ 2) := by
  have h1 : (3 / 2 + s) ≠ 0 := by positivity
  have h2 : k ≠ 0 := hk.ne'
  have h3 : s ≠ 0 := hs.ne'
  field_simp
  have : k ^ 2 = s ^ 2 - 25 / 4 := by linarith
  rw [this]; ring

theorem riccati_residual_pos (k s : ℝ) (hk : 0 < k) (hs : 0 < s)
    (hks : s ^ 2 = k ^ 2 + 25 / 4) :
    0 < 1 - (k / (3 / 2 + s)) ^ 2 - 3 * (k / (3 / 2 + s)) / k
        - (3 * s / 2 + 25 / 4) / (s * (3 / 2 + s) ^ 2) := by
  rw [riccati_residual k s hk hs hks]
  have : 5 / 2 < s := by nlinarith
  apply div_pos (by nlinarith) (by positivity)

/-- The derivative of `A_lo(k) = k/(3/2 + √(k²+25/4))`. -/
theorem hasDerivAt_Alo (k : ℝ) (hk : 0 < k) :
    HasDerivAt (fun x : ℝ => x / (3 / 2 + Real.sqrt (x ^ 2 + 25 / 4)))
      ((3 * Real.sqrt (k ^ 2 + 25 / 4) / 2 + 25 / 4) /
        (Real.sqrt (k ^ 2 + 25 / 4) * (3 / 2 + Real.sqrt (k ^ 2 + 25 / 4)) ^ 2)) k := by
  have hpos : 0 < k ^ 2 + 25 / 4 := by positivity
  have hq : HasDerivAt (fun x : ℝ => x ^ 2 + 25 / 4) (2 * k) k := by
    simpa using (hasDerivAt_pow 2 k).add_const (25 / 4 : ℝ)
  have hden : HasDerivAt (fun x : ℝ => 3 / 2 + Real.sqrt (x ^ 2 + 25 / 4))
      (2 * k / (2 * Real.sqrt (k ^ 2 + 25 / 4))) k := by
    simpa using (hq.sqrt hpos.ne').const_add (3 / 2 : ℝ)
  have hs : 0 < Real.sqrt (k ^ 2 + 25 / 4) := Real.sqrt_pos.mpr hpos
  have hne : (3 / 2 + Real.sqrt (k ^ 2 + 25 / 4)) ≠ 0 := by positivity
  refine ((hasDerivAt_id' k).div hden hne).congr_deriv ?_
  obtain ⟨s, hsdef⟩ : ∃ s, s = Real.sqrt (k ^ 2 + 25 / 4) := ⟨_, rfl⟩
  have hss : s ^ 2 = k ^ 2 + 25 / 4 := by rw [hsdef]; exact Real.sq_sqrt hpos.le
  rw [← hsdef] at hs ⊢
  have : s ≠ 0 := hs.ne'
  have : 3 / 2 + s ≠ 0 := by positivity
  field_simp
  rw [show k ^ 2 = s ^ 2 - 25 / 4 by linarith]
  ring

/-- The quadratic part `T(A) = A²/4 - A/(2k) - 2A²/k²` is decreasing on `[0, k/4]`
when `k² ≤ 12`. -/
theorem T_antitone (k A B : ℝ) (hk : 0 < k) (hk12 : k ^ 2 ≤ 12)
    (hB : 0 ≤ B) (hBA : B ≤ A) (hA : A ≤ k / 4) :
    A ^ 2 / 4 - A / (2 * k) - 2 * A ^ 2 / k ^ 2 ≤ B ^ 2 / 4 - B / (2 * k) - 2 * B ^ 2 / k ^ 2 := by
  have hk2 : 0 < k ^ 2 := by positivity
  rw [← sub_nonneg]
  have e : (B ^ 2 / 4 - B / (2 * k) - 2 * B ^ 2 / k ^ 2) - (A ^ 2 / 4 - A / (2 * k) - 2 * A ^ 2 / k ^ 2)
      = (A - B) * (2 * k - (A + B) * (k ^ 2 - 8)) / (4 * k ^ 2) := by
    field_simp; ring
  rw [e]
  apply div_nonneg _ (by positivity)
  apply mul_nonneg (by linarith)
  rcases le_or_gt (k ^ 2) 8 with h8 | h8
  · nlinarith
  · have : (A + B) * (k ^ 2 - 8) ≤ (k / 2) * (k ^ 2 - 8) :=
      mul_le_mul_of_nonneg_right (by linarith) (by linarith)
    nlinarith

/-- The final step of the zonal bound: with `A_lo ≤ A ≤ k/4` and `A' = 1 - A² - 3A/k ≥ 0`,
`4A'² + T(A) < 0` for `0 < k`, `9k² ≤ 88`. -/
theorem zonal_final (k A : ℝ) (hk : 0 < k) (hk88 : 9 * k ^ 2 ≤ 88)
    (hlo : k / (3 / 2 + Real.sqrt (k ^ 2 + 25 / 4)) ≤ A) (hhi : A ≤ k / 4)
    (hAp : 0 ≤ 1 - A ^ 2 - 3 * A / k) :
    4 * (1 - A ^ 2 - 3 * A / k) ^ 2 + (A ^ 2 / 4 - A / (2 * k) - 2 * A ^ 2 / k ^ 2) < 0 := by
  set s := Real.sqrt (k ^ 2 + 25 / 4) with hsdef
  have hpos : 0 < k ^ 2 + 25 / 4 := by positivity
  have hs : 0 < s := Real.sqrt_pos.mpr hpos
  have hss : s ^ 2 = k ^ 2 + 25 / 4 := Real.sq_sqrt hpos.le
  have hs52 : 5 / 2 < s := by nlinarith
  set Alo := k / (3 / 2 + s) with hAlo
  have hAlo0 : 0 ≤ Alo := by positivity
  -- A' decreases in A
  have hmono : 1 - A ^ 2 - 3 * A / k ≤ 1 - Alo ^ 2 - 3 * Alo / k := by
    have : 3 * Alo / k ≤ 3 * A / k := by
      apply div_le_div_of_nonneg_right _ hk.le; linarith
    nlinarith
  have hy : Alo / k = 2 / (2 * s + 3) := by
    rw [hAlo]; field_simp; ring
  have hAloy : Alo = k * (2 / (2 * s + 3)) := by
    rw [← hy]; field_simp
  have halpha : 1 - Alo ^ 2 - 3 * Alo / k = 16 / (2 * s + 3) ^ 2 := by
    rw [← alpha_identity s hs, hAloy]
    have : k ≠ 0 := hk.ne'
    have hk2 : k ^ 2 = s ^ 2 - 25 / 4 := by linarith
    field_simp
    rw [hk2]; ring
  have hT := T_antitone k A Alo hk (by linarith) hAlo0 hlo hhi
  have hTlo : Alo ^ 2 / 4 - Alo / (2 * k) - 2 * Alo ^ 2 / k ^ 2
      = (s ^ 2 - 25 / 4) * (2 / (2 * s + 3)) ^ 2 / 4 - (2 / (2 * s + 3)) / 2
          - 2 * (2 / (2 * s + 3)) ^ 2 := by
    rw [hAloy]
    have : k ≠ 0 := hk.ne'
    have hk2 : k ^ 2 = s ^ 2 - 25 / 4 := by linarith
    field_simp
    rw [hk2]; ring
  have hsq : (1 - A ^ 2 - 3 * A / k) ^ 2 ≤ (16 / (2 * s + 3) ^ 2) ^ 2 := by
    rw [← halpha]; exact pow_le_pow_left₀ hAp hmono 2
  have hphi := phi_neg s hs52 (by nlinarith)
  rw [← phi_identity s hs] at hphi
  linarith
end YangMills.StrongCoupling
