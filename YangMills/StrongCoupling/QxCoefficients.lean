import Mathlib

/-!
# The transverse half of the sharp quarter: coefficient positivity

The transverse part of the chordal-Lipschitz bound for the von Mises–Fisher law
on `S³` reduces to the positivity of
`D(k) = (8 + 2k²) I₁(k)² - k² I₀(k)² - k I₁(2k)`,
a power series in `k²`. Writing `I₀`, `I₁` as their series, the coefficients of
`I₀²` and `I₁²` are Cauchy products (`a`, `b`), and that of `k I₁(2k)` is `e`.

* `a_closed`, `b_closed` : the Cauchy products in closed form, by Vandermonde;
* `Dcoef_closed` : the coefficient of `k^{2(m+2)}` in `D` times a positive factor
  equals `C(2m+2,m+1)(4(2m+3) + m(m+2)(m+3)) - 4^{m+1}(m+2)(m+3)`;
* `coeff_ineq` : that number is `≥ 0` for every `m`, from
  `16^m ≤ 4m C(2m,m)²` (`sixteen_pow_le`);
* `Dcoef_nonneg` : every coefficient of `D` is nonnegative.
-/

open Nat Finset

namespace YangMills.StrongCoupling

theorem sixteen_pow_le (m : ℕ) (hm : 1 ≤ m) :
    16 ^ m ≤ 4 * m * centralBinom m ^ 2 := by
  induction m, hm using Nat.le_induction with
  | base => simp [Nat.centralBinom, Nat.choose]
  | succ m hm ih =>
    have h := Nat.succ_mul_centralBinom_succ m
    have key : 16 ^ m * (m + 1) ≤ (2 * m + 1) ^ 2 * centralBinom m ^ 2 := by
      calc 16 ^ m * (m + 1) ≤ 4 * m * centralBinom m ^ 2 * (m + 1) :=
            Nat.mul_le_mul_right _ ih
        _ ≤ (2 * m + 1) ^ 2 * centralBinom m ^ 2 := by nlinarith [Nat.zero_le (centralBinom m ^ 2)]
    have hpos : 0 < m + 1 := by omega
    have : 16 ^ (m + 1) * (m + 1) ≤ 4 * (m + 1) * centralBinom (m + 1) ^ 2 * (m + 1) := by
      have e : 4 * (m + 1) * centralBinom (m + 1) ^ 2 * (m + 1)
          = 4 * ((m + 1) * centralBinom (m + 1)) ^ 2 := by ring
      rw [e, h, pow_succ]
      nlinarith [key]
    exact Nat.le_of_mul_le_mul_right this hpos

theorem central_mul_ge (m : ℕ) (hm : 6 ≤ m) :
    4 ^ m ≤ centralBinom m * (m - 1) := by
  have h1 := sixteen_pow_le m (by omega)
  have h2 : 4 * m ≤ (m - 1) ^ 2 := by
    obtain ⟨k, rfl⟩ : ∃ k, m = k + 6 := ⟨m - 6, by omega⟩
    simp only [show k + 6 - 1 = k + 5 by omega]; nlinarith
  have : (4 ^ m) ^ 2 ≤ (centralBinom m * (m - 1)) ^ 2 := by
    calc (4 ^ m) ^ 2 = 16 ^ m := by rw [← pow_mul, mul_comm, pow_mul]; norm_num
      _ ≤ 4 * m * centralBinom m ^ 2 := h1
      _ ≤ (m - 1) ^ 2 * centralBinom m ^ 2 := Nat.mul_le_mul_right _ h2
      _ = (centralBinom m * (m - 1)) ^ 2 := by ring
  exact (Nat.pow_le_pow_iff_left (by norm_num)).1 this

/-- The coefficient inequality: for every `m`,
`4^m (m+1)(m+2) ≤ C(2m,m) (4(2m+1) + (m-1)(m+1)(m+2))`. -/
theorem coeff_ineq (m : ℕ) :
    (4 : ℤ) ^ m * ((m + 1) * (m + 2)) ≤
      (centralBinom m : ℤ) * (4 * (2 * m + 1) + ((m : ℤ) - 1) * (m + 1) * (m + 2)) := by
  rcases lt_or_ge m 6 with h | h
  · interval_cases m <;> simp [Nat.centralBinom, Nat.choose]
  · have h3 := central_mul_ge m h
    have h3' : (4 : ℤ) ^ m ≤ (centralBinom m : ℤ) * ((m : ℤ) - 1) := by
      have : ((m - 1 : ℕ) : ℤ) = (m : ℤ) - 1 := by omega
      rw [← this]; exact_mod_cast h3
    have hc : (0 : ℤ) ≤ centralBinom m := by positivity
    have hp : (0 : ℤ) < (m + 1) * (m + 2) := by positivity
    nlinarith [mul_le_mul_of_nonneg_right h3' hp.le, mul_nonneg hc (show (0:ℤ) ≤ 4 * (2 * m + 1) by positivity)]

/-- Coefficient of `k^{2m}` in `I₀(k)²`, as a Cauchy product of the series
`I₀(k) = Σ (k/2)^{2j}/(j!)²`. -/
noncomputable def a (m : ℕ) : ℚ :=
  (1 / 4 ^ m) * ∑ j ∈ range (m + 1), 1 / ((j ! : ℚ) ^ 2 * ((m - j)! : ℚ) ^ 2)

/-- Coefficient of `k^{2m+2}` in `I₁(k)²`, from `I₁(k) = Σ (k/2)^{2j+1}/(j!(j+1)!)`. -/
noncomputable def b (m : ℕ) : ℚ :=
  (1 / 4 ^ (m + 1)) *
    ∑ j ∈ range (m + 1), 1 / ((j ! : ℚ) * (j + 1)! * (m - j)! * (m - j + 1)!)

theorem a_closed (m : ℕ) : a m = (centralBinom m : ℚ) / (4 ^ m * (m ! : ℚ) ^ 2) := by
  have hs : ∑ j ∈ range (m + 1), 1 / ((j ! : ℚ) ^ 2 * ((m - j)! : ℚ) ^ 2)
      = ∑ j ∈ range (m + 1), ((m.choose j : ℕ) : ℚ) ^ 2 / (m ! : ℚ) ^ 2 := by
    refine sum_congr rfl fun j hj => ?_
    have hj : j ≤ m := Nat.lt_succ_iff.mp (mem_range.mp hj)
    rw [Nat.cast_choose ℚ hj]
    have h1 : (j ! : ℚ) ≠ 0 := by positivity
    have h2 : ((m - j)! : ℚ) ≠ 0 := by positivity
    have h3 : (m ! : ℚ) ≠ 0 := by positivity
    field_simp
  have hc : (∑ j ∈ range (m + 1), ((m.choose j : ℕ) : ℚ) ^ 2) = centralBinom m := by
    rw [Nat.centralBinom_eq_two_mul_choose, ← Nat.sum_range_choose_sq]; push_cast; rfl
  rw [a, hs, ← sum_div, hc]
  field_simp

theorem vandermonde_shift (m : ℕ) :
    ∑ j ∈ range (m + 1), (m.choose j) * ((m + 2).choose (j + 1)) = centralBinom (m + 1) := by
  have hv : centralBinom (m + 1) =
      ∑ k ∈ range (m + 1 + 1), m.choose k * (m + 2).choose (m + 1 - k) := by
    rw [Nat.centralBinom_eq_two_mul_choose, show 2 * (m + 1) = m + (m + 2) by ring,
      Nat.add_choose_eq, Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk]
  rw [hv, Finset.sum_range_succ _ (m + 1), Nat.choose_succ_self, zero_mul, add_zero]
  refine sum_congr rfl fun j hj => ?_
  have hj : j ≤ m := Nat.lt_succ_iff.mp (mem_range.mp hj)
  congr 1
  rw [← Nat.choose_symm (show j + 1 ≤ m + 2 by omega)]
  congr 1; omega

theorem b_closed (m : ℕ) :
    b m = (centralBinom (m + 1) : ℚ) / (4 ^ (m + 1) * (m ! : ℚ) * ((m + 2)! : ℚ)) := by
  have hs : ∑ j ∈ range (m + 1), 1 / ((j ! : ℚ) * (j + 1)! * (m - j)! * (m - j + 1)!)
      = ∑ j ∈ range (m + 1),
          ((m.choose j : ℕ) : ℚ) * (((m + 2).choose (j + 1) : ℕ) : ℚ) / ((m ! : ℚ) * (m + 2)!) := by
    refine sum_congr rfl fun j hj => ?_
    have hj : j ≤ m := Nat.lt_succ_iff.mp (mem_range.mp hj)
    rw [Nat.cast_choose ℚ hj, Nat.cast_choose ℚ (show j + 1 ≤ m + 2 by omega),
      show m + 2 - (j + 1) = m - j + 1 by omega]
    have h1 : (j ! : ℚ) ≠ 0 := by positivity
    have h1' : ((j + 1)! : ℚ) ≠ 0 := by positivity
    have h2 : ((m - j)! : ℚ) ≠ 0 := by positivity
    have h2' : ((m - j + 1)! : ℚ) ≠ 0 := by positivity
    have h3 : (m ! : ℚ) ≠ 0 := by positivity
    have h4 : ((m + 2)! : ℚ) ≠ 0 := by positivity
    field_simp
  have hc : (∑ j ∈ range (m + 1), ((m.choose j : ℕ) : ℚ) * (((m + 2).choose (j + 1) : ℕ) : ℚ))
      = centralBinom (m + 1) := by
    rw [← vandermonde_shift]; push_cast; rfl
  rw [b, hs, ← sum_div, hc]
  field_simp

/-- Coefficient of `k^{2m+2}` in `k I₁(2k) = Σ k^{2m+2}/(m!(m+1)!)`. -/
noncomputable def e (m : ℕ) : ℚ := 1 / ((m ! : ℚ) * (m + 1)!)

/-- Coefficient of `k^{2n}` in `D(k) = (8 + 2k²) I₁(k)² - k² I₀(k)² - k I₁(2k)`. -/
noncomputable def Dcoef : ℕ → ℚ
  | 0 => 0
  | 1 => 8 * b 0 - a 0 - e 0
  | m + 2 => 8 * b (m + 1) + 2 * b m - a (m + 1) - e (m + 1)

theorem Dcoef_one : Dcoef 1 = 0 := by
  simp [Dcoef, a, b, e]; norm_num

theorem Dcoef_closed (m : ℕ) :
    Dcoef (m + 2) * (4 ^ (m + 1) * ((m + 1)! : ℚ) * ((m + 2)! : ℚ) * ((m + 2) * (m + 3)))
      = (centralBinom (m + 1) : ℚ) * (4 * (2 * m + 3) + m * (m + 2) * (m + 3))
        - 4 ^ (m + 1) * ((m + 2) * (m + 3)) := by
  have hc := Nat.succ_mul_centralBinom_succ (m + 1)
  have hc' : ((m + 2 : ℕ) : ℚ) * (centralBinom (m + 2) : ℚ)
      = 2 * (2 * ((m + 1 : ℕ) : ℚ) + 1) * (centralBinom (m + 1) : ℚ) := by exact_mod_cast hc
  push_cast at hc'
  have f1 : ((m + 1)! : ℚ) = (m + 1) * (m ! : ℚ) := by push_cast [Nat.factorial_succ]; ring
  have f2 : ((m + 2)! : ℚ) = (m + 2) * (m + 1) * (m ! : ℚ) := by
    push_cast [Nat.factorial_succ]; ring
  have f3 : ((m + 1 + 2)! : ℚ) = (m + 3) * (m + 2) * (m + 1) * (m ! : ℚ) := by
    push_cast [Nat.factorial_succ]; ring
  have hf : (m ! : ℚ) ≠ 0 := by positivity
  have hC : (centralBinom (m + 2) : ℚ) =
      2 * (2 * (m : ℚ) + 3) * (centralBinom (m + 1) : ℚ) / (m + 2) := by
    rw [eq_div_iff (by positivity)]; linarith
  simp only [Dcoef, a_closed, b_closed, e]
  rw [f3, f2, f1, hC, show m + 1 + 1 = m + 2 by ring]
  field_simp
  ring

theorem Dcoef_nonneg (n : ℕ) : 0 ≤ Dcoef n := by
  match n with
  | 0 => simp [Dcoef]
  | 1 => rw [Dcoef_one]
  | m + 2 =>
    have hK : (0 : ℚ) < 4 ^ (m + 1) * ((m + 1)! : ℚ) * ((m + 2)! : ℚ) * ((m + 2) * (m + 3)) := by
      positivity
    have h := coeff_ineq (m + 1)
    have h' : (0 : ℚ) ≤ (centralBinom (m + 1) : ℚ) * (4 * (2 * m + 3) + m * (m + 2) * (m + 3))
        - 4 ^ (m + 1) * ((m + 2) * (m + 3)) := by
      have : (0 : ℤ) ≤ (centralBinom (m + 1) : ℤ) *
          (4 * (2 * ((m + 1 : ℕ) : ℤ) + 1) + (((m + 1 : ℕ) : ℤ) - 1) * ((m + 1 : ℕ) + 1) * ((m + 1 : ℕ) + 2))
          - 4 ^ (m + 1) * ((((m + 1 : ℕ) : ℤ) + 1) * (((m + 1 : ℕ) : ℤ) + 2)) := by linarith
      have hq : (0 : ℚ) ≤ ((_ : ℤ) : ℚ) := Int.cast_nonneg this
      push_cast at hq
      convert hq using 2 <;> ring
    rw [← Dcoef_closed m] at h'
    exact nonneg_of_mul_nonneg_left h' hK

end YangMills.StrongCoupling
