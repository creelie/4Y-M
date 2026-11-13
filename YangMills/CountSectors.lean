import Mathlib

/-!
# The count of states and the centre sectors

Finite and elementary parts of Section `ssec:count` of the paper
(Proposition `prop:count`, Theorem `thm:plane`, Lemma `lem:transport`,
Proposition `prop:countscale`, Lemma `lem:flux`, Theorem `thm:sectors`,
Proposition `prop:centreblind`): the bound `T/(T-λ) ≤ 2`, the criterion
`2P e^{-mT/2} ≤ d`, the trace inequality of the plane model, the minimum of
`x³/(x-a)`, the constant `926`, the charged-defect identity, the sector
arithmetic `1/642 · 641/640 · 642/641 = 1/640`, the constants `0.02492`
and `0.1505`, and the numbers `3.443`, `5.947` and `6.01` of the
centre-blind example.
-/

namespace YangMills.CountSectors

/-- `T/(T-λ) ≤ 2` when `0 ≤ λ ≤ T/2`. -/
theorem ratio_le_two (T l : ℝ) (hT : 0 < T) (hl : l ≤ T / 2) :
    T / (T - l) ≤ 2 := by
  have h : 0 < T - l := by linarith
  rw [div_le_iff₀ h]; linarith

/-- `2P e^{-mT/2} ≤ d` once `mT ≥ 2 log(2P/d)`. -/
theorem count_criterion (P d mT : ℝ) (hP : 0 < P) (hd : 0 < d)
    (h : 2 * Real.log (2 * P / d) ≤ mT) : 2 * P * Real.exp (-(mT / 2)) ≤ d := by
  have hq : 0 < 2 * P / d := by positivity
  have h1 : Real.exp (-(mT / 2)) ≤ Real.exp (-Real.log (2 * P / d)) :=
    Real.exp_le_exp.mpr (by linarith)
  have h2 : Real.exp (-Real.log (2 * P / d)) = d / (2 * P) := by
    rw [Real.exp_neg, Real.exp_log hq, inv_div]
  rw [h2] at h1
  calc 2 * P * Real.exp (-(mT / 2)) ≤ 2 * P * (d / (2 * P)) :=
        mul_le_mul_of_nonneg_left h1 (by positivity)
    _ = d := by field_simp

/-- The plane model: `q^ℓ ≤ q^T + q^ℓ - 1` for `q ≥ 1`. -/
theorem plane_trace (q : ℝ) (hq : 1 ≤ q) (T ℓ : ℕ) :
    q ^ ℓ ≤ q ^ T + q ^ ℓ - 1 := by
  have := one_le_pow₀ (n := T) hq
  linarith

/-- The plane model: `ℓ ≥ 2T` gives `2(ℓ - T) ≥ ℓ`. -/
theorem plane_half (T ℓ : ℕ) (h : 2 * T ≤ ℓ) : ℓ ≤ 2 * (ℓ - T) := by
  omega

/-- The plane model: three factors of `(V/ℓ)(ℓ/2)` give `3V/2`. -/
theorem plane_sum (V l1 l2 l3 : ℝ) (h1 : 0 < l1) (h2 : 0 < l2) (h3 : 0 < l3) :
    V / l1 * (l1 / 2) + V / l2 * (l2 / 2) + V / l3 * (l3 / 2) = 3 / 2 * V := by
  field_simp; ring

/-- `x³ - (27/4) a² (x - a) = (x - 3a/2)² (x + 3a)`. -/
theorem cube_identity (x a : ℝ) :
    x ^ 3 - 27 / 4 * a ^ 2 * (x - a) = (x - 3 / 2 * a) ^ 2 * (x + 3 * a) := by
  ring

/-- On `x > a > 0`, `x³/(x - a) ≥ (27/4) a²`. -/
theorem cube_div_min (x a : ℝ) (ha : 0 < a) (hx : a < x) :
    27 / 4 * a ^ 2 ≤ x ^ 3 / (x - a) := by
  have h : 0 < x - a := by linarith
  rw [le_div_iff₀ h]
  have := cube_identity x a
  have : 0 ≤ (x - 3 / 2 * a) ^ 2 * (x + 3 * a) := by
    apply mul_nonneg (sq_nonneg _); linarith
  linarith

/-- `4·96³ = 3 538 944` and `27·96³ = 23 887 872`. -/
theorem scale_constants : (4 : ℝ) * 96 ^ 3 = 3538944 ∧ (27 : ℝ) * 96 ^ 3 = 23887872 := by
  constructor <;> norm_num

/-- With `a ≥ 5.768` (`a = log 320`): `926³ ≤ 27·96³·a²`. -/
theorem scale_926 (a : ℝ) (ha : 5.768 ≤ a) : (926 : ℝ) ^ 3 ≤ 27 * 96 ^ 3 * a ^ 2 := by
  have h : (5.768 : ℝ) ^ 2 ≤ a ^ 2 := by nlinarith
  nlinarith

/-- The charged defect: `D = W(1 + D₀ + D)` with `W < 1` gives `D = W(1+D₀)/(1-W)`. -/
theorem charged_defect (W D0 D : ℝ) (hW : W < 1) (h : D = W * (1 + D0 + D)) :
    D = W * (1 + D0) / (1 - W) := by
  have : 1 - W ≠ 0 := by linarith
  field_simp
  linarith

/-- Sector arithmetic of Theorem `thm:sectors`(a). -/
theorem sector_product : (1 / 642 : ℝ) * (641 / 640) * (642 / 641) = 1 / 640 := by
  norm_num

/-- `v/640 + v/640 = v/320`. -/
theorem sector_sum (v : ℝ) : v / 640 + v / 640 = v / 320 := by ring

/-- The constant `κ₂ = 16/642 ≤ 0.02493`. -/
theorem kappa_two : (16 / 642 : ℝ) ≤ 0.02493 ∧ 0.02492 ≤ (16 / 642 : ℝ) := by
  constructor <;> norm_num

/-- `3 e^{1/320}/20 ≤ 0.1505`, from `e^{1/320} ≤ 320/319`. -/
theorem twist_extensive (e : ℝ) (he : e ≤ 320 / 319) : 3 * e / 20 ≤ 0.1505 := by
  have : 3 * e / 20 ≤ 3 * (320 / 319) / 20 := by linarith
  have h2 : 3 * (320 / 319 : ℝ) / 20 ≤ 0.1505 := by norm_num
  linarith

/-- `e^{x} ≤ 1/(1-x)` gives `e^{1/320} ≤ 320/319`. -/
theorem exp_small : Real.exp (1 / 320) ≤ 320 / 319 := by
  have h := Real.add_one_le_exp (-(1 / 320 : ℝ))
  have hexp : Real.exp (1 / 320) * Real.exp (-(1 / 320)) = 1 := by
    rw [← Real.exp_add]; simp
  have hpos := Real.exp_pos (1 / 320 : ℝ)
  nlinarith

/-- `(1+a/p)(1+b/p)(1+c/p) - 1` expanded. -/
theorem mode_expand (a b c p : ℝ) (hp : p ≠ 0) :
    (1 + a / p) * (1 + b / p) * (1 + c / p) - 1 =
      (a + b + c) / p + (a * b + a * c + b * c) / p ^ 2 + a * b * c / p ^ 3 := by
  field_simp; ring

/-- With every `ρᵢ ≥ ρ₀ > 0`: `ρ₁ ≤ ρ₁ρ₂ρ₃/ρ₀²`. -/
theorem rho_one (a b c r : ℝ) (hr : 0 < r) (ha : 0 ≤ a) (hb : r ≤ b) (hc : r ≤ c) :
    a ≤ a * b * c / r ^ 2 := by
  rw [le_div_iff₀ (by positivity)]
  have : r ^ 2 ≤ b * c := by nlinarith
  nlinarith

/-- With every `ρᵢ ≥ ρ₀ > 0`: `ρ₁ρ₂ ≤ ρ₁ρ₂ρ₃/ρ₀`. -/
theorem rho_two (a b c r : ℝ) (hr : 0 < r) (ha : 0 ≤ a) (hb : 0 ≤ b) (hc : r ≤ c) :
    a * b ≤ a * b * c / r := by
  rw [le_div_iff₀ hr]
  have : 0 ≤ a * b := mul_nonneg ha hb
  nlinarith

/-- `3ζ(2)/2.399² + 3ζ(3)/2.399 + ζ(4) ≤ 3.443` with the decimal bounds of `ζ`. -/
theorem zeta_combination (z2 z3 z4 : ℝ) (h2 : z2 ≤ 1.64494) (h3 : z3 ≤ 1.20206)
    (h4 : z4 ≤ 1.08233) : 3 * z2 / 2.399 ^ 2 + 3 * z3 / 2.399 + z4 ≤ 3.443 := by
  have e2 : 3 * z2 / 2.399 ^ 2 ≤ 3 * 1.64494 / 2.399 ^ 2 := by
    apply div_le_div_of_nonneg_right _ (by norm_num); linarith
  have e3 : 3 * z3 / 2.399 ≤ 3 * 1.20206 / 2.399 := by
    apply div_le_div_of_nonneg_right _ (by norm_num); linarith
  have : 3 * 1.64494 / (2.399 : ℝ) ^ 2 + 3 * 1.20206 / 2.399 + 1.08233 ≤ 3.443 := by
    norm_num
  linarith

/-- `3.443 · c³ ≤ 5.947` for `c ≤ 1.19979`, and `5.947 + 0.0574 ≤ 6.01`. -/
theorem blind_total (c : ℝ) (hc0 : 0 ≤ c) (hc : c ≤ 1.19979) :
    3.443 * c ^ 3 ≤ 5.947 ∧ (5.947 : ℝ) + 0.0574 ≤ 6.01 := by
  constructor
  · have : c ^ 3 ≤ 1.19979 ^ 3 := pow_le_pow_left₀ hc0 hc 3
    have h2 : (3.443 : ℝ) * 1.19979 ^ 3 ≤ 5.947 := by norm_num
    nlinarith
  · norm_num

/-- The zero mode: `0.4588/8 ≤ 0.0574`. -/
theorem zero_mode : (0.4588 : ℝ) / 8 ≤ 0.0574 := by norm_num

end YangMills.CountSectors
