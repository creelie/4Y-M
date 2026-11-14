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
and `0.1505`, the numbers `3.443`, `5.947` and `6.01` of the
centre-blind example, the constants of the finite gauge groups
(Proposition `prop:finiteG`) and of the Maxwell model
(Proposition `prop:maxwell`).
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

/-! ### Finite gauge groups (Proposition `prop:finiteG`) -/

/-- The window bound at `T₀ = 3`: `(6·48³·3⁴/400)/(1 - 1/100) · 10⁻³⁶ ≤ 1.36·10⁻³¹`. -/
theorem finite_window_S :
    (6 * 48 ^ 3 * 3 ^ 4 / 400 : ℝ) / (1 - 1 / 100) * (1 / 10 ^ 36) ≤ 136 / 10 ^ 33 := by
  norm_num

/-- The window bound decreases: `(T+1)⁴ · 10^{4T²} ≤ T⁴ · 10^{4(T+1)²}` for `T ≥ 1`. -/
theorem finite_window_decreasing (T : ℕ) (hT : 1 ≤ T) :
    (T + 1) ^ 4 * 10 ^ (4 * T ^ 2) ≤ T ^ 4 * 10 ^ (4 * (T + 1) ^ 2) := by
  have hsplit : 4 * (T + 1) ^ 2 = 4 * T ^ 2 + (8 * T + 4) := by ring
  rw [hsplit, pow_add]
  have h1 : (T + 1) ^ 4 ≤ 16 * T ^ 4 := by
    have : T + 1 ≤ 2 * T := by omega
    calc (T + 1) ^ 4 ≤ (2 * T) ^ 4 := Nat.pow_le_pow_left this 4
      _ = 16 * T ^ 4 := by ring
  have h2 : 16 ≤ 10 ^ (8 * T + 4) := by
    calc 16 ≤ 10 ^ 4 := by norm_num
      _ ≤ 10 ^ (8 * T + 4) := Nat.pow_le_pow_right (by norm_num) (by omega)
  calc (T + 1) ^ 4 * 10 ^ (4 * T ^ 2) ≤ 16 * T ^ 4 * 10 ^ (4 * T ^ 2) :=
        Nat.mul_le_mul_right _ h1
    _ ≤ 10 ^ (8 * T + 4) * T ^ 4 * 10 ^ (4 * T ^ 2) := by
        apply Nat.mul_le_mul_right; exact Nat.mul_le_mul_right _ h2
    _ = T ^ 4 * (10 ^ (4 * T ^ 2) * 10 ^ (8 * T + 4)) := by ring

/-- `e^S - 1 ≤ S/(1-S)` for `0 ≤ S < 1`. -/
theorem exp_sub_one_le (S : ℝ) (_h0 : 0 ≤ S) (h1 : S < 1) :
    Real.exp S - 1 ≤ S / (1 - S) := by
  have hpos : 0 < 1 - S := by linarith
  have h : 1 - S ≤ Real.exp (-S) := by linarith [Real.add_one_le_exp (-S)]
  have he : Real.exp S * Real.exp (-S) = 1 := by rw [← Real.exp_add]; simp
  have hS : Real.exp S * (1 - S) ≤ 1 := by
    calc Real.exp S * (1 - S) ≤ Real.exp S * Real.exp (-S) :=
          mul_le_mul_of_nonneg_left h (Real.exp_pos S).le
      _ = 1 := he
  rw [le_div_iff₀ hpos]; nlinarith

/-- The twisted ratio on the window: `e^S - 1 ≤ 1.4·10⁻³¹ < 1 - 6.3·10⁻³`. -/
theorem finite_twist_small (S : ℝ) (h0 : 0 ≤ S) (h : S ≤ 136 / 10 ^ 33) :
    Real.exp S - 1 ≤ 14 / 10 ^ 32 ∧ (14 / 10 ^ 32 : ℝ) < 1 - 63 / 10 ^ 4 := by
  refine ⟨?_, by norm_num⟩
  have h1 : S < 1 := by linarith [show (136 / 10 ^ 33 : ℝ) < 1 by norm_num]
  have hpos : 0 < 1 - S := by linarith
  calc Real.exp S - 1 ≤ S / (1 - S) := exp_sub_one_le S h0 h1
    _ ≤ 14 / 10 ^ 32 := by
        rw [div_le_iff₀ hpos]
        nlinarith

/-- The twist costs at least `70`: `e^{70} ≤ 1/(1.4·10⁻³¹)`. -/
theorem finite_twist_cost : Real.exp 70 ≤ 1 / (14 / 10 ^ 32) := by
  have h : Real.exp 70 = Real.exp 1 ^ 70 := by
    rw [← Real.exp_nat_mul]; norm_num
  rw [h]
  have h1 : Real.exp 1 ^ 70 ≤ (2.7182818286 : ℝ) ^ 70 :=
    pow_le_pow_left₀ (Real.exp_pos 1).le Real.exp_one_lt_d9.le 70
  have h2 : (2.7182818286 : ℝ) ^ 70 ≤ 1 / (14 / 10 ^ 32) := by norm_num
  linarith

/-- The cube of side `n = 3` at `T = n⁴`: `6·3⁷/396 · 10⁻¹⁸ ≤ 3.4·10⁻¹⁷`. -/
theorem finite_cube_S : (6 * 3 ^ 7 / 396 : ℝ) * (1 / 10 ^ 18) ≤ 34 / 10 ^ 18 := by
  norm_num

/-- `ε ≤ 3.4·10⁻¹⁷` and `(1-D)/(1+D) ≤ ε` give `D ≥ 1/2`. -/
theorem finite_cube_defect (D ε : ℝ) (hD : 0 ≤ D) (he : ε ≤ 34 / 10 ^ 18)
    (h : (1 - D) / (1 + D) ≤ ε) : 1 / 2 ≤ D := by
  have hp : 0 < 1 + D := by linarith
  rw [div_le_iff₀ hp] at h
  nlinarith

/-- `|G|/(4|G|-3) ≤ 2/5` for `|G| ≥ 2`. -/
theorem finite_hom_ratio (n : ℝ) (hn : 2 ≤ n) : n / (4 * n - 3) ≤ 2 / 5 := by
  have hp : 0 < 4 * n - 3 := by linarith
  rw [div_le_iff₀ hp]; linarith

/-- `2/5 + 1/10 = 1/2 < e^{-1/160}`. -/
theorem finite_purity : (2 / 5 + 1 / 10 : ℝ) < Real.exp (-(1 / 160)) := by
  have h := Real.add_one_le_exp (-(1 / 160 : ℝ))
  norm_num at h ⊢; linarith

/-! ### The Maxwell model (Proposition `prop:maxwell`) -/

/-- `3ζ(2)/2.27797² + 3ζ(3)/2.27797 + ζ(4) ≤ 3.6164` with decimal bounds of `ζ`. -/
theorem maxwell_zeta (z2 z3 z4 : ℝ) (h2 : z2 ≤ 1.64494) (h3 : z3 ≤ 1.20206)
    (h4 : z4 ≤ 1.08233) : 3 * z2 / 2.27797 ^ 2 + 3 * z3 / 2.27797 + z4 ≤ 3.6164 := by
  have e2 : 3 * z2 / 2.27797 ^ 2 ≤ 3 * 1.64494 / 2.27797 ^ 2 := by
    apply div_le_div_of_nonneg_right _ (by norm_num); linarith
  have e3 : 3 * z3 / 2.27797 ≤ 3 * 1.20206 / 2.27797 := by
    apply div_le_div_of_nonneg_right _ (by norm_num); linarith
  have : 3 * 1.64494 / (2.27797 : ℝ) ^ 2 + 3 * 1.20206 / 2.27797 + 1.08233 ≤ 3.6164 := by
    norm_num
  linarith

/-- The photon constant: `2 · 3.6164 · c³ ≤ 10.69` for `c³ ≤ 1.47761`. -/
theorem maxwell_photon (c3 : ℝ) (_h0 : 0 ≤ c3) (h : c3 ≤ 1.47761) :
    2 * (3.6164 * c3) ≤ 10.69 := by nlinarith

/-- `C(1) ≤ 13.37` from `log(1+√(2304π)) ≤ 4.4553` and `log(1+√(576/π)) ≤ 2.677`. -/
theorem maxwell_count_one (L1 L2 : ℝ) (h1 : L1 ≤ 4.4553) (h2 : L2 ≤ 2.677) :
    10.69 + 3 / 8 * (L1 + L2) ≤ 13.37 := by linarith

/-! ### The twist cost for SU(N) (Proposition `prop:twistcost`) -/

/-- `(1 - x) ≤ 2 W^{1/4}` in the form used: `(1-x)⁴/16 ≤ W` and `W ≤ D` give `(1-x)⁴ ≤ 16 D`. -/
theorem twist_fourth_power (x W D : ℝ) (h1 : (1 - x) ^ 4 / 16 ≤ W) (h2 : W ≤ D) :
    (1 - x) ^ 4 ≤ 16 * D := by linarith

/-- `(0.04967/2)⁴ ≥ 3.804·10⁻⁷`, so `(1-x)⁴ ≤ 16·3.804·10⁻⁷` forces `1 - x ≤ 0.04967`. -/
theorem twist_root_general (u : ℝ) (_hu : 0 ≤ u) (h : u ^ 4 ≤ 16 * (3804 / 10 ^ 10)) :
    u ≤ 0.04967 := by
  by_contra hc
  push Not at hc
  have : (0.04967 : ℝ) ^ 4 < u ^ 4 := by
    exact pow_lt_pow_left₀ hc (by norm_num) (by norm_num)
  have h2 : (16 : ℝ) * (3804 / 10 ^ 10) ≤ 0.04967 ^ 4 := by norm_num
  linarith

/-- The even case: `(1-x)⁴ ≤ 16·2.261·10⁻⁷` forces `1 - x ≤ 0.04362`. -/
theorem twist_root_even (u : ℝ) (_hu : 0 ≤ u) (h : u ^ 4 ≤ 16 * (2261 / 10 ^ 10)) :
    u ≤ 0.04362 := by
  by_contra hc
  push Not at hc
  have : (0.04362 : ℝ) ^ 4 < u ^ 4 := by
    exact pow_lt_pow_left₀ hc (by norm_num) (by norm_num)
  have h2 : (16 : ℝ) * (2261 / 10 ^ 10) ≤ 0.04362 ^ 4 := by norm_num
  linarith

/-- The constants of the base scale: `(1 - 0.04967)/4 ≥ 0.2375` and `4(1 - 0.04362) ≥ 3.825`. -/
theorem twist_base_constants :
    (0.2375 : ℝ) ≤ (1 - 0.04967) / 4 ∧ (3.825 : ℝ) ≤ 4 * (1 - 0.04362) := by
  constructor <;> norm_num

/-! ### Finite gauge groups: the exponential base scale (Proposition `prop:finiteG`(b)) -/

/-- The link ordering leaves `|Λ|+3` free links: `(|Λ'|+2) + (n-1)|Λ'| + 1 = n|Λ'| + 3`. -/
theorem finiteG_free_links (V' n : ℕ) (hn : 1 ≤ n) :
    (V' + 2) + (n - 1) * V' + 1 = n * V' + 3 := by
  obtain ⟨k, rfl⟩ : ∃ k, n = k + 1 := ⟨n - 1, by omega⟩
  simp only [Nat.add_sub_cancel]
  ring

/-- `Tr ρ² ≤ a(1 + |G|³ a (e^x - 1))` with `a ≤ 1/2`: the right side increases in `a`. -/
theorem finiteG_trace_mono (a G3 E : ℝ) (ha : 0 ≤ a) (ha2 : a ≤ 1 / 2) (hG : 0 ≤ G3) (hE : 0 ≤ E) :
    a * (1 + G3 * a * E) ≤ 1 / 2 * (1 + G3 * (1 / 2) * E) := by
  have h1 : G3 * a * E ≤ G3 * (1 / 2) * E := by
    apply mul_le_mul_of_nonneg_right _ hE
    exact mul_le_mul_of_nonneg_left ha2 hG
  have h2 : 0 ≤ 1 + G3 * a * E := by positivity
  nlinarith

/-- With `t ≥ 1.97507/|G|³` and `|G|³ ≥ 8`: `t/(1+t) ≥ 1.584/|G|³`, written for `c = t |G|³`. -/
theorem finiteG_t_ratio (c G3 : ℝ) (hc : 1.97507 ≤ c) (hG : 8 ≤ G3) :
    1.584 * (G3 + c) ≤ c * G3 := by nlinarith

/-- `96/1.584 < 60.61` and `60.61 < 2.791⁴`. -/
theorem finiteG_exp_constants : (96 : ℝ) / 1.584 < 60.61 ∧ (60.61 : ℝ) < 2.791 ^ 4 := by
  constructor <;> norm_num

/-! ### Two dimensions, SU(2) (Proposition `prop:twod`) -/

/-- `Σ_{n≥1} qⁿ = q/(1-q) ≤ d` iff `q ≤ d/(1+d)`. -/
theorem twod_geometric (q d : ℝ) (hq1 : q < 1)
    (h : q * (1 + d) ≤ d) : q / (1 - q) ≤ d := by
  rw [div_le_iff₀ (by linarith)]
  linarith

/-- The base scale in string units: `1.8528² ≤ 3.4332` and `4.5784 ≤ 2.1398²`. -/
theorem twod_string_units : (1.8528 : ℝ) ^ 2 ≤ 3.4332 ∧ (4.5784 : ℝ) ≤ 2.1398 ^ 2 := by
  constructor <;> norm_num

end YangMills.CountSectors
