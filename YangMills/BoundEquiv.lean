import Mathlib

/-!
# What the one bound is

Finite and elementary parts of Section `ssec:whatB` of the paper
(Theorem `thm:Bequiv`, Propositions `prop:defect` and `prop:twist`): the
normalisation `γ = 1/(320 (96T)³)` of (B*) at `T₀ = 2T`, the constant
`2·96³ = 1 769 472`, the monotonicity of `(e^{αv} − 1)/v`, the constant
`4.6·10⁻⁷` of the Stefan–Boltzmann bound, the defect thresholds `d_*`,
`e^{1/320} − 1`, and the bounds `6.3·10⁻³` and `4.6·10⁻⁷` on the centre twist.
-/

namespace YangMills.BoundEquiv

/-- `2·96³ = 1 769 472`. -/
theorem two_mul_cube : (2 : ℝ) * 96 ^ 3 = 1769472 := by norm_num

/-- At `T₀ = 2T` the constant `γ = 1/(320 (96T)³)` has `2γ(96 T₀)³ = 1/20`. -/
theorem gamma_normalisation (T : ℝ) (hT : 0 < T) :
    2 * (1 / (320 * (96 * T) ^ 3)) * (96 * (2 * T)) ^ 3 = 1 / 20 := by
  field_simp
  ring

/-- `e^{W - mT} ≤ 1/320` whenever `mT ≥ W + log 320`. -/
theorem exp_le_of_margin (W mT : ℝ) (h : W + Real.log 320 ≤ mT) :
    Real.exp (W - mT) ≤ 1 / 320 := by
  have h320 : Real.exp (-Real.log 320) = 1 / 320 := by
    rw [Real.exp_neg, Real.exp_log (by norm_num)]
    norm_num
  rw [← h320]
  exact Real.exp_le_exp.mpr (by linarith)

/-- `(e^{αv} − 1)/v` is nondecreasing on `v > 0` (for every real `α`). -/
theorem expm1_div_mono (α v w : ℝ) (hv : 0 < v) (hvw : v ≤ w) :
    (Real.exp (α * v) - 1) / v ≤ (Real.exp (α * w) - 1) / w := by
  have hw : 0 < w := lt_of_lt_of_le hv hvw
  rw [div_le_div_iff₀ hv hw]
  have h1 : 0 ≤ 1 - v / w := by
    rw [sub_nonneg, div_le_one hw]; exact hvw
  have h2 : 0 ≤ v / w := div_nonneg hv.le hw.le
  have hc := convexOn_exp.2 (Set.mem_univ (0 : ℝ)) (Set.mem_univ (α * w)) h1 h2 (by ring)
  have harg : (1 - v / w) • (0 : ℝ) + (v / w) • (α * w) = α * v := by
    simp only [smul_eq_mul, mul_zero, zero_add]; field_simp
  rw [harg] at hc
  have key : Real.exp (α * v) ≤ (1 - v / w) * Real.exp (α * 0) + (v / w) * Real.exp (α * w) := by
    simpa [smul_eq_mul] using hc
  simp only [mul_zero, Real.exp_zero, mul_one] at key
  have : (Real.exp (α * v) - 1) * w ≤ (Real.exp (α * w) - 1) * v := by
    have hw' : w ≠ 0 := hw.ne'
    have e : (1 - v / w) + v / w * Real.exp (α * w) - 1 = v / w * (Real.exp (α * w) - 1) := by
      ring
    have k2 : Real.exp (α * v) - 1 ≤ v / w * (Real.exp (α * w) - 1) := by linarith
    calc (Real.exp (α * v) - 1) * w ≤ v / w * (Real.exp (α * w) - 1) * w :=
          mul_le_mul_of_nonneg_right k2 hw.le
      _ = (Real.exp (α * w) - 1) * v := by field_simp
  linarith

/-- The Stefan–Boltzmann constant of Theorem `thm:Bequiv`(b): `8 y₀ / 96³ ≤ 4.6·10⁻⁷`. -/
theorem sb_constant (y₀ : ℝ) (h1 : y₀ ≤ 1 / 20) :
    8 * y₀ / 96 ^ 3 ≤ 4.6e-7 := by
  have : 8 * y₀ / 96 ^ 3 ≤ 8 * (1 / 20) / 96 ^ 3 := by
    apply div_le_div_of_nonneg_right _ (by norm_num); linarith
  have h2 : 8 * (1 / 20 : ℝ) / 96 ^ 3 ≤ 4.6e-7 := by norm_num
  linarith

/-- `T_i ≤ T < 2 T_i` gives `T_i⁻³ < 8 T⁻³`. -/
theorem cube_ratio (Ti T : ℝ) (hTi : 0 < Ti) (h1 : Ti ≤ T) (h2 : T < 2 * Ti) :
    1 / Ti ^ 3 < 8 / T ^ 3 := by
  have hT : 0 < T := lt_of_lt_of_le hTi h1
  rw [div_lt_div_iff₀ (by positivity) (by positivity)]
  have : T ^ 3 < (2 * Ti) ^ 3 := by
    exact pow_lt_pow_left₀ h2 hT.le (by norm_num)
  nlinarith

/-- On the window the volume factor is `γ (48 T₀)³ = y₀/16 ≤ 1/320`. -/
theorem window_volume (T₀ y₀ : ℝ) (hT : 0 < T₀) (h1 : y₀ ≤ 1 / 20) :
    y₀ / (2 * (96 * T₀) ^ 3) * (48 * T₀) ^ 3 ≤ 1 / 320 := by
  have e : y₀ / (2 * (96 * T₀) ^ 3) * (48 * T₀) ^ 3 = y₀ / 16 := by
    field_simp
    ring
  rw [e]; linarith

/-- The smallest window volume: `γ (2 T₀)³ = 4 y₀ / 96³`. -/
theorem window_min_volume (T₀ y₀ : ℝ) (hT : 0 < T₀) :
    y₀ / (2 * (96 * T₀) ^ 3) * (2 * T₀) ^ 3 = 4 * y₀ / 96 ^ 3 := by
  field_simp
  ring

/-- `e^x − 1 ≤ x/(1 − x)` for `x < 1`. -/
theorem exp_sub_one_le_div (x : ℝ) (h1 : x < 1) :
    Real.exp x - 1 ≤ x / (1 - x) := by
  have hpos : 0 < 1 - x := by linarith
  have h := Real.add_one_le_exp (-x)
  have hexp : Real.exp x * Real.exp (-x) = 1 := by
    rw [← Real.exp_add]; simp
  have hle : Real.exp x * (1 - x) ≤ 1 := by
    calc Real.exp x * (1 - x) ≤ Real.exp x * Real.exp (-x) :=
          mul_le_mul_of_nonneg_left (by linarith) (Real.exp_pos x).le
      _ = 1 := hexp
  rw [le_div_iff₀ hpos]
  nlinarith

/-- Proposition `prop:defect`(c): `e^{1/320} − 1 ≤ 1/319 < 3.2·10⁻³`. -/
theorem window_defect : Real.exp (1 / 320) - 1 ≤ 1 / 319 ∧ (1 / 319 : ℝ) < 3.2e-3 := by
  refine ⟨?_, by norm_num⟩
  have := exp_sub_one_le_div (1 / 320) (by norm_num)
  norm_num at this ⊢
  linarith

/-- The threshold `d_* = e^{4 y₀/96³} − 1` at `y₀ = 1/20` lies in `[2.26·10⁻⁷, 2.27·10⁻⁷]`. -/
theorem dstar_bounds :
    2.26e-7 ≤ Real.exp (4 * (1 / 20) / 96 ^ 3) - 1 ∧
      Real.exp (4 * (1 / 20) / 96 ^ 3) - 1 ≤ 2.27e-7 := by
  constructor
  · have h := Real.add_one_le_exp (4 * (1 / 20 : ℝ) / 96 ^ 3)
    have h2 : (2.26e-7 : ℝ) ≤ 4 * (1 / 20) / 96 ^ 3 := by norm_num
    linarith
  · have h := exp_sub_one_le_div (4 * (1 / 20 : ℝ) / 96 ^ 3) (by norm_num)
    have h2 : (4 * (1 / 20 : ℝ) / 96 ^ 3) / (1 - 4 * (1 / 20) / 96 ^ 3) ≤ 2.27e-7 := by
      norm_num
    linarith

/-- `(1 − D)/(1 + D) ≥ 1 − 2D` for `D ≥ 0`. -/
theorem twist_ratio (D : ℝ) (hD : 0 ≤ D) : 1 - 2 * D ≤ (1 - D) / (1 + D) := by
  rw [le_div_iff₀ (by linarith)]
  nlinarith

/-- Proposition `prop:twist`: the two numerical bounds `6.3·10⁻³` and `4.6·10⁻⁷`. -/
theorem twist_constants :
    2 * (1 / 319 : ℝ) ≤ 6.3e-3 ∧ 2 * (2.27e-7 : ℝ) ≤ 4.6e-7 := by
  constructor <;> norm_num

/-- The criterion of Proposition `prop:defect`(b): a gap `e^{-m}` over the time
`T₀ − T'` turns `D(T') ≤ B` into `D(T₀) ≤ d_*` once `m (T₀ − T') ≥ log(B/d_*)`. -/
theorem defect_criterion (m gap B dstar : ℝ) (hB : 0 < B) (hd : 0 < dstar)
    (h : Real.log (B / dstar) ≤ m * gap) : Real.exp (-(m * gap)) * B ≤ dstar := by
  have h1 : Real.exp (-(m * gap)) ≤ Real.exp (-Real.log (B / dstar)) :=
    Real.exp_le_exp.mpr (by linarith)
  have h2 : Real.exp (-Real.log (B / dstar)) = dstar / B := by
    rw [Real.exp_neg, Real.exp_log (div_pos hB hd), inv_div]
  rw [h2] at h1
  calc Real.exp (-(m * gap)) * B ≤ dstar / B * B := mul_le_mul_of_nonneg_right h1 hB.le
    _ = dstar := by field_simp

end YangMills.BoundEquiv
