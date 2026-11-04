import Mathlib

/-!
# Finite arithmetic of the layer theorem and the block-map growth

These are the finite inequalities behind `thm:layerbound`, `cor:twoislands`
and `prop:massgrowth` of the paper.

* `schur_row`, `gram_norm_inv` : a row of `𝒢 = BB*` has norm at most
  `1/4 + 681κ`, and `1/(1/4 + 681κ) = 4/(1 + 2724κ)`.
* `inf_add_ge` : the infimum of a sum over a constraint set is at least the
  sum of the infima.
* `calm_background`, `chart_inclusion` : `η + 2·14r + 2r + 32r = η + 62r`, and
  the path links move by at most `31.5r < 32r` at `ρ₁ = 5r/4`.
* `hessian_lower`, `hessian_weighted`, `rate_admissible` : the lower bound
  `0.0185 λ ≥ λ/60`, the weighted norm `2068.4 λ`, and `α_G = 2.5·10⁻⁴`
  within the admissible rate.
* `covariance_constant`, `cross_ratio_constant` : `𝒜 = 2/m = 120/λ =
  480/β_W` and `4𝒜 = 1920/β_W`.
* `radius_threshold`, `coupling_threshold` : `2758 ρ ≤ 6.0·10⁻¹³` at
  `ρ = 2.17·10⁻¹⁶`, and `4·240·c₆ ≤ 8.93·10¹⁷` for `c₆ ≤ 9.3·10¹⁴`.
-/

namespace YangMills.StrongCoupling.IslandLayer

/-- A row of `𝒢`: the diagonal block and two further blocks. -/
theorem schur_row (κ : ℝ) :
    (3 / 16 + 227 * κ) + 2 * (1 / 32 + 227 * κ) = 1 / 4 + 681 * κ := by ring

theorem gram_norm_inv (κ : ℝ) (hκ : 0 ≤ κ) :
    1 / (1 / 4 + 681 * κ) = 4 / (1 + 2724 * κ) := by
  have h1 : (0 : ℝ) < 1 / 4 + 681 * κ := by positivity
  have h2 : (0 : ℝ) < 1 + 2724 * κ := by positivity
  field_simp
  ring

/-- The gain of the infimum in `prop:massgrowth`: `ε|w|²/(1/4 + 681κ)`
equals `4ε|w|²/(1 + 2724κ)`. -/
theorem gain_eq (κ ε q : ℝ) (hκ : 0 ≤ κ) :
    ε * (q / (1 / 4 + 681 * κ)) = 4 * ε / (1 + 2724 * κ) * q := by
  have h1 : (0 : ℝ) < 1 / 4 + 681 * κ := by positivity
  have h2 : (0 : ℝ) < 1 + 2724 * κ := by positivity
  field_simp
  ring

/-- The infimum of a sum is at least the sum of the infima. -/
theorem inf_add_ge {X : Type*} (S : Set X) (N T : X → ℝ)
    (hN : BddBelow (N '' S)) (hT : BddBelow (T '' S)) (ε : ℝ) (hε : 0 ≤ ε) :
    ∀ x ∈ S, sInf (N '' S) + ε * sInf (T '' S) ≤ N x + ε * T x := by
  intro x hx
  have h1 : sInf (N '' S) ≤ N x := csInf_le hN ⟨x, hx, rfl⟩
  have h2 : sInf (T '' S) ≤ T x := csInf_le hT ⟨x, hx, rfl⟩
  nlinarith [mul_le_mul_of_nonneg_left h2 hε]

/-- Calmness of the background: at most two path links per plaquette at
`14r`, two further links at `r`, and the chart radius `32r`. -/
theorem calm_background (η r : ℝ) : η + 2 * (14 * r) + 2 * r + 32 * r = η + 62 * r := by
  ring

theorem normalisation : (92 : ℝ) + 62 * 43 = 2758 := by norm_num

/-- `ρ₁ = 5r/4 ≤ 10⁻⁴` for `r ≤ 8·10⁻⁵`, and the path links move by
`14r + 14ρ₁ = 31.5r < 32r`. -/
theorem chart_inclusion (r : ℝ) (hr : 0 < r) (hr' : r ≤ 8e-5) :
    5 / 4 * r ≤ 1e-4 ∧ 14 * r + 14 * (5 / 4 * r) = 31.5 * r ∧
      31.5 * r < 32 * r := by
  refine ⟨by linarith, by ring, by linarith⟩

/-- The lower bound on the Hessian in the label coordinates. -/
theorem hessian_lower :
    (1 / 30 - 1 / 240 : ℝ) / (1 + 3.6e-5) ^ 2 - 4.35e-3 - 1.512 / 240 ≥ 0.0185 ∧
      (0.0185 : ℝ) ≥ 1 / 60 := by
  constructor <;> norm_num

/-- The weighted norm of the Hessian in the label coordinates. -/
theorem hessian_weighted :
    (2068 : ℝ) * (1 + 7.0e-5) ^ 2 + 9.41e-3 + 5.88 / 240 ≤ 2068.4 := by
  norm_num

/-- `α_G = 2.5·10⁻⁴` lies below the admissible rate `0.197(240·2068.4)^{-1/2}`,
written without the square root. -/
theorem rate_admissible : (2.5e-4 : ℝ) ^ 2 * (240 * 2068.4) ≤ 0.197 ^ 2 := by
  norm_num

/-- `𝒜 = 2/m` with `m = λ/60` is `120/λ`, and `480/β_W` at `λ = β_W/4`. -/
theorem covariance_constant (β : ℝ) (hβ : 0 < β) :
    2 / ((β / 4) / 60) = 480 / β ∧ 2 / ((β / 4) / 60) = 120 / (β / 4) := by
  constructor <;> field_simp <;> ring

theorem cross_ratio_constant (β : ℝ) : 4 * (480 / β) = 1920 / β := by ring

theorem radius_threshold : (2758 : ℝ) * 2.17e-16 ≤ 6.0e-13 := by norm_num

theorem coupling_threshold (c₆ : ℝ) (h : c₆ ≤ 9.3e14) :
    4 * 240 * c₆ ≤ 8.93e17 := by linarith

end YangMills.StrongCoupling.IslandLayer
