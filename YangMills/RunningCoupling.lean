import Mathlib

/-!
# One-loop running coupling: the intermediate-coupling window is always entered

With `b₀ = 11 N / (48 π²)`, the one-loop solution of `μ dg/dμ = -b₀ g³` is
`g²(μ) = 1 / (2 b₀ log (μ / Λ))` for `μ > Λ`. It tends to `0` as `μ → ∞`
(asymptotic freedom) and to `+∞` as `μ → Λ⁺`, and it takes every value
`g_*² > 0` at the scale `μ_* = Λ exp (1 / (2 b₀ g_*²))`, which does not depend
on the lattice spacing. Hence along any continuum trajectory with
`β(a) = 2N / g²(1/a)`, the effective coupling at the fixed physical scale
`μ_*` equals `g_*²` for every spacing `a < 1/μ_*`: the bare coupling leaves
the intermediate window, but the effective coupling on physical scales does
not.
-/

open Real Filter Topology

namespace YangMills

/-- The one-loop coefficient `b₀ = 11 N / (48 π²)`. -/
noncomputable def b0 (N : ℕ) : ℝ := 11 * N / (48 * π ^ 2)

theorem b0_pos {N : ℕ} (hN : 1 ≤ N) : 0 < b0 N := by
  unfold b0
  have : (1 : ℝ) ≤ N := by exact_mod_cast hN
  positivity

/-- One-loop running coupling `g²(μ)`. -/
noncomputable def gsq (N : ℕ) (Λ μ : ℝ) : ℝ := 1 / (2 * b0 N * Real.log (μ / Λ))

/-- The scale at which the coupling equals `g`. -/
noncomputable def scaleOf (N : ℕ) (Λ g : ℝ) : ℝ := Λ * Real.exp (1 / (2 * b0 N * g))

theorem scaleOf_gt {N : ℕ} (hN : 1 ≤ N) {Λ g : ℝ} (hΛ : 0 < Λ) (hg : 0 < g) :
    Λ < scaleOf N Λ g := by
  unfold scaleOf
  have : 0 < 1 / (2 * b0 N * g) := by have := b0_pos hN; positivity
  have h1 : 1 < Real.exp (1 / (2 * b0 N * g)) := Real.one_lt_exp_iff.mpr this
  nlinarith

/-- **Every coupling is attained.** -/
theorem gsq_scaleOf {N : ℕ} (hN : 1 ≤ N) {Λ g : ℝ} (hΛ : 0 < Λ) (hg : 0 < g) :
    gsq N Λ (scaleOf N Λ g) = g := by
  have hb := b0_pos hN
  unfold gsq scaleOf
  rw [mul_div_cancel_left₀ _ hΛ.ne', Real.log_exp]
  field_simp

/-- **Asymptotic freedom.** -/
theorem gsq_tendsto_zero {N : ℕ} (hN : 1 ≤ N) {Λ : ℝ} (hΛ : 0 < Λ) :
    Tendsto (gsq N Λ) atTop (𝓝 0) := by
  have hb := b0_pos hN
  have h1 : Tendsto (fun μ : ℝ => μ / Λ) atTop atTop :=
    tendsto_id.atTop_div_const hΛ
  have h2 : Tendsto (fun μ : ℝ => 2 * b0 N * Real.log (μ / Λ)) atTop atTop :=
    (Real.tendsto_log_atTop.comp h1).const_mul_atTop (by positivity)
  have h3 := h2.inv_tendsto_atTop
  refine h3.congr fun μ => ?_
  simp [gsq]

/-- The coupling at the spacing-independent scale `μ_*` equals `g` for every
lattice spacing `a < 1/μ_*`, and that scale lies strictly between `Λ` and the
cutoff `1/a`. -/
theorem window_entered {N : ℕ} (hN : 1 ≤ N) {Λ g : ℝ} (hΛ : 0 < Λ) (hg : 0 < g)
    {a : ℝ} (ha : 0 < a) (hsmall : a < 1 / scaleOf N Λ g) :
    Λ < scaleOf N Λ g ∧ scaleOf N Λ g < 1 / a ∧ gsq N Λ (scaleOf N Λ g) = g := by
  have hs := scaleOf_gt hN hΛ hg
  have hspos : 0 < scaleOf N Λ g := hΛ.trans hs
  refine ⟨hs, ?_, gsq_scaleOf hN hΛ hg⟩
  rw [lt_div_iff₀ ha]
  rw [lt_div_iff₀ hspos] at hsmall
  linarith

end YangMills
