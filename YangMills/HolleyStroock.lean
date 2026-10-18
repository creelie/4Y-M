import Mathlib

/-!
# The Holley–Stroock perturbation lemma on a finite state space

Let `ν` be a probability vector on a finite set `Ω`, `V : Ω → ℝ` with
`m ≤ V ≤ M`, and `μ = ν e^{-V} / Z` the perturbed probability vector. For a
nonnegative energy density `Γ` (for instance a carré du champ `|∇f|²`, which
does not depend on the measure):

* `holleyStroock_poincare` : `Var_ν f ≤ C Σ ν Γ` implies
  `Var_μ f ≤ C e^{M - m} Σ μ Γ`;
* `holleyStroock_entropy` : `Ent_ν g ≤ C Σ ν Γ` implies
  `Ent_μ g ≤ C e^{M - m} Σ μ Γ` for `g ≥ 0`; with `g = f²` this is the
  log-Sobolev form.

The factor `e^{M - m} = e^{osc V}` is the whole content: for `V = β S` with the
Wilson action, `osc V` grows linearly with the number of plaquettes
(`YangMills.wilsonAction_osc_le`).
-/

open Real Finset

namespace YangMills

set_option linter.unusedSectionVars false

variable {Ω : Type*} [Fintype Ω]

/-- Mean of `f` under the weight vector `w`. -/
def wmean (w f : Ω → ℝ) : ℝ := ∑ ω, w ω * f ω

/-- Variance of `f` under `w`. -/
def wvar (w f : Ω → ℝ) : ℝ := ∑ ω, w ω * (f ω - wmean w f) ^ 2

/-- Entropy of `g ≥ 0` under `w`. -/
noncomputable def went (w g : Ω → ℝ) : ℝ :=
  ∑ ω, w ω * (g ω * Real.log (g ω)) - wmean w g * Real.log (wmean w g)

/-- The perturbed probability vector `ν e^{-V} / Z`. -/
noncomputable def gibbs (ν V : Ω → ℝ) (ω : Ω) : ℝ :=
  ν ω * Real.exp (-V ω) / ∑ ω', ν ω' * Real.exp (-V ω')

/-- The function whose `w`-average, minimised over `c > 0`, is the entropy. -/
noncomputable def entKernel (x c : ℝ) : ℝ := x * Real.log x - x * Real.log c - x + c

lemma entKernel_nonneg {x c : ℝ} (hx : 0 ≤ x) (hc : 0 < c) : 0 ≤ entKernel x c := by
  unfold entKernel
  rcases hx.eq_or_lt with rfl | hx
  · simp; exact hc.le
  · have h := Real.log_le_sub_one_of_pos (div_pos hc hx)
    rw [Real.log_div hc.ne' hx.ne'] at h
    have h2 : x * (Real.log c - Real.log x) ≤ x * (c / x - 1) :=
      mul_le_mul_of_nonneg_left h hx.le
    have h3 : x * (c / x - 1) = c - x := by field_simp
    nlinarith

section

variable {w : Ω → ℝ} (hw0 : ∀ ω, 0 ≤ w ω) (hw1 : ∑ ω, w ω = 1)
include hw0 hw1

theorem wvar_le (f : Ω → ℝ) (c : ℝ) : wvar w f ≤ ∑ ω, w ω * (f ω - c) ^ 2 := by
  set m := wmean w f with hm
  have key : ∑ ω, w ω * (f ω - c) ^ 2 - ∑ ω, w ω * (f ω - m) ^ 2 = (m - c) ^ 2 := by
    rw [← Finset.sum_sub_distrib]
    have h : ∀ ω, w ω * (f ω - c) ^ 2 - w ω * (f ω - m) ^ 2
        = (m - c) * (2 * (w ω * f ω)) - (m - c) * ((c + m) * w ω) := by intro; ring
    rw [Finset.sum_congr rfl (fun ω _ => h ω), Finset.sum_sub_distrib,
      ← Finset.mul_sum, ← Finset.mul_sum, ← Finset.mul_sum, ← Finset.mul_sum, hw1]
    have : ∑ ω, w ω * f ω = m := rfl
    rw [this]; ring
  have : 0 ≤ (m - c) ^ 2 := sq_nonneg _
  unfold wvar; linarith

theorem went_le {g : Ω → ℝ} (hg : ∀ ω, 0 ≤ g ω) {c : ℝ} (hc : 0 < c) :
    went w g ≤ ∑ ω, w ω * entKernel (g ω) c := by
  set m := wmean w g with hm
  have hm0 : 0 ≤ m := Finset.sum_nonneg fun ω _ => mul_nonneg (hw0 ω) (hg ω)
  have hsum : ∑ ω, w ω * entKernel (g ω) c
      = ∑ ω, w ω * (g ω * Real.log (g ω)) - m * Real.log c - m + c := by
    unfold entKernel
    have h : ∀ ω, w ω * (g ω * Real.log (g ω) - g ω * Real.log c - g ω + c)
        = w ω * (g ω * Real.log (g ω)) - Real.log c * (w ω * g ω) - w ω * g ω + c * w ω := by
      intro; ring
    rw [Finset.sum_congr rfl (fun ω _ => h ω), Finset.sum_add_distrib, Finset.sum_sub_distrib,
      Finset.sum_sub_distrib, ← Finset.mul_sum, ← Finset.mul_sum, hw1]
    have : ∑ ω, w ω * g ω = m := rfl
    rw [this]; ring
  rw [hsum]; unfold went
  rw [← hm]
  suffices h : -(m * Real.log m) ≤ -(m * Real.log c) - m + c by linarith
  rcases hm0.eq_or_lt with h0 | hpos
  · rw [← h0]; simp; exact hc.le
  · have := entKernel_nonneg hpos.le hc
    unfold entKernel at this; linarith

theorem went_eq_sum_entKernel (g : Ω → ℝ) :
    went w g = ∑ ω, w ω * entKernel (g ω) (wmean w g) := by
  unfold entKernel went
  have h : ∀ ω, w ω * (g ω * Real.log (g ω) - g ω * Real.log (wmean w g) - g ω + wmean w g)
      = w ω * (g ω * Real.log (g ω)) - Real.log (wmean w g) * (w ω * g ω) - w ω * g ω
        + wmean w g * w ω := by intro; ring
  rw [Finset.sum_congr rfl (fun ω _ => h ω), Finset.sum_add_distrib, Finset.sum_sub_distrib,
    Finset.sum_sub_distrib, ← Finset.mul_sum, ← Finset.mul_sum, hw1]
  have : ∑ ω, w ω * g ω = wmean w g := rfl
  rw [this]; ring

end

section Gibbs

variable {ν V : Ω → ℝ} (hν0 : ∀ ω, 0 ≤ ν ω) (hν1 : ∑ ω, ν ω = 1)
  {m M : ℝ} (hm : ∀ ω, m ≤ V ω) (hM : ∀ ω, V ω ≤ M)
include hν0 hν1 hM

lemma partition_pos : 0 < ∑ ω, ν ω * Real.exp (-V ω) := by
  have : ∑ ω, ν ω * Real.exp (-M) ≤ ∑ ω, ν ω * Real.exp (-V ω) :=
    Finset.sum_le_sum fun ω _ =>
      mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr (by linarith [hM ω])) (hν0 ω)
  rw [← Finset.sum_mul, hν1, one_mul] at this
  exact (Real.exp_pos _).trans_le this

lemma gibbs_nonneg (ω : Ω) : 0 ≤ gibbs ν V ω :=
  div_nonneg (mul_nonneg (hν0 ω) (Real.exp_pos _).le) (partition_pos hν0 hν1 hM).le

lemma gibbs_sum : ∑ ω, gibbs ν V ω = 1 := by
  unfold gibbs
  rw [← Finset.sum_div]
  exact div_self (partition_pos hν0 hν1 hM).ne'

include hm

/-- Comparison of `μ`-averages of a nonnegative function with `ν`-averages. -/
lemma gibbs_sum_le {h : Ω → ℝ} (hh : ∀ ω, 0 ≤ h ω) :
    ∑ ω, gibbs ν V ω * h ω
      ≤ Real.exp (-m) / (∑ ω, ν ω * Real.exp (-V ω)) * ∑ ω, ν ω * h ω := by
  have hZ := partition_pos hν0 hν1 hM
  set Z := ∑ ω, ν ω * Real.exp (-V ω) with hZdef
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum fun ω _ => ?_
  have e : gibbs ν V ω * h ω = Real.exp (-V ω) / Z * (ν ω * h ω) := by
    unfold gibbs; rw [← hZdef]; ring
  rw [e]
  have : Real.exp (-V ω) ≤ Real.exp (-m) := Real.exp_le_exp.mpr (by linarith [hm ω])
  exact mul_le_mul_of_nonneg_right (div_le_div_of_nonneg_right this hZ.le)
    (mul_nonneg (hν0 ω) (hh ω))

/-- The reverse comparison, for the energy term. -/
lemma sum_le_gibbs_sum {h : Ω → ℝ} (hh : ∀ ω, 0 ≤ h ω) :
    ∑ ω, ν ω * h ω
      ≤ Real.exp M * (∑ ω, ν ω * Real.exp (-V ω)) * ∑ ω, gibbs ν V ω * h ω := by
  have hZ := partition_pos hν0 hν1 hM
  set Z := ∑ ω, ν ω * Real.exp (-V ω) with hZdef
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum fun ω _ => ?_
  have e : Real.exp M * Z * (gibbs ν V ω * h ω)
      = Real.exp M * Real.exp (-V ω) * (ν ω * h ω) := by
    unfold gibbs; rw [← hZdef]; field_simp
  rw [e]
  have : 1 ≤ Real.exp M * Real.exp (-V ω) := by
    rw [← Real.exp_add]; exact Real.one_le_exp (by linarith [hM ω])
  have h1 : 0 ≤ ν ω * h ω := mul_nonneg (hν0 ω) (hh ω)
  nlinarith

/-- **Holley–Stroock for the Poincaré inequality.** -/
theorem holleyStroock_poincare {Γ : Ω → ℝ} (hΓ : ∀ ω, 0 ≤ Γ ω) {f : Ω → ℝ} {C : ℝ}
    (hC : 0 ≤ C) (hPI : wvar ν f ≤ C * ∑ ω, ν ω * Γ ω) :
    wvar (gibbs ν V) f ≤ C * Real.exp (M - m) * ∑ ω, gibbs ν V ω * Γ ω := by
  have hZ := partition_pos hν0 hν1 hM
  set Z := ∑ ω, ν ω * Real.exp (-V ω)
  have hK : 0 ≤ Real.exp (-m) / Z := div_nonneg (Real.exp_pos _).le hZ.le
  calc wvar (gibbs ν V) f
      ≤ ∑ ω, gibbs ν V ω * (f ω - wmean ν f) ^ 2 :=
        wvar_le (gibbs_nonneg hν0 hν1 hM) (gibbs_sum hν0 hν1 hM) f _
    _ ≤ Real.exp (-m) / Z * wvar ν f :=
        gibbs_sum_le hν0 hν1 hm hM fun ω => sq_nonneg _
    _ ≤ Real.exp (-m) / Z * (C * ∑ ω, ν ω * Γ ω) := mul_le_mul_of_nonneg_left hPI hK
    _ ≤ Real.exp (-m) / Z * (C * (Real.exp M * Z * ∑ ω, gibbs ν V ω * Γ ω)) :=
        mul_le_mul_of_nonneg_left
          (mul_le_mul_of_nonneg_left (sum_le_gibbs_sum hν0 hν1 hm hM hΓ) hC) hK
    _ = C * Real.exp (M - m) * ∑ ω, gibbs ν V ω * Γ ω := by
        rw [Real.exp_sub, Real.exp_neg]; field_simp

/-- **Holley–Stroock for the entropy inequality** (log-Sobolev form with
`g = f²`). -/
theorem holleyStroock_entropy {Γ : Ω → ℝ} (hΓ : ∀ ω, 0 ≤ Γ ω) {g : Ω → ℝ}
    (hg : ∀ ω, 0 ≤ g ω) {C : ℝ} (hC : 0 ≤ C) (hLS : went ν g ≤ C * ∑ ω, ν ω * Γ ω) :
    went (gibbs ν V) g ≤ C * Real.exp (M - m) * ∑ ω, gibbs ν V ω * Γ ω := by
  have hZ := partition_pos hν0 hν1 hM
  set Z := ∑ ω, ν ω * Real.exp (-V ω)
  have hK : 0 ≤ Real.exp (-m) / Z := div_nonneg (Real.exp_pos _).le hZ.le
  have hRHS : 0 ≤ C * Real.exp (M - m) * ∑ ω, gibbs ν V ω * Γ ω :=
    mul_nonneg (mul_nonneg hC (Real.exp_pos _).le)
      (Finset.sum_nonneg fun ω _ => mul_nonneg (gibbs_nonneg hν0 hν1 hM ω) (hΓ ω))
  have hm0 : 0 ≤ wmean ν g := Finset.sum_nonneg fun ω _ => mul_nonneg (hν0 ω) (hg ω)
  rcases hm0.eq_or_lt with h0 | hpos
  · -- `g = 0` on the support of `ν`, hence on that of `μ`: both entropies vanish
    have hzero : ∀ ω, ν ω * g ω = 0 := by
      intro ω
      have := (Finset.sum_eq_zero_iff_of_nonneg
        (fun ω _ => mul_nonneg (hν0 ω) (hg ω))).mp h0.symm ω (Finset.mem_univ _)
      exact this
    have hμ : ∀ ω, gibbs ν V ω * g ω = 0 := by
      intro ω; unfold gibbs
      rw [div_mul_eq_mul_div, mul_assoc, mul_comm (Real.exp _), ← mul_assoc, hzero ω]; simp
    have hμlog : ∀ ω, gibbs ν V ω * (g ω * Real.log (g ω)) = 0 := by
      intro ω; rw [← mul_assoc, hμ ω, zero_mul]
    have hent : went (gibbs ν V) g = 0 := by
      unfold went wmean
      simp [hμ, hμlog]
    rw [hent]; exact hRHS
  · calc went (gibbs ν V) g
        ≤ ∑ ω, gibbs ν V ω * entKernel (g ω) (wmean ν g) :=
          went_le (gibbs_nonneg hν0 hν1 hM) (gibbs_sum hν0 hν1 hM) hg hpos
      _ ≤ Real.exp (-m) / Z * ∑ ω, ν ω * entKernel (g ω) (wmean ν g) :=
          gibbs_sum_le hν0 hν1 hm hM fun ω => entKernel_nonneg (hg ω) hpos
      _ = Real.exp (-m) / Z * went ν g := by rw [went_eq_sum_entKernel hν0 hν1 g]
      _ ≤ Real.exp (-m) / Z * (C * ∑ ω, ν ω * Γ ω) := mul_le_mul_of_nonneg_left hLS hK
      _ ≤ Real.exp (-m) / Z * (C * (Real.exp M * Z * ∑ ω, gibbs ν V ω * Γ ω)) :=
          mul_le_mul_of_nonneg_left
            (mul_le_mul_of_nonneg_left (sum_le_gibbs_sum hν0 hν1 hm hM hΓ) hC) hK
      _ = C * Real.exp (M - m) * ∑ ω, gibbs ν V ω * Γ ω := by
          rw [Real.exp_sub, Real.exp_neg]; field_simp

end Gibbs

end YangMills
