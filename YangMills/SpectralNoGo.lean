import Mathlib

/-!
# The spacing-uniform gap no-go theorem

Reflection positivity gives, for a gauge-invariant observable localised in a
time slice, the spectral representation of its truncated temporal two-point
function in lattice units,
`S_a(t) = ∫ exp (-(E t / a)) dρ_a(E)`,
with `ρ_a` a finite positive measure on the spectrum of `-log T_a` above the
vacuum. A lattice gap `δ` uniform in the spacing means `ρ_a (Iio δ) = 0` for
every `a`.

This file proves, for arbitrary positive renormalisations `Z_a`:

* `twoPoint_ratio_le` : `S_a(t₂) ≤ exp (-(δ (t₂ - t₁) / a)) S_a(t₁)`;
* `nogo_later_limit_zero` : if `Z_a S_a(t₁)` and `Z_a S_a(t₂)` converge as
  `a → 0⁺` (along any filter, so along any sequence of spacings), the later
  limit is `0`;
* `kallenLehmann_pos` : a nonzero finite spectral measure supported in
  `[Δ, ∞)` with `Δ > 0` gives a two-point function strictly positive at every
  `t > 0`;
* `nogo_no_gapped_limit` : no family with a spacing-uniform lattice gap has
  renormalised two-point functions converging to such a continuum two-point
  function.
-/

open MeasureTheory Filter Topology Real Set

namespace YangMills

/-- The temporal two-point function at spacing `a` and separation `t`, from a
spectral measure `ρ` on lattice-unit energies. On `t ∈ aℤ` this is the lattice
two-point function; the formula defines it for every real `t`. -/
noncomputable def twoPoint (ρ : Measure ℝ) (a t : ℝ) : ℝ :=
  ∫ E, Real.exp (-(E * t / a)) ∂ρ

lemma ae_ge_of_Iio_null {ρ : Measure ℝ} {δ : ℝ} (h : ρ (Iio δ) = 0) :
    ∀ᵐ E ∂ρ, δ ≤ E := by
  rw [ae_iff]
  have : {E : ℝ | ¬δ ≤ E} = Iio δ := by ext E; simp
  rw [this]; exact h

lemma integrable_twoPoint {ρ : Measure ℝ} [IsFiniteMeasure ρ] {δ a t : ℝ}
    (hδ : 0 ≤ δ) (ha : 0 < a) (ht : 0 ≤ t) (h : ρ (Iio δ) = 0) :
    Integrable (fun E => Real.exp (-(E * t / a))) ρ := by
  refine Integrable.of_bound (C := 1) ?_ ?_
  · exact (Real.continuous_exp.comp
      ((continuous_id.mul continuous_const).div_const a).neg).aestronglyMeasurable
  · filter_upwards [ae_ge_of_Iio_null h] with E hE
    rw [Real.norm_eq_abs, abs_of_pos (Real.exp_pos _)]
    apply Real.exp_le_one_iff.mpr
    have : 0 ≤ E * t / a := div_nonneg (mul_nonneg (hδ.trans hE) ht) ha.le
    linarith

lemma twoPoint_nonneg (ρ : Measure ℝ) (a t : ℝ) : 0 ≤ twoPoint ρ a t :=
  integral_nonneg fun _ => (Real.exp_pos _).le

/-- The ratio bound. Every renormalisation constant cancels in it. -/
theorem twoPoint_ratio_le {ρ : Measure ℝ} [IsFiniteMeasure ρ] {δ a t₁ t₂ : ℝ}
    (hδ : 0 ≤ δ) (ha : 0 < a) (ht₁ : 0 ≤ t₁) (ht : t₁ ≤ t₂)
    (h : ρ (Iio δ) = 0) :
    twoPoint ρ a t₂ ≤ Real.exp (-(δ * (t₂ - t₁) / a)) * twoPoint ρ a t₁ := by
  unfold twoPoint
  rw [← integral_const_mul]
  apply integral_mono_ae (integrable_twoPoint hδ ha (ht₁.trans ht) h)
    ((integrable_twoPoint hδ ha ht₁ h).const_mul _)
  filter_upwards [ae_ge_of_Iio_null h] with E hE
  rw [← Real.exp_add]
  apply Real.exp_le_exp.mpr
  have hdt : 0 ≤ t₂ - t₁ := sub_nonneg.mpr ht
  have key : δ * (t₂ - t₁) ≤ E * (t₂ - t₁) := mul_le_mul_of_nonneg_right hE hdt
  have : (δ * (t₂ - t₁)) / a ≤ (E * (t₂ - t₁)) / a := div_le_div_of_nonneg_right key ha.le
  have e1 : E * t₂ / a = E * (t₂ - t₁) / a + E * t₁ / a := by ring
  rw [e1]; linarith

/-- `exp (-(c / a)) → 0` as `a → 0⁺` for `c > 0`. -/
lemma tendsto_exp_neg_div {c : ℝ} (hc : 0 < c) :
    Tendsto (fun a : ℝ => Real.exp (-(c / a))) (𝓝[>] 0) (𝓝 0) := by
  have h1 : Tendsto (fun a : ℝ => c * a⁻¹) (𝓝[>] 0) atTop :=
    tendsto_inv_nhdsGT_zero.const_mul_atTop hc
  have h2 := Real.tendsto_exp_neg_atTop_nhds_zero.comp h1
  refine h2.congr fun a => ?_
  simp [Function.comp, div_eq_mul_inv]

/-- **No-go theorem.** Let a family of lattice theories, indexed along any
filter `l` on which the spacing `a i` tends to `0⁺`, have a lattice gap `δ > 0`
uniform in the spacing. For any positive renormalisations `Z i`, if the
renormalised two-point functions converge at two separations `t₁ < t₂`, then
the limit at `t₂` is zero. -/
theorem nogo_later_limit_zero {ι : Type*} {l : Filter ι} [l.NeBot]
    {a : ι → ℝ} (ha : Tendsto a l (𝓝[>] 0))
    {ρ : ι → Measure ℝ} [∀ i, IsFiniteMeasure (ρ i)] {δ : ℝ} (hδ : 0 < δ)
    (hsupp : ∀ i, ρ i (Iio δ) = 0) {Z : ι → ℝ} (hZ : ∀ i, 0 < Z i)
    {t₁ t₂ : ℝ} (ht₁ : 0 ≤ t₁) (ht : t₁ < t₂) {L₁ L₂ : ℝ}
    (h₁ : Tendsto (fun i => Z i * twoPoint (ρ i) (a i) t₁) l (𝓝 L₁))
    (h₂ : Tendsto (fun i => Z i * twoPoint (ρ i) (a i) t₂) l (𝓝 L₂)) :
    L₂ = 0 := by
  have hc : 0 < δ * (t₂ - t₁) := mul_pos hδ (sub_pos.mpr ht)
  have hexp : Tendsto (fun i => Real.exp (-(δ * (t₂ - t₁) / a i))) l (𝓝 0) :=
    (tendsto_exp_neg_div hc).comp ha
  have hupper : Tendsto
      (fun i => Real.exp (-(δ * (t₂ - t₁) / a i)) * (Z i * twoPoint (ρ i) (a i) t₁))
      l (𝓝 0) := by
    simpa using hexp.mul h₁
  have hpos : ∀ᶠ i in l, 0 < a i := ha (self_mem_nhdsWithin)
  have hzero : Tendsto (fun i => Z i * twoPoint (ρ i) (a i) t₂) l (𝓝 0) := by
    refine tendsto_of_tendsto_of_tendsto_of_le_of_le' tendsto_const_nhds hupper ?_ ?_
    · exact Eventually.of_forall fun i => mul_nonneg (hZ i).le (twoPoint_nonneg _ _ _)
    · filter_upwards [hpos] with i hi
      have := twoPoint_ratio_le (ρ := ρ i) hδ.le hi ht₁ ht.le (hsupp i)
      calc Z i * twoPoint (ρ i) (a i) t₂
          ≤ Z i * (Real.exp (-(δ * (t₂ - t₁) / a i)) * twoPoint (ρ i) (a i) t₁) :=
            mul_le_mul_of_nonneg_left this (hZ i).le
        _ = _ := by ring
  exact tendsto_nhds_unique h₂ hzero

/-- **Källén–Lehmann positivity.** A nonzero finite spectral measure supported
in `[Δ, ∞)` with `Δ > 0` gives a two-point function `∫ exp (-E t) dσ`
strictly positive at every `t > 0`. -/
theorem kallenLehmann_pos {σ : Measure ℝ} [IsFiniteMeasure σ] (hσ : σ ≠ 0)
    {Δ t : ℝ} (hΔ : 0 < Δ) (ht : 0 < t) (hsupp : σ (Iio Δ) = 0) :
    0 < ∫ E, Real.exp (-(E * t)) ∂σ := by
  have hint : Integrable (fun E => Real.exp (-(E * t))) σ := by
    have := integrable_twoPoint (ρ := σ) hΔ.le one_pos ht.le hsupp
    simpa using this
  rw [integral_pos_iff_support_of_nonneg (fun _ => (Real.exp_pos _).le) hint]
  have : Function.support (fun E : ℝ => Real.exp (-(E * t))) = univ := by
    ext E; simp [(Real.exp_pos _).ne']
  rw [this]
  exact Measure.measure_univ_pos.mpr hσ

/-- **No gapped continuum limit from a spacing-uniform lattice gap.** If the
renormalised two-point functions of a family with a spacing-uniform lattice gap
converged at every `t > 0` to the two-point function of a continuum theory with
gap `Δ > 0` and nonzero spectral measure, we would have a contradiction. -/
theorem nogo_no_gapped_limit {ι : Type*} {l : Filter ι} [l.NeBot]
    {a : ι → ℝ} (ha : Tendsto a l (𝓝[>] 0))
    {ρ : ι → Measure ℝ} [∀ i, IsFiniteMeasure (ρ i)] {δ : ℝ} (hδ : 0 < δ)
    (hsupp : ∀ i, ρ i (Iio δ) = 0) {Z : ι → ℝ} (hZ : ∀ i, 0 < Z i)
    {σ : Measure ℝ} [IsFiniteMeasure σ] (hσ : σ ≠ 0) {Δ : ℝ} (hΔ : 0 < Δ)
    (hσsupp : σ (Iio Δ) = 0) :
    ¬ ∀ t > 0, Tendsto (fun i => Z i * twoPoint (ρ i) (a i) t) l
        (𝓝 (∫ E, Real.exp (-(E * t)) ∂σ)) := by
  intro h
  have h0 := nogo_later_limit_zero ha hδ hsupp hZ zero_le_one one_lt_two
    (h 1 one_pos) (h 2 two_pos)
  exact (kallenLehmann_pos hσ hΔ two_pos hσsupp).ne' h0

/-- The renormalised target is incompatible with a spacing-uniform gap: if
`Δ_lat(a) / a → m` as `a → 0⁺`, then `Δ_lat(a) → 0`, so no `δ > 0` bounds it
below for all small `a`. -/
theorem target_excludes_uniform_gap {Δlat : ℝ → ℝ} {m : ℝ}
    (h : Tendsto (fun a => Δlat a / a) (𝓝[>] 0) (𝓝 m)) :
    ¬ ∃ δ > 0, ∀ a > 0, δ ≤ Δlat a := by
  rintro ⟨δ, hδ, hbound⟩
  have hid : Tendsto (fun a : ℝ => a) (𝓝[>] 0) (𝓝 0) :=
    tendsto_nhdsWithin_of_tendsto_nhds tendsto_id
  have hprod : Tendsto (fun a => Δlat a / a * a) (𝓝[>] 0) (𝓝 0) := by
    simpa using h.mul hid
  have heq : ∀ᶠ a in 𝓝[>] (0:ℝ), Δlat a / a * a = Δlat a := by
    filter_upwards [self_mem_nhdsWithin] with a (ha : 0 < a)
    field_simp
  have hΔ : Tendsto Δlat (𝓝[>] 0) (𝓝 0) := (tendsto_congr' heq).mp hprod
  have hge : ∀ᶠ a in 𝓝[>] (0:ℝ), δ ≤ Δlat a := by
    filter_upwards [self_mem_nhdsWithin] with a (ha : 0 < a)
    exact hbound a ha
  have := ge_of_tendsto hΔ hge
  linarith

end YangMills
