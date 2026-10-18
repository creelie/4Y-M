import YangMills.SpectralNoGo

/-!
# The critical scaling of the lattice gap

The no-go theorem of `SpectralNoGo.lean` excludes a lattice gap bounded below
in lattice units. This file records the bound that is *not* excluded and that
the Clay statement needs: a lattice gap `Δ_lat(a) ≥ m_a · a` with `m_a → m`.
Under that bound the limiting two-point functions satisfy the ratio inequality
`L(t₂) ≤ exp (-(m (t₂ - t₁))) L(t₁)`, which is the statement that the limit has
spectral gap at least `m`. Nothing here proves such a bound for Yang–Mills; it
only shows what the residual lower bound has to be.
-/

open MeasureTheory Filter Topology Real Set

namespace YangMills

/-- If the lattice gaps scale as `m i * a i` with `m i → m`, and the renormalised
two-point functions converge at `t₁ ≤ t₂`, the limits obey the continuum ratio
bound with rate `m`. -/
theorem limit_ratio_le {ι : Type*} {l : Filter ι} [l.NeBot]
    {a : ι → ℝ} (ha : ∀ᶠ i in l, 0 < a i)
    {ρ : ι → Measure ℝ} [∀ i, IsFiniteMeasure (ρ i)]
    {mseq : ι → ℝ} {m : ℝ} (hm : Tendsto mseq l (𝓝 m)) (hm0 : ∀ i, 0 ≤ mseq i)
    (hsupp : ∀ i, ρ i (Iio (mseq i * a i)) = 0)
    {Z : ι → ℝ} (hZ : ∀ i, 0 < Z i)
    {t₁ t₂ : ℝ} (ht₁ : 0 ≤ t₁) (ht : t₁ ≤ t₂) {L₁ L₂ : ℝ}
    (h₁ : Tendsto (fun i => Z i * twoPoint (ρ i) (a i) t₁) l (𝓝 L₁))
    (h₂ : Tendsto (fun i => Z i * twoPoint (ρ i) (a i) t₂) l (𝓝 L₂)) :
    L₂ ≤ Real.exp (-(m * (t₂ - t₁))) * L₁ := by
  have hrhs : Tendsto
      (fun i => Real.exp (-(mseq i * (t₂ - t₁))) * (Z i * twoPoint (ρ i) (a i) t₁))
      l (𝓝 (Real.exp (-(m * (t₂ - t₁))) * L₁)) := by
    refine Tendsto.mul ?_ h₁
    exact (Real.continuous_exp.tendsto _).comp ((hm.mul_const _).neg)
  refine le_of_tendsto_of_tendsto h₂ hrhs ?_
  filter_upwards [ha] with i hi
  · have hδ : 0 ≤ mseq i * a i := mul_nonneg (hm0 i) hi.le
    have := twoPoint_ratio_le (ρ := ρ i) hδ hi ht₁ ht (hsupp i)
    have e : mseq i * a i * (t₂ - t₁) / a i = mseq i * (t₂ - t₁) := by
      field_simp
    rw [e] at this
    calc Z i * twoPoint (ρ i) (a i) t₂
        ≤ Z i * (Real.exp (-(mseq i * (t₂ - t₁))) * twoPoint (ρ i) (a i) t₁) :=
          mul_le_mul_of_nonneg_left this (hZ i).le
      _ = _ := by ring

end YangMills
