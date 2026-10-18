import Mathlib

/-!
# Bounds on the Wilson action and the volume obstruction

For a unitary `N × N` matrix `U`, `|Re tr U| ≤ N`. Hence the Wilson plaquette
term `N - Re tr U_P` lies in `[0, 2N]`, the total action over a finite set `P`
of plaquettes lies in `[0, 2N |P|]`, and its oscillation is at most `2N |P|`.
The Holley–Stroock factor `exp (-(β · osc S))` therefore tends to `0` as the
number of plaquettes grows, at every fixed `β > 0`.

The plaquette holonomies are an arbitrary map `hol : P → Ω → U(N)` from
configurations to unitary matrices; any product of link variables around a
plaquette is of this form, since the unitary group is closed under products.
-/

open Matrix Filter Topology Real

namespace YangMills

variable {N : ℕ}

theorem norm_trace_le_of_unitary {U : Matrix (Fin N) (Fin N) ℂ}
    (hU : U ∈ Matrix.unitaryGroup (Fin N) ℂ) : ‖U.trace‖ ≤ N := by
  unfold Matrix.trace
  calc ‖∑ i, U.diag i‖ ≤ ∑ i, ‖U.diag i‖ := norm_sum_le _ _
    _ ≤ ∑ _i : Fin N, (1 : ℝ) :=
        Finset.sum_le_sum fun i _ => entry_norm_bound_of_unitary hU i i
    _ = N := by simp

theorem abs_re_trace_le_of_unitary {U : Matrix (Fin N) (Fin N) ℂ}
    (hU : U ∈ Matrix.unitaryGroup (Fin N) ℂ) : |U.trace.re| ≤ N :=
  (Complex.abs_re_le_norm _).trans (norm_trace_le_of_unitary hU)

/-- The Wilson plaquette term `N - Re tr U`. -/
noncomputable def plaquetteTerm (U : Matrix (Fin N) (Fin N) ℂ) : ℝ := N - U.trace.re

theorem plaquetteTerm_nonneg {U : Matrix (Fin N) (Fin N) ℂ}
    (hU : U ∈ Matrix.unitaryGroup (Fin N) ℂ) : 0 ≤ plaquetteTerm U := by
  have := abs_re_trace_le_of_unitary hU
  unfold plaquetteTerm; linarith [le_abs_self U.trace.re]

theorem plaquetteTerm_le {U : Matrix (Fin N) (Fin N) ℂ}
    (hU : U ∈ Matrix.unitaryGroup (Fin N) ℂ) : plaquetteTerm U ≤ 2 * N := by
  have := abs_re_trace_le_of_unitary hU
  unfold plaquetteTerm; linarith [neg_abs_le U.trace.re]

variable {Ω P : Type*} [Fintype P]

/-- The Wilson action `S(ω) = Σ_P (N - Re tr U_P(ω))`. -/
noncomputable def wilsonAction (hol : P → Ω → Matrix (Fin N) (Fin N) ℂ) (ω : Ω) : ℝ :=
  ∑ p, plaquetteTerm (hol p ω)

theorem wilsonAction_nonneg {hol : P → Ω → Matrix (Fin N) (Fin N) ℂ}
    (hU : ∀ p ω, hol p ω ∈ Matrix.unitaryGroup (Fin N) ℂ) (ω : Ω) :
    0 ≤ wilsonAction hol ω :=
  Finset.sum_nonneg fun p _ => plaquetteTerm_nonneg (hU p ω)

theorem wilsonAction_le {hol : P → Ω → Matrix (Fin N) (Fin N) ℂ}
    (hU : ∀ p ω, hol p ω ∈ Matrix.unitaryGroup (Fin N) ℂ) (ω : Ω) :
    wilsonAction hol ω ≤ 2 * N * Fintype.card P := by
  unfold wilsonAction
  calc ∑ p, plaquetteTerm (hol p ω) ≤ ∑ _p : P, (2 * N : ℝ) :=
        Finset.sum_le_sum fun p _ => plaquetteTerm_le (hU p ω)
    _ = 2 * N * Fintype.card P := by simp [mul_comm]

/-- The oscillation of `β S` is at most `2 N β |P|`. -/
theorem wilsonAction_osc_le {hol : P → Ω → Matrix (Fin N) (Fin N) ℂ}
    (hU : ∀ p ω, hol p ω ∈ Matrix.unitaryGroup (Fin N) ℂ) {β : ℝ} (hβ : 0 ≤ β)
    (ω ω' : Ω) :
    β * wilsonAction hol ω - β * wilsonAction hol ω' ≤ 2 * N * β * Fintype.card P := by
  have h1 := wilsonAction_le hU ω
  have h2 := wilsonAction_nonneg hU ω'
  have : β * wilsonAction hol ω ≤ β * (2 * N * Fintype.card P) :=
    mul_le_mul_of_nonneg_left h1 hβ
  nlinarith

/-- **Volume obstruction.** The Holley–Stroock factor obtained by perturbing the
product Haar measure by the full Wilson weight vanishes as the number of
plaquettes `m` tends to infinity, at every fixed `β > 0` and `N ≥ 1`. -/
theorem holleyStroock_factor_tendsto_zero {β : ℝ} (hβ : 0 < β) (hN : 1 ≤ N) :
    Tendsto (fun m : ℕ => Real.exp (-(2 * N * β * m))) atTop (𝓝 0) := by
  have hc : 0 < 2 * (N : ℝ) * β := by
    have : (1 : ℝ) ≤ N := by exact_mod_cast hN
    positivity
  have h1 : Tendsto (fun m : ℕ => 2 * (N : ℝ) * β * m) atTop atTop :=
    tendsto_natCast_atTop_atTop.const_mul_atTop hc
  exact Real.tendsto_exp_neg_atTop_nhds_zero.comp h1

/-! ### Plaquette variables are not independent in dimension three or more

On a unit cube with abelian link variables `θ_e` (the maximal torus, or the
gauge group `U(1)`), the six oriented plaquette angles sum to zero identically.
So the plaquette variables of a lattice of dimension `d ≥ 3` satisfy exact
linear relations and cannot be independent under the product Haar measure on
links. The labels are: edges along `x` at heights `(y,z)`, along `y` at
`(x,z)`, along `z` at `(x,y)`, each coordinate in `{0,1}`. -/

/-- The oriented plaquette angle `θ₁ + θ₂ - θ₃ - θ₄`. -/
def plaq {G : Type*} [AddCommGroup G] (θ₁ θ₂ θ₃ θ₄ : G) : G := θ₁ + θ₂ - θ₃ - θ₄

theorem cube_bianchi {G : Type*} [AddCommGroup G]
    (x y z : Fin 2 → Fin 2 → G) :
    -- faces z = 0 and z = 1 (normal z)
    plaq (x 0 0) (y 1 0) (x 1 0) (y 0 0) - plaq (x 0 1) (y 1 1) (x 1 1) (y 0 1)
    -- faces y = 0 and y = 1 (normal y)
    - plaq (x 0 0) (z 1 0) (x 0 1) (z 0 0) + plaq (x 1 0) (z 1 1) (x 1 1) (z 0 1)
    -- faces x = 0 and x = 1 (normal x)
    + plaq (y 0 0) (z 0 1) (y 0 1) (z 0 0) - plaq (y 1 0) (z 1 1) (y 1 1) (z 1 0)
    = 0 := by
  unfold plaq; abel

end YangMills
