import Mathlib
import YangMills.HolleyStroock

/-!
# Bounded nonnegative local terms do not give a volume-uniform cumulant bound

The second derivative of the finite-volume free energy
`f_Λ(β) = -|Λ|⁻¹ log ∫ e^{-β S}` is `-Var_β(S) / |Λ|`. A bound on it uniform in
`|Λ|` therefore needs `Var(S) = O(|Λ|)`, a decorrelation property. Boundedness
and nonnegativity of each local term do not give it: below, `n` local terms
each in `[0, 2N]`, all equal to `2N Y` with `Y` a fair coin, give
`Var(S_n) / n = N² n`, unbounded in the volume `n`. Independence of the local
terms is what would exclude this; plaquette variables are not independent in
dimension `d ≥ 3` (`YangMills.cube_bianchi`).
-/

open Real

namespace YangMills

/-- The fair coin on `Bool`. -/
noncomputable def coin : Bool → ℝ := fun _ => 1 / 2

/-- One local term, equal to `2N Y` with `Y ∈ {0,1}`. -/
def localTerm (N : ℕ) (b : Bool) : ℝ := if b then 2 * N else 0

theorem localTerm_mem (N : ℕ) (b : Bool) : 0 ≤ localTerm N b ∧ localTerm N b ≤ 2 * N := by
  cases b <;> simp [localTerm]

/-- The action of volume `n`: the sum of `n` copies of the local term. -/
def corrAction (N n : ℕ) (b : Bool) : ℝ := ∑ _i : Fin n, localTerm N b

theorem corrAction_var (N n : ℕ) :
    wvar coin (corrAction N n) = (N : ℝ) ^ 2 * n ^ 2 := by
  simp [wvar, wmean, coin, corrAction, localTerm]
  ring

/-- The variance per unit volume is `N² n`, unbounded in `n` for `N ≥ 1`. -/
theorem corrAction_var_per_volume_unbounded {N : ℕ} (hN : 1 ≤ N) (K : ℝ) :
    ∃ n : ℕ, 0 < n ∧ K < wvar coin (corrAction N n) / n := by
  obtain ⟨n, hn⟩ := exists_nat_gt (max K 0)
  have hnpos : (0 : ℝ) < n := lt_of_le_of_lt (le_max_right _ _) hn
  refine ⟨n, by exact_mod_cast hnpos, ?_⟩
  rw [corrAction_var, show (N : ℝ) ^ 2 * n ^ 2 / n = N ^ 2 * n by field_simp]
  have hN1 : (1 : ℝ) ≤ (N : ℝ) ^ 2 := by
    have : (1 : ℝ) ≤ N := by exact_mod_cast hN
    nlinarith
  have : (n : ℝ) ≤ N ^ 2 * n := by nlinarith
  linarith [le_max_left K 0]

end YangMills
