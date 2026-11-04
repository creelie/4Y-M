import Mathlib

/-!
# Finite and algebraic parts of the averaged cube condition and of the ends

This file checks the counting, the SU(2) trace identity and the arithmetic of
the constants on which the averaged cube condition `BP(n, ϖ)` and the bounds
on the ends of the pair factors rest. The analytic steps (Haar measure,
couplings, the spectral argument) are proved in the paper and are not
formalised here.

* `card_box`, `card_box_inter_shift` : the cube `Q = [0,n)^4` and its
  intersection with a translate;
* `count_unit`, `count_diag` : the points `z ∉ Q` with `z + Δ ∈ Q`, for a unit
  vector and for `e_ν - e_μ`; `partner_count` : `18 n³ - 3 n²`;
* `trace_commutator_two` : `tr(A B adj A adj B) = det B · x² + det A · y² + w²
  - x y w - 2 det A det B`, the Fricke identity for `SL(2)`;
* `commutator_bound` : the bound `M > 0.64` of the almost-commuting lemma;
* `bp_ranges`, `bp_fixed_constants`, `bpflat_constant` : the constants of the
  strong-coupling ranges and of the fixed-cube obstruction;
* `ends_constants` : the decimals of the derivative bounds in the ends.
-/

open Finset

namespace YangMills.Averaged

/-! ## Cubes in `ℤ⁴` -/

/-- The cube `[0,n)^4`. -/
noncomputable def box (n : ℕ) : Finset (Fin 4 → ℤ) :=
  Fintype.piFinset fun _ => Ico (0 : ℤ) n

/-- The set of `z` with `z + Δ ∈ [0,n)^4`. -/
noncomputable def shiftBox (n : ℕ) (Δ : Fin 4 → ℤ) : Finset (Fin 4 → ℤ) :=
  Fintype.piFinset fun i => Ico (-Δ i) (n - Δ i)

theorem card_box (n : ℕ) : (box n).card = n ^ 4 := by
  simp [box, Fintype.card_piFinset]

theorem card_shiftBox (n : ℕ) (Δ : Fin 4 → ℤ) : (shiftBox n Δ).card = n ^ 4 := by
  simp [shiftBox, Fintype.card_piFinset]

theorem mem_shiftBox (n : ℕ) (Δ z : Fin 4 → ℤ) :
    z ∈ shiftBox n Δ ↔ z + Δ ∈ box n := by
  simp only [shiftBox, box, Fintype.mem_piFinset, mem_Ico, Pi.add_apply]
  constructor <;> intro h i <;> have := h i <;> constructor <;> omega

theorem box_inter_shift (n : ℕ) (Δ : Fin 4 → ℤ) :
    box n ∩ shiftBox n Δ =
      Fintype.piFinset fun i => Ico (max 0 (-Δ i)) (min (n : ℤ) (n - Δ i)) := by
  ext z
  simp only [box, shiftBox, mem_inter, Fintype.mem_piFinset, mem_Ico]
  constructor
  · rintro ⟨h1, h2⟩ i
    have := h1 i; have := h2 i
    constructor
    · exact max_le (by omega) (by omega)
    · exact lt_min (by omega) (by omega)
  · intro h
    refine ⟨fun i => ?_, fun i => ?_⟩ <;> have := h i <;>
      constructor <;> omega

theorem card_box_inter_shift (n : ℕ) (Δ : Fin 4 → ℤ) :
    (box n ∩ shiftBox n Δ).card =
      ∏ i, (min (n : ℤ) (n - Δ i) - max 0 (-Δ i)).toNat := by
  rw [box_inter_shift, Fintype.card_piFinset]
  simp [Int.card_Ico]

/-- Points `z ∉ Q` with `z + Δ ∈ Q`: their number is `n⁴ - |Q ∩ (Q - Δ)|`. -/
theorem card_outside (n : ℕ) (Δ : Fin 4 → ℤ) :
    (shiftBox n Δ \ box n).card + (box n ∩ shiftBox n Δ).card = n ^ 4 := by
  rw [← card_shiftBox n Δ, inter_comm, card_sdiff_add_card_inter]

theorem factor_zero (n : ℕ) :
    (min (n : ℤ) (n - 0) - max 0 (-0)).toNat = n := by simp

theorem factor_plus (n : ℕ) (hn : 1 ≤ n) :
    (min (n : ℤ) (n - 1) - max 0 (-1)).toNat = n - 1 := by
  rw [min_eq_right (by omega), max_eq_left (by omega)]; omega

theorem factor_minus (n : ℕ) (hn : 1 ≤ n) :
    (min (n : ℤ) (n - (-1)) - max 0 (-(-1))).toNat = n - 1 := by
  rw [min_eq_left (by omega), max_eq_right (by omega)]; omega

/-- For `Δ = ±e_j` the cube meets its translate in `n³(n-1)` points, so
`n³` points `z ∉ Q` have `z + Δ ∈ Q`. -/
theorem count_unit (n : ℕ) (hn : 1 ≤ n) (j : Fin 4) (s : ℤ) (hs : s = 1 ∨ s = -1) :
    (shiftBox n (Pi.single j s) \ box n).card = n ^ 3 := by
  have h := card_outside n (Pi.single j s)
  rw [card_box_inter_shift] at h
  have hprod : ∏ i, (min (n : ℤ) (n - (Pi.single j s : Fin 4 → ℤ) i) - max 0 (-(Pi.single j s : Fin 4 → ℤ) i)).toNat
      = n ^ 3 * (n - 1) := by
    rw [← Finset.mul_prod_erase _ _ (mem_univ j)]
    have hj : (min (n : ℤ) (n - (Pi.single j s : Fin 4 → ℤ) j) - max 0 (-(Pi.single j s : Fin 4 → ℤ) j)).toNat = n - 1 := by
      rcases hs with rfl | rfl
      · simp
      · simp
    have hrest : ∀ i ∈ (univ : Finset (Fin 4)).erase j,
        (min (n : ℤ) (n - (Pi.single j s : Fin 4 → ℤ) i) - max 0 (-(Pi.single j s : Fin 4 → ℤ) i)).toNat = n := by
      intro i hi
      rw [Pi.single_eq_of_ne (ne_of_mem_erase hi)]
      exact factor_zero n
    rw [hj, prod_congr rfl hrest, prod_const, card_erase_of_mem (mem_univ j)]
    simp [mul_comm]
  rw [hprod] at h
  have e : n ^ 3 * (n - 1) + n ^ 3 = n ^ 4 := by
    obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_le hn
    simp only [Nat.add_sub_cancel_left]; ring
  exact Nat.add_right_cancel (h.trans (e.symm.trans (Nat.add_comm _ _)))

/-- For `Δ = e_ν - e_μ` the cube meets its translate in `n²(n-1)²` points, so
`2n³ - n²` points `z ∉ Q` have `z + Δ ∈ Q`. -/
theorem count_diag (n : ℕ) (hn : 1 ≤ n) (μ ν : Fin 4) (hμν : μ ≠ ν) :
    (shiftBox n (Pi.single ν 1 - Pi.single μ 1) \ box n).card = 2 * n ^ 3 - n ^ 2 := by
  have h := card_outside n (Pi.single ν 1 - Pi.single μ 1)
  rw [card_box_inter_shift] at h
  set F : Fin 4 → ℕ := fun i =>
    (min (n : ℤ) (n - (Pi.single ν 1 - Pi.single μ 1 : Fin 4 → ℤ) i)
      - max 0 (-(Pi.single ν 1 - Pi.single μ 1 : Fin 4 → ℤ) i)).toNat
  have hν : F ν = n - 1 := by
    simp only [F, Pi.sub_apply, Pi.single_eq_same, Pi.single_eq_of_ne hμν.symm, sub_zero]
    exact factor_plus n hn
  have hμ : F μ = n - 1 := by
    simp only [F, Pi.sub_apply, Pi.single_eq_same, Pi.single_eq_of_ne hμν, zero_sub]
    exact factor_minus n hn
  have hrest : ∀ i, i ≠ μ → i ≠ ν → F i = n := by
    intro i hiμ hiν
    simp only [F, Pi.sub_apply, Pi.single_eq_of_ne hiμ, Pi.single_eq_of_ne hiν, sub_zero]
    exact factor_zero n
  have hprod : ∏ i, F i = n ^ 2 * (n - 1) ^ 2 := by
    rw [← Finset.mul_prod_erase _ _ (mem_univ ν),
      ← Finset.mul_prod_erase _ _ (mem_erase.mpr ⟨hμν, mem_univ μ⟩)]
    have hr : ∀ i ∈ ((univ : Finset (Fin 4)).erase ν).erase μ, F i = n := by
      intro i hi
      exact hrest i (ne_of_mem_erase hi) (ne_of_mem_erase (mem_of_mem_erase hi))
    rw [prod_congr rfl hr, prod_const, card_erase_of_mem
      (mem_erase.mpr ⟨hμν, mem_univ μ⟩), card_erase_of_mem (mem_univ ν), hν, hμ]
    simp; ring
  change (shiftBox n (Pi.single ν 1 - Pi.single μ 1) \ box n).card + ∏ i, F i = n ^ 4 at h
  rw [hprod] at h
  obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_le hn
  simp only [Nat.add_sub_cancel_left] at h ⊢
  have e1 : (1 + k) ^ 4 = (1 + k) ^ 2 * k ^ 2 + (2 * (1 + k) ^ 3 - (1 + k) ^ 2) := by
    have : (1 + k) ^ 2 ≤ 2 * (1 + k) ^ 3 := by nlinarith [Nat.zero_le k]
    zify [this]; ring
  exact Nat.add_right_cancel (h.trans (e1.trans (Nat.add_comm _ _)))

/-- The partners of a link: for each of the three directions `μ ≠ ν` the
translates `e_ν, e_μ, 0, -e_μ, -e_μ, e_ν - e_μ` contribute `6n³ - n²`, in all
`18n³ - 3n²` (Lemma `lem:bpcount`(i)). -/
theorem partner_count (n : ℕ) (hn : 1 ≤ n) (μ ν : Fin 4) (hμν : μ ≠ ν) :
    (shiftBox n (Pi.single ν 1) \ box n).card
      + (shiftBox n (Pi.single μ 1) \ box n).card
      + (shiftBox n 0 \ box n).card
      + 2 * (shiftBox n (Pi.single μ (-1)) \ box n).card
      + (shiftBox n (Pi.single ν 1 - Pi.single μ 1) \ box n).card
      = 6 * n ^ 3 - n ^ 2 := by
  have h0 : shiftBox n 0 \ box n = ∅ := by
    rw [sdiff_eq_empty_iff_subset]
    intro z hz
    simpa using (mem_shiftBox n 0 z).1 hz
  rw [count_unit n hn ν 1 (Or.inl rfl), count_unit n hn μ 1 (Or.inl rfl), h0,
    count_unit n hn μ (-1) (Or.inr rfl), count_diag n hn μ ν hμν, card_empty]
  have : n ^ 2 ≤ 2 * n ^ 3 := by nlinarith [Nat.one_le_pow 2 n hn]
  omega

theorem partner_total (n : ℕ) (hn : 1 ≤ n) :
    3 * (6 * n ^ 3 - n ^ 2) = 18 * n ^ 3 - 3 * n ^ 2 := by
  have : n ^ 2 ≤ 6 * n ^ 3 := by nlinarith [Nat.one_le_pow 2 n hn]
  omega


/-! ## The trace of a commutator in `SL(2)` -/

/-- The Fricke identity, for every `2 × 2` matrix over a commutative ring:
`tr(A B adj A adj B) = det B · (tr A)² + det A · (tr B)² + (tr AB)²
 - tr A · tr B · tr AB - 2 det A det B`. -/
theorem trace_commutator_two {R : Type*} [CommRing R] (A B : Matrix (Fin 2) (Fin 2) R) :
    Matrix.trace (A * B * A.adjugate * B.adjugate)
      = B.det * (Matrix.trace A) ^ 2 + A.det * (Matrix.trace B) ^ 2
        + (Matrix.trace (A * B)) ^ 2
        - Matrix.trace A * Matrix.trace B * Matrix.trace (A * B)
        - 2 * A.det * B.det := by
  simp [Matrix.trace_fin_two, Matrix.adjugate_fin_two, Matrix.det_fin_two,
    Matrix.mul_apply, Fin.sum_univ_two]
  ring

/-- For `det A = det B = 1`, where `adj A = A⁻¹`: `tr(ABA⁻¹B⁻¹) = x² + y² + w² - xyw - 2`. -/
theorem trace_commutator_sl2 {R : Type*} [CommRing R] (A B : Matrix (Fin 2) (Fin 2) R)
    (hA : A.det = 1) (hB : B.det = 1) :
    Matrix.trace (A * B * A.adjugate * B.adjugate)
      = (Matrix.trace A) ^ 2 + (Matrix.trace B) ^ 2 + (Matrix.trace (A * B)) ^ 2
        - Matrix.trace A * Matrix.trace B * Matrix.trace (A * B) - 2 := by
  rw [trace_commutator_two, hA, hB]; ring

/-- The bound of Lemma `lem:bpcommute` for `N = 2`: if the real traces satisfy
`x² + y² + w² - xyw ≥ 4 - √2`, then `M = x² + y² + w² > 0.64`. -/
theorem commutator_bound (x y w : ℝ) (h : 4 - Real.sqrt 2 ≤ x ^ 2 + y ^ 2 + w ^ 2 - x * y * w) :
    0.64 < x ^ 2 + y ^ 2 + w ^ 2 := by
  by_contra hc
  push Not at hc
  have hs : Real.sqrt 2 < 1.5 := by
    rw [Real.sqrt_lt' (by norm_num)]; norm_num
  have bx : |x| ≤ 0.8 := by rw [abs_le]; constructor <;> nlinarith [sq_nonneg y, sq_nonneg w]
  have by' : |y| ≤ 0.8 := by rw [abs_le]; constructor <;> nlinarith [sq_nonneg x, sq_nonneg w]
  have bw : |w| ≤ 0.8 := by rw [abs_le]; constructor <;> nlinarith [sq_nonneg x, sq_nonneg y]
  have hp : |x * y * w| ≤ 0.512 := by
    rw [abs_mul, abs_mul]
    have h1 : |x| * |y| ≤ 0.8 * 0.8 := mul_le_mul bx by' (abs_nonneg _) (by norm_num)
    calc |x| * |y| * |w| ≤ (0.8 * 0.8) * 0.8 :=
          mul_le_mul h1 bw (abs_nonneg _) (by norm_num)
      _ = 0.512 := by norm_num
  have := neg_abs_le (x * y * w)
  linarith

/-! ## Constants of the averaged condition -/

/-- The strong-coupling ranges of Corollary `cor:bpranges`(i) for `SU(2)`:
`α_n = β(18n-3)/((4-18β)n²)`. -/
noncomputable def alphaSU2 (β : ℝ) (n : ℝ) : ℝ := β * (18 * n - 3) / ((4 - 18 * β) * n ^ 2)

theorem alphaSU2_two_eighth : alphaSU2 (1 / 8) 2 = 33 / 56 := by
  norm_num [alphaSU2]

theorem alphaSU2_nine_fifth : alphaSU2 (1 / 5) 9 = 159 / 162 := by
  norm_num [alphaSU2]

/-- `α_n < 1` exactly when `β < 4n²/(18n² + 18n - 3)`, for `0 ≤ β < 2/9`. -/
theorem alphaSU2_lt_one_iff (β n : ℝ) (_hβ : 0 ≤ β) (hβ' : β < 2 / 9) (hn : 1 ≤ n) :
    alphaSU2 β n < 1 ↔ β < 4 * n ^ 2 / (18 * n ^ 2 + 18 * n - 3) := by
  have hd : 0 < (4 - 18 * β) * n ^ 2 := by
    have : 0 < 4 - 18 * β := by linarith
    positivity
  have hd' : 0 < 18 * n ^ 2 + 18 * n - 3 := by nlinarith
  rw [alphaSU2, div_lt_one hd, lt_div_iff₀ hd']
  constructor <;> intro h <;> nlinarith

theorem threshold_two : (4 : ℝ) * 2 ^ 2 / (18 * 2 ^ 2 + 18 * 2 - 3) = 16 / 105 := by
  norm_num

/-- The decimals of Corollary `cor:bpranges`(ii) for `SU(3)` at `β_W = 0.35`. -/
theorem su3_constants :
    (0.578568 : ℝ) < 0.76064 ^ 2 ∧
    (2.1 : ℝ) * 0.76064 * (640 / 1191) ≤ 0.8584 ∧
    (0.8584 : ℝ) / (18 * (1 - 0.8584)) ≤ 0.3368 ∧
    (0.3368 : ℝ) * (18 * 6 - 3) / 6 ^ 2 ≤ 0.9824 := by
  norm_num

/-- The constants of Theorem `thm:bpfixed`: for `N = 2` the Dirichlet bound
`8 N² (N-1)² L³ n ε` against the variance `0.64 (1 - 1/2)/3` gives `300`, and
`L ≤ 4n` gives `19200 n⁴ ε`; the other case gives `2048 ≤ 19200`; and
`ε ≤ 2 (log β + 6.11)/β` gives `38400`. -/
theorem bp_fixed_constants :
    (8 : ℝ) * 2 ^ 2 * (2 - 1) ^ 2 / (0.64 * (1 - 1 / 2) / 3) = 300 ∧
    (300 : ℝ) * 4 ^ 3 = 19200 ∧ (2 : ℝ) * (4 * 4 ^ 4) ≤ 19200 ∧
    (8 : ℝ) * 4 ^ 3 * 4 = 2048 ∧ (19200 : ℝ) * 2 = 38400 := by
  norm_num

/-- The constant of Lemma `lem:bpflat`: `log 32 + 1 + (2/3) log(3π³/8) ≤ 6.11`. -/
theorem bpflat_constant :
    Real.log 32 + 1 + 2 / 3 * Real.log (3 * Real.pi ^ 3 / 8) ≤ 6.11 := by
  have hpi := Real.pi_lt_d6
  have hpi0 := Real.pi_pos
  have h32 : Real.log 32 = 5 * Real.log 2 := by
    rw [show (32 : ℝ) = 2 ^ 5 by norm_num, Real.log_pow]; norm_num
  have hl2 := Real.log_two_lt_d9
  have he := Real.exp_one_gt_d9
  have hq : 1 + 0.4663 + 0.4663 ^ 2 / 2 ≤ Real.exp 0.4663 :=
    Real.quadratic_le_exp_of_nonneg (by norm_num)
  have hexp : Real.exp 2.4663 = Real.exp 1 ^ 2 * Real.exp 0.4663 := by
    rw [sq, ← Real.exp_add, ← Real.exp_add]; norm_num
  have hlt : 3 * Real.pi ^ 3 / 8 < Real.exp 2.4663 := by
    rw [hexp]
    have h1 : 3 * Real.pi ^ 3 / 8 < 3 * 3.141593 ^ 3 / 8 := by
      have : Real.pi ^ 3 < 3.141593 ^ 3 := by
        apply pow_lt_pow_left₀ hpi hpi0.le (by norm_num)
      linarith
    have h2 : (2.7182818283 : ℝ) ^ 2 < Real.exp 1 ^ 2 :=
      pow_lt_pow_left₀ he (by norm_num) (by norm_num)
    have h3 : (0 : ℝ) < Real.exp 0.4663 := Real.exp_pos _
    nlinarith
  have hlog : Real.log (3 * Real.pi ^ 3 / 8) < 2.4663 := by
    rw [Real.log_lt_iff_lt_exp (by positivity)]; exact hlt
  rw [h32]; linarith

/-! ## Decimals of the derivative bounds in the ends -/

theorem sqrt_two_bounds : (1.41421 : ℝ) < Real.sqrt 2 ∧ Real.sqrt 2 < 1.414214 := by
  constructor
  · rw [Real.lt_sqrt (by norm_num)]; norm_num
  · rw [Real.sqrt_lt' (by norm_num)]; norm_num

/-- Lemma `lem:endboundary`: `8√2 π > 35.54`, `7/35.54 < 0.19697`, and
`0.19697 · 1.706·10⁻⁷ < 3.4·10⁻⁸`. -/
theorem endboundary_constants :
    (35.54 : ℝ) < 8 * Real.sqrt 2 * Real.pi ∧ (7 : ℝ) / 35.54 < 0.19697 ∧
    (0.19697 : ℝ) * 1.706e-7 < 3.4e-8 := by
  obtain ⟨h1, _⟩ := sqrt_two_bounds
  have := Real.pi_gt_d6
  refine ⟨by nlinarith, by norm_num, by norm_num⟩

/-- Lemma `lem:endchart` and Corollary `cor:endchart`. -/
theorem endchart_constants :
    (1.828 : ℝ) * 70.06 < 128.07 ∧ (1.828 : ℝ) ^ 2 < 3.3416 ∧
    (3.3416 : ℝ) * 53.96 < 180.4 ∧ (3.3416 : ℝ) * 15.134 < 50.6 ∧
    128.07 + 3 * Real.sqrt 2 < 132.32 := by
  obtain ⟨_, h2⟩ := sqrt_two_bounds
  refine ⟨by norm_num, by norm_num, by norm_num, by norm_num, by linarith⟩

/-- Corollary `cor:tieareab`: `14/(2√2 π) < 1.57556` and `50.6 · 1.5756 < 79.73`. -/
theorem tiearea_constants :
    (14 : ℝ) / (2 * Real.sqrt 2 * Real.pi) < 1.57556 ∧ (50.6 : ℝ) * 1.5756 < 79.73 := by
  obtain ⟨h1, _⟩ := sqrt_two_bounds
  have hpi := Real.pi_gt_d6
  refine ⟨?_, by norm_num⟩
  have hpos : 0 < 2 * Real.sqrt 2 * Real.pi := by positivity
  rw [div_lt_iff₀ hpos]
  nlinarith

/-- Proposition `prop:jweighted` and Theorem `thm:endlinear`. -/
theorem endlinear_constants :
    (1.828 : ℝ) * 20.113 < 36.77 ∧
    (77 : ℝ) * 0.6034 * 36.77 * 1.1825 < 2020.2 ∧
    (11.87 : ℝ) + 2020.2 * 1.9414 < 3934 ∧ (2020.2 : ℝ) * 120.09 < 242610 ∧
    (7 : ℝ) * 1.9811 * 1.0913 < 15.134 ∧ (25.592 : ℝ) * 38.4385 < 983.72 ∧
    (1.828 : ℝ) * 3934 < 7192 ∧ 128.07 + (1.828 : ℝ) * 242610 < 443620 ∧
    (1.828 : ℝ) * 43.59 < 79.69 ∧ (1.828 : ℝ) * 14888 < 27216 ∧
    (1 : ℝ) / 0.515099 < 1.9414 ∧ 3 * Real.sqrt 2 * 18.12 < 76.877 ∧
    6 * (2.064 : ℝ) < 12.39 ∧ 6 * Real.sqrt 2 < 8.486 := by
  obtain ⟨_, h2⟩ := sqrt_two_bounds
  refine ⟨by norm_num, by norm_num, by norm_num, by norm_num, by norm_num, by norm_num,
    by norm_num, by norm_num, by norm_num, by norm_num, by norm_num, by nlinarith,
    by norm_num, by nlinarith⟩

end YangMills.Averaged
