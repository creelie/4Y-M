import Mathlib

/-!
# The thermal hypothesis at one scale: the finite and algebraic parts

The paper reduces hypothesis (Z) to an entropy bound at one physical scale
and a Casimir condition on the halving of the spatial torus. The steps of
that reduction which are finite arithmetic or elementary inequalities are
checked here.

* `exp_mul_le_chord`, `weighted_exp_le_cosh` : the convexity bound
  `e^{xX} ≤ ((1+X)/2)e^x + ((1-X)/2)e^{-x}` and `c_0 ≤ cosh x` for a centred
  variable with values in `[-1,1]` (closed polymers, part (i));
* `one_sub_div_le` : `1 - e^{-x}/c_0 ≤ e^x - 1` when `1 ≤ c_0 ≤ cosh x`;
* `subtree_three`, `one_add_two_div_pow_le`, `choose_le_exp_pow`,
  `subtree_count_bound` : the count of rooted subtrees,
  `t_k = Δ/(k-1) · C((Δ-1)k, k-2) ≤ (eΔ)^{k-1}` for `k ≥ 4` (part (iii));
* `kp_v`, `kp_constant` : `v = 5/6` and `(6/(5e)) v^6/(1-v) ≤ 0.89`;
* `su2_range`, `su2_range_wilson` : `24 e² (e^{4β} - 1) ≤ 1` for
  `β ≤ 1.405 × 10⁻³`, that is `β_W ≤ 5.62 × 10⁻³`;
* `wrap_constant` : `(6M + 3/2) L³ ≤ 6 L³ (M+1)` for `M ≥ 4`;
* `exp_sub_one_le_two_mul` : `e^x - 1 ≤ 2x` on `[0,1]` (no fixed rate);
* `sum_pow_add_le` : `Σ λ_i^{a+b} ≤ (Σ λ_i^a)(Σ λ_i^b)` for `λ_i ≥ 0`
  (cutting a side);
* `log_one_add_le_of_prod` : `log(1+D) ≤ Mκ + Σ D_s(M')D_s(M'')` from the
  product inequality (halving);
* `exp_sub_one_le_contraction`, `contraction` : the contraction lemma with
  the constants `0.0943`, `1.049`, `0.84`, `1.05`;
* `geom_sum_le`, `weighted_geom_split`, `onescale_constants` : the geometric
  sums of the proof of (Z) from one scale and the constants `6.6`, `6.7`;
* `zzero` : `e^{c - Δ s_0/2 + Δ/(2N)} ≤ 1/10` for large `N` when `s_0 > (2/Δ)(c + log 10)`;
* `chain_corr_sq`, `chain_corr_gt_half` : the pinned Ising chain with
  weights `8, 1/4, 2` has `corr(s_1,s_3)² = 1143318969/2580030436 > 1/4`.
-/

open Finset Real

namespace YangMills.OneScale

/-- Convexity of the exponential on `[-1,1]`. -/
theorem exp_mul_le_chord (x X : ℝ) (h1 : -1 ≤ X) (h2 : X ≤ 1) :
    exp (x * X) ≤ (1 + X) / 2 * exp x + (1 - X) / 2 * exp (-x) := by
  have ha : 0 ≤ (1 + X) / 2 := by linarith
  have hb : 0 ≤ (1 - X) / 2 := by linarith
  have hab : (1 + X) / 2 + (1 - X) / 2 = 1 := by ring
  have h := convexOn_exp.2 (Set.mem_univ x) (Set.mem_univ (-x)) ha hb hab
  simp only [smul_eq_mul] at h
  have hx : (1 + X) / 2 * x + (1 - X) / 2 * (-x) = x * X := by ring
  rwa [hx] at h

/-- A centred average of `e^{xX}` over values `X ∈ [-1,1]` is at most `cosh x`. -/
theorem weighted_exp_le_cosh {ι : Type*} (s : Finset ι) (w X : ι → ℝ) (x : ℝ)
    (hw : ∀ i ∈ s, 0 ≤ w i) (hsum : ∑ i ∈ s, w i = 1)
    (hmean : ∑ i ∈ s, w i * X i = 0) (hX : ∀ i ∈ s, -1 ≤ X i ∧ X i ≤ 1) :
    ∑ i ∈ s, w i * exp (x * X i) ≤ cosh x := by
  calc ∑ i ∈ s, w i * exp (x * X i)
      ≤ ∑ i ∈ s, w i * ((1 + X i) / 2 * exp x + (1 - X i) / 2 * exp (-x)) := by
        apply sum_le_sum; intro i hi
        exact mul_le_mul_of_nonneg_left
          (exp_mul_le_chord x (X i) (hX i hi).1 (hX i hi).2) (hw i hi)
    _ = ∑ i ∈ s, (w i * ((exp x + exp (-x)) / 2)
          + (w i * X i) * ((exp x - exp (-x)) / 2)) := by
        apply sum_congr rfl; intro i _; ring
    _ = (∑ i ∈ s, w i) * ((exp x + exp (-x)) / 2)
          + (∑ i ∈ s, w i * X i) * ((exp x - exp (-x)) / 2) := by
        rw [sum_add_distrib, sum_mul, sum_mul]
    _ = cosh x := by rw [hsum, hmean, cosh_eq]; ring

/-- The negative side of the activity bound: `1 - e^{-x}/c_0 ≤ e^x - 1`. -/
theorem one_sub_div_le (x c : ℝ) (hx : 0 ≤ x) (hc : 0 < c) (hcc : c ≤ cosh x) :
    1 - exp (-x) / c ≤ exp x - 1 := by
  have h1 : exp (-x) / cosh x ≤ exp (-x) / c :=
    div_le_div_of_nonneg_left (exp_pos _).le hc hcc
  suffices h : 1 - exp (-x) / cosh x ≤ exp x - 1 by linarith
  have hv0 : 0 < exp x := exp_pos x
  have hv1 : 1 ≤ exp x := one_le_exp hx
  have hcosh : exp (-x) / cosh x = 2 / (exp x ^ 2 + 1) := by
    rw [cosh_eq, exp_neg]; field_simp; try ring
  rw [hcosh]
  have hpos : 0 < exp x ^ 2 + 1 := by positivity
  have hkey : 2 / (exp x ^ 2 + 1) ≥ 2 - exp x := by
    rw [ge_iff_le, le_div_iff₀ hpos]
    nlinarith [mul_nonneg hv0.le (sq_nonneg (exp x - 1))]
  linarith

/-- `t_3`: the Lagrange formula at `k = 3` agrees with the direct count. -/
theorem subtree_three (Δ : ℝ) : Δ / 2 * (3 * (Δ - 1)) = Δ * (Δ - 1) / 2 + Δ * (Δ - 1) := by
  ring

theorem one_add_two_div_pow_le (n : ℕ) (hn : 0 < n) : (1 + 2 / (n : ℝ)) ^ n ≤ exp 2 := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have h : 1 + 2 / (n : ℝ) ≤ exp (2 / n) := by
    have := add_one_le_exp (2 / (n : ℝ)); linarith
  calc (1 + 2 / (n : ℝ)) ^ n ≤ exp (2 / n) ^ n := pow_le_pow_left₀ (by positivity) h n
    _ = exp 2 := by rw [← exp_nat_mul]; congr 1; field_simp

/-- `C(n,m) ≤ (en/m)^m` for `m ≥ 1`. -/
theorem choose_le_exp_pow (n m : ℕ) (hm : 0 < m) :
    (n.choose m : ℝ) ≤ (exp 1 * n / m) ^ m := by
  have hm' : (0 : ℝ) < m := by exact_mod_cast hm
  have hfac : (0 : ℝ) < (m.factorial : ℝ) := by exact_mod_cast m.factorial_pos
  have h1 : (n.choose m : ℝ) ≤ (n : ℝ) ^ m / m.factorial := Nat.choose_le_pow_div m n
  have h2 : (m : ℝ) ^ m / m.factorial ≤ exp m := pow_div_factorial_le_exp (m : ℝ) hm'.le m
  have h3 : exp (m : ℝ) = exp 1 ^ m := by rw [← exp_nat_mul, mul_one]
  have h4 : (1 : ℝ) / m.factorial ≤ exp 1 ^ m / (m : ℝ) ^ m := by
    rw [h3] at h2
    rw [div_le_div_iff₀ hfac (by positivity)]
    rw [div_le_iff₀ hfac] at h2
    linarith
  calc (n.choose m : ℝ) ≤ (n : ℝ) ^ m / m.factorial := h1
    _ = (n : ℝ) ^ m * (1 / m.factorial) := by ring
    _ ≤ (n : ℝ) ^ m * (exp 1 ^ m / (m : ℝ) ^ m) :=
        mul_le_mul_of_nonneg_left h4 (by positivity)
    _ = (exp 1 * n / m) ^ m := by rw [div_pow, mul_pow]; ring

/-- The count of subtrees with `k ≥ 4` vertices through the root,
`Δ/(k-1) · C((Δ-1)k, k-2)`, is at most `(eΔ)^{k-1}`. -/
theorem subtree_count_bound (Δ k : ℕ) (hΔ : 1 ≤ Δ) (hk : 4 ≤ k) :
    (Δ : ℝ) / ((k : ℝ) - 1) * (((Δ - 1) * k).choose (k - 2) : ℝ)
      ≤ (exp 1 * Δ) ^ (k - 1) := by
  obtain ⟨m, rfl⟩ : ∃ m, k = m + 2 := ⟨k - 2, by omega⟩
  obtain ⟨δ, rfl⟩ : ∃ δ, Δ = δ + 1 := ⟨Δ - 1, by omega⟩
  have hm : 2 ≤ m := by omega
  have hm0 : (0 : ℝ) < m := by exact_mod_cast (show 0 < m by omega)
  simp only [add_tsub_cancel_right]
  have hc := choose_le_exp_pow (δ * (m + 2)) m (by omega)
  have hpow := one_add_two_div_pow_le m (by omega)
  have he : 0 < exp 1 := exp_pos 1
  have hkey : (exp 1 * ((δ * (m + 2) : ℕ) : ℝ) / m) ^ m
      = (exp 1 * δ) ^ m * (1 + 2 / (m : ℝ)) ^ m := by
    rw [← mul_pow]; congr 1; push_cast; field_simp; try ring
  have hA : ((δ * (m + 2)).choose m : ℝ) ≤ (exp 1 * δ) ^ m * exp 2 := by
    calc ((δ * (m + 2)).choose m : ℝ) ≤ _ := hc
      _ = (exp 1 * δ) ^ m * (1 + 2 / (m : ℝ)) ^ m := hkey
      _ ≤ (exp 1 * δ) ^ m * exp 2 := mul_le_mul_of_nonneg_left hpow (by positivity)
  have hcast : ((m + 2 : ℕ) : ℝ) - 1 = (m : ℝ) + 1 := by push_cast; ring
  rw [hcast]
  have hpos : (0 : ℝ) < (m : ℝ) + 1 := by positivity
  have hδpow : (exp 1 * δ) ^ m ≤ (exp 1 * ((δ + 1 : ℕ) : ℝ)) ^ m := by
    apply pow_le_pow_left₀ (by positivity); push_cast; nlinarith
  have he2 : exp 2 = exp 1 * exp 1 := by rw [← exp_add]; norm_num
  have hratio : exp 1 / ((m : ℝ) + 1) ≤ 1 := by
    rw [div_le_one hpos]
    have := exp_one_lt_d9
    have : (2 : ℝ) ≤ m := by exact_mod_cast hm
    linarith
  calc ((δ + 1 : ℕ) : ℝ) / ((m : ℝ) + 1) * ((δ * (m + 2)).choose m : ℝ)
      ≤ ((δ + 1 : ℕ) : ℝ) / ((m : ℝ) + 1) * ((exp 1 * δ) ^ m * exp 2) :=
        mul_le_mul_of_nonneg_left hA (by positivity)
    _ ≤ ((δ + 1 : ℕ) : ℝ) / ((m : ℝ) + 1) * ((exp 1 * ((δ + 1 : ℕ) : ℝ)) ^ m * exp 2) :=
        mul_le_mul_of_nonneg_left
          (mul_le_mul_of_nonneg_right hδpow (exp_pos 2).le) (by positivity)
    _ = (exp 1 * ((δ + 1 : ℕ) : ℝ)) ^ (m + 1) * (exp 1 / ((m : ℝ) + 1)) := by
        rw [he2, pow_succ]; ring
    _ ≤ (exp 1 * ((δ + 1 : ℕ) : ℝ)) ^ (m + 1) * 1 :=
        mul_le_mul_of_nonneg_left hratio (by positivity)
    _ = (exp 1 * ((δ + 1 : ℕ) : ℝ)) ^ (m + 2 - 1) := by rw [mul_one]; rfl

/-- `v = 20e · ε' e^{1+d} = 5/6` when `e^d = 1/(24 e² ε')`. -/
theorem kp_v : 20 * exp 1 * (1 / (24 * exp 1)) = (5 / 6 : ℝ) := by
  have := exp_pos 1; field_simp; ring

/-- The Kotecký–Preiss sum `(6/(5e)) v^6/(1-v) ≤ 0.89 ≤ 1` at `v = 5/6`. -/
theorem kp_constant : 6 / (5 * exp 1) * ((5 / 6 : ℝ) ^ 6 / (1 - 5 / 6)) ≤ 0.89 := by
  have h := exp_one_gt_d9
  have hpos : 0 < 5 * exp 1 := by positivity
  rw [div_mul_eq_mul_div, div_le_iff₀ hpos]
  norm_num
  linarith

/-- The strong-coupling range for `SU(2)`: `24 e² (e^{4β} - 1) ≤ 1` for `β ≤ 1.405 × 10⁻³`. -/
theorem su2_range (β : ℝ) (h0 : 0 ≤ β) (hβ : β ≤ 1.405e-3) :
    24 * exp 1 ^ 2 * (exp (4 * β) - 1) ≤ 1 := by
  have he := exp_one_lt_d9
  have he0 := exp_pos 1
  have h1 : exp (4 * β) ≤ exp (0.00562) := exp_le_exp.mpr (by linarith)
  have h2 := exp_bound' (x := 0.00562) (by norm_num) (by norm_num) (n := 3) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at h2
  norm_num at h2
  have hsq : exp 1 ^ 2 ≤ 2.7182818286 ^ 2 := pow_le_pow_left₀ he0.le he.le 2
  have hA : exp (4 * β) - 1 ≤ 0.00563584 := by linarith
  have hB : 0 ≤ exp (4 * β) - 1 := by linarith [one_le_exp (by linarith : 0 ≤ 4 * β)]
  calc 24 * exp 1 ^ 2 * (exp (4 * β) - 1) ≤ 24 * 2.7182818286 ^ 2 * 0.00563584 := by
        apply mul_le_mul (mul_le_mul_of_nonneg_left hsq (by norm_num)) hA hB (by positivity)
    _ ≤ 1 := by norm_num

/-- The same range in the Wilson coupling `β_W = 4β ≤ 5.62 × 10⁻³`. -/
theorem su2_range_wilson (βW : ℝ) (h0 : 0 ≤ βW) (hβ : βW ≤ 5.62e-3) :
    24 * exp 1 ^ 2 * (exp βW - 1) ≤ 1 := by
  have := su2_range (βW / 4) (by linarith) (by linarith)
  rwa [show 4 * (βW / 4) = βW by ring] at this

/-- `(6M + 3/2) L³ ≤ 6 L³ (M + 1)` for `M ≥ 4`. -/
theorem wrap_constant (L M : ℝ) (hL : 0 ≤ L) (_hM : 4 ≤ M) :
    (6 * M + 3 / 2) * L ^ 3 ≤ 6 * L ^ 3 * (M + 1) := by
  have : 0 ≤ L ^ 3 := by positivity
  nlinarith

/-- `e^x - 1 ≤ 2x` on `[0,1]` (no fixed rate along a tuning). -/
theorem exp_sub_one_le_two_mul (x : ℝ) (h0 : 0 ≤ x) (h1 : x ≤ 1) : exp x - 1 ≤ 2 * x := by
  have h := exp_bound' h0 h1 (n := 1) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at h
  norm_num at h
  linarith

/-- Cutting a side: `Σ λ_i^{a+b} ≤ (Σ λ_i^a)(Σ λ_i^b)` for nonnegative `λ_i`. -/
theorem sum_pow_add_le {ι : Type*} (s : Finset ι) (lam : ι → ℝ)
    (h : ∀ i ∈ s, 0 ≤ lam i) (a b : ℕ) :
    ∑ i ∈ s, lam i ^ (a + b) ≤ (∑ i ∈ s, lam i ^ a) * (∑ i ∈ s, lam i ^ b) := by
  rw [sum_mul]
  apply sum_le_sum
  intro i hi
  rw [pow_add]
  apply mul_le_mul_of_nonneg_left _ (pow_nonneg (h i hi) a)
  exact single_le_sum (fun j hj => pow_nonneg (h j hj) b) hi

/-- Halving: from `1 + D ≤ e^{c} ∏ (1 + D_s)` follows `log(1 + D) ≤ c + Σ D_s`. -/
theorem log_one_add_le_of_prod {ι : Type*} (s : Finset ι) (D c : ℝ) (Ds : ι → ℝ)
    (hD : 0 ≤ D) (hDs : ∀ i ∈ s, 0 ≤ Ds i)
    (h : 1 + D ≤ exp c * ∏ i ∈ s, (1 + Ds i)) :
    log (1 + D) ≤ c + ∑ i ∈ s, Ds i := by
  have hprod : 0 < ∏ i ∈ s, (1 + Ds i) := prod_pos (fun i hi => by linarith [hDs i hi])
  calc log (1 + D) ≤ log (exp c * ∏ i ∈ s, (1 + Ds i)) := log_le_log (by linarith) h
    _ = c + ∑ i ∈ s, log (1 + Ds i) := by
        rw [log_mul (exp_pos c).ne' hprod.ne', log_exp,
          log_prod (fun i hi => by linarith [hDs i hi])]
    _ ≤ c + ∑ i ∈ s, Ds i := by
        gcongr with i hi
        have := log_le_sub_one_of_pos (show 0 < 1 + Ds i by linarith [hDs i hi])
        linarith

/-- `e^y - 1 ≤ 1.049 y` for `0 ≤ y ≤ 0.0943`. -/
theorem exp_sub_one_le_contraction (y : ℝ) (h0 : 0 ≤ y) (h1 : y ≤ 0.0943) :
    exp y - 1 ≤ 1.049 * y := by
  have h := exp_bound' h0 (by linarith) (n := 4) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at h
  norm_num at h
  have h2 : y ^ 2 ≤ 0.0943 * y := by nlinarith
  have h3 : y ^ 3 ≤ 0.0943 ^ 2 * y := by nlinarith [pow_nonneg h0 2]
  have h4 : y ^ 4 ≤ 0.0943 ^ 3 * y := by nlinarith [pow_nonneg h0 3]
  nlinarith

/-- The contraction lemma. -/
theorem contraction (d η : ℝ) (hd0 : 0 ≤ d) (hd : d ≤ 1 / 10) (hη0 : 0 ≤ η) (hη : η ≤ 1 / 70) :
    exp (η + 8 * d ^ 2) - 1 ≤ 1 / 10 ∧
      exp (η + 8 * d ^ 2) - 1 ≤ 0.84 * d + 1.05 * η := by
  have hd2 : 8 * d ^ 2 ≤ 0.8 * d := by nlinarith
  have hy0 : 0 ≤ η + 8 * d ^ 2 := by positivity
  have hy1 : η + 8 * d ^ 2 ≤ 0.0943 := by nlinarith
  have h := exp_sub_one_le_contraction _ hy0 hy1
  constructor <;> nlinarith

/-- `Σ_{j<n} 0.84^j ≤ 1/0.16`. -/
theorem geom_sum_le (n : ℕ) : ∑ j ∈ range n, (0.84 : ℝ) ^ j ≤ 1 / 0.16 := by
  rw [geom_sum_eq (by norm_num) n]
  have h0 : 0 ≤ (0.84 : ℝ) ^ n := by positivity
  have hid : ((0.84 : ℝ) ^ n - 1) / (0.84 - 1) = (1 - (0.84 : ℝ) ^ n) / 0.16 := by ring
  rw [hid]
  exact div_le_div_of_nonneg_right (by linarith) (by norm_num)

/-- The split of `Σ_{j<J} 0.84^j η_j` at `j_0`: the terms before `j_0` are at most `ε`,
all terms at most `1/70`. -/
theorem weighted_geom_split (J j0 : ℕ) (η : ℕ → ℝ) (ε : ℝ) (hε : 0 ≤ ε)
    (_hη0 : ∀ j, 0 ≤ η j) (hsmall : ∀ j < j0, η j ≤ ε) (hall : ∀ j, η j ≤ 1 / 70) :
    ∑ j ∈ range J, (0.84 : ℝ) ^ j * η j ≤ ε / 0.16 + (0.84 : ℝ) ^ j0 / (70 * 0.16) := by
  rcases le_or_gt j0 J with hJ | hJ
  · rw [← sum_range_add_sum_Ico _ hJ]
    have hA : ∑ j ∈ range j0, (0.84 : ℝ) ^ j * η j ≤ ε / 0.16 := by
      calc ∑ j ∈ range j0, (0.84 : ℝ) ^ j * η j ≤ ∑ j ∈ range j0, (0.84 : ℝ) ^ j * ε := by
            apply sum_le_sum; intro j hj
            exact mul_le_mul_of_nonneg_left (hsmall j (mem_range.mp hj)) (by positivity)
        _ = (∑ j ∈ range j0, (0.84 : ℝ) ^ j) * ε := by rw [sum_mul]
        _ ≤ 1 / 0.16 * ε := mul_le_mul_of_nonneg_right (geom_sum_le j0) hε
        _ = ε / 0.16 := by ring
    have hB : ∑ j ∈ Ico j0 J, (0.84 : ℝ) ^ j * η j ≤ (0.84 : ℝ) ^ j0 / (70 * 0.16) := by
      calc ∑ j ∈ Ico j0 J, (0.84 : ℝ) ^ j * η j
          ≤ ∑ j ∈ Ico j0 J, (0.84 : ℝ) ^ j * (1 / 70) := by
            apply sum_le_sum; intro j _
            exact mul_le_mul_of_nonneg_left (hall j) (by positivity)
        _ = (0.84 : ℝ) ^ j0 * (∑ k ∈ range (J - j0), (0.84 : ℝ) ^ k) * (1 / 70) := by
            rw [sum_Ico_eq_sum_range, ← sum_mul, mul_sum]
            congr 1; apply sum_congr rfl; intro k _; rw [pow_add]
        _ ≤ (0.84 : ℝ) ^ j0 * (1 / 0.16) * (1 / 70) := by
            apply mul_le_mul_of_nonneg_right _ (by norm_num)
            exact mul_le_mul_of_nonneg_left (geom_sum_le _) (by positivity)
        _ = (0.84 : ℝ) ^ j0 / (70 * 0.16) := by ring
    linarith
  · have hA : ∑ j ∈ range J, (0.84 : ℝ) ^ j * η j ≤ ε / 0.16 := by
      calc ∑ j ∈ range J, (0.84 : ℝ) ^ j * η j ≤ ∑ j ∈ range J, (0.84 : ℝ) ^ j * ε := by
            apply sum_le_sum; intro j hj
            exact mul_le_mul_of_nonneg_left
              (hsmall j (lt_trans (mem_range.mp hj) hJ)) (by positivity)
        _ = (∑ j ∈ range J, (0.84 : ℝ) ^ j) * ε := by rw [sum_mul]
        _ ≤ 1 / 0.16 * ε := mul_le_mul_of_nonneg_right (geom_sum_le J) hε
        _ = ε / 0.16 := by ring
    have : 0 ≤ (0.84 : ℝ) ^ j0 / (70 * 0.16) := by positivity
    linarith

/-- The constants of the proof of (Z) from one scale: `1.05/0.16 ≤ 6.6`,
`1.05/(70 · 0.16) ≤ 1/10`, and `6.6 ε + ε/10 ≤ 6.7 ε`. -/
theorem onescale_constants (ε : ℝ) (_hε : 0 ≤ ε) :
    (1.05 : ℝ) / 0.16 ≤ 6.6 ∧ (1.05 : ℝ) / (70 * 0.16) ≤ 1 / 10 ∧
      6.6 * ε + ε / 10 ≤ 6.7 * ε := by
  refine ⟨by norm_num, by norm_num, by linarith⟩

/-- `(Z_0)` from a gap and a pressure bound: with sides in `[⌊Ns_0⌋, 2Ns_0+1]` the
exponent is `c - Δ s_0/2 + Δ/(2N)`, which is at most `-log 10` once
`N ≥ Δ/(Δ s_0 - 2c - 2 log 10)`, the denominator being positive when
`s_0 > (2/Δ)(c + log 10)`. -/
theorem zzero (c Δ s0 N : ℝ) (hΔ : 0 < Δ) (hs : 2 / Δ * (c + log 10) < s0)
    (hN : Δ / (Δ * s0 - 2 * c - 2 * log 10) ≤ N) :
    exp (c - Δ * s0 / 2 + Δ / (2 * N)) ≤ 1 / 10 := by
  have h1 : 2 * (c + log 10) < Δ * s0 := by
    have := mul_lt_mul_of_pos_left hs hΔ
    rwa [← mul_assoc, mul_div_cancel₀ _ hΔ.ne'] at this
  have hg : 0 < Δ * s0 - 2 * c - 2 * log 10 := by linarith
  have hN0 : 0 < N := lt_of_lt_of_le (div_pos hΔ hg) hN
  have h2 : Δ / (2 * N) ≤ (Δ * s0 - 2 * c - 2 * log 10) / 2 := by
    rw [div_le_div_iff₀ (by positivity) (by norm_num)]
    rw [div_le_iff₀ hg] at hN
    nlinarith
  have h3 : c - Δ * s0 / 2 + Δ / (2 * N) ≤ -log 10 := by linarith
  calc exp (c - Δ * s0 / 2 + Δ / (2 * N)) ≤ exp (-log 10) := exp_le_exp.mpr h3
    _ = 1 / 10 := by rw [exp_neg, exp_log (by norm_num)]; norm_num

/-! ### The pinned Ising chain -/

/-- The weights `w(+,+) = 8`, `w(+,-) = 1/4`, `w(-,-) = 2`; `true` is `+`. -/
def w : Bool → Bool → ℚ
  | true, true => 8
  | true, false => 1 / 4
  | false, true => 1 / 4
  | false, false => 2

def sgn : Bool → ℚ
  | true => 1
  | false => -1

/-- The partition function of five spins pinned to `-` at both ends. -/
def Z : ℚ := ∑ a : Bool, ∑ b : Bool, ∑ c : Bool, ∑ d : Bool, ∑ e : Bool,
  w false a * w a b * w b c * w c d * w d e * w e false

def S1 : ℚ := ∑ a : Bool, ∑ b : Bool, ∑ c : Bool, ∑ d : Bool, ∑ e : Bool,
  w false a * w a b * w b c * w c d * w d e * w e false * sgn b

def S3 : ℚ := ∑ a : Bool, ∑ b : Bool, ∑ c : Bool, ∑ d : Bool, ∑ e : Bool,
  w false a * w a b * w b c * w c d * w d e * w e false * sgn d

def S13 : ℚ := ∑ a : Bool, ∑ b : Bool, ∑ c : Bool, ∑ d : Bool, ∑ e : Bool,
  w false a * w a b * w b c * w c d * w d e * w e false * (sgn b * sgn d)

/-- `corr(s_1, s_3)² = 1143318969/2580030436` exactly. -/
theorem chain_corr_sq :
    (S13 * Z - S1 * S3) ^ 2 / ((Z ^ 2 - S1 ^ 2) * (Z ^ 2 - S3 ^ 2))
      = 1143318969 / 2580030436 := by
  simp only [Z, S1, S3, S13, Fintype.sum_bool, w, sgn]
  norm_num

/-- So the correlation exceeds `1/2`. -/
theorem chain_corr_gt_half : (1 / 2 : ℚ) ^ 2 < 1143318969 / 2580030436 := by norm_num

/-! ### Doubling defects, symmetrised rings and the finite-size bound -/

/-- Lemma doubling: for `0 ≤ r_i ≤ 1`, `1 ≤ 𝒵(2M) ≤ 𝒵(M)` as finite sums. -/
theorem doubling_sums {ι : Type*} (s : Finset ι) (r : ι → ℝ) (M : ℕ)
    (h0 : ∀ i ∈ s, 0 ≤ r i) (h1 : ∀ i ∈ s, r i ≤ 1) :
    1 ≤ 1 + ∑ i ∈ s, r i ^ (2 * M) ∧
      1 + ∑ i ∈ s, r i ^ (2 * M) ≤ 1 + ∑ i ∈ s, r i ^ M := by
  refine ⟨?_, ?_⟩
  · have : 0 ≤ ∑ i ∈ s, r i ^ (2 * M) :=
      Finset.sum_nonneg fun i hi => pow_nonneg (h0 i hi) _
    linarith
  · have : ∑ i ∈ s, r i ^ (2 * M) ≤ ∑ i ∈ s, r i ^ M :=
      Finset.sum_le_sum fun i hi => pow_le_pow_of_le_one (h0 i hi) (h1 i hi) (by omega)
    linarith

/-- `𝒵(M) ≤ Z(M)²/Z(2M) ≤ 𝒵(M)²` from `1 ≤ 𝒵(2M) ≤ 𝒵(M)`. -/
theorem doubling_sandwich {a b : ℝ} (hb : 1 ≤ b) (hba : b ≤ a) :
    a ≤ a ^ 2 / b ∧ a ^ 2 / b ≤ a ^ 2 := by
  have hbpos : 0 < b := by linarith
  constructor
  · rw [le_div_iff₀ hbpos]; nlinarith
  · rw [div_le_iff₀ hbpos]; nlinarith

/-- The expansion of `16F(n) − F(2n)` into four doubling brackets. -/
theorem four_direction_identity (F0 F1 F2 F3 F4 : ℝ) :
    16 * F0 - F4 =
      8 * (2 * F0 - F1) + 4 * (2 * F1 - F2) + 2 * (2 * F2 - F3) + (2 * F3 - F4) := by
  ring

/-- `(Z₀)` from the doubling defect: `1 + D ≤ e^ψ` and `ψ ≤ log(11/10)` give `D ≤ 1/10`. -/
theorem zzero_of_psi {D ψ : ℝ} (hD : 1 + D ≤ exp ψ) (hψ : ψ ≤ log (11 / 10)) :
    D ≤ 1 / 10 := by
  have h := Real.exp_le_exp.mpr hψ
  rw [Real.exp_log (by norm_num)] at h
  linarith

/-- `log(11/10) ≥ 1/11`. -/
theorem log_eleven_tenths : (1 / 11 : ℝ) ≤ log (11 / 10) := by
  have h := Real.one_sub_inv_le_log_of_pos (by norm_num : (0 : ℝ) < 11 / 10)
  have e : (1 : ℝ) - (11 / 10)⁻¹ = 1 / 11 := by norm_num
  rw [e] at h
  exact h

/-- Bernoulli: `k x ≤ (1+x)^k − 1` for `x ≥ 0` (the rings' `D₀`, and the Hilbert–Schmidt
bound for a product of `k` pairs). -/
theorem bernoulli_lower (x : ℝ) (hx : 0 ≤ x) (k : ℕ) : (k : ℝ) * x ≤ (1 + x) ^ k - 1 := by
  have := one_add_mul_le_pow (by linarith : (-2 : ℝ) ≤ x) k
  linarith

/-- For `0 ≤ r ≤ 1`: `log(1+r) − ½ log(1+r²) ≥ ½ log(1+r)`. -/
theorem ring_difference (r : ℝ) (h0 : 0 ≤ r) (h1 : r ≤ 1) :
    log (1 + r) / 2 ≤ log (1 + r) - log (1 + r ^ 2) / 2 := by
  have hr : r ^ 2 ≤ r := by nlinarith
  have := Real.log_le_log (by positivity : (0 : ℝ) < 1 + r ^ 2) (by linarith : 1 + r ^ 2 ≤ 1 + r)
  linarith

/-- `x ↦ x⁴ e^{−m x/2}` decreases on `x ≥ 8/m`. -/
theorem quartic_exp_decreasing {m x y : ℝ} (_hm : 0 < m) (hy : 0 < y) (h8 : 8 ≤ m * y)
    (hxy : y ≤ x) :
    x ^ 4 * exp (-(m * x) / 2) ≤ y ^ 4 * exp (-(m * y) / 2) := by
  have h1 : (x / y) ^ 4 ≤ exp (m * (x - y) / 2) := by
    have e1 : x / y = (x - y) / y + 1 := by field_simp; ring
    have h2 : ((x - y) / y + 1) ^ 4 ≤ (exp ((x - y) / y)) ^ 4 :=
      pow_le_pow_left₀ (by positivity) (Real.add_one_le_exp _) 4
    have h3 : (exp ((x - y) / y)) ^ 4 = exp (4 * ((x - y) / y)) := by
      rw [← Real.exp_nat_mul]; norm_num
    have h4 : 4 * ((x - y) / y) ≤ m * (x - y) / 2 := by
      have e : 4 * ((x - y) / y) = 4 * (x - y) / y := by ring
      rw [e, div_le_iff₀ hy]
      nlinarith [mul_nonneg (sub_nonneg.mpr hxy) (sub_nonneg.mpr h8)]
    calc (x / y) ^ 4 = ((x - y) / y + 1) ^ 4 := by rw [e1]
      _ ≤ (exp ((x - y) / y)) ^ 4 := h2
      _ = exp (4 * ((x - y) / y)) := h3
      _ ≤ exp (m * (x - y) / 2) := Real.exp_le_exp.mpr h4
  have e2 : x ^ 4 = (x / y) ^ 4 * y ^ 4 := by field_simp
  have e3 : m * (x - y) / 2 + -(m * x) / 2 = -(m * y) / 2 := by ring
  calc x ^ 4 * exp (-(m * x) / 2) = (x / y) ^ 4 * (y ^ 4 * exp (-(m * x) / 2)) := by
        rw [e2]; ring
    _ ≤ exp (m * (x - y) / 2) * (y ^ 4 * exp (-(m * x) / 2)) := by gcongr
    _ = y ^ 4 * exp (-(m * y) / 2) := by
        rw [mul_left_comm, ← Real.exp_add, e3]

/-- The prefactor of the finite-size bound: `e^{m/N} ≤ e` for `N ≥ m`. -/
theorem fsc_prefactor {m N K : ℝ} (hm : 0 < m) (hN : m ≤ N) (hK : 0 ≤ K) :
    6 * K * exp (m / N) ≤ 6 * exp 1 * K := by
  have hNpos : 0 < N := lt_of_lt_of_le hm hN
  have h1 : m / N ≤ 1 := by rw [div_le_one hNpos]; exact hN
  have h2 := Real.exp_le_exp.mpr h1
  nlinarith [Real.exp_pos (m / N)]

/-- Two slab sums `σ_* ≤ 3E` give the factor `6`, and `2σ_* + E = 7E`. -/
theorem fsc_factors {V s1 s2 E : ℝ} (hV : 0 ≤ V) (h1 : s1 ≤ 3 * E) (h2 : s2 ≤ 3 * E) :
    V * s1 + V * s2 ≤ 6 * V * E ∧ 2 * (3 * E) + E = 7 * E := by
  constructor
  · nlinarith
  · ring

/-- The volumes of the eight halves add up to the volume of the torus. -/
theorem halves_volume (a a' b b' c c' : ℝ) :
    (a + a') * (b + b') * (c + c') =
      a * b * c + a * b * c' + a * b' * c + a * b' * c' +
        a' * b * c + a' * b * c' + a' * b' * c + a' * b' * c' := by
  ring

/-- `X` uniform on `{−1, 0, 1}` and `Y = X² − 2/3`: `Y` is centred and uncorrelated with `X`,
while `Y` is a function of `X`. -/
theorem uncorrelated_example :
    ((-1 : ℚ) ^ 2 - 2 / 3) + ((0 : ℚ) ^ 2 - 2 / 3) + ((1 : ℚ) ^ 2 - 2 / 3) = 0 ∧
    (-1 : ℚ) * ((-1 : ℚ) ^ 2 - 2 / 3) + 0 * ((0 : ℚ) ^ 2 - 2 / 3)
      + 1 * ((1 : ℚ) ^ 2 - 2 / 3) = 0 := by
  norm_num

/-! ### The gap against the Casimir defect (Lemma `lem:dual`, Theorem `thm:gapcasimir`,
Proposition `prop:dichotomy`, Theorem `thm:Crate`). -/

/-- The single-cut defect: `log 𝒵_u(a) + log 𝒵_u(b) − log 𝒵_u(a+b) ≤ log 𝒵_u(a) + log 𝒵_u(b)`
because `𝒵 ≥ 1`; and the identity of `lem:dual`(i) in the form of logarithms. -/
theorem single_cut_defect {za zb zab : ℝ} (hab : 1 ≤ zab) :
    log za + log zb - log zab ≤ log za + log zb := by
  have := Real.log_nonneg hab
  linarith

/-- `lem:dual`(i): `M κ_p = Φ + log 𝒵_t(M) − Σ_s log 𝒵_s(M)` from `log Z = −M E₀ + log 𝒵`
for the torus and its halves (written for eight halves). -/
theorem dual_identity (M Et Zt : ℝ) (Es Zs : Fin 8 → ℝ) :
    M * (Et - ∑ i, Es i) =
      ((∑ i, (-M * Es i + Zs i)) - (-M * Et + Zt)) + Zt - ∑ i, Zs i := by
  have h : ∑ i, (-M * Es i + Zs i) = -(M * ∑ i, Es i) + ∑ i, Zs i := by
    rw [Finset.sum_add_distrib, Finset.mul_sum, ← Finset.sum_neg_distrib]
    simp only [neg_mul]
  rw [h]; ring

/-- `lem:dual`(iii): `15 X (m+1)/m' ≤ 30 (1 + 1/m') X` when `m + 1 ≤ 2 m' + 2`. -/
theorem dual_coefficient {m m' X : ℝ} (hm' : 0 < m') (hX : 0 ≤ X) (h : m + 1 ≤ 2 * m' + 2) :
    15 * X * ((m + 1) / m') ≤ 30 * (1 + 1 / m') * X := by
  have : (m + 1) / m' ≤ 2 * (1 + 1 / m') := by
    rw [div_le_iff₀ hm']; field_simp; nlinarith
  nlinarith

/-- The recursion of `thm:gapcasimir`: the gain of the gap over half the period against the
cost of the cut, under (K). -/
theorem gap_beats_casimir {Δ m N κ θ : ℝ} (hN : 0 < N)
    (hK : (m + 1) * κ ≤ (1 - θ) * Δ * m / N) :
    1.1 ^ 8 * exp (-(Δ * m) / (2 * N) + (m + 1) * κ / 2) ≤ 1.1 ^ 8 * exp (-(θ * Δ * m) / (2 * N)) := by
  have hpos : (0 : ℝ) < 1.1 ^ 8 := by positivity
  refine mul_le_mul_of_nonneg_left (Real.exp_le_exp.mpr ?_) hpos.le
  have h1 : (m + 1) * κ / 2 ≤ (1 - θ) * Δ * m / N / 2 := by linarith
  have h2 : -(Δ * m) / (2 * N) + (1 - θ) * Δ * m / N / 2 = -(θ * Δ * m) / (2 * N) := by
    field_simp; ring
  linarith

/-- The base of the induction: `1.1^8 e^{−x} ≤ 1/10` when `x ≥ log(10 · 1.1^8)`. -/
theorem base_tenth {x : ℝ} (hx : log (10 * 1.1 ^ 8) ≤ x) : 1.1 ^ 8 * exp (-x) ≤ 1 / 10 := by
  have hpos : (0 : ℝ) < 10 * 1.1 ^ 8 := by positivity
  have h1 : exp (-x) ≤ exp (-log (10 * 1.1 ^ 8)) := Real.exp_le_exp.mpr (by linarith)
  have h2 : exp (-log (10 * 1.1 ^ 8)) = 1 / (10 * 1.1 ^ 8) := by
    rw [Real.exp_neg, Real.exp_log hpos, one_div]
  calc 1.1 ^ 8 * exp (-x) ≤ 1.1 ^ 8 * (1 / (10 * 1.1 ^ 8)) := by rw [← h2]; gcongr
    _ = 1 / 10 := by field_simp

/-- The halving positions: for `M ∈ {m, m+1}`, `M − ⌊M/2⌋ ≥ m/2` and `⌊M/2⌋ ≤ (m+1)/2`. -/
theorem halving_positions (m : ℕ) :
    (m : ℝ) / 2 ≤ ((m : ℝ) - (m / 2 : ℕ)) ∧ (m / 2 : ℕ) ≤ ((m : ℝ) + 1) / 2 ∧
    (m : ℝ) / 2 ≤ (((m + 1 : ℕ) : ℝ) - ((m + 1) / 2 : ℕ)) ∧
    (((m + 1) / 2 : ℕ) : ℝ) ≤ ((m : ℝ) + 1) / 2 := by
  have h1 : (m / 2 : ℕ) * 2 ≤ m := Nat.div_mul_le_self m 2
  have h2 : ((m + 1) / 2 : ℕ) * 2 ≤ m + 1 := Nat.div_mul_le_self (m + 1) 2
  have h1' : ((m / 2 : ℕ) : ℝ) * 2 ≤ m := by exact_mod_cast h1
  have h2' : (((m + 1) / 2 : ℕ) : ℝ) * 2 ≤ (m : ℝ) + 1 := by exact_mod_cast h2
  push_cast
  refine ⟨by linarith, by linarith, by linarith, by linarith⟩

/-- `lem:convexalt`: if `x ≤ A(B e^{bx} − 1)` and `A B b e^{bx} ≤ 1/2`, then `x ≤ 2A(B−1)`. -/
theorem convex_alternative {A B b x : ℝ} (hA : 0 < A) (hB : 1 ≤ B) (_hb : 0 < b) (hx : 0 ≤ x)
    (h : x ≤ A * (B * exp (b * x) - 1)) (hstar : A * B * b * exp (b * x) ≤ 1 / 2) :
    x ≤ 2 * A * (B - 1) := by
  -- `e^{bx} − 1 ≤ b x e^{bx}` from `1 − e^{−bx} ≤ bx`
  have h1 : 1 - b * x ≤ exp (-(b * x)) := by linarith [Real.add_one_le_exp (-(b * x))]
  have hE : 0 < exp (b * x) := Real.exp_pos _
  have h2 : exp (b * x) - 1 ≤ b * x * exp (b * x) := by
    have : exp (-(b * x)) * exp (b * x) = 1 := by rw [← Real.exp_add]; simp
    nlinarith
  have hAB : 0 ≤ A * B := by positivity
  nlinarith [mul_le_mul_of_nonneg_left h2 hAB]

/-- The threshold `x_*` of `lem:convexalt`: `x ≤ log(1/(2AbB))/b` gives `A B b e^{bx} ≤ 1/2`. -/
theorem convex_threshold {A B b x : ℝ} (hA : 0 < A) (hB : 0 < B) (hb : 0 < b)
    (hx : x ≤ log (1 / (2 * A * b * B)) / b) : A * B * b * exp (b * x) ≤ 1 / 2 := by
  have hp : 0 < 1 / (2 * A * b * B) := by positivity
  have h1 : b * x ≤ log (1 / (2 * A * b * B)) := by
    rwa [le_div_iff₀ hb, mul_comm] at hx
  have h2 : exp (b * x) ≤ 1 / (2 * A * b * B) := by
    rw [← Real.exp_log hp]; exact Real.exp_le_exp.mpr h1
  calc A * B * b * exp (b * x) ≤ A * B * b * (1 / (2 * A * b * B)) := by gcongr
    _ = 1 / 2 := by field_simp

/-- `prop:dichotomy`: the positions `τ − ⌊τ/2⌋ ≥ (m−1)/4` and `⌊τ/2⌋ ≤ (m+2)/4` for
`τ ∈ {m', m'+1}`, `m' = ⌊m/2⌋`, and the coefficient `b_j ≤ 7.8` once `m' ≥ 50`, `m ≥ 2m'`. -/
theorem dichotomy_positions {m m' τ τ2 : ℝ} (hm'l : (m - 1) / 2 ≤ m') (hm'u : m' ≤ m / 2)
    (hτ : m' ≤ τ) (hτ' : τ ≤ m' + 1) (hτ2 : τ - τ2 ≥ τ / 2) (_hτ2' : τ2 ≤ τ / 2) :
    (m - 1) / 4 ≤ τ - τ2 ∧ τ2 ≤ (m + 2) / 4 := by
  constructor <;> linarith

theorem dichotomy_coefficient {m m' : ℝ} (hm' : 50 ≤ m') (hm : 2 * m' ≤ m) :
    30 * (m + 2) * (m' + 1) / (4 * (m + 1) * m') ≤ 7.8 := by
  have hm'pos : 0 < m' := by linarith
  have hmpos : 0 < m + 1 := by linarith
  rw [div_le_iff₀ (by positivity)]
  nlinarith

/-- `thm:Crate`: `(m+2)(m−4) ≤ (m+1)(m−2)`, and the constants `e^{0.8}` and `30 e^{0.8}(1 + 1/99)`. -/
theorem crate_arith {m : ℝ} (hm : 0 ≤ m) : (m + 2) * (m - 4) ≤ (m + 1) * (m - 2) := by nlinarith

theorem crate_constants : (8 : ℝ) * (1 / 10) = 0.8 ∧ (30 : ℝ) * (1 + 1 / 99) * 2.2256 < 67.5 := by
  constructor <;> norm_num

/-- The exponent of the dyadic bridge `cor:gapbridge`(ii): with `c_{j+1} = c_j/2` and
`N κ ≤ (1−θ)Δ`, `−Δ c/(2N) + (c/2) κ ≤ −θ Δ c/(2N)`. -/
theorem dyadic_exponent {Δ c N κ θ : ℝ} (hN : 0 < N) (hc : 0 ≤ c) (hK : N * κ ≤ (1 - θ) * Δ) :
    -(Δ * c) / (2 * N) + c / 2 * κ ≤ -(θ * Δ * c) / (2 * N) := by
  have h : c / 2 * κ ≤ c / 2 * ((1 - θ) * Δ / N) := by
    apply mul_le_mul_of_nonneg_left _ (by linarith)
    rw [le_div_iff₀ hN]; linarith
  have e : -(Δ * c) / (2 * N) + c / 2 * ((1 - θ) * Δ / N) = -(θ * Δ * c) / (2 * N) := by
    field_simp; ring
  linarith

/-- The product `15 > 1` of the two coefficients, which is why positivity alone does not
contract: `30 · (1/2) = 15`. -/
theorem positivity_product : (30 : ℝ) * (1 / 2) = 15 ∧ (1 : ℝ) < 15 := by norm_num

/-! ### The pinned two-state chain (Lemma `lem:pinnedchain`, Proposition `prop:fieldrings`). -/

/-- The partition function, one-spin and two-spin sums of the pinned chain, as polynomials in
`σ = sin α`, `γ = cos α` and `Q`, and the two identities `WZ − M² = 4Qσ²γ²(1−Q)²(Qγ²+σ²)²`
(a polynomial identity) and `Z² − M² = 4σ²γ²(1−Q)²(1+Q)(Qγ²+σ²)(Q²γ²+σ²)` (modulo `σ²+γ²=1`). -/
theorem pinned_chain_identities (σ γ Q : ℝ) (h : σ ^ 2 + γ ^ 2 = 1) :
    let Z := σ ^ 2 + Q ^ 3 * γ ^ 2
    let M := -Q ^ 3 * γ ^ 4 + (Q ^ 3 - 2 * Q ^ 2 - 2 * Q + 1) * γ ^ 2 * σ ^ 2 - σ ^ 4
    let W := Q ^ 3 * γ ^ 6 - 2 * Q ^ 3 * γ ^ 4 * σ ^ 2 + Q ^ 3 * γ ^ 2 * σ ^ 4 + 8 * Q ^ 2 * γ ^ 4 * σ ^ 2
      - 4 * Q ^ 2 * γ ^ 2 * σ ^ 4 - 4 * Q * γ ^ 4 * σ ^ 2 + 8 * Q * γ ^ 2 * σ ^ 4 + γ ^ 4 * σ ^ 2
      - 2 * γ ^ 2 * σ ^ 4 + σ ^ 6
    W * Z - M ^ 2 = 4 * Q * σ ^ 2 * γ ^ 2 * (1 - Q) ^ 2 * (Q * γ ^ 2 + σ ^ 2) ^ 2 ∧
    Z ^ 2 - M ^ 2 = 4 * σ ^ 2 * γ ^ 2 * (1 - Q) ^ 2 * (1 + Q) * (Q * γ ^ 2 + σ ^ 2) * (Q ^ 2 * γ ^ 2 + σ ^ 2) := by
  intro Z M W
  constructor
  · ring
  · linear_combination (-(σ ^ 2 + Q ^ 3 * γ ^ 2) ^ 2 * (γ ^ 2 + σ ^ 2 + 1)) * h

/-- The matrix entries: `Z = aᵀD(Q³)a`, `M = aᵀD(Q)SD(Q²)a`, `W = aᵀD(Q)SD(Q)SD(Q)a` with
`a = (σ, γ)`, `D(x) = diag(1, x)` and `S` the spin operator in the eigenbasis, expand to the
polynomials above. -/
theorem pinned_chain_entries (σ γ Q : ℝ) (h : σ ^ 2 + γ ^ 2 = 1) :
    let s2 := 2 * σ * γ            -- sin 2α
    let c2 := γ ^ 2 - σ ^ 2        -- cos 2α
    -- S = [[c2, -s2], [-s2, -c2]], D(x) a = (σ, x γ)
    let SDa := fun x : ℝ => (c2 * σ - s2 * (x * γ), -s2 * σ - c2 * (x * γ))
    (σ * σ + γ * (Q ^ 3 * γ) = σ ^ 2 + Q ^ 3 * γ ^ 2) ∧
    (σ * (SDa (Q ^ 2)).1 + Q * γ * (SDa (Q ^ 2)).2 =
      -Q ^ 3 * γ ^ 4 + (Q ^ 3 - 2 * Q ^ 2 - 2 * Q + 1) * γ ^ 2 * σ ^ 2 - σ ^ 4) ∧
    (σ * (c2 * (SDa Q).1 - s2 * (Q * (SDa Q).2)) + Q * γ * (-s2 * (SDa Q).1 - c2 * (Q * (SDa Q).2)) =
      Q ^ 3 * γ ^ 6 - 2 * Q ^ 3 * γ ^ 4 * σ ^ 2 + Q ^ 3 * γ ^ 2 * σ ^ 4 + 8 * Q ^ 2 * γ ^ 4 * σ ^ 2
      - 4 * Q ^ 2 * γ ^ 2 * σ ^ 4 - 4 * Q * γ ^ 4 * σ ^ 2 + 8 * Q * γ ^ 2 * σ ^ 4 + γ ^ 4 * σ ^ 2
      - 2 * γ ^ 2 * σ ^ 4 + σ ^ 6) := by
  intro s2 c2 SDa
  refine ⟨by ring, ?_, ?_⟩
  · simp only [SDa, s2, c2]; linear_combination (0 : ℝ) * h
  · simp only [SDa, s2, c2]; linear_combination (0 : ℝ) * h

/-- The correlation `C(t,Q) = Q(Q+t)/((1+Q)(Q²+t))` exceeds `1/2` exactly when `t < Q²`
(for `0 < Q < 1`, `t > 0`). -/
theorem pinned_chain_half {t Q : ℝ} (ht : 0 < t) (hQ : 0 < Q) (hQ1 : Q < 1) :
    1 / 2 < Q * (Q + t) / ((1 + Q) * (Q ^ 2 + t)) ↔ t < Q ^ 2 := by
  have hpos : 0 < (1 + Q) * (Q ^ 2 + t) := by positivity
  rw [lt_div_iff₀ hpos]
  constructor <;> intro h <;> nlinarith

/-- `t ↦ (Q+t)/(Q²+t)` is decreasing: its numerator derivative sign `Q² − Q < 0`; in the form
`t < t' → (Q+t')/(Q²+t') < (Q+t)/(Q²+t)`. -/
theorem pinned_chain_monotone {t t' Q : ℝ} (ht : 0 ≤ t) (htt : t < t') (hQ : 0 < Q) (hQ1 : Q < 1) :
    (Q + t') / (Q ^ 2 + t') < (Q + t) / (Q ^ 2 + t) := by
  have ht' : 0 < t' := lt_of_le_of_lt ht htt
  rw [div_lt_div_iff₀ (add_pos (by positivity) ht') (add_pos_of_pos_of_nonneg (by positivity) ht)]
  nlinarith [mul_pos hQ (sub_pos.mpr hQ1)]

/-- The scaling of `prop:fieldrings`: `(√(1+r²) − r)(√(1+r²) + r) = 1`, so `t_∞ < Q₀²` iff
`√(1+r²) + r > 1/Q₀`; and the numbers of the example, `(√401 − 20)² < 6.25·10⁻⁴ < 1/400`. -/
theorem field_rings_scaling (r : ℝ) :
    (Real.sqrt (1 + r ^ 2) - r) * (Real.sqrt (1 + r ^ 2) + r) = 1 := by
  have h : Real.sqrt (1 + r ^ 2) ^ 2 = 1 + r ^ 2 := Real.sq_sqrt (by positivity)
  nlinarith [h]

theorem field_rings_example :
    (Real.sqrt 401 - 20) ^ 2 < 6.25e-4 ∧ (6.25e-4 : ℝ) < (1 / 20) ^ 2 ∧ 148 / 0.9 ^ 4 < (227 : ℝ) - 1 := by
  have h1 : (20 : ℝ) < Real.sqrt 401 := by
    rw [Real.lt_sqrt (by norm_num)]; norm_num
  have h2 : Real.sqrt 401 < 20.025 := by
    rw [Real.sqrt_lt' (by norm_num)]; norm_num
  refine ⟨?_, by norm_num, by norm_num⟩
  nlinarith

/-! ### The Casimir defect of the free field (`ssec:freecasimir`) -/

/-- The transfer kernel of one mode, `lem:freevacuum`: with `c + 1/c = 2 + a`, the Gaussian
`e^{-bx²/2}` with `b = c − 1 − a/2` is reproduced, since the exponent `−1/2 + 1/(2c) − a/4`
equals `−b/2`; and then `b = (c − 1/c)/2 = sinh ω > 0`. -/
theorem free_mode_eigenfunction {a c : ℝ} (hc : 0 < c) (h : c + 1 / c = 2 + a) :
    -1 / 2 + 1 / (2 * c) - a / 4 = -(c - 1 - a / 2) / 2 ∧ c - 1 - a / 2 = (c - 1 / c) / 2 := by
  have hc' : c ≠ 0 := ne_of_gt hc
  constructor
  · field_simp
    have : a = c + 1 / c - 2 := by linarith
    subst this; field_simp; ring
  · have : a = c + 1 / c - 2 := by linarith
    subst this; ring

/-- `c + 1/c = 2 + a` is the quadratic `c² − (2+a)c + 1 = 0`. -/
theorem free_mode_quadratic {a c : ℝ} (hc : 0 < c) :
    c + 1 / c = 2 + a ↔ c ^ 2 - (2 + a) * c + 1 = 0 := by
  have hc' : c ≠ 0 := ne_of_gt hc
  constructor
  · intro h; field_simp at h; nlinarith [h]
  · intro h; field_simp; nlinarith [h]

/-- The prefactor of the Bessel series, `prop:cfunction`(ii):
`(4π)^{-1/2} · 8 (4π)^{-3/2} = 1/(2π²)`. -/
theorem casimir_prefactor :
    1 / Real.sqrt (4 * Real.pi) * (8 / Real.sqrt (4 * Real.pi) ^ 3) = 1 / (2 * Real.pi ^ 2) := by
  have hpi := Real.pi_pos
  have hs : Real.sqrt (4 * Real.pi) ^ 2 = 4 * Real.pi := Real.sq_sqrt (by positivity)
  have hs0 : Real.sqrt (4 * Real.pi) ≠ 0 := by positivity
  have h4 : Real.sqrt (4 * Real.pi) ^ 4 = 16 * Real.pi ^ 2 := by
    calc Real.sqrt (4 * Real.pi) ^ 4 = (Real.sqrt (4 * Real.pi) ^ 2) ^ 2 := by ring
      _ = (4 * Real.pi) ^ 2 := by rw [hs]
      _ = 16 * Real.pi ^ 2 := by ring
  rw [div_mul_div_comm, one_mul, ← pow_succ', h4]
  field_simp
  ring

/-- `prop:cfunction`(ii): each `w` contributes `16 − 1 = 15` times `|w|^{-4}` at `x = 0`, and
with `Z₄ > 9` and `π < 3.15`, `c(0) = 15 Z₄/(2π²) > 6.8`. -/
theorem casimir_zero_lower {Z : ℝ} (hZ : 9 < Z) : 6.8 < 15 * Z / (2 * Real.pi ^ 2) := by
  have hpi := Real.pi_pos
  have hlt : Real.pi < 3.15 := Real.pi_lt_d2
  have hp2 : Real.pi ^ 2 < 9.9225 := by nlinarith
  rw [lt_div_iff₀ (by positivity)]
  nlinarith

theorem casimir_sixteen_minus_one : (16 : ℝ) - 1 = 15 ∧ (1 : ℝ) - 8 * 2 = -15 := by norm_num

/-- Lagrange's identity in three dimensions, `prop:freegluons`: the curl form is
`|q|²|A|² − ⟨q,A⟩²` on real vectors. -/
theorem lagrange_identity_three (a b c u v w : ℝ) :
    (a * v - b * u) ^ 2 + (a * w - c * u) ^ 2 + (b * w - c * v) ^ 2
      = (a ^ 2 + b ^ 2 + c ^ 2) * (u ^ 2 + v ^ 2 + w ^ 2) - (a * u + b * v + c * w) ^ 2 := by
  ring

/-- `lem:freelimit`(ii): the derivative of `y²/(2(1+y))` is `(1 − (1+y)⁻²)/2`. -/
theorem psi_lower_derivative {y : ℝ} (hy : 0 ≤ y) :
    (2 * y * (1 + y) - y ^ 2) / (2 * (1 + y) ^ 2) = (1 - 1 / (1 + y) ^ 2) / 2 := by
  have h : (1 + y) ≠ 0 := by positivity
  field_simp
  ring

/-- `lem:freelimit`(ii) with `x = Lk`, `L ≥ m/3`, `t = m²τ`: `2t · y²/(2(1+y))` with `y = x/2t`
is `x²/(4t + 2x)`. -/
theorem psi_exponent {x t : ℝ} (ht : 0 < t) (hx : 0 ≤ x) :
    2 * t * ((x / (2 * t)) ^ 2 / (2 * (1 + x / (2 * t)))) = x ^ 2 / (4 * t + 2 * x) := by
  have h1 : (2 * t) ≠ 0 := by positivity
  have h2 : (4 * t + 2 * x) ≠ 0 := by positivity
  field_simp
  ring

/-- The exponent in the middle region of `thm:freelimit`: for `k ≥ 1`, `0 < τ ≤ 1` and
`m ≥ τ^{-1/2}`, `min(k²/(36τ), mk/6)/2 ≥ k/(72√τ)`. -/
theorem middle_region_exponent {k τ m : ℝ} (hk : 1 ≤ k) (hτ : 0 < τ) (hτ1 : τ ≤ 1)
    (hm : 1 / Real.sqrt τ ≤ m) :
    k / (72 * Real.sqrt τ) ≤ min (k ^ 2 / (36 * τ)) (m * k / 6) / 2 := by
  have hs : 0 < Real.sqrt τ := Real.sqrt_pos.mpr hτ
  have hs1 : Real.sqrt τ ≤ 1 := by
    rw [Real.sqrt_le_one]; exact hτ1
  have hss : Real.sqrt τ * Real.sqrt τ = τ := Real.mul_self_sqrt hτ.le
  have hk0 : 0 < k := by linarith
  have h1 : k / (36 * Real.sqrt τ) ≤ k ^ 2 / (36 * τ) := by
    rw [div_le_div_iff₀ (by positivity) (by positivity)]
    nlinarith [mul_pos hk0 hs, hss, hs1, mul_le_mul_of_nonneg_left hs1 hk0.le]
  have h2 : k / (36 * Real.sqrt τ) ≤ m * k / 6 := by
    have hm' : 1 ≤ m * Real.sqrt τ := by
      rw [div_le_iff₀ hs] at hm; linarith
    rw [div_le_div_iff₀ (by positivity) (by norm_num)]
    nlinarith [mul_le_mul_of_nonneg_left hm' hk0.le, mul_pos hk0 hs]
  have h3 : k / (36 * Real.sqrt τ) ≤ min (k ^ 2 / (36 * τ)) (m * k / 6) := le_min h1 h2
  have h4 : k / (72 * Real.sqrt τ) = k / (36 * Real.sqrt τ) / 2 := by
    field_simp; ring
  rw [h4]
  exact div_le_div_of_nonneg_right h3 (by norm_num)

/-- The tail estimates of `thm:freelimit`: `Σ_{k≥1} δ^k ≤ 2δ` for `δ ≤ 1/2` is the geometric
bound `δ/(1−δ) ≤ 2δ`, and `(p + 4δ)³ − p³ ≤ 108 δ` for `0 ≤ p ≤ 1`, `0 ≤ δ ≤ 1/2`. -/
theorem tail_region_bounds {p δ : ℝ} (hp0 : 0 ≤ p) (hp : p ≤ 1) (hδ0 : 0 ≤ δ) (hδ : δ ≤ 1 / 2) :
    δ / (1 - δ) ≤ 2 * δ ∧ (p + 4 * δ) ^ 3 - p ^ 3 ≤ 108 * δ := by
  constructor
  · rw [div_le_iff₀ (by linarith)]; nlinarith
  · nlinarith [mul_nonneg hp0 hδ0, mul_nonneg (mul_nonneg hp0 hp0) hδ0, mul_nonneg hδ0 hδ0,
      mul_nonneg (mul_nonneg hδ0 hδ0) hδ0, mul_nonneg hp0 (mul_nonneg hδ0 hδ0)]

/-- The free gluons, `prop:freegluons`: `2 dim G · c(0)` for `SU(2)` and `SU(3)` with
`c(0) = 15 Z₄/(2π²)` is `45 Z₄/π²` and `120 Z₄/π²`. -/
theorem free_gluon_multiplicity (Z : ℝ) :
    2 * 3 * (15 * Z / (2 * Real.pi ^ 2)) = 45 * Z / Real.pi ^ 2 ∧
    2 * 8 * (15 * Z / (2 * Real.pi ^ 2)) = 120 * Z / Real.pi ^ 2 := by
  have hpi := Real.pi_pos
  constructor <;> (field_simp; ring)


/-! ### The pressure bound at one scale (Section `ssec:pressure`) -/

/-- Proposition `prop:exchange`(i): `ζ(h) ≤ δ = 2ζ(h) - ζ(2h) ≤ 2ζ(h)`. -/
theorem exchange_delta_bounds (ζh ζ2h : ℝ) (h0 : 0 ≤ ζ2h) (h1 : ζ2h ≤ ζh) :
    ζh ≤ 2 * ζh - ζ2h ∧ 2 * ζh - ζ2h ≤ 2 * ζh := by
  constructor <;> linarith

/-- Proposition `prop:exchange`(ii): the rearrangement giving
`ζ(τ) ≤ log 𝒵_{(b,c,ℓ)}(τ)/(ℓ-τ)`. -/
theorem tube_rearrangement (ζτ ζℓ A ℓ τ : ℝ) (hτ : 0 < τ) (hℓ : τ < ℓ)
    (hmono : ζℓ ≤ ζτ) (hA : ζτ ≤ A / ℓ + (τ / ℓ) * ζℓ) : ζτ ≤ A / (ℓ - τ) := by
  have hℓ0 : 0 < ℓ := by linarith
  have hd : 0 < ℓ - τ := by linarith
  rw [le_div_iff₀ hd]
  have h1 : ζτ * ℓ ≤ A + τ * ζℓ := by
    have := mul_le_mul_of_nonneg_right hA hℓ0.le
    rw [add_mul, div_mul_cancel₀ _ hℓ0.ne', mul_assoc, mul_comm ζℓ ℓ, ← mul_assoc,
      div_mul_cancel₀ _ hℓ0.ne'] at this
    linarith
  nlinarith

/-- `y/2 ≤ log(1+y) ≤ y` for `0 ≤ y ≤ 1`, used in Proposition `prop:pressrings`. -/
theorem log_one_add_two_sided (y : ℝ) (h0 : 0 ≤ y) (h1 : y ≤ 1) :
    y / 2 ≤ Real.log (1 + y) ∧ Real.log (1 + y) ≤ y := by
  have hpos : 0 < 1 + y := by linarith
  constructor
  · have := Real.one_sub_inv_le_log_of_pos hpos
    have hinv : (1 + y)⁻¹ * (1 + y) = 1 := inv_mul_cancel₀ hpos.ne'
    have key : (1 - y / 2 - (1 + y)⁻¹) * (1 + y) = y * (1 - y) / 2 := by
      rw [sub_mul, hinv]; ring
    have hk : 0 ≤ (1 - y / 2 - (1 + y)⁻¹) * (1 + y) := by
      rw [key]; exact div_nonneg (mul_nonneg h0 (by linarith)) (by norm_num)
    have h5 : (0 : ℝ) ≤ 1 - y / 2 - (1 + y)⁻¹ :=
      le_of_mul_le_mul_right (by rw [zero_mul]; exact hk) hpos
    linarith
  · have := Real.log_le_sub_one_of_pos hpos; linarith

/-- Proposition `prop:pressrings`(ii): the Casimir defect of the symmetrised rings
under halving of the cube of side `2h`, from the vacuum energy
`E₀ = -4V log cosh K - Σ (V/ℓ_i) φ(ℓ_i)`. -/
theorem rings_cube_defect (h lc φh φ2h : ℝ) (hh : h ≠ 0) :
    (-(4 * (8 * h ^ 3)) * lc - 3 * ((8 * h ^ 3) / (2 * h)) * φ2h)
      - 8 * (-(4 * h ^ 3) * lc - 3 * (h ^ 3 / h) * φh)
      = 12 * h ^ 2 * (2 * φh - φ2h) := by
  field_simp; ring

/-- `2k + 4 ≤ 4·2^k`, the exponent count in Proposition `prop:pressrings`(iii). -/
theorem kdy_exponent (k : ℕ) : 2 * k + 4 ≤ 4 * 2 ^ k := by
  induction k with
  | zero => norm_num
  | succ n ih => rw [pow_succ]; omega

/-- Proposition `prop:pressrings`(iii): `4^k N^4 ≤ N^(4·2^k)` for `N ≥ 2`, that is
`4^k N^(3-4·2^k) ≤ 1/N`. -/
theorem kdy_power (N k : ℕ) (hN : 2 ≤ N) : 4 ^ k * N ^ 4 ≤ N ^ (4 * 2 ^ k) := by
  have h1 : (4 : ℕ) ^ k = 2 ^ (2 * k) := by rw [pow_mul]; norm_num
  have h2 : 2 ^ (2 * k) ≤ N ^ (2 * k) := Nat.pow_le_pow_left hN _
  calc 4 ^ k * N ^ 4 = 2 ^ (2 * k) * N ^ 4 := by rw [h1]
    _ ≤ N ^ (2 * k) * N ^ 4 := Nat.mul_le_mul_right _ h2
    _ = N ^ (2 * k + 4) := by rw [← pow_add]
    _ ≤ N ^ (4 * 2 ^ k) := Nat.pow_le_pow_right (by omega) (kdy_exponent k)

/-- Corollary `cor:pressuv`: `t - t⁴/4 ≤ 3/4` for `t ≥ 0`, by the factorisation
`3/4 - t + t⁴/4 = (t-1)²(t²+2t+3)/4`. -/
theorem uv_quartic_bound (t : ℝ) (ht : 0 ≤ t) : t - t ^ 4 / 4 ≤ 3 / 4 := by
  have hid : 3 / 4 - t + t ^ 4 / 4 = (t - 1) ^ 2 * (t ^ 2 + 2 * t + 3) / 4 := by ring
  have hnn : 0 ≤ (t - 1) ^ 2 * (t ^ 2 + 2 * t + 3) / 4 := by positivity
  linarith

/-- Corollary `cor:pressuv`: if `2L ≤ (3/4)Δ s_*` with `8E s_*³ = Δ`, then
`4096 L³ E ≤ 27 Δ⁴`. -/
theorem uv_threshold (Δ E L s : ℝ) (hE : 0 < E) (hL : 0 ≤ L) (hs : 0 ≤ s)
    (hcube : 8 * E * s ^ 3 = Δ) (h : 2 * L ≤ 3 / 4 * Δ * s) :
    4096 * L ^ 3 * E ≤ 27 * Δ ^ 4 := by
  have hΔ : 0 ≤ Δ := by rw [← hcube]; positivity
  have h3 : (2 * L) ^ 3 ≤ (3 / 4 * Δ * s) ^ 3 :=
    pow_le_pow_left₀ (by positivity) h 3
  have h4 : 512 * E * (2 * L) ^ 3 ≤ 512 * E * (3 / 4 * Δ * s) ^ 3 :=
    mul_le_mul_of_nonneg_left h3 (by positivity)
  calc 4096 * L ^ 3 * E = 512 * E * (2 * L) ^ 3 := by ring
    _ ≤ 512 * E * (3 / 4 * Δ * s) ^ 3 := h4
    _ = 27 * Δ ^ 3 * (8 * E * s ^ 3) := by ring
    _ = 27 * Δ ^ 4 := by rw [hcube]; ring

/-- Lemma `lem:freemodes`: completing the square in the one-mode kernel, with
`c = 1 + a/2 + b` and `c² - (2+a)c + 1 = 0`. -/
theorem free_mode_square (x y a b c : ℝ) (hc : c = 1 + a / 2 + b) (hc0 : c ≠ 0)
    (hq : c ^ 2 - (2 + a) * c + 1 = 0) :
    -(1 / 2) * (x - y) ^ 2 - a / 4 * (x ^ 2 + y ^ 2) - b / 2 * y ^ 2
      = -(c / 2) * (y - x / c) ^ 2 - b / 2 * x ^ 2 := by
  have hb : b = c - 1 - a / 2 := by linarith
  have hq' : 1 + 2 * b * c - c ^ 2 = 0 := by rw [hb]; linear_combination hq
  have hdiff : (-(1 / 2) * (x - y) ^ 2 - a / 4 * (x ^ 2 + y ^ 2) - b / 2 * y ^ 2)
      - (-(c / 2) * (y - x / c) ^ 2 - b / 2 * x ^ 2)
      = x ^ 2 * (1 + 2 * b * c - c ^ 2) / (2 * c) := by
    rw [hb]; field_simp; ring
  rw [hq', mul_zero, zero_div, sub_eq_zero] at hdiff
  exact hdiff

/-- Theorem `thm:pressfree`(ii): Jordan's inequality `sin²(k/2) ≥ (k/π)²` on `[0, π]`. -/
theorem jordan_half_square (k : ℝ) (h0 : 0 ≤ k) (h1 : k ≤ Real.pi) :
    (k / Real.pi) ^ 2 ≤ Real.sin (k / 2) ^ 2 := by
  have hk : 0 ≤ k / 2 := by linarith
  have hk' : k / 2 ≤ Real.pi / 2 := by linarith
  have hs := Real.mul_le_sin hk hk'
  have heq : 2 / Real.pi * (k / 2) = k / Real.pi := by
    rw [div_mul_div_comm, mul_comm Real.pi 2, mul_div_mul_left _ _ (two_ne_zero)]
  rw [heq] at hs
  have hnn : 0 ≤ k / Real.pi := div_nonneg h0 Real.pi_pos.le
  exact pow_le_pow_left₀ hnn hs 2

/-- Theorem `thm:pressfree`(iv): the lattice gap `2N arsinh(μ/2N)` is at most `μ`. -/
theorem free_gap_le_mass (μ N : ℝ) (hμ : 0 ≤ μ) (hN : 0 < N) :
    2 * N * Real.arsinh (μ / (2 * N)) ≤ μ := by
  have hy : 0 ≤ μ / (2 * N) := by positivity
  have h1 : Real.arsinh (μ / (2 * N)) ≤ Real.sinh (Real.arsinh (μ / (2 * N))) :=
    Real.self_le_sinh_iff.mpr (Real.arsinh_nonneg_iff.mpr hy)
  rw [Real.sinh_arsinh] at h1
  have h2 := mul_le_mul_of_nonneg_left h1 (by positivity : (0:ℝ) ≤ 2 * N)
  rwa [mul_div_cancel₀ _ (by positivity : (2 * N : ℝ) ≠ 0)] at h2

/-! ### Extensivity and the ascent in time (Section ssec:ascent) -/

/-- Lemma `lem:cutcost`: the cutting coefficient `12T/P` is at most `1/2` when `P ≥ 24T`. -/
theorem cut_coefficient_half (T P : ℝ) (hT : 0 < T) (hP : 24 * T ≤ P) :
    12 * T / P ≤ 1 / 2 := by
  have hP0 : 0 < P := by linarith
  rw [div_le_iff₀ hP0]; linarith

/-- Proposition `prop:extensive`: a finite `f` with `f ≤ a f + g`, `0 ≤ a < 1`, satisfies
`f ≤ g/(1-a)`. -/
theorem extensive_fixed_point (f g a : ℝ) (_ha0 : 0 ≤ a) (ha1 : a < 1) (h : f ≤ a * f + g) :
    f ≤ g / (1 - a) := by
  have h1 : 0 < 1 - a := by linarith
  rw [le_div_iff₀ h1]; nlinarith

/-- Theorem `thm:ascent`: `e^y − 1 ≤ (e − 1) y` for `0 ≤ y ≤ 1`, by the convexity of the
exponential between `0` and `1`. -/
theorem exp_sub_one_le_chord (y : ℝ) (h0 : 0 ≤ y) (h1 : y ≤ 1) :
    Real.exp y - 1 ≤ (Real.exp 1 - 1) * y := by
  have hc := convexOn_exp.2 (Set.mem_univ (0:ℝ)) (Set.mem_univ (1:ℝ))
    (by linarith : 0 ≤ 1 - y) h0 (by ring : (1 - y) + y = 1)
  simp only [smul_eq_mul, mul_zero, zero_add, mul_one, Real.exp_zero] at hc
  linarith

/-- Theorem `thm:ascent`: the recursion `y_{i+1} = 16(e−1) e^{−δT} y_i` contracts when
`δT ≥ log(16(e−1))`. -/
theorem ascent_contraction (δ T y : ℝ) (hy : 0 ≤ y)
    (h : Real.log (16 * (Real.exp 1 - 1)) ≤ δ * T) :
    16 * (Real.exp 1 - 1) * Real.exp (-(δ * T)) * y ≤ y := by
  have hpos : 0 < 16 * (Real.exp 1 - 1) := by
    have := Real.add_one_lt_exp (one_ne_zero : (1:ℝ) ≠ 0); linarith
  have h2 : 16 * (Real.exp 1 - 1) ≤ Real.exp (δ * T) := by
    rw [← Real.exp_log hpos]; exact Real.exp_le_exp.mpr h
  have h3 : 16 * (Real.exp 1 - 1) * Real.exp (-(δ * T)) ≤ 1 := by
    rw [Real.exp_neg, mul_inv_le_iff₀ (Real.exp_pos _), one_mul]; exact h2
  calc 16 * (Real.exp 1 - 1) * Real.exp (-(δ * T)) * y ≤ 1 * y :=
        mul_le_mul_of_nonneg_right h3 hy
    _ = y := one_mul y

/-- Theorem `thm:ascentZ` and Corollary `cor:branch`: the volume ratios `64/96³ < 7.234·10⁻⁵`,
`512/96³ < 5.79·10⁻⁴`, the constant `30(1+1/6)·5.8·10⁻⁴ < 0.021`, and the small-branch
threshold `log 15.6 + 0.8 + 5.8·10⁻⁴·7.8 < 3.56`. -/
theorem ascent_constants :
    (64 : ℝ) / 96 ^ 3 < 7.234e-5 ∧ (512 : ℝ) / 96 ^ 3 < 5.79e-4 ∧
    (30 : ℝ) * (1 + 1 / 6) * 5.8e-4 < 0.021 := by
  refine ⟨by norm_num, by norm_num, by norm_num⟩

theorem small_branch_threshold :
    Real.log 15.6 + 0.8 + 0.00058 * 7.8 < (3.56 : ℝ) := by
  have h : Real.log 15.6 < 2.75 := by
    rw [Real.log_lt_iff_lt_exp (by norm_num)]
    have h1 : Real.exp 2.75 = Real.exp 1 ^ 2 * Real.exp 0.09375 ^ 8 := by
      rw [← Real.exp_nat_mul, ← Real.exp_nat_mul, ← Real.exp_add]; norm_num
    have hq := Real.quadratic_le_exp_of_nonneg (by norm_num : (0:ℝ) ≤ 0.09375)
    have h4 : (1 + 0.09375 + 0.09375 ^ 2 / 2 : ℝ) ^ 8 ≤ Real.exp 0.09375 ^ 8 :=
      pow_le_pow_left₀ (by norm_num) hq 8
    have he2 : (2.7182818283 : ℝ) ^ 2 ≤ Real.exp 1 ^ 2 :=
      pow_le_pow_left₀ (by norm_num) Real.exp_one_gt_d9.le 2
    have hm := mul_le_mul he2 h4 (by positivity) (by positivity)
    rw [h1]; refine lt_of_lt_of_le ?_ hm; norm_num
  generalize Real.log 15.6 = L at h ⊢
  norm_num at h ⊢
  linarith

/-- Proposition `prop:fieldringsall`(i): from `R² − sinh²h = ε²` with `ε > 0` and
`R + sinh h > 0`, `(R − sinh h)/ε = ε/(R + sinh h)`. -/
theorem field_ring_t_identity (R s ε : ℝ) (hε : 0 < ε) (hRs : 0 < R + s)
    (h : R ^ 2 - s ^ 2 = ε ^ 2) : (R - s) / ε = ε / (R + s) := by
  rw [div_eq_div_iff hε.ne' hRs.ne']; nlinarith [h]

/-- Proposition `prop:fieldringsall`(ii): for a `2×2` probability table the covariance of the
two signs is `4(p₁₁p₂₂ − p₁₂p₂₁)`. -/
theorem table_covariance (p11 p12 p21 p22 : ℝ) (h : p11 + p12 + p21 + p22 = 1) :
    (p11 - p12 - p21 + p22) - (p11 + p12 - p21 - p22) * (p11 - p12 + p21 - p22)
      = 4 * (p11 * p22 - p12 * p21) := by
  linear_combination (-(p11 - p12 - p21 + p22)) * h

/-- Proposition `prop:fieldringsall`(iii): `M_c/n ≤ 4` for `M_c = 2n+4`, `n ≥ 2`. -/
theorem tunnel_window_ratio (n : ℝ) (hn : 2 ≤ n) : (2 * n + 4) / n ≤ 4 := by
  rw [div_le_iff₀ (by linarith)]; linarith

/-- Proposition `prop:fieldringsall`(iv): `C(0,Q) = 1/(1+Q) > 1/2` for `0 < Q < 1`. -/
theorem field_ring_limit (Q : ℝ) (hQ : 0 < Q) (hQ1 : Q < 1) :
    Q * (Q + 0) / ((1 + Q) * (Q ^ 2 + 0)) = 1 / (1 + Q) ∧ 1 / 2 < 1 / (1 + Q) := by
  constructor
  · field_simp; ring
  · rw [div_lt_div_iff₀ (by norm_num) (by linarith)]; linarith

/-- The strong-coupling example: `θ_{1/10} = 3/11`, `c = 9/121`, and `m − 1 = 282` satisfies
`282 ≥ 148/(1−2c)⁴ > 281`. -/
theorem strong_coupling_example :
    (3 : ℝ) / 10 / (2 - 9 / 10) = 3 / 11 ∧
    (281 : ℝ) < 148 / (1 - 2 * (3 / 11) ^ 2) ^ 4 ∧ 148 / (1 - 2 * (3 / 11) ^ 2) ^ 4 ≤ (282 : ℝ) := by
  refine ⟨by norm_num, by norm_num, by norm_num⟩

end YangMills.OneScale
