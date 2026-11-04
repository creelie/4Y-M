# Stratified Gauss–Green for the moving fibre, and Conjecture B for local F

Labels: **proved** (modulo the cited standard theorems), **sketch**, **open**.
Setting: SU(2), ρ_G = 1, one type-(b) inner link j of the average k = (y, μ), private links Pr (8 links),
constraint map c_t = (V♯ of all averages) with the inner link moved to e^{tξ}j, |ξ| = 1. Q does not depend on
j, since j is in no plaquette meeting A_W ([CD] §4). Ambient space M = SU(2)^n with Haar measure dU.

This note answers the two items `AVERAGED_SHELL_MASS.md` leaves open: it **identifies the signed flux** of the
moving constraint, and **bounds its total variation**. It works with fibre integrals averaged over c, which is
the only form used downstream (Lemma `mixture`, `CALM_LINK.md` §2).

## 1. Pieces (proved)

Split M, up to a null set, into finitely many open semialgebraic pieces Ω_i:
- Ω_calm = {m < 1/2}; Ω_large = {m > 1};
- for each path x, Ω_x = {1/2 < m < 1, δ_x > δ_{x'} for every x' not identified with x}.

SU(2) = S³ ⊂ ℝ⁴ is algebraic, and δ_x² is a polynomial in the link entries, so every piece is semialgebraic. Its
topological boundary is a compact semialgebraic set of dimension ≤ n − 1. Such a set has finite H^{n−1} measure,
so each Ω_i has finite perimeter, and ∂*Ω_i ⊂ ∂Ω_i. [Federer 3.4.8 / Łojasiewicz; for semialgebraic sets,
e.g. Bochnak–Coste–Roy and Yomdin–Comte, *Tame geometry*, ch. 4.]

On each piece the branch formula V♯_i (χ affine or constant, logs principal, the maximiser fixed) is real-analytic
on a neighbourhood of the closure of Ω_i. Indeed δ_x ≤ 1 < √2 keeps the eigen-angles of v_x below π, so the logs
are analytic there. The proof of `privatesubmersion` gives J_i J_i* ≥ σ₀² on Ω_i, hence on its closure by
continuity, hence on a neighbourhood. So

  Y_i = −J_i*(J_i J_i*)⁻¹ ∂_t c_i,   Y_large = 0,

is real-analytic, hence Lipschitz, on a neighbourhood of each closure Ω̄_i, with |Y_i| ≤ (Λ/16 + 2π/3)/σ₀ ≤ 69.7
(`SHELL_TRANSPORT.md` §1; numerically ≤ 0.84) and |div Y_i| ≤ C_G (compact closure, finitely many branches).
J_i Y_i = −∂_t c_i holds identically.

## 2. The identity (proved)

**Lemma GG.** For bounded Borel φ on the constraint space and g Lipschitz on M,

  ∫ g φ(c_t) dU − ∫ g φ(c_0) dU = ∫₀ᵗ ∫ φ(c_s) dD_s(g) ds,

where the divergence measure is

  D_s(g) = Σ_i 1_{Ω_i} div(g Y_i) dU  −  Σ_{Σ_{ik}} g [Y·n]_{ik} dH^{n−1},

and [Y·n]_{ik} = Y_i·n_i + Y_k·n_k is the jump of the normal component across the common reduced boundary
Σ_{ik} = ∂*Ω_i ∩ ∂*Ω_k.

*Proof.* Take φ ∈ C¹ first. The map t ↦ c_t(U) is Lipschitz for a.e. U and differentiable off ∪∂Ω_i. So
d/dt ∫gφ(c_t) = ∫ g Dφ(c_t)[∂_t c_t] by dominated convergence. On Ω_i, ∂_t c = −J_i Y_i, so the integrand is
−g ∇_{Y_i}(φ∘c_t) = −div(g φ(c_t) Y_i) + φ(c_t) div(g Y_i). Apply the Gauss–Green theorem for sets of finite
perimeter with a Lipschitz field [Federer 4.5.6; Evans–Gariepy 5.8] on each Ω_i to the field g φ(c_t) Y_i,
which is Lipschitz on Ω̄_i. M is compact without boundary, so the boundary terms pair up on the common interfaces
and give the jump measure. This proves the derivative form for C¹ φ. Integrate in t. Both sides are integrals of
φ against finite Borel measures (the push-forwards under c_t of the measures above), so equality extends from
C¹ to bounded Borel φ by the monotone-class theorem. ∎

The weight e^{−λQ} and the localisation 1 − θ_t are absorbed into g. θ_t depends on t through ψ(·, E_t), which adds
the bounded term ∫ φ(c_t) g ∂_tθ_t dU.

## 3. Localisation fix (proved)

[CD]'s θ uses ψ = Σ_x δ_x², and ψ ≥ ρ²/64 only gives m ≥ ρ/32 (sixteen paths). Replace it with
ψ₄ = (Σ_x δ_x⁸)^{1/4}, which satisfies m² ≤ ψ₄ ≤ 2m². Take θ = 1 on [0, ρ²/16] and θ = 0 on [ρ²/4, ∞). Then:
- supp θ ⊂ {m < ρ/2} = R;
- supp(1 − θ) ⊂ {m ≥ ρ/(4√2)} ⊂ {m ≥ ρ/8};
- θ(ψ₄) is smooth (ψ₄ is smooth where it is nonzero, and θ ≡ 1 near 0);
- θ ≡ 1 on the chart (paths within O(κ + s) ≪ ρ/8 of P_{2y}).

θ is a proof device in [CD], not part of the model, so every [CD] statement about θ goes through unchanged.

## 4. Conjecture B for bounded local F (proved; C_G finite but not computed)

Average Ξ_j(F) over the fibre value c with any weight φ ≥ 0. Write f(c) = E_{π_0(c)}F, which is bounded and Borel.
Then

  ∫ φ Z_0 Ξ_j(F) dc = d/dt|₀ [ ∫ F(1−θ_t) e^{−λQ} φ(c_t) dU − ∫ (1−θ_t) e^{−λQ} (φf)(c_t) dU ].

Apply Lemma GG with g = F(1−θ)e^{−λQ} and g = (1−θ)e^{−λQ}. With |F| ≤ 1 and F Lipschitz (constant L_F), the
absolutely continuous part is bounded pointwise by

  e^{−λQ} 1_{m ≥ ρ/8} ( L_F|Y| + λ|DQ[Y]| + |div Y| + |∇θ||Y| + |∂_tθ| ) ≤ C_G (1 + L_F + λ) e^{−λQ} 1_{m ≥ ρ/8},

using |DQ| ≤ C per link (compactness). The jump part is supported on ∪Σ_{ik} ⊂ {m = 1/2} ∪ {m = 1} ∪ T, with T the
tie set of the maximiser inside the shell. Its total variation is ≤ 2·69.7 × the e^{−λQ}-weighted H^{n−1} mass of
these surfaces:
- {m = r}, r = 1/2 or 1: bounded by (4)–(5) of `AVERAGED_SHELL_MASS.md` by C(K + 1) μ(m ≥ r), with K = O(β_W).
- T: **proved.** Fix a pair a ≠ b of non-identified paths in the shell. Each path is a product of one inner link
  (u or w) and one private link, and **path a has an inner link e_a that path b does not contain**. If that link is the moved link j (held
  fixed, not integrated), swap a and b: distinct non-identified paths have distinct inner links, so at most one of
  them is j. Then
  h = δ_a − δ_b depends on e_a only through δ_a = d_HS(e_a, g), where the centre g is a function of the other
  links.
  - *Conditional law.* Q = Σ_p δ(U_p, I)² = Σ_p (4 − 2 Re tr U_p) (main.tex l. 6522) is affine in the entries of
    each single link. So with all other links frozen, the density of e_a under Haar·e^{−λQ} is ∝ exp(κ a·u) on S³,
    with κ ≤ 2λ·Σ|staples| ≤ 24√2 λ (six plaquettes; κ = 0 if e_a is in no plaquette of Q). This is exactly the
    hypothesis of (1) in `SU2_RADIAL_HAZARD.md`, uniformly in the frozen links. (Concentration κ here, not the
    plaquette threshold K of `CALM_LINK.md`.)
  - *Surface measure by slicing.* |∇h| ≤ 8 (speed bound (3) of `AVERAGED_SHELL_MASS.md`, twice), and
    |∂_{e_a} h| = dδ/dR ≥ 1 on R ≤ π/2 (and the shell has δ ≤ 1, i.e. R < 0.73). The coarea/slicing inequality
    then gives H^{n−1}⌊{h = 0} ≤ 8 · (H² of the sphere {d(e_a, g) = δ_b} in the e_a variable) ⊗ (the other links).
  - *Integrate.* The e^{−λQ}-weighted H² mass of that sphere is the conditional radial density f_g at δ_b. By (1)
    it is ≤ e sec²(1)(κ+1)·P(R_a ≥ δ_b | frozen) ≤ e sec²(1)(κ+1)·P(m ≥ 1/2 | frozen) on the shell, pointwise in
    the frozen links. Integrating over them gives
    weighted H^{n−1}(T_{ab}) ≤ 8 e sec²(1)(κ+1) μ(m ≥ 1/2).
  The same slicing argument, with h = δ_a − r and e_a the inner link of the maximising path, proves the hazard
  hypothesis that (4) of `AVERAGED_SHELL_MASS.md` assumes for the levels {m = r}. So that bound is now
  unconditional too, with κ ≤ 24√2 λ in place of its K.
  There are fewer than 16² pairs, and |[Y·n]| ≤ 2·69.7. Every bound is pointwise in the level, so no
  a.e.-level/continuity issue arises.

Dividing by ∫ φ Z_0 dc and using |f| ≤ 1:

  |⟨Ξ_j(F)⟩_φ| ≤ C_G (1 + L_F + λ + K) ⟨π(m ≥ ρ/8)⟩_φ.

With κ, K = O(λ), this is **Conjecture B of [CD] for bounded Lipschitz F, in φ-averaged form**. It is
volume-independent, and the averaged form is what Lemma `mixture` consumes. It is **proved**, with C_G a finite constant
(compactness, finitely many analytic branches) that has not been evaluated.

## 5. What this does not do (open)

The bound has no decay in the distance from j to supp F. The cube coefficient sums Ξ over all j, so it needs
|Ξ_j(F)| ≲ e^{−α d(j, F)}. That is Conjecture L, whose large-field part is (R1) / D7. See `SHELL_TRANSPORT.md` §3
and `CALM_LINK.md` §4. This note does not move that wall.
