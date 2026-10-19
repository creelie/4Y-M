# Four-dimensional SU(3) in the forest gauge: β_W ≤ 45/64

Status: proved by hand; constants certified in ball arithmetic (`checks/su3/certify_su3.py` for Lemmas B–E,
`checks/su3/certify_tilted.py` for Lemmas G–J, `checks/su3/certify_blend.py` for Lemmas D′–H′, python-flint 0.9.0),
and the evaluation of Section 5 at every staple sum by rigorous floating-point error bounds and ball
arithmetic (`checks/su3/exact/`); not yet read by an outside referee. Ranges known before for four-dimensional SU(3): β_W < 9/32 = 0.281
(Theorem B of the paper in `../prd/`, sharp Hessian 4d), β_W < 9/29 = 0.310 (gauge-fixed curvature,
`curvature.md`); Shen–Zhu–Zhu (arXiv 2204.12737) has a smaller explicit range. Search of 19 Oct 2026 (arXiv
and web, "SU(3) lattice Yang–Mills strong coupling uniqueness Dobrushin explicit beta"): no larger
explicit SU(3) range found; arXiv 2506.00284 claims a full SU(3) mass gap and is not used.

Normalisation of the paper: weight exp(Nβ Σ_p Re Tr Q_p), N = 3, β_W = 9β, β ≥ 0. On SU(3) the
inner product of the Lie algebra is ⟨X,Y⟩ = −Tr(XY), so Ric = 3/2, the Casimir of the fundamental
representation is λ_F = 8/3 and that of the adjoint is 6. The distance used for Wasserstein bounds is
d(U,V) = ‖U − V‖_HS, which is at most the geodesic distance. A(M) = (M − M*)/2 and P removes the trace;
K = conv(SU(3)).

## 1. Reduction to one link

Use the line field F of `forest.md` and its specification γ^F (links of F set to 1). The conditional
law of an active link e is

    ν_b(dU) ∝ exp(Re Tr(U b)) dU,   b = Nβ Σ_{p ∋ e} S_p,

where S_p is the staple of e in p (a product of SU(3) matrices and of the identity on F-links, or its
adjoint), so b ∈ rK with r = 6Nβ = 2β_W. Changing a partner e′ changes one staple by a matrix of norm
d(U_{e′}, U′_{e′}). Along the segment b_t from b to b′ (inside rK by convexity),
d/dt ν_{b_t}(g) = Cov_{ν_{b_t}}(g, Re Tr(U(b′ − b))), so

    W_1(ν_b, ν_{b′}) ≤ Nβ Λ(r) d(U_{e′}, U′_{e′}),   Λ(r) = sup_{b ∈ rK} sup_{‖h‖_HS = 1} L_h(b),
    L_h(b) = sup { Cov_{ν_b}(g, φ_h) : g 1-Lipschitz for d },   φ_h(U) = Re Tr(U h).

Every active link has 12 active partners (Lean `twelve_neighbours`), so the Dobrushin row and column
sums are α = 12 Nβ Λ(r) = 2r Λ(r). If α < 1, the argument of `forest.md` applies word for word with
SU(2) replaced by SU(3) and the chordal metric by d: the comparison lemma of the paper (stated there for
SU(2), whose proof uses only that δ is a bounded metric and that the C_{ee′} bound the Wasserstein
influences), the free-boundary limit ω of the Wilson measures for every gauge-invariant local
observable, translation and hypercubic invariance, reflection positivity, decay at every rate below
−log α, and Prop. gapdecay.

## 2. The single-link estimate

**Lemma A (gradients).** In the left trivialisation, ∇φ_c(U) = −P A(cU), since
d/dt Re Tr(U e^{tX} c) = Re Tr(cUX) = −⟨P A(cU), X⟩ for X ∈ su(3). Hence

1. |∇φ_h|² ≤ ‖A(hU)‖²_HS ≤ ‖hU‖²_HS = ‖h‖²_HS, because ‖M‖² = ‖A(M)‖² + ‖(M + M*)/2‖²;
2. for V = φ_b with b ∈ rK, |∇V| ≤ g r with g² = max_{W ∈ SU(3)} ‖P A(W)‖²_HS: bU ∈ rK, and the convex
   function M ↦ ‖P A(M)‖² attains its maximum over K at an extreme point, which lies in SU(3);
3. Hess V(X,X) = Re Tr(U X² b) ≤ ‖X²‖_{S_1} ‖b‖_op ≤ r |X|², since K lies in the unit ball of the
   operator norm (on a group with bi-invariant metric, ∇_X X = 0 for left-invariant X);
4. Δφ_c = −λ_F φ_c, and the Haar mean of φ_h² is ‖h‖²/6, so the Haar mean of |∇φ_h|² is
   λ_F/6 = 4/9 for ‖h‖ = 1.

**Lemma B (the constant g).** For W with eigenvalues e^{iφ_k}, φ_1 + φ_2 + φ_3 = 0,
‖P A(W)‖² = Σ sin²φ_k − (Σ sin φ_k)²/3. Its maximum is 2.065471…, attained at φ = (−2t, t, t),
t = 0.935929…; a branch and bound over [−π, π]² in ball arithmetic proves g² ≤ 2.06548.

**Lemma C (rough Laplacian on 1-forms).** The bottom of the spectrum of ∇*∇ on L² 1-forms of SU(3)
is μ_0 = 7/6.

Proof. Trivialise 1-forms on the left: ω = Σ_a ω_a θ^a. With X_a an orthonormal basis acting as
left-invariant vector fields and ∇_{X_a} Y = ½[X_a, Y], ∇*∇ = −Σ_a (ρ(X_a) + ½ ad X_a)². On the
isotypic part V_π ⊗ V_π* ⊗ su(3) of L² ⊗ su(3), with A_a = ρ_π(X_a), B_a = ad X_a and T_a = A_a + B_a
(the tensor product action), A_a + ½B_a = ½(A_a + T_a) and A_aT_a + T_aA_a = A_a² + T_a² − B_a², so on
the component σ ⊂ π ⊗ ad

    ∇*∇ = (C_π + C_σ)/2 − C_ad/4 = (C_π + C_σ)/2 − 3/2,

with C(p,q) = 2(p² + q² + pq + 3p + 3q)/3. The components of (p,q) ⊗ (1,1) are (p+1,q+1),
(p+2,q−1), (p−1,q+2), (p+1,q−2), (p−2,q+1), (p−1,q−1) and (p,q) with multiplicity
2, 1 or 0 as both, one or none of p, q are positive (labels must be nonnegative). The four smallest
Casimirs among them give

    (p,q):       2(p²+q²+pq+3p+3q)/3 − 3/2,
    (p−1,q−1):   (2p²+2q²+2pq+3p+3q−3)/3 − 3/2,
    (p−2,q+1):   (2p²+2q²+2pq+3p+6q)/3 − 3/2,  and (p+1,q−2) symmetrically,

each at least (2(p²+q²+pq) + 3(p+q) − 3)/3 − 3/2 ≥ 3/2 when p + q ≥ 2. For π = (0,0), σ = (1,1)
gives 3/2; for π = (1,0) the components are (2,1), (0,2), (1,0), and σ = (1,0) gives
8/3 − 3/2 = 7/6, the minimum; (0,1) is the conjugate. The enumeration in `checks/su3/certify_su3.py` (part A)
repeats this for p, q ≤ 40 and checks the dimensions of every decomposition. ∎

**Lemma D (spectral gap).** Let gap(r) = (√μ_0 − g r/2)² + 3/2 − r for 0 ≤ g r ≤ 2√μ_0. For b ∈ rK
the spectral gap of L = Δ + ⟨∇V, ∇·⟩ on L²(ν_b) is at least gap(r).

Proof. Let −Lf = λf, λ > 0, ω = df. Bochner's formula integrated against ν_b, with ν_b⟨ω, ∇Lf⟩ =
−λ ν_b|ω|², gives

    λ ν_b|ω|² = ν_b|∇ω|² + ν_b((Ric − Hess V)(ω,ω)) ≥ ν_b|∇ω|² + (3/2 − r) ν_b|ω|²

by Lemma A(3). Put η = e^{V/2} ω. Then e^{V/2}∇ω = ∇η − ½ dV ⊗ η, so in L²(dU), by Minkowski,
Lemma A(2) and Lemma C, ‖e^{V/2}∇ω‖ ≥ ‖∇η‖ − ½ g r ‖η‖ ≥ (√μ_0 − g r/2)‖η‖. Since
ν_b|∇ω|² = ‖e^{V/2}∇ω‖²/Z and ν_b|ω|² = ‖η‖²/Z, λ ≥ (√μ_0 − g r/2)² + 3/2 − r. ∎

At b = 0 the bound is 7/6 + 3/2 = 8/3, the first eigenvalue of SU(3), so it is sharp there.

**Lemma E (the gradient moment).** Let J(r) = ∫_0^1 g r / √gap(t r) dt and
m(r) = sin²(min(arcsin(2/3) + J(r)/2, π/2)). For b ∈ rK and ‖h‖ = 1, ν_b|∇φ_h|² ≤ m(r).

Proof. u(t) = ν_{tb}(F), F = |∇φ_h|², is smooth, and 0 < u < 1 because 0 ≤ F ≤ 1 and F < 1 on an
open set. u′ = Cov_{ν_{tb}}(F, V), Var F ≤ u(1 − u), and Var_{ν_{tb}} V ≤ ν_{tb}|∇V|²/gap(tr) ≤
g²r²/gap(tr) by Lemma D (tb ∈ trK). So u′ ≤ √(u(1−u)) g r/√gap(tr), that is
(arcsin √u)′ ≤ g r/(2√gap(tr)), and u(0) = 4/9 by Lemma A(4). ∎

**Lemma F (corrector).** For b ∈ rK and ‖h‖ = 1, with Γ = ⟨∇V, ∇φ_h⟩,

    L_h(b) ≤ (3/8) √(ν_b|∇φ_h|²) + (3/8) √(Var_{ν_b}(Γ)/gap(r)).

Proof. Let −Lψ = φ_h − ν_bφ_h. For g Lipschitz, Cov(g, φ_h) = ν_b⟨∇g, ∇ψ⟩ ≤ ν_b|∇ψ|, as |∇g| ≤ 1
a.e. Write ψ = φ_h/λ_F + χ. Since −L(φ_h/λ_F) = φ_h − Γ/λ_F, χ solves −Lχ = R := (Γ − ν_bΓ)/λ_F.
Then ν_b|∇ψ| ≤ √(ν_b|∇φ_h|²)/λ_F + √(ν_b|∇χ|²), and ν_b|∇χ|² = ν_b(Rχ) ≤ √(Var R) √(Var χ) ≤
√(Var R) √(ν_b|∇χ|²/gap(r)) by Lemma D, so ν_b|∇χ|² ≤ Var(Γ)/(λ_F² gap(r)). With 1/λ_F = 3/8 this is
the claim. ∎

Since |Γ| ≤ |∇V||∇φ_h| ≤ g r |∇φ_h|, Var(Γ) ≤ ν_b(Γ²) ≤ g²r² ν_b|∇φ_h|², and Lemma E gives the first bound

    Λ(r) ≤ Λ̄(r) := (3/8) √m(r) (1 + g r/√gap(r)).

At r = 0, Λ̄(0) = 1/4, which is the Cauchy–Schwarz value 1/√(2(N² − 1)); the true value of L at b = 0
lies between 1/6 and 1/4.

## 3. Tilted moments

The first bound loses twice: Lemma E follows the moment along a path with the crude covariance bound,
and Lemma F then uses |Γ| ≤ g r|∇φ_h|. Both quantities can be bounded by moments of the Haar measure
tilted by exp(r Re Tr U), which are computed exactly. For a matrix c write
G_{c,h} = ⟨∇φ_c, ∇φ_h⟩; by Lemma A, with ⟨PA(M), PA(M′)⟩ = ½ReTr(MM′*) − ½ReTr(MM′) −
⅓ ImTr M ImTr M′,

    G_{c,h}(U) = ½ ReTr(c h*) − ½ ReTr(cUhU) − ⅓ ImTr(cU) ImTr(hU),

real bilinear in (c, h), with Haar mean (4/9) ReTr(c h*) (Lemma A(4) and polarisation). Thus
|∇φ_h|² = G_{h,h} and Γ = G_{b,h}.

**Lemma G (tilted gradient moment).** For b ∈ rK and ‖h‖ = 1,

    ν_b|∇φ_h|² ≤ m̃(r) := max_{‖h‖ = 1} E_H[ |∇φ_h|² e^{r Re Tr U} ].

Proof. Z(b) = E_H e^{φ_b} ≥ e^{E_H φ_b} = 1 by Jensen, so ν_b(F) ≤ E_H[F e^{φ_b}] for F ≥ 0. For F =
|∇φ_h|² the right side is convex in b, so its maximum over rK = conv(r SU(3)) is attained at some
b = rW, W ∈ SU(3). The substitution U ↦ UW⁻¹ (Haar invariance) turns e^{Re Tr(rUW)} into
e^{r Re Tr U} and |∇φ_h|²(UW⁻¹) = ‖PA(hUW⁻¹)‖² into ‖PA(W⁻¹hU)‖² = |∇φ_{W⁻¹h}|²(U), because
conjugation by W preserves A, P and the norm; and ‖W⁻¹h‖ = 1. ∎

**Lemma H (tilted variance).** For b ∈ rK and ‖h‖ = 1, Var_{ν_b}(Γ) ≤ r² κ̃(r), where

    κ̃(r) := max { E_H[ (G_{W,h} − (4/9) ReTr(W h*))² e^{r Re Tr U} ] : W ∈ T, ‖h‖ = 1 },

and T is the diagonal torus diag(e^{iθ_1}, e^{iθ_2}, e^{−i(θ_1+θ_2)}) of SU(3).

Proof. Write b = r a with a ∈ K and g_a = G_{a,h} − (4/9)ReTr(a h*), which is linear in a. For any
constant c, Var(Γ) = r² Var(G_{a,h}) ≤ r² ν_b(G_{a,h} − c)², so with Z(b) ≥ 1,
Var(Γ) ≤ r² Φ(a, b), Φ(a, b) := E_H[g_a² e^{φ_b}]. Bound Φ(a, b) by its maximum over K × rK. It is
convex in a for fixed b and convex in b for fixed a, so the maximum is attained at a = W ∈ SU(3),
b = rW′ ∈ r SU(3). The substitution U ↦ UW′⁻¹ replaces (W, h) by (W′⁻¹W, W′⁻¹h) and the tilt by
e^{r Re Tr U}, as in Lemma G. The tilted Haar measure is invariant under U ↦ VUV⁻¹, V ∈ SU(3), which
replaces (W, h) by (V⁻¹WV, V⁻¹hV); choosing V that diagonalises W puts W in T. ∎

Both maxima are maxima of largest eigenvalues of quadratic forms in h ∈ M_3(C) = R^18 (real coordinates in the
basis E_ij, iE_ij, orthonormal for ReTr(XY*)): Q_pq = E_H[G_{e_p,e_q} e^{r Re Tr U}] for m̃, and
N(θ)_pq = E_H[g_{W,e_p} g_{W,e_q} e^{r Re Tr U}] for κ̃. Writing W = Σ_a w_a d_a with d_a the six real
and imaginary diagonal units and w(θ) = (cos θ_1, sin θ_1, cos θ_2, sin θ_2, cos θ_3, sin θ_3),
θ_3 = −θ_1 − θ_2, gives N(θ) = Σ_{a,b} w_a w_b K_ab with fixed symmetric 18 × 18 matrices
K_ab = K_ba. Conjugation by the diagonal matrices of SU(3) commutes with every W ∈ T and multiplies
h_ij by e^{i(δ_j − δ_i)}, so by Schur's lemma each K_ab is block diagonal for the real isotypic
decomposition of M_3(C): the six coordinates of the diagonal entries, and for each pair i < j the
four coordinates of h_ij, h_ji.

**Lemma I (exact tilted moments).** Every entry of Q and of the K_ab is an exact rational combination
of the numbers S(n) = Σ_{k ∈ Z} I_{k−n_1}(r) I_{k−n_2}(r) I_{k−n_3}(r), n ∈ Z³, |n_i| ≤ 12.

Proof. The entries are tilted Haar moments of products of at most four entries of U and U*. For
X_1, …, X_n ∈ {U, U*}, n ≤ 4, the tensor T_X = E_H[X_1 ⊗ ⋯ ⊗ X_n e^{r Re Tr U}] commutes with V^{⊗n}
for every V ∈ U(3), since U ↦ VUV* preserves Haar measure and Re Tr U. By Schur–Weyl duality
T_X = Σ_σ m_σ P_σ over the permutation operators of (C³)^{⊗n}. Pairing with P_τ gives the linear
system Σ_σ 3^{c(τ⁻¹σ)} m_σ = t_τ, where c counts cycles and t_τ = Tr(P_τ* T_X) is the tilted mean of a
product over the cycles of τ of traces Tr(U^k), k the signed number of factors U in the cycle (U* = U⁻¹
collapses the products). For n ≤ 3 the matrix is invertible; for n = 4 its kernel is spanned by the
sign vector s, because the antisymmetriser vanishes on (C³)^{⊗4}, and t ⊥ s for the same reason, so
the solution of the system with matrix G + ssᵀ also solves the original one; any two solutions define
the same operator.
The t_τ are tilted means of class functions. By the Weyl integration formula they equal
(1/6)(2π)^{-2} ∫∫ f |Δ|² e^{r Σ cos θ_j} dθ_1 dθ_2 on T, and with e^{r cos θ} = Σ_k I_k(r) e^{ikθ}
the torus integral of z_1^{n_1} z_2^{n_2} z_3^{n_3} e^{r Σ cos θ_j} is S(n), since z_1^{a_1} z_2^{a_2} z_3^{a_3}
integrates to 1 if a_1 = a_2 = a_3 and to 0 otherwise. The degrees give |n_i| ≤ 12. ∎

For 0 ≤ r ≤ R ≤ 2, I_m(r) ≤ (R/2)^m e^{R²/4}/m!, which decreases in m, so the terms |k| > 40 of S(n) add
at most 4[(R/2)^{29} e^{R²/4}/29!]³. At r = 0 the computation returns Q = (4/9)·Id, the value of Lemma A(4).

**Lemma J (the torus).** Let θ range over a square of half-width ρ about θ_c, and put
N_c = N(θ_c), D_i = 2 Σ_{a,b} w_a(θ_c) ∂_i w_b(θ_c) K_ab. Then, block by block,

    λ_max N(θ) ≤ max_{±,±} λ_max(N_c ± ρ D_1 ± ρ D_2) + 8 ρ² Σ_{a,b} ‖K_ab‖_F.

Proof. Each w_a is cos or sin of θ_1, θ_2 or θ_3, so along the segment θ_c + t s, |s_i| ≤ ρ, the second
derivative of w_a w_b is at most (|ℓ_a| + |ℓ_b|)² ≤ 16ρ² with ℓ_a ∈ {s_1, s_2, −s_1 − s_2}; Taylor's
formula bounds the remainder after the linear term by 8ρ² Σ‖K_ab‖_op ≤ 8ρ² Σ‖K_ab‖_F. The function s ↦ λ_max(N_c + s_1D_1
+ s_2D_2) is convex, so its maximum over the square is at a corner. ∎

`checks/su3/certify_tilted.py` builds Q and the K_ab exactly (it checks that the imaginary parts and the
off-block entries vanish identically), evaluates them in ball arithmetic at r = 0.98, proves
λ_max(Q) < 0.7123 by an LDLᵀ factorisation of 0.7123·Id − Q with certified positive pivots, covers
[0, 2π)² by 8464 dyadic squares (total area checked exactly) and proves the bound of Lemma J below
0.1830 on each by the same factorisation. Hence m̃(0.98) < 0.7123 (the float value is 0.712207) and
κ̃(0.98) < 0.1830. Lemmas F, G, H give the second bound

    Λ(r) ≤ Λ̄₂(r) := (3/8) √m̃(r) + (3/8) r √(κ̃(r)/gap(r)).

## 4. A sharper gap, a normalised moment and a blended corrector

Three further steps remove part of each loss of Section 3; what remains is measured in Section 7.

**Lemma D′ (sharper gap).** For b ∈ rK the spectral gap of ν_b is at least

    gap₂(r) = (39/64) μ_0 + 3/2 − r/2 = 283/128 − r/2.

Proof. Keep λ, ω = df and η = e^{V/2}ω from Lemma D, and take norms in L²(dU). Bochner's formula reads

    λ‖η‖² = ‖∇η − ½ dV ⊗ η‖² + ∫ (Ric − Hess V)(η, η).                                    (1)

Put A = ∇η and B = ½ dV ⊗ η. Then ⟨A, B⟩ = ½ ∫ Σ_a ⟨∇_{X_a}η, η⟩ X_aV = ¼ ∫ ⟨∇V, ∇|η|²⟩ = (λ_F/4) ∫ V|η|²,
because V = φ_b is a matrix coefficient of the fundamental representation, so −ΔV = λ_F V. Dropping
‖(5/8)A − B‖² ≥ 0 from ‖A − B‖² leaves

    ‖A − B‖² ≥ (39/64)‖A‖² − (3/4)⟨A, B⟩ = (39/64)‖∇η‖² − ½ ∫ V|η|² ≥ (39/64) μ_0‖η‖² − ½ ∫ V|η|²

by Lemma C. With Ric = 3/2, (1) gives λ‖η‖² ≥ ((39/64)μ_0 + 3/2)‖η‖² − ∫ (Hess V(η, η) + ½ V|η|²), and
the claim follows from the pointwise bound

    Hess V(X, X) + ½ V|X|² = Re Tr(M X²) + ½ Re Tr M ≤ (r/2)|X|²,   M = bU ∈ rK.              (2)

The left side of (2), maximised over |X| = 1, is convex in M, so it suffices to prove (2) at M = rW,
W ∈ SU(3), and after a conjugation at W = diag(e^{iθ_k}) with θ_1 + θ_2 + θ_3 = 0; put c_k = cos θ_k.
Write X = iH with H Hermitian, traceless and Tr H² = 1, with diagonal d and off-diagonal entries H_jk.
Then

    Re Tr(W X²) = −Σ_k c_k d_k² − Σ_{j<k} (c_j + c_k)|H_jk|²,   Σ_k d_k² + 2 Σ_{j<k} |H_jk|² = 1,

so the maximum over X is the larger of q_off = max_{j<k} −(c_j + c_k)/2 and
q_diag = max{−Σ_k c_k d_k² : Σ_k d_k = 0, Σ_k d_k² = 1}. If l is the third index,
q_off + ½Σ_k c_k = c_l/2 ≤ ½. q_diag is minus the smaller eigenvalue of diag(c) compressed to the plane Σ_k d_k = 0; the eigenvalues μ of
the compression solve Σ_k 1/(c_k − μ) = 0, that is 3μ² − 2e_1μ + e_2 = 0 with e_1, e_2 the elementary
symmetric functions of the c_k, so q_diag = (√(e_1² − 3e_2) − e_1)/3. Then q_diag + e_1/2 ≤ ½ means
√(e_1² − 3e_2) ≤ (3 − e_1)/2, and since both sides are nonnegative this is e_1² − 4e_2 + 2e_1 − 3 ≤ 0.
When θ_1 + θ_2 + θ_3 = 0, Σ_k c_k² = 1 + 2c_1c_2c_3, and with this identity the left side equals
−2(1 − c_1)(1 − c_2)(1 − c_3) ≤ 0. ∎

The weight 5/8 is the one for which the coefficient of V is ½, and (2) is an equality at W = I. Lemma D
is the better bound for r < 0.24 (at b = 0 it gives the exact value 8/3, Lemma D′ gives 2.21); at
r = 1.14, gap₂ = 1.64094 against 0.428 for Lemma D. `checks/su3/gap_lemma.py` checks (2) and the closed
form of the maximum over X on a grid of the torus and at random W, and compares gap₂ with Galerkin
values (Section 7).

**Lemma G′ (normalised gradient moment).** For b ∈ rK and ‖h‖ = 1,

    ν_b|∇φ_h|² ≤ m̂(r) := m̃(r)/(m̃(r) + (5/9) e^{−r/5}).

Proof. Let F = |∇φ_h|², so 0 ≤ F ≤ 1 by Lemma A(1) and E_H F = 4/9 by Lemma A(4). Then ν_b F = N/(N + N′)
with N = E_H[F e^{φ_b}] ≤ m̃(r) (proof of Lemma G) and N′ = E_H[(1 − F) e^{φ_b}], and N/(N + N′) increases
in N and decreases in N′. Jensen's inequality for the probability measure (9/5)(1 − F) dU and E_H φ_b = 0
give

    N′ ≥ (5/9) exp((9/5) E_H[(1 − F)φ_b]) = (5/9) exp(−(9/5) Re Tr(E_H[F U] b)).

Here E_H[F U] = adj(h)/9. Indeed, Haar measure is invariant under U ↦ ωU, ω³ = 1, so a moment of degree
j in the entries of U and k in those of Ū vanishes unless j ≡ k mod 3. In F·U = (½‖h‖² − ½ReTr(hUhU) − ⅓(ImTr hU)²)U
only −¼Tr(hUhU)U and (1/12)Tr(hU)²U survive, and E_H[U_ia U_jb U_kc] = ε_ijk ε_abc/6 gives
E_H[Tr(hU)² U] = adj(h)/3 and E_H[Tr(hUhU) U] = −adj(h)/3, so E_H[FU] = adj(h)/12 + adj(h)/36. (Part (A)
of `checks/su3/certify_blend.py` checks the identity exactly, coefficient by coefficient.) The singular
values of adj(h) are s_2s_3, s_1s_3, s_1s_2 when those of h are s_1, s_2, s_3, so its trace norm is
((Σ s_k)² − Σ s_k²)/2 ≤ (3 − 1)/2 = 1 and |Re Tr(adj(h) b)| ≤ r. So N′ ≥ (5/9) e^{−r/5}. ∎

**Lemma F′ (blended corrector).** For b ∈ rK, ‖h‖ = 1 and s ≥ 0, with Γ as in Lemma F and
X_s = (1 − sλ_F)φ_h + sΓ,

    L_h(b) ≤ s √(ν_b|∇φ_h|²) + √(Var_{ν_b}(X_s)/gap₂(r)).

Proof. As for Lemma F, with ψ = sφ_h + χ: since −L(sφ_h) = sλ_Fφ_h − sΓ, χ solves −Lχ = X_s − ν_bX_s, and
ν_b|∇χ|² ≤ Var(X_s)/gap₂ by Lemma D′. The choice s = 1/λ_F = 3/8 is Lemma F. ∎

**Lemma H′ (blended variance).** Put c_1 = 1 − sλ_F and c_2 = s r. For b ∈ rK and ‖h‖ = 1,

    Var_{ν_b}(X_s) ≤ κ_s(r) := max{ E_H[(c_1 φ_h + c_2 g_{W,h})² e^{r Re Tr U}] : W ∈ T, ‖h‖ = 1 },

with g_{W,h} = G_{W,h} − (4/9) ReTr(W h*) as in Lemma H.

Proof. As for Lemma H, with b = r a: X_s minus the constant c_2 (4/9) ReTr(a h*) is c_1φ_h + c_2 g_{a,h},
which is affine in a, so Var(X_s) ≤ E_H[(c_1φ_h + c_2 g_{a,h})² e^{φ_b}], which is convex in a and in b
separately. The substitutions U ↦ UW′⁻¹ and U ↦ VUV⁻¹ turn φ_h into φ_{W′⁻¹h} and φ_{V⁻¹hV}, the same
changes of h as in g. ∎

The quadratic form in h is N_s(θ) = c_1² P + c_1c_2 Σ_a w_a (R_a + R_aᵀ) + c_2² Σ_{a,b} w_a w_b K_ab with
P_pq = E_H[φ_{e_p} φ_{e_q} e^{r Re Tr U}] and (R_a)_pq = E_H[φ_{e_p} g_{d_a,e_q} e^{r Re Tr U}]. Diagonal
conjugation preserves the whole integrand, so P and the R_a have the block structure of the K_ab, and
every entry is again a rational combination of the S(n) of Lemma I. In Lemma J the linear part adds
c_1c_2 Σ_a ∂_i w_a (R_a + R_aᵀ) to D_i; the second derivative of w_a along a segment is at most 4ρ², so
the remainder becomes ρ² (2|c_1c_2| Σ_a ‖R_a + R_aᵀ‖_F + 8 c_2² Σ_{a,b} ‖K_ab‖_F).

`checks/su3/certify_blend.py` checks E_H[FU] = adj(h)/9 exactly (every polarised coefficient, with the
exact Haar values S(n) = 1 if n_1 = n_2 = n_3 and 0 otherwise), proves m̃(1.14) < 0.8051 (float value
0.805063), so m̂(1.14) < 0.64543, builds P, the R_a and the K_ab exactly, and at r = 1.14, s = 0.322
covers the torus by 9778 dyadic squares on each of which the bound of Lemma J is below 0.0508. The
float maximum over a 90 × 90 grid is 0.049821, at the central elements ωI, ω²I of T (0.04950 at I).
Lemmas D′, F′, G′, H′ give the third bound

    Λ(r) ≤ Λ̄₃(r) := s √m̂(r) + √(κ_s(r)/gap₂(r)),

and the run takes about 25 seconds (output in `checks/su3/certify_blend.txt`).

## 5. Lemma F′ with exact moments at every staple sum

Sections 3 and 4 bound the moments in Lemma F′ by numerators at extreme points of rK and the
normalisation Z ≥ 1. This section evaluates Lemma F′ itself at every b ∈ rK, with the moments of ν_b
computed from the character expansion of the Haar measure and bounded rigorously, at r = 45/32.

**Lemma K (orbits).** (i) For V, W ∈ SU(3) and every h, L_h(VbW) = L_{V⁻¹hW⁻¹}(b) and L_h(b̄) = L_{h̄}(b),
so sup_{‖h‖=1} L_h is constant on the orbits of b ↦ VbW and b ↦ b̄. (ii) Every b ∈ C^{3×3} lies in the orbit
of a matrix diag(a_1, a_2, z_3) with a_1 ≥ a_2 ≥ |z_3| and Im z_3 ≥ 0. (iii) A diagonal matrix lies in K if
and only if it lies in conv T, where T is the diagonal torus of SU(3).

Proof. (i) If U has law ν_{VbW}, then U′ = WUV has law ν_b, because Re Tr(UVbW) = Re Tr(WUVb) and Haar
measure is invariant on both sides. The map U ↦ WUV preserves d, so it maps 1-Lipschitz functions to
1-Lipschitz functions, and φ_h(U) = Re Tr(W⁻¹U′V⁻¹h) = φ_{V⁻¹hW⁻¹}(U′) with ‖V⁻¹hW⁻¹‖ = ‖h‖. For b̄ use
U ↦ Ū, which preserves Haar measure and d, with Re Tr(Ūb̄) = Re Tr(Ub) and φ_h(Ū) = φ_{h̄}(U).
(ii) Write b = V_1 Σ W_1 with V_1, W_1 unitary and Σ = diag(σ_1, σ_2, σ_3), σ_1 ≥ σ_2 ≥ σ_3 ≥ 0. With
D = diag(1, 1, det V_1) and D′ = diag(1, 1, det W_1), V = V_1D⁻¹ and W = D′⁻¹W_1 lie in SU(3) and
b = V diag(σ_1, σ_2, σ_3 det V_1 det W_1) W; a complex conjugation makes Im z_3 ≥ 0.
(iii) conv T ⊂ K. Let c be diagonal with nonzero entries of distinct moduli, and let U maximise
Re Tr(cU) over SU(3). For X ∈ su(3), Re Tr(cUX) = d/dt Re Tr(cUe^{tX}) at t = 0 vanishes, so the
anti-Hermitian part of cU is a multiple of the identity: cU = H + iλI with H Hermitian and λ real. Then
H² + λ²I = (cU)(cU)* = cc* is diagonal with distinct entries, so H, which commutes with it, is diagonal,
and U = c⁻¹(H + iλI) is diagonal, that is U ∈ T. Hence max_{SU(3)} Re Tr(cU) = max_T Re Tr(ct) for such c,
and by continuity for every diagonal c. If a diagonal d is not in conv T, separation in the real space
of diagonal matrices gives a diagonal c with Re Tr(cd) > max_T Re Tr(ct) = max_{SU(3)} Re Tr(cU), which is
max_K Re Tr(c ·) because the functional is linear; so d ∉ K. ∎

Since every b ∈ rK has operator norm at most r,

    Λ(r) = sup { sup_{‖h‖=1} L_h(b) : b ∈ 𝓡_r },
    𝓡_r = { diag(a_1, a_2, x_3 + iy_3) ∈ r conv T : a_1 ≥ a_2 ≥ |x_3 + iy_3|, y_3 ≥ 0 },

a subset of the box [0, r]² × [−r, r] × [0, r] in the coordinates (a_1, a_2, x_3, y_3), and b is linear in
these coordinates: b = Σ_k b_k F_k with F = (e_11, e_22, e_33, i e_33).

**Lemma L (moments at a diagonal staple sum).** For D, D′ ∈ C^{3×3} put
Z(D, D′) = E_H exp(½Tr(UD) + ½Tr(ŪD′)). Then

    Z(D, D′) = Σ_{p ≥ q ≥ 0} s_{(p,q)}(DD′ᵀ) G_{p,q}(det D) G_{p,q}(det D′) / d_{p,q},
    G_{p,q}(z) = Σ_{k ≥ 0} c_{(p+k, q+k, k)} z^k,   c_λ = 2^{−|λ|}/∏ hooks(λ),
    d_{p,q} = (p − q + 1)(q + 1)(p + 2)/2,

with s_{(p,q)} the Schur polynomial, and the moments of degree j in U and k in Ū are

    E_H[u^α ū^β e^{Re Tr(Ub)}] = 2^{|α|+|β|} ∂_D^α ∂_{D′}^β Z(b, b̄),

the entry u_{ij} pairing with D_{ji} and ū_{ij} with D′_{ji}.

Proof. e^{½Tr X} = Σ_n Tr(X)^n/(2^n n!) and Tr(X)^n = Σ_{λ ⊢ n} f^λ s_λ(X) with f^λ = n!/∏ hooks(λ), so
e^{½Tr X} = Σ_λ c_λ s_λ(X) over partitions with at most three parts. For U ∈ SU(3),
s_{λ+(k,k,k)}(UD) = det(D)^k s_λ(UD), so e^{½Tr(UD)} = Σ_{p≥q} G_{p,q}(det D) χ_{p,q}(UD), where χ_{p,q} is the
character of the irreducible representation with highest weight (p, q, 0); these are pairwise
inequivalent on SU(3) and have dimension d_{p,q}. Since Ū = (U⁻¹)ᵀ and characters are invariant under
transposition, χ_{p,q}(ŪD′) = χ_{p,q}(U⁻¹D′ᵀ), and the convolution formula
E_H[χ_ρ(UA) χ_{ρ′}(U⁻¹B)] = δ_{ρρ′} χ_ρ(AB)/d_ρ gives the expansion. The bound below shows absolute
convergence, so the sums and E_H commute; the moments follow by differentiating
exp(½Tr(UD) + ½Tr(ŪD′)) under E_H, since Re Tr(Ub) = ½Tr(Ub) + ½Tr(Ūb̄). ∎

*Truncation.* The computation keeps p ≤ P = 24 and k ≤ K = 8 in both factors G. The truncated sum is
E_H[A_t(UD) B_t(ŪD′)], where A_t is the part of Σ_λ c_λ s_λ(UD) over the partitions with λ_1 − λ_3 ≤ P and
λ_3 ≤ K, and B_t likewise. Every omitted partition has |λ| > P, because 3(K + 1) > P. As |s_λ(X)| ≤ dim(λ)
‖X‖^{|λ|} (operator norm) and Σ_{λ ⊢ n} f^λ dim(λ) = 3^n (Schur–Weyl), the omitted part is at most
T_P(‖D‖) = Σ_{n > P} (3‖D‖/2)^n/n! in modulus. On the polydisc |D_{ij} − b_{ij}| ≤ ρ, ‖D‖ ≤ ‖b‖ + 3ρ, and
Cauchy's estimate bounds the α-th derivative of the omitted part at b by α! ρ^{−|α|} T_P(‖b‖ + 3ρ); the
derivatives of the full exponential at b are at most 2^{−|α|} e^{V/2}, V = Σ_k |b_k|, because
|Tr(Ub)| ≤ Σ_k |b_k| for diagonal b. With ε_α = 2^{|α|} α! ρ^{−|α|} T_P(‖b‖ + 3ρ), the truncated moment
differs from the true one by at most ε_α e^{V/2} + (e^{V/2} + ε_α) ε_β; ρ is chosen to minimise ρ^{−4}T_P.

*Rounding.* The series is a Taylor polynomial of total degree 4 in the 36 real directions of
(D − b, D′ − b̄), computed in double precision. Every operation is repeated on absolute values and
carries an operation count n, so that the computed coefficient differs from the exact one by at most
γ_n = nu/(1 − nu), u = 2⁻⁵³, times the absolute-value coefficient (running error analysis). The
polynomials φ_h, ∇φ_h·∇φ_{h′} and G_{F_k, h} are exact Gaussian rationals; normalisation by Z, the
quadratic forms and their error matrices are formed with the same bookkeeping. At b = diag(1.36, 1.02,
1.02) the total error of Z is 1.3·10⁻¹², and of the 18 × 18 forms below 5·10⁻¹⁰.

*Block structure.* For diagonal b, ν_b is invariant under U ↦ tUt⁻¹, t ∈ T, so every quadratic form
below is invariant under h ↦ t⁻¹ht and splits into four blocks: the diagonal entries of h (six real
coordinates) and the three pairs {h_jk, h_kj}. Means of φ_h and of G_{F_k,h} vanish off the diagonal
block. The off-block entries of the exact forms are zero, so the computation keeps the blocks with
their error bounds and drops the rest.

**Lemma M (interpolation on a box).** Let Q = ∏_k [l_k, l_k + w_k] in the coordinates (a_1, a_2, x_3, y_3),
with corners b_i, and σ = w_1² + w_2² + max(w_3², w_4²). For b ∈ Q let ω_i(b) ≥ 0 be the multilinear
interpolation weights, Σ_i ω_i(b) = 1. Then for every b ∈ Q:

(i) Z(b) ≥ e^{−σ/8} Σ_i ω_i(b) Z(b_i);
(ii) ν_b(f) ≤ e^{σ/8} max_i ν_{b_i}(f) for every f ≥ 0 that does not depend on b;
(iii) if Y_b = (A + Σ_k (b_k − β_k) B_k)² with functions A, B_k of U and a fixed point β, then

    ν_b(Y_b) ≤ e^{σ/8} ( max_i ν_{b_i}(Y_{b_i}) + Σ_k (w_k²/4) max_i Z(b_i) ν_{b_i}(B_k²) / min_i Z(b_i) ).

Proof. E(b) = e^{Re Tr(Ub)} = ∏_k e^{c_k b_k} with c = (Re u_11, Re u_22, Re u_33, −Im u_33), so |c_k| ≤ 1
and c_3² + c_4² ≤ 1. In each coordinate the chord of e^{c_k x} over [l_k, l_k + w_k] lies above the
function and, by Hoeffding's lemma for a Bernoulli variable, below e^{c_k x} e^{c_k² w_k²/8}. Multiplying,

    E(b) ≤ Σ_i ω_i(b) E(b_i) ≤ E(b) e^{Σ_k c_k² w_k²/8} ≤ E(b) e^{σ/8}.                       (3)

Integrating the right inequality gives (i). For (ii), ∫ f E(b) ≤ Σ_i ω_i ∫ f E(b_i) ≤ max_i ν_{b_i}(f)
Σ_i ω_i Z(b_i), and (i). For (iii) let F(b) = Y_b E(b) and let I_k be linear interpolation in the k-th
coordinate, so that the multilinear interpolant is I = I_1I_2I_3I_4 and
F − IF = Σ_k I_1⋯I_{k−1}(F − I_kF). Along the k-th coordinate F = L² e^{c_k x} times a positive factor, with L
affine of slope B_k, and ∂_k²F = ((c_k L + 2B_k)² − 2B_k²) E ≥ −2B_k² E. With the Green kernel G ≥ 0 of the
interval, ∫ G = t(1 − t)w_k²/2 ≤ w_k²/8,

    F − I_kF = ∫ G (−∂_k²F) ≤ 2B_k² ∫ G E ≤ (w_k²/4) B_k² Σ_i ω_i^{(k)} E(b_i)

by the left inequality of (3), with weights ω_i^{(k)} ≥ 0 of sum one that do not depend on U. The operators
I_j have nonnegative weights of sum one, so F(b) ≤ Σ_i ω_i F(b_i) + Σ_k (w_k²/4) B_k² Σ_i ω_i^{(k)} E(b_i).
Integrating over U and dividing by Z(b), bounded below by (i), gives (iii). ∎

**The bound on a box.** Fix s = 3j/8192 with an integer j, so that c_1 = 1 − sλ_F = 1 − j/1024. With
G_{b,h} = Σ_k b_k G_{F_k,h} and any constant κ, Var_{ν_b}(X_s) ≤ ν_b[(X_s − κ)²]. Take κ = κ(b)·h with κ(b)
affine in b, κ(b) = κ_0 + Σ_k (b_k − β_k) g_k, its coefficients on the dyadic grid 2⁻³⁰ℤ fitted to the means
at the corners (β the centre of Q). Then X_s − κ = A + Σ_k (b_k − β_k)B_k with A = c_1φ_h + sG_{β,h} − κ_0·h
and B_k = sG_{F_k,h} − g_k·h, and Lemmas F′ and M give, for all b ∈ Q,

    sup_h L_h(b) ≤ B_Q := s √(e^{σ/8} max_i λ_M(b_i)) + √( e^{σ/8} (max_i λ_V(b_i) + e_N/min_i Z(b_i)) / gap₂(r) ),
    e_N = Σ_k (w_k²/4) max_i Z(b_i) λ_k(b_i),

where λ_M(b_i), λ_V(b_i), λ_k(b_i) are the largest eigenvalues of the forms h ↦ ν_{b_i}|∇φ_h|²,
h ↦ ν_{b_i}[(c_1φ_h + sG_{b_i,h} − κ(b_i)·h)²] and h ↦ ν_{b_i}[B_k²]. Each is bounded above in ball
arithmetic: an LDLᵀ factorisation of tI − A with certified positive pivots, A the computed symmetric
block, plus the Frobenius norm of its error matrix (Weyl's inequality). The integer j is chosen per box
by a floating-point search; the choice does not affect validity.

**Discarded boxes.** A box is dropped when it has no point with a_1 ≥ a_2 ≥ |x_3 + iy_3|, or when a linear
functional ℓ(b) = Re Tr(c*b) separates it from r conv T: min_Q ℓ > r max_T ℓ, where max_T ℓ is bounded
above by its maximum over a 512 × 512 grid of the torus plus ½H(δ²/2), δ = 2π/512,
H = |c_1| + |c_2| + 2|c_3| a bound for the Hessian on the torus. The functional c is the dual solution of a
linear program for the gauge of the centre with respect to a 72 × 72 grid of T; only the final inequality
is used, and it is checked with a relative margin 10⁻¹².

**Result at r = 45/32.** Boxes of side r/8 cover [0, r]² × [−r, r] × [0, r]: 8192 boxes, of which 6762 have no
point with a_1 ≥ a_2 ≥ |z_3| and 589 are separated from r conv T. Of the other 841, 797 satisfy
2rB_Q < 1; the remaining 44 are cut into sixteen boxes of side r/16 each, and of these 704 boxes 68 survive
the two tests and all satisfy 2rB_Q < 1. The 865 certified boxes have j between 776 and 1020; the largest
value of 2rB_Q is 0.999657, on the box [1.0547, 1.2305] × [0.7031, 0.8789] × [0.3516, 0.5273] ×
[0.7031, 0.8789] in (a_1, a_2, x_3, y_3), and the median is 0.966. The rounding and truncation errors of
the moments at every corner enter the bound B_Q. On a 5⁴ grid of the worst box, 55 points lie in rK, where 2r sup_h
of the floating-point bound of Lemma F′ is at most 0.972; the largest value on the grid, 0.992, is at a
corner outside rK (`checks/su3/exact/worstbox.py`, output `worstbox.txt`). `checks/su3/exact/cover.py` runs in about half an hour on four cores and writes
`cover_log.txt` and `cover_boxes.txt` (every box with its j and bound). Hence

    2rΛ(r) ≤ 0.999657 at r = 45/32.


## 6. The theorem

Λ̄ is increasing on [0, 2√μ_0/g): gap decreases, J and m increase. For r′ ≤ r, r′K ⊂ rK, because K is
convex and contains 0 (the mean of I, ωI, ω²I with ω³ = 1), so Λ is nondecreasing and the bound at
one r holds below it. In ball arithmetic (g² ≤ 2.06548; for Λ̄ upper Riemann sums for J with 4000 steps):

| bound | r = 2β_W | ingredients | 2r Λ̄ | rate −log(2rΛ̄) |
|---|---|---|---|---|
| Λ̄ (Lemmas E, F), `certify_su3.py` | 0.730 | gap 1.07864 | 0.994249 | 0.00577 |
| Λ̄ (Lemmas E, F), `certify_su3.py` | 0.732 | gap 1.07504 | 0.999841 | 0.00016 |
| Λ̄₂ (Lemmas F–J), `certify_tilted.py` | 0.98 | m̃ < 0.7123, κ̃ < 0.1830, gap 0.661306 | 0.999236 | 0.00076 |
| Λ̄₃ (Lemmas D′–H′), `certify_blend.py` | 1.14 | s = 0.322, m̂ < 0.64543, κ_s < 0.0508, gap₂ 1.64094 | 0.990976 | 0.00907 |
| Lemma F′ at every b (Section 5), `exact/cover.py` | 45/32 | 865 boxes, gap₂ 193/128 | 0.999657 | 0.00034 |

So Λ(45/32) ≤ 0.999657/(45/16) < 0.355435, and for β_W ≤ 45/64, α ≤ 4β_W Λ(45/32) < 1.42174 β_W ≤ 0.999657.

**Theorem.** Four-dimensional SU(3) lattice Yang–Mills with 0 ≤ β_W ≤ 45/64 = 0.703125: the forest-gauge
specification γ^F has a unique Gibbs measure; the free-boundary limit ω of the Wilson measures exists
for every gauge-invariant local observable along every exhausting sequence of boxes, is translation
and hypercubic invariant and reflection positive, has exponentially decaying truncated correlations,
and has a transfer-matrix gap m ≥ −log(1.4218 β_W) > 0. For β_W ≤ 0.366 also
m ≥ −log(4β_W Λ̄(2β_W)), and since Λ̄ → 1/4, m ≥ −log β_W + o(1) as β_W → 0.

This passes 9/29 = 0.310 (gauge-fixed curvature) and 9/32 = 0.281 (Theorem B).

## 7. Where the estimate loses

Section 4 left three losses, measured at r = 1.14, s = 0.322 with the exact moments and a Galerkin
computation (`checks/su3/losses.py`, `checks/su3/gap_lemma.py`); Section 5 removes the first two by
evaluating Lemma F′ with the exact moments of ν_b at b itself.

1. Normalisation. At b = rI, Z = 1.4632; Lemma G′ gives m̂ = 0.6454 (m̃ = 0.8051), against the normalised
   moment max_h ν_{rI}|∇φ_h|² = 0.5502.
2. The variance in Lemma H′. The decoupled form has its maximum 0.04982 at the central elements ωI, ω²I
   of T (0.04950 at W = I), while the largest variance of X_s under ν_{rI} is 0.01825. The factor 2.7
   comes from Z ≥ 1 and from bounding a variance by a second moment about the Haar mean; shifting
   X_s by a multiple of Re Tr(W h*), the only shift compatible with the reductions, does not lower the
   maximum.
3. Lemma D′. gap₂(1.14) = 1.641, against Galerkin values 2.29 at b = rI, 2.31 at b = r diag(0, 0, 1)
   and 2.30 at b = (r/2)e^{iπ/3} I, upper bounds for the true gap (polynomials of degree two; the same
   computation gives 2.58 at b = 0, where the gap is 8/3).

The route bounds sup_h √(ν_b|∇ψ|²), the H^{-1}(ν_b) norm of φ_h − ν_bφ_h, and needs it at every
b ∈ rK. Its Galerkin values (polynomials of degree two, `checks/su3/ceiling.py`) are not largest at the
extreme points r SU(3): along ρI they rise to 0.317 at ρ ≈ 1.1 and then fall (0.305 at ρ = 1.55), while
interior points of rK such as r diag(0, 0, 1) = r(I + diag(−1, −1, 1))/2 and (r/2)e^{iπ/3} I = r(I + ωI)/2
keep growing (0.3245 and 0.3255 at r = 1.55). Twice r times the largest value found is 0.79 at r = 1.25,
0.89 at r = 1.4, 0.97 at r = 1.5 and 1.009 at r = 1.55. These are lower bounds at a few matrices only, so
single-link Dobrushin estimated through this norm stops at or below β_W ≈ 0.77 even with exact constants.

With the moments exact, what is left is Lemma D′ and the limit of Lemma F′ itself. In floating point
(`checks/su3/exact/scan.py`, `face.py`, outputs `scan_1.4.txt`, `face.txt`) the largest values of
2r sup_h of the bound of Lemma F′ found on 𝓡_r lie on the face diag(r, ρ, ρ) with ρ between 0.75r and 0.8r:
0.9458 at r = 1.36, 0.9767 at r = 1.4 and 1.0080 at r = 1.44, so Lemma F′ with Lemma D′ stops near
r = 1.43, β_W ≈ 0.715. Between that value and 0.77 the room is in the gap bound (item 3); past 0.77 one
needs blocks of several links. `checks/su3/tilted_mc.py` and `checks/su3/exact/mc_check.py` compare the exact
moments of Sections 3 and 5 with Haar Monte Carlo estimates (agreement within two standard errors).
