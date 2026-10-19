# Four-dimensional SU(3) in the forest gauge: β_W ≤ 0.49

Status: proved by hand; constants certified in ball arithmetic (`checks/su3/certify_su3.py` for Lemmas B–E,
`checks/su3/certify_tilted.py` for Lemmas G–J, python-flint 0.9.0); not yet read by an outside referee. Ranges known before for four-dimensional SU(3): β_W < 9/32 = 0.281
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

For 0 ≤ r ≤ 1, I_m(r) ≤ (r/2)^m e^{r²/4}/m!, which decreases in m, so the terms |k| > 40 of S(n) add at
most 4[(1/2)^{29} e^{1/4}/29!]³. At r = 0 the computation returns Q = (4/9)·Id, the value of Lemma A(4).

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

## 4. The theorem

Λ̄ is increasing on [0, 2√μ_0/g): gap decreases, J and m increase. For r′ ≤ r, r′K ⊂ rK, because K is
convex and contains 0 (the mean of I, ωI, ω²I with ω³ = 1), so Λ is nondecreasing and the bound at
one r holds below it. In ball arithmetic (g² ≤ 2.06548; for Λ̄ upper Riemann sums for J with 4000 steps):

| bound | r = 2β_W | ingredients | 2r Λ̄ | rate −log(2rΛ̄) |
|---|---|---|---|---|
| Λ̄ (Lemmas E, F), `certify_su3.py` | 0.730 | gap 1.07864 | 0.994249 | 0.00577 |
| Λ̄ (Lemmas E, F), `certify_su3.py` | 0.732 | gap 1.07504 | 0.999841 | 0.00016 |
| Λ̄₂ (Lemmas F–J), `certify_tilted.py` | 0.98 | m̃ < 0.7123, κ̃ < 0.1830, gap 0.661306 | 0.999236 | 0.00076 |

So Λ(0.98) ≤ Λ̄₂(0.98) < 0.50982, and α = 2rΛ(r) ≤ 4β_W · 0.50982 < 1 for β_W ≤ 0.49.

**Theorem.** Four-dimensional SU(3) lattice Yang–Mills with 0 ≤ β_W ≤ 0.49: the forest-gauge
specification γ^F has a unique Gibbs measure; the free-boundary limit ω of the Wilson measures exists
for every gauge-invariant local observable along every exhausting sequence of boxes, is translation
and hypercubic invariant and reflection positive, has exponentially decaying truncated correlations,
and has a transfer-matrix gap m ≥ −log(2.0393 β_W) > 0. For β_W ≤ 0.366 also
m ≥ −log(4β_W Λ̄(2β_W)), and since Λ̄ → 1/4, m ≥ −log β_W + o(1) as β_W → 0.

This passes 9/29 = 0.310 (gauge-fixed curvature) and 9/32 = 0.281 (Theorem B).

## 5. Where the estimate loses

Three losses remain, measured at r = 0.98 with the exact moments and a Galerkin computation
(`checks/su3/losses.py`):

1. Z(b) ≥ 1 in Lemmas G and H. At b = rI, Z = 1.3185, and the normalised moment ν_b|∇φ_h|² is at most
   0.540 against m̃ = 0.712.
2. Decoupling a = b/r from b in Lemma H. The maximum, about 0.183, of λ_max N(θ) is attained away from W = I; at
   W = I, the alignment of the extreme points b = rW′, the bound is 0.122 and the true variance of Γ/r
   under ν_{rI} is 0.084.
3. Lemma D. Its value 0.661 compares with the Galerkin value 2.29 at b = rI, an upper bound for the true
   gap (polynomials of degree two; the same computation gives 2.58 at b = 0, where the gap is 8/3).

With the true moments in Lemma F and the gap of Lemma D, a Monte Carlo evaluation gives 2rΛ < 1 up to
r ≈ 1.07 (β_W ≈ 0.53, `checks/su3/route1.py`); a certified gap bound near the true one would go further.
`checks/su3/tilted_mc.py` compares the exact tilted moments with a Haar Monte Carlo estimate.
