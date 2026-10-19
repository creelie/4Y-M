# Four-dimensional SU(3) in the forest gauge: β_W ≤ 0.366

Status: proved by hand; constants certified in ball arithmetic (`checks/su3/certify_su3.py`, python-flint 0.9.0);
not yet read by an outside referee. Ranges known before for four-dimensional SU(3): β_W < 9/32 = 0.281
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

**Lemma F (corrector).** For b ∈ rK and ‖h‖ = 1, L_h(b) ≤ (3/8) √(ν_b|∇φ_h|²) (1 + g r/√gap(r)).

Proof. Let −Lψ = φ_h − ν_bφ_h. For g Lipschitz, Cov(g, φ_h) = ν_b⟨∇g, ∇ψ⟩ ≤ ν_b|∇ψ|, as |∇g| ≤ 1
a.e. Write ψ = φ_h/λ_F + χ. Since −L(φ_h/λ_F) = φ_h − Γ/λ_F with Γ = ⟨∇V, ∇φ_h⟩, χ solves
−Lχ = R := Γ/λ_F − ν_b(Γ)/λ_F. Then ν_b|∇ψ| ≤ √(ν_b|∇φ_h|²)/λ_F + √(ν_b|∇χ|²), and
ν_b|∇χ|² = ν_b(Rχ) ≤ Var(R)/gap(r) ≤ ν_b(Γ²)/(λ_F² gap(r)) ≤ g²r² ν_b|∇φ_h|²/(λ_F² gap(r)), because
|Γ| ≤ |∇V||∇φ_h| ≤ g r |∇φ_h|. With 1/λ_F = 3/8 this is the claim. ∎

Hence

    Λ(r) ≤ Λ̄(r) := (3/8) √m(r) (1 + g r/√gap(r)).

At r = 0, Λ̄(0) = 1/4, which is the Cauchy–Schwarz value 1/√(2(N² − 1)); the true value of L at b = 0
lies between 1/6 and 1/4.

## 3. The theorem

Λ̄ is increasing on [0, 2√μ_0/g): gap decreases, J and m increase. In ball arithmetic
(`checks/su3/certify_su3.py`, part C, upper Riemann sums for J with 4000 steps, g² ≤ 2.06548):

| r = 2β_W | gap(r) | 2r Λ̄(r) | rate −log(2rΛ̄) |
|---|---|---|---|
| 0.730 | 1.07864 | 0.994249 | 0.00577 |
| 0.732 | 1.07504 | 0.999841 | 0.00016 |

**Theorem.** Four-dimensional SU(3) lattice Yang–Mills with 0 ≤ β_W ≤ 0.366: the forest-gauge
specification γ^F has a unique Gibbs measure; the free-boundary limit ω of the Wilson measures exists
for every gauge-invariant local observable along every exhausting sequence of boxes, is translation
and hypercubic invariant and reflection positive, has exponentially decaying truncated correlations,
and has a transfer-matrix gap m ≥ −log(4β_W Λ̄(2β_W)) > 0. For small β_W, Λ̄ → 1/4 and
m ≥ −log β_W + o(1).

This passes 9/29 = 0.310 (gauge-fixed curvature) and 9/32 = 0.281 (Theorem B).

## 4. Where the estimate loses

Monte Carlo at r = 0.732 (10^5 Haar points, `checks/su3/sanity.py`): the true gap of ν_b is about 2.32 for
b = rI (Galerkin upper bound) against 1.075 in Lemma D, and ν_b|∇φ_h|² is at most about 0.52 against
m = 0.82 in Lemma E. With the true moments in Lemma F and the gap of Lemma D, a Monte Carlo
evaluation gives 2rΛ < 1 up to r ≈ 1.07 (β_W ≈ 0.53, `checks/su3/route1.py`). Certifying those moments uniformly
over b ∈ rK is open.
