# Forest gauge fixing in four dimensions

Three results that go past the ranges of the paper in `../prd/` by fixing the
gauge on one link of every plaquette. All three are proved by hand and have
not yet been read by an outside referee. They are written up as the paper in
`../forest/`; the notes below carry the full computations and the scripts.

## The line field

Read each site x of Z^4 modulo 2 and set σ(x) by the table (index
x₀ + 2x₁ + 4x₂ + 8x₃)

    [0, 0, 2, 3, 3, 1, 2, 1, 1, 2, 1, 3, 3, 2, 0, 0].

The fixed set F consists of the forward links (x, σ(x)). It has three
properties, each checked by `decide` in
`YangMills/StrongCoupling/GaugeFix.lean` (seven theorems, axioms `propext` at
most):

1. every plaquette of Z^4 contains exactly one link of F
   (`one_fixed_per_plaquette`);
2. σ(x + e_σ(x)) = σ(x), so F is a disjoint union of straight lines, one
   through every site; in a box with free boundary it is a finite union of
   segments, hence a forest (`straight_lines`, `one_per_site`);
3. a link outside F lies in six plaquettes and has exactly twelve partners
   outside F in them, against eighteen without gauge fixing
   (`twelve_neighbours`); in two of the six plaquettes the fixed link is the
   one opposite to it (`two_opposite`).

The seventh theorem, `budget`, checks the arithmetic 4·21/10 + 2·5/4 +
4·9/10 = 29/2 of the per-link bound in `curvature.md`, and `threshold` the
counts 6·3 = 18 and 6·2 = 12.

## Results

| Note | Statement | Status |
|---|---|---|
| [`forest.md`](forest.md) | Four-dimensional SU(2): the free-boundary limit of the Wilson measure exists for β_W < 1/3, is translation and hypercubic invariant and reflection positive, has exponential decay and a transfer-matrix gap m ≥ −log(3β_W). The value 1/3 is the ceiling of single-link Dobrushin after any forest gauge fixing. | proved, not refereed |
| [`curvature.md`](curvature.md) | On the gauge-fixed directions ±Hess S₀ ≤ 29/2 at every field (against 16 for the ungauged measure), so Bakry–Émery, uniqueness, decay and, for β ≥ 0, a mass gap hold for \|β\| < 1/29 for every SU(N) and \|β\| < (N−2)/(58N) for SO(N); for SU(3) this is β_W < 9/29. The best constant h_F lies between 9 + √17 and 29/2; numerically h_F = 13.9740139454, attained by an abelian field of period two. | proved, not refereed; h_F open |
| [`su3.md`](su3.md) | Four-dimensional SU(3): single-link Dobrushin with twelve partners. The one-link measure ν_b ∝ exp(Re Tr(Ub)) has spectral gap at least 283/128 − r/2 (Bochner, the bottom 7/6 of the rough Laplacian on 1-forms and the pointwise bound Hess V + V/2 ≤ r/2), which bounds a corrector for the Lipschitz rate (blended, then of degree two); the gradient moment and the variance of the corrector are bounded by moments of the Haar measure tilted by exp(r Re Tr U), computed exactly (Schur–Weyl, Weyl integration, Bessel sums), and the moment is normalised through E_H[\|∇φ_h\|² U] = adj(h)/9. The corrector is then evaluated at every staple sum: the staple sums reduce to diagonal matrices in r conv T, where the moments of the one-link measure come from the character expansion of the Haar measure with rigorous truncation and rounding bounds, and 865 boxes cover the diagonal staple sums at β_W = 45/64. A corrector of degree two in the matrix entries, with exact moments of degree six, then gives a cover by 2341 boxes at β_W = 3/4. Uniqueness, decay, reflection positivity and a gap m ≥ −log(1.3332 β_W) follow for β_W ≤ 3/4, past 9/29 = 0.310. | proved, not refereed; constants in ball arithmetic |

## Checks

The notes do not rely on the scripts; they confirm the finite steps.

| Script | Checks |
|---|---|
| `checks/verify_gauge.py` | on a box of side 6 with no periodicity: one fixed link per plaquette, F a forest (union–find), twelve partners (eighteen ungauged, fourteen in temporal gauge) |
| `checks/onefixed.py`, `checks/enum2.py` | SAT search for a period-two set with one fixed link per plaquette that lifts to a forest, and enumeration of all such sets |
| `checks/curvature/exact4.py`, `minors_table.py` | the local lemma of `curvature.md` in exact arithmetic (sympy) for every position of the frozen link |
| `checks/curvature/local_un.py` | the local lemma under adversarial search over U(2), U(3), U(4) |
| `checks/curvature/ident_exact.py`, `ident.py` | the compressed curl–curl operator at the identity field on tori of side 2 and 4 |
| `checks/curvature/verify_global.py` | the global identity and the bound 29/2 on tori, SU(N) and U(N) |
| `checks/curvature/cert.py`, `sym.py`, `pair.py` | the local certificates and the semidefinite search for their weights |
| `checks/curvature/thess.py`, `fast.py`, `analyse_best.py`, `enum_sig.py` | the gradient search for the extremal field (best value 13.9740139454, field in `extremal_field_L4.pt`) |
| `checks/su3/certify_su3.py` | the bottom 7/6 of the rough Laplacian (Casimirs and dimensions for p, q ≤ 40), the constant g² ≤ 2.06548 by branch and bound, and 2rΛ̄(r) < 1 at r = 0.730, 0.732, all in ball arithmetic (python-flint) |
| `checks/su3/tilted_moments.py`, `certify_tilted.py` | exact tilted Haar moments of degree four; m̃(0.98) < 0.7123, κ̃(0.98) < 0.1830 over 8464 squares of the torus, and 2rΛ̄₂ = 0.99924 < 1 at r = 0.98, in ball arithmetic (output in `certify_tilted.txt`, about 20 s) |
| `checks/su3/certify_blend.py` | the identity E_H[\|∇φ_h\|² U] = adj(h)/9 exactly; m̃(1.14) < 0.8051, the blended variance κ_s(1.14) < 0.0508 at s = 0.322 over 9778 squares of the torus, gap₂(1.14) = 283/128 − 0.57 and 2rΛ̄₃ = 0.99098 < 1 at r = 1.14, in ball arithmetic (output in `certify_blend.txt`, about 25 s) |
| `checks/su3/gap_lemma.py` | the pointwise bound Re Tr(WX²) + Re Tr(W)/2 ≤ 1/2 of Lemma D′ and its closed form on the torus and at random W; Galerkin gaps at r = 1.14 (`gap_lemma.txt`) |
| `checks/su3/tilted_mc.py`, `losses.py` | Monte Carlo comparison of the exact moments; the three remaining losses at r = 1.14 (`losses.txt`) |
| `checks/su3/exact/series.py`, `rseries.py`, `epoly.py`, `corner.py` | moments of the one-link measure at diagonal staple sums by the character expansion of the Haar measure, with truncation (Cauchy) and rounding (running error analysis) bounds; the exact polynomials φ_h, ∇φ_h·∇φ_{h′} and G_{F_k,h}; the normalised quadratic forms at a corner with their error matrices |
| `checks/su3/exact/box.py`, `cover.py` | the bound of Section 5 of `su3.md` on a box of diagonal staple sums (multilinear interpolation, Hoeffding, LDLᵀ in ball arithmetic) and the cover at r = 45/32: 865 boxes, largest bound 0.999657 (`cover_log.txt`, `cover_boxes.txt`, about half an hour on four cores) |
| `checks/su3/deg2/ring6.py`, `rs6z.py`, `tsmul.c`, `cmul.py` | moments of degree six of the one-link measure at diagonal staple sums (character expansion with truncation λ_1 − λ_3 ≤ 34, λ_3 ≤ 11, computed in the subring of torus weight zero; C kernel for the products), with truncation and rounding bounds |
| `checks/su3/deg2/polys.py`, `box2.py`, `cover2.py` | the exact polynomials of the corrector of degree two (Gaussian-integer coefficients, partial sums below 2⁵⁰), the bound of Section 6 of `su3.md` on a box, and the cover at r = 3/2: 2341 boxes, largest bound 0.999887 (`r3_2/cover_log.txt`, `r3_2/cover_boxes.txt`, about an hour on four cores) |
| `checks/su3/deg2/face2.py` | floating-point values of both correctors on the face diag(r, ρ, ρ) at r = 1.44, 1.50, 1.54 (`face2.txt`) |
| `checks/su3/exact/fast.py`, `scan.py`, `face.py`, `worstbox.py`, `mc_check.py` | floating-point values of the bound of Lemma F′ with exact moments over the staple sums (`scan_1.4.txt`, `face.txt`, `worstbox.txt`), and a Monte Carlo check of the moments (`mc_check.txt`) |
| `checks/su3/ceiling.py` | Galerkin lower bounds for the H^{-1} norm that the route bounds, at matrices of rK for r = 1.25 to 1.55 (`ceiling.txt`) |
| `checks/su3/sanity.py`, `lip.py` | Monte Carlo comparison of the gap and gradient-moment bounds with Galerkin values |
| `checks/su3/route1.py`, `scan2.py` | Monte Carlo estimates of the true single-link rate (`route1.txt`, `scan2.txt`) |
| `checks/groups/ricci.py` | κ_G = h^∨/(2ℓ_ρ) for every compact simple group in its least faithful representation, from the Killing form and the Dynkin index, and the ranges N\|β\| < 2κ_G/29 in β_W (all closed forms exact) |
| `checks/curvature/abelian.py` | the same maximum over abelian fields: 64 angles in one U(1) subgroup, value 13.9740139454 (`abelian.txt`, angles in `abelian_extremal_L4.pt`) |

The four open problems of the paper, with what is known of each, are listed
in [`../OPEN.md`](../OPEN.md).
