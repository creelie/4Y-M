# Relative-chart audit and a local curvature witness (November 2026)

**Latest constructive result (proved deductions; conditional domain):**
`CANONICAL_CHART.md` gives measurable background selection, an explicit
density and weighted covariance decay for a specified product graph law
around a calm R(c). This does not match the actual physical label domain.
`COARSE_LABELS.md` adds coarse-curvature witnesses and sharper overlap losses.
(I1), D7 and D8 remain open.

**New obstruction (proved):** `ROUGH_BACKGROUND.md` gives a four-link
negative constrained-action Hessian direction at R(c) with every relative
physical link label good. Coarse-background calmness is an additional
hypothesis; relative-link goodness alone cannot imply input (a).

**Follow-up (proved deductions / open sketch):** `TOUCHING_CHART.md` constructs
a local graph solving all touching constraints around an assumed calm completion,
with finite-range solved links and an explicit disintegration density. It also
computes a conservative finite-intersection loss for input (b):
[6*2401*z(r/50)]^(|S|/2401), when 6*2401*z(r/50)<1.
This does not prove a Hessian bound or match the actual label-conditioned
domain. Relative decay and (I1) remain open; no island threshold is certified.

This continues `LARGE_FIELD.md` §3a, then input (b). **Proved** below means a mathematical deduction from the
stated lemmas of `main.tex`, not a formal verification. Nothing here proves (I1), D7 or (D).

## 1. Input (b) with a canonical, measurable reference (proved)

Use the tree gauge and the Hilbert–Schmidt chordal distance δ of the paper. Let d be the bi-invariant geodesic
distance. Write c = V^sharp(U), and let R(c) be the reconstruction defined before Theorem `pathcoercive`:
inside-block links are I and links leaving B(y) in direction μ are c_(y,μ). This is a continuous function of c.
Every two-link path P_x of R(c) is exactly c_(y,μ), so

    V^sharp(R(c)) = c.

Thus R(c) is an explicitly selected point of the full fibre for every c, without a minimisation or a choice of
gauge depending on U. It need not be calm.

**Lemma 1 (local witness).** Fix y. Let P_y be all plaquettes with a vertex in a block of {y}^+ (the collar used
in Lemma `smallframes`). Suppose δ(U_p,I) ≤ κ ≤ 10^(-8)ρ_G on P_y. Then every active link with base in B(y) obeys

    d(U_e, R_e(c)) < 50κ.

At κ = 0 read the conclusion as equality to zero.

*Proof.* For an inside-block link, Lemma `smallframes`(i) gives δ(U_e,I) ≤ 3κ, and Lemma `logbound` gives
d(U_e,I) ≤ πκ. For a link leaving the block in direction μ, put v = U_(2y+e_μ,μ). Lemma `smallframes`(i),(iv) gives

    δ(U_e,v) ≤ 21κ,
    c_(y,μ) = v exp M,       ||M|| ≤ 8πκ.

Consequently δ(U_e,c_(y,μ)) ≤ (21+8π)κ. This is below ρ_G ≤ 1 under the stated hypothesis, so the principal log
lies in the Lie algebra and supplies a geodesic of length at most (π/3)(21+8π)κ < 50κ. The distance is no larger
than that length. The tree links agree identically. ∎

**Corollary 1 (input (b), corrected constant).** For 0 < r ≤ 50·10^(-8)ρ_G, define a bad block by

    b_y(U) = {some active e based in B(y) has d(U_e,R_e(c)) > r}.

Then

    b_y(U) ⊂ union over p in P_y of {δ(U_p,I) > r/50}.

*Proof.* Apply Lemma 1 at κ = r/50 and take its contrapositive. ∎

This is a fixed-threshold, fixed-size local witness, independent of the number of blocks. Although R(c) depends
on U through c before disintegration, the centres are fixed on each full fibre. The statement holds pointwise,
so it remains true after averaging over c. The witness event is a function of local plaquette defects.

**Corollary 2 (two calm completions).** If U and U_0 have the same averages c and their plaquette defects on P_y
are bounded by κ and κ_0 respectively, with each below 10^(-8)ρ_G, then

    d(U_e,(U_0)_e) ≤ 50(κ+κ_0)

for every active e based in B(y). This follows by the triangle inequality through R_e(c). In particular, a
reference of curvature κ_0 costs an additive error; it cannot be silently omitted from input (b).

**Sketch (finite intersections; constants not computed).** Choose blocks from a specified set S greedily so
their plaquette collars P_y are disjoint. Since these collars have bounded radius, at least |S|/D blocks can be
chosen for a fixed lattice constant D. For each chosen block take the union over at most M = sup_y |P_y|
witness plaquettes. If a torus chessboard estimate gives a bound p(ρ)^n for n distinct witnesses, the result is
at most [M p(r/50)]^(|S|/D), when M p(r/50) < 1. Any restrictions on witness orientations/reflection classes
must also be included in D. This is not the single-witness ledger q(r/4); the extraction and union losses must
be computed before quoting a new island threshold.

## 2. What the relative chart actually requires

**Proved (reading of the definitions).** The product coordinates in `chartdecay` are the independent entries
η_e of a tangent vector in ker DC, for e in L'. Physical links are

    U_e(η) = (U_0)_e exp((Sη)_e + h_e(Sη)).

Thus Φ(Ω_r) is a product in η, not in physical links. Lemma `kernels` and the proof of `chartdecay` control the
nonlocal correction h by decaying kernels; they do not identify physical link balls with coordinate balls.
In particular, imposing physical link balls pulls back to nonlinear constraints involving h and the section S.
Convexity of that pullback is not proved by near-isometry.

**Proved (limitation of a sandwich).** Inclusions D_- ⊂ D ⊂ D_+ alone do not transfer a Poincaré inequality or
covariance decay from D_- and D_+ to D. Already in one dimension take D_- = (-1,1), D_+ = (-3,3), and
D = (-1,1) union (2,5/2). A positive smooth density, even exp(-βx²), conditioned on D has positive variance for
a smooth function constant with different values on the two components and zero derivative on D. Hence it
has no finite gradient Poincaré constant. Both bounding intervals are convex. A smooth ambient function can
make its transition in the omitted interval. This example does not show that the particular gauge-theory
pullback is disconnected; it proves that the proposed sandwich argument is insufficient.

**Proved (h is not a projection onto another fibre).** In `fibrechart`(i), h is constructed after choosing a
background U_0 in the target fibre. Its constraint equation is C(U_0 exp X) = 0, where C is centred at the
averages of U_0. Therefore that theorem cannot select a point with different prescribed averages by applying
h to an arbitrary reference. Selecting a background from a conditioned exterior requires a separate
inhomogeneous constraint argument, or an explicit extension that already has all the required averages.

**Proved (a conditional-law check).** For the full disintegration given all averages c, freezing links off A_W
retains every constraint V^sharp_k = c_k whose support meets A_W, including pairs based outside W. By contrast,
the fibre F of `main.tex`, subsection “Small fields on the fibres of the cut-off average”, fixes only averages
based in W. A pair based in y-e_μ can use inside-block links of B(y) through its two-link paths. Thus the
conditional law is not automatically the law of that theorem. Adding private links for straddling pairs is
the repair proposed in `ISLANDS.md` §1, but its chart geometry and constants must then be checked. A spatial
collar alone does not remove incoming constraints at the new boundary.

**Sketch (a possible replacement).** R(c) supplies fixed centres for a physical-link split and proves input (b).
One could parameterise a window by physical non-private links, solve all touching averages for private links,
and prove a Hessian bound for the resulting density (including the constraint Jacobian). To use product
convexity, restrictions on solved private links must be either redundant on that product or separately handled.
The Hessian and decay estimates for this parametrisation are not the estimates already proved for Φ.

Alternatively define the split in a fixed tangent product chart. Then each scalar coordinate restriction is
fixed there, but the action and link observables have the decaying nonlocal tails of h; their bad-side
coefficients must be re-established in that metric. Recentring charts from the exterior changes the support
under exterior variation and requires boundary transport estimates.

## 3. Two further checks on the proposed route

**Proved (concentration is centred).** Strong convexity gives concentration about a mean or minimiser, not about
the identity under every exterior. For example on [-r,r] the density proportional to exp[-β(x-2r)²/2] has
Hessian β, but P(x > r/4) tends to one as β tends to infinity. To see this, the potential decreases throughout
the interval and its minimum is uniquely at r; the energy separation of any fixed subinterval bounded away
from r is positive. Thus the absolute plaquette rarity claimed from Herbst in `LARGE_FIELD.md` §3 needs
additional control of the conditional mean. Torus rarity via Corollary 1 does not use that claim.

**Proved (scope of the fixed-support derivative).** For an unconstrained link on a fixed positive-Haar-mass
set C, differentiation gives -Cov(f,β∂S), as in `LARGE_FIELD.md`. On a fibre, the constraint surface and its
coarea density also enter. A single physical link is generally not a free coordinate. The fixed-support
calculation alone therefore does not prove the bad-side coefficient for the constrained law μ^ω.

**Sketch (environment term remains).** Lemma 3 of `ANNEALED.md` proves a scalar recursion *if that recursion
holds*. Neither it nor torus witness rarity establishes the recursion for Cov_ω(μ^ω F,μ^ω H). Locality under
changing labels, regularity of hard indicators, and a coefficient uniform in the conditioning still need proof.

## 4. Status and next work

- **Proved:** canonical full-fibre selection R(c); the local witness at r/50; comparison of two calm
  completions with the additive reference-curvature error.
- **Proved:** the proposed domain sandwich does not by itself justify covariance transfer; the paper's h
  does not solve inhomogeneous constraints; incoming averages must be retained under full disintegration.
- **Sketch:** a physical non-private-link chart with all touching constraints is a concrete next target.
  Its deterministic Hessian and domain estimates remain open, as do constrained bad-side coefficients and the
  environment-covariance bootstrap hypothesis.
- **Computation:** the r/4 ledger in `LARGE_FIELD.md` is conditional arithmetic, not an established threshold.
  Even replacing 4 by 50 omits the finite-intersection extraction losses. No new threshold is certified here.

Input (b) has a rigorous version for a specified reference and a conservative constant. The relative chart and
the passage to (I1) remain open; this audit does not reduce the single-scale large-field problem to two completed
deterministic checks. D8 is unchanged.
