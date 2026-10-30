# The four open problems of the strong-coupling paper

The table `tab:gaps` of `prd/main.tex` lists four problems that separate the
results of the paper from the Clay mass gap. Normalisation as in the paper:
weight exp(Nβ Σ_p Re Tr Q_p), β_W = N²β. A problem counts as solved only when
every case below it is proved.

## 1. Intermediate fixed coupling

Uniqueness of the Gibbs measure, decay of correlations and a transfer-matrix
gap at a fixed lattice spacing, for every coupling up to the region where
Monte Carlo puts the crossover (β_W of order 2.2 for SU(2)).

| Range | Status | Where |
|---|---|---|
| SU(2), β_W < 2/9 | proved; the ceiling of every single-link Wasserstein–Dobrushin argument on the ungauged measure | `prd/`, Theorems E and F |
| every SU(N), \|β\| < 1/(8d) | proved by curvature; no curvature bound uniform over configurations goes further on the ungauged measure | `prd/`, Theorems B and C |
| SU(2), 2/9 ≤ β_W < 1/3 | proved by forest gauge fixing, not yet refereed; 1/3 is the ceiling of single-link Dobrushin after any forest gauge fixing | `gauge-fixed/forest.md` |
| every SU(N), \|β\| < 1/29 (SU(3): β_W < 9/29) | proved by the gauge-fixed Hessian bound 29/2, not yet refereed | `gauge-fixed/curvature.md` |
| SU(3), β_W ≤ 3/4 | proved by single-link Dobrushin in the forest gauge with an analytic spectral-gap bound for the one-link measure and a corrector of degree two in the matrix entries, evaluated with the exact moments of degree six of the one-link measure at every staple sum (reduction to diagonal staple sums, character expansion of the Haar measure with rigorous truncation and rounding bounds, a cover of the diagonal staple sums by 2341 boxes, ball arithmetic), not yet refereed; the blended corrector of the earlier version stopped near β_W ≈ 0.715, the corrector of degree two reaches 1.0018 only at r = 1.54 in floating point, close to the Galerkin values at matrices of the allowed set rK put the limit of single-link Dobrushin at or below β_W ≈ 0.77, beyond which blocks of several links are needed | `gauge-fixed/su3.md` |
| SU(3), β_W ≤ 0.76 | proved by plaquette blocks in the forest gauge: α₃(0.76) ≤ 0.99950 from a second cover of the staple sums at r = 38/25 (1063 boxes, 2rΛ(r) ≤ 1.029967) and the Hessian and third cumulant of log Z on 784 boxes, ball arithmetic; not yet refereed; floating point puts the end of plaquette blocks near β_W ≈ 0.78 | `gauge-fixed/blocks/su3/` |
| every compact simple group, N\|β\| < 2κ_G/29 | proved by the gauge-fixed Hessian bound 29/2 and the Ricci curvature κ_G = h^∨/(2ℓ_ρ), with gap m_G(β) = arccosh(1 + (κ_G/(N\|β\|) − 29/2)/12); not yet refereed | `forest/`, Theorem on every compact simple group |
| the best gauge-fixed Hessian constant h_F | 9 + √17 ≤ h_F ≤ 14.420001 (Romanov's computer-assisted certificate, github.com/romanov360/yang-mills-4d-forest-certificate, checked with our own programs; 29/2 by hand), so \|β\| < 1/28.840002 (SU(2) β_W < 0.1387, SU(3) β_W < 0.3121); numerically 13.9740139454, attained by an abelian field; certificates on one plaquette and one vertex stop at 14.3839 | open |
| SU(2) from β_W = 1/3 to about 2.2 | no method; block (Dobrushin–Shlosman) criteria are the natural tool, and the paper's Theorem G gives one that is stable in every 't Hooft sector, but no box has been certified at these couplings | open |

## 2. The crossover to Bałaban's regime

Bałaban's renormalisation group controls the effective actions only at weak
coupling and in finite volume. Joining it to a strong-coupling statement
requires a criterion that one finite box can verify at a coupling where the
renormalisation group has already brought the effective action. Theorem G of
the paper is such a criterion at strong coupling; no criterion has been
verified in the window, and the hierarchical (Migdal–Kadanoff) models, where
confinement at every coupling is known, have no proved transfer to Z⁴.
Open.

## 3. The continuum limit with the Osterwalder–Schrader axioms

In two dimensions the continuum theory is known (Driver, Sengupta, Lévy). In
three and four dimensions no non-trivial limit with the axioms has been
constructed; lattice holonomies without smoothing converge to a trivial
(Haar white-noise) limit, so a smoothing that keeps reflection positivity
and a proof that the limit is not Gaussian are both needed. Open.

## 4. The gap in physical units

A positive gap at fixed spacing, as in problem 1, is a gap in lattice units.
The Clay problem asks for Δ_lat(a)/a ≥ m > 0 as a → 0 along the
renormalisation-group trajectory, together with the limit of problem 3
(Hypothesis R of `tex/main.tex`). The no-go theorem of `tex/main.tex` shows
that a gap bounded below in lattice units uniformly in the spacing forces
every renormalised two-point function to vanish, so the bound has to scale
with a. Open.

`forest/` (Section "The continuum step at one physical scale") reduces
problem 4 to one estimate: if the averaged-cube condition D(n, α) holds with
n a ≤ ℓ for a fixed physical length ℓ, uniformly along the continuum tuning,
then every continuum limit has a gap ≥ log(1/α)/(2ℓ). The lattice gap this
gives is of order a/ℓ in lattice units, so the no-go theorem does not apply.
The cube estimate (D) and the existence of the limit (problem 3) are open;
a linear-response count does not exclude (D) for ℓ a multiple of the
correlation length, but in gauges built on straight lines the gauge-fixed
links may fluctuate on all scales, in which case (D) has to be posed for
block variables.
The subsection "What the cube estimate requires" proves that D(n, α) holds
at strong coupling, that it forces the forest-gauge links of a cube with the
identity outside to stay at mean distance at least (1−α)D/3 from the identity
(at least 1−α for SU(2), by four colours and tetrahedral quaternions, and
(1−α)D/2 for every compact simple group with a representation metric, where
the cube keeps the fraction 1−α of the Haar disorder; hence it fails at every fixed n once β is large), and that bounds on the
conditional law of single lines, uniform in the exterior, cannot give it.
None of this decides (D) along the tuning.
The subsection "Block-constrained cubes" poses (D) for block variables, as
above: in the block-rooted axial gauge with full-block averages fixed, the
quadratic approximation has a gap 1/52 in every volume and satisfies the
constrained cube condition at every coupling once n ≥ 2·10¹⁵/α. One step of
composition is proved (Corollary "compose": cube coefficients at most the
constrained ones plus those of the averages; Theorem "effectiveupper": the
effective forms are equivalent to the coarse forms, Q_k <= S_k <= 256 Q_k,
uniformly in k (Whitney extension with nodes at the block centres); Proposition "effective": the
averages carry at least the coarse quadratic form at every scale; Theorem
"effectivelocal": in Feynman gauge S_k is a sum over the six planes of
convolution forms in the coarse curvature with kernels |K_P(y)| <= C e^{-η|y|_1},
C and η independent of k; Theorem "constrainedscales": the constrained cube
condition D_c(n,α) for the Gaussian with form S_k holds for n ≥ 10^827/α at
every scale k, and Corollary "compose" holds at every scale, constants
independent of k). Turning the cube condition for the averages into that of
the law one scale up, the cube condition for the coarsest averages, the large
fields and the crossover are open. For the group-valued links the straight path between
block roots gives a coercivity bound D_v² ≤ 224 Q for all fields (subsection
"Path blocks"), but that map does not compose across scales (thin paths),
and for every compact connected group no continuous covariant block map treats the transverse
positions of a block alike (Proposition "nomean"). A cut-off of Bałaban's
averaging operation is continuous, covariant and coercive for every field and
agrees with the block averages to first order (Theorem "cutoff"), but its
constants are those of one step. On torus-valued fields of small curvature at
every scale it composes, with nonincreasing energy (Theorem "torusscales"),
and for every compact group it composes on fields whose curvatures at the
successive scales have a small sum (Theorem "smallscales"). Under the SU(2)
Wilson measure the plaquettes farther than ρ from I form clusters whose
physical size tends to zero in probability, by reflection positivity and the
chessboard estimate, and the cut-off average differs from the path average
only next to them (Theorem "defects", Corollaries "defectclusters",
"activelinks"); such plaquettes are present along the one-loop tuning
(Proposition "manydefects"), and under a frustrated exterior a plaquette is bad
with probability tending to one. The composition is local and the defects
spread by at most three spacings of the current scale per step (Theorem
"defectscales"), so the cut-off average passes the clusters up to a coarse
spacing of order (β_W/log β_W)^{1/4} a (Corollary "defectscales"). In a box
of side l whose boundary links lie in a maximal torus with plaquette angles
at most ρ_s < π/3, the boundary field extends inside with every plaquette at
most Cρ_s, C < 5.46·10⁵ (Lemma "boxextension": facet-by-facet axial lift with
McShane extensions, then a radial pullback of Whitney forms), and conditionally
on such an exterior the fraction ε of bad unit cubes satisfies
−log ν(E_ε) ≥ β_W ρ² ε l⁴/16 − 6(l+1)⁴(log β_W + C_0 + 1/2 + C²β_W ρ_s²)
(Proposition "conddensity"; check gauge-fixed/checks/cube/box_extension.py).
Such exteriors are a null set. For arbitrary SU(2) boundary links with
plaquette angles at most ρ_s and l²ρ_s ≤ ε_0 = 1/(100 C_8²), C_8 = 85 743 550,
the boundary field extends inside with every plaquette at most Kρ_s,
K = 42 C_8 + 20 (Lemma "nonabelianfilling": coordinate-path gauges on the
facets glued through principal logarithms and McShane extensions, the radial
cubical extension of the logarithms, and transports bounded by subdividing
each plaquette; check gauge-fixed/checks/cube/nonabelian_filling.py), so the
conditional estimate holds with κ = Kρ_s on these open sets of exteriors.
At the level of Corollary "defectscales" this allows boxes up to a constant
times (β_W/log β_W)^{1/4} spacings; along the tuning their physical size
tends to zero, and for calm exteriors on cubes of fixed physical size the
conditional estimate is open (Remark "boxnonabelian").
On the fibres of the cut-off average the constrained gap survives the curvature
of the fibres: at any background, however far its links are from I, whose
plaquettes near the constrained blocks are at most κ ≤ 10⁻⁸ρ_G, every tangent
vector ξ of the fibre has |ξ|² ≤ 52(1 + 1.4·10⁶κ) Q_U(ξ), Q_U the covariant
curvature form, and the Hessian of the SU(2) Wilson action along the
constrained directions is at least (β/2)(1/52 − 2.7·10⁴κ); the second
fundamental form of the fibre costs a further O(κ) (Lemma "smallframes",
Theorem "smallfibre"; check gauge-fixed/checks/cube/small_field_fibre.py).
On a chart Φ(η) = U exp(η + h(η)) of the fibre over the convex set
K_s = {η tangent, max_e |η_e| ≤ s}, s ≤ 1/(6000 c_3), the Hessian of Q∘Φ is
at least (1/26 − 4·10⁵c_3(κ+s)), the Haar, coarea and area factors e^{−ψ}
have |ψ″| ≤ c_6, and for λ ≥ 60c_6 the constrained law restricted to the chart
satisfies a Brascamp–Lieb inequality with constant 60/λ (240/β for SU(2));
c_3, c_6 depend only on G, by compactness (Theorem "fibrechart"; check
gauge-fixed/checks/cube/fibre_chart.py). The decay of the constrained
influences (Helffer–Sjöstrand on the chart, which has a boundary, with the
non-local but exponentially decaying Hessian) and the mass off the chart
are open (Remark "fibrechartopen").
A block
map that is coercive for all nonabelian fields and composes across the
remaining scales is open.
For SU(2) the same subsection makes the size explicit: D(n, α) needs
n⁴ ≥ 2(1−α) β_W/(5 L(β_W)), L(b) = 1 + log 24 + (3/2) log 9b (Proposition
"cube size"), a bound met along the tuning; the n² lost in its proof is of the
right order (a stationary flux configuration). For every SU(N) the same bound
holds with L_N(b) = 1 + ((N²−1)/2) log(81π²b) + N(N−1) log(12/11)
(Proposition "cube size for SU(N)"). The continuum step also holds
for observables smoothed by the spatial gradient flow (Proposition "smoothed
observables"). Renormalised by one condition at a single physical
separation, the flowed fields have continuum limits along subsequences for
every tuning (Theorem "renormalised flowed observables"); a nontrivial limit
exists if and only if limsup C(σ)/C(0) > 0 for the flowed two-point function
at a fixed physical separation σ, and that inequality is not proved. In a
periodic spatial box the assumption C(0) > 0 is proved for product test
functions (strictly positive transfer matrix); the ratio is untouched. C(0)
is the reflected norm of the smeared packet, not a point coincidence (for the
free massless field the point ratio tends to zero). A lower bound on the ratio
gives tight spectral measures after renormalising at σ (Proposition "shifted
renormalisation"), and ⟨Ψ, (1−T)Ψ⟩ ≤ Ca for one slab would give the ratio. For
product test functions g ⊗ h the packet ratio lies between ratios of the
correlations D(r) of the spatially smeared slice field at two fixed physical
separations, so the open estimate is liminf D(r′)/D(r) > 0 for one pair
r < r′ with r ≤ 2t₀ (Proposition "slice ratio"; a bound at a later pair only
covers a packet shifted forward in time); it holds for the free field (limit
0.6337 at the parameters of Remark "free slice") and for two-dimensional
SU(N) lattice gauge theory on a circle of length ℓ with β = 2N/(g²a²), where
the energies −(ℓ/a²) log ĉ_λ(β) converge to g²ℓC₂(λ)/2 (Remark "two
dimensions"); in two dimensions there is no ultraviolet problem, and the
spatial flow alone does not give the inequality (a free field of mass M_k with
aM_k fixed has ratio → 0 under the same flow). For the free field,
conditioning an energy packet on the averages over blocks of side two leaves
at least the fraction 2^(−d)d^(−2) of its variance, at every lattice spacing
(Proposition "block averages of the free field"), so a route that replaces
the packet by conditional expectations on block averages has to measure its
errors in the reflected norm C(0), not in L². Every bound proved on such
ratios is an upper bound. Tree bounds on the truncated correlations of the
flowed fields (T1)–(T3) give uniformly bounded higher correlations, all
n-point Schwinger functions as reflection-positive distributions along a
subsequence, and a Gaussian limit when the physical volume of the tree weight
tends to zero (Theorem "higher correlations"). With amplitude scales r_k in
the tree bounds, the same hypotheses also give a time semigroup e^{−tH},
the ratio bound C(σ)/C(0) ≥ 1/2 (so they are at least as strong as the open
ratio), and, with (D), the mass gap for every polynomial in the smeared
fields. (T1)–(T3) are not proved along the continuum tuning. Their order-two
cases alone force the lattice gap m_k ≤ a_k log 2/σ and put a fixed fraction
γ of the tree weight at time distance at least 2t₀/a_k − 1, so they fail at
every fixed coupling in the strong coupling range (where C(σ)/C(0) ≤
e^{−m_G σ/a_k} → 0) and cannot come from any expansion whose convergence
gives a lattice gap bounded below (Proposition "tree bounds are critical").
For a packet smooth in time they give the ratio with an explicit rate,
C(σ)/C(0) ≥ e^{−2√B₁σ}, tight spectral measures and a limit packet of finite
energy (Proposition "energy control from the tree bounds"); this is still
conditional on (T1), (T2). The estimates uniform in the coupling (SU(N) force
moments, the Ward identity, the SU(2) plaquette variance, bare blocking) are
local and do not bear on the ratio, the tree bounds or (D). A three-state reflection-positive chain with a fixed gap
and a cutoff-independent two-point function has divergent four-point
functions, so the two-point function alone cannot replace them.

Monte Carlo (not used by any proof; `gauge-fixed/checks/cube/eps_mc.py`,
heat bath, four chains, 300 + 300 sweeps): the mean distance ε_V of the
forest-gauge links of a cube with the identity outside, which D(n, α) needs to
be at least 1−α for SU(2), grows with n at every coupling tried:

| β_W | n = 2 | n = 4 | n = 6 | n = 8 | n = 10 | n = 12 |
|---|---|---|---|---|---|---|
| 2.0 | 0.557 | 0.700 | 0.827 | 0.948 | 1.047 | 1.110 |
| 2.3 | 0.519 | 0.641 | 0.750 | 0.852 | 0.937 | 1.009 |
| 2.6 | 0.485 | 0.596 | 0.700 | 0.787 | 0.867 | 0.936 |
| 3.0 | 0.450 | 0.550 | 0.642 | 0.725 | 0.797 | 0.864 |

It passes 1 near n = 9 and n = 12 for β_W = 2.0, 2.3 and has not reached it at
n = 12 for β_W = 2.6, 3.0, far more slowly
than the correlation length grows, so this necessary condition excludes nothing
at a fixed physical scale. A coupled-chain estimate of the influence
coefficients themselves (`influence_mc.py`) is uninformative in this range:
the two chains decouple, and it returns only the trivial upper estimate.

Proposition "Gaussian benchmark" computes the cube coefficient of a Gaussian
field exactly: K_c = N_n at every n and every stiffness when the field is
massless, so D(n, α) fails for every α < 1, while a mass m gives D(n, α) at
a fixed physical scale. A proof of (D) along the tuning therefore has to
produce in the forest-gauge variables the loss that a mass supplies there,
from the Wilson action alone.

A further route that stops: joint-cube transport. Moving one boundary link j
along a geodesic changes the conditional cube measure μ_t ∝ e^{−H_t}, and
d/dt E_t f = −Cov_t(f, Ḣ_t) exactly, for all noncommuting links of the cube.
Solving the weighted Poisson equation L_t ψ = Ḣ_t − E_t Ḣ_t on G^V (L_t the
Langevin generator of μ_t) gives the transport field ∇ψ, and W1 is at most
∫ E_t Σ_e |∇_e ψ| dt. Two losses stop it. A Poincaré constant λ controls only
the ℓ² norm, E|∇ψ|² ≤ λ^{−1}‖Ḣ‖², so passing to the ℓ¹ distance d_V costs
|V|^{1/2}; avoiding that loss needs the decay of ∇_e ψ away from j, which is
the content of (D) itself. And λ cannot be bounded uniformly by curvature:
at the legal configuration U_j = −1 with every other link 1, the six
plaquettes through j have holonomy −1 and the weighted curvature in the
direction of j is 2 − 6β_W (unit S³), negative for β_W > 1/3 at every cube
size (check (j) of `cube_requirements.py`). The Bakry–Émery range β_W < 4/29
of the paper lies below that threshold. For the Gaussian benchmark the
Poincaré constant is of order β/n², while this route needs order β²n², so
even the quadratic approximation is out of its reach.

The weak-coupling approximation of the forest gauge itself is no better than
the massless benchmark. Linearise U = exp(iA) with A = 0 on F. The fixed links
form disjoint monotone lines (every site has one incoming and one outgoing
fixed link), and an active link joins two different lines. A gauge
transformation φ equal to 1 on one line ℓ and 0 elsewhere keeps A = 0 on F and
leaves the action invariant, so it translates the conditional law of the cube
by dφ restricted to V; for a Gaussian the cube coefficients then give
‖dφ|_V‖₁ ≤ Σ_j C_{V,j} |dφ(j)|. Summing over all lines, each link has
exactly two lines through its endpoints, so Σ_c K_c ≥ |V| = 48 N_n: some class
has K_c ≥ N_n and D(n, α) fails for every α < 1 at every n. The exact values
(`gauge-fixed/checks/cube/gaussian_forest.py`, sparse solve, independent of β)
are far larger and grow with n:

| n | max K_c/N_n | mean | min |
|---|---|---|---|
| 2 | 6.03 | 3.21 | 1.57 |
| 4 | 7.88 | 4.71 | 2.59 |
| 6 | 9.38 | 6.14 | 3.40 |
| 8 | 11.24 | 7.52 | 3.83 |

The lines are not what makes this work (Proposition "every forest gauge",
v4.13.0). For any 2Z⁴-invariant set F of links with no closed path whose trees
are infinite and cover the sites, the 2Z⁴-periodic flat fields
ψ(x,μ) = c_μ + h(x) − h(x+e_μ) form a space of dimension 19, F has at most 16
orbits, so some ψ ≠ 0 vanishes on F. It costs no action, so the mean of the
cube law given ψ outside is ψ on V, and Σ_c |ψ_c| K_c ≥ N_n Σ_c |ψ_c|: D(n, α)
fails for every α < 1 at every n and every stiffness, in every such gauge
(`gauge-fixed/checks/cube/gauss_every_gauge.py`, 24 checks; dimension 7 for
the forest field, 10 for the axial field). A mass, or the block averages held
fixed as in D_c, removes these fields; within the quadratic approximation the
loss has to come from the Wilson action beyond its quadratic part.

As in an axial gauge, a change of one boundary link is carried along the fixed
lines into the whole cube. Below the correlation length, where the theory is
weakly coupled, the cube coefficients of the forest gauge therefore follow this
growth, and (D) along the tuning requires the nonperturbative regime to reduce
them by a factor of order the cube size.
