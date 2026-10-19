# Gauge-fixed Dobrushin bound for four-dimensional SU(2): beta_W < 1/3

Status: proved by hand; not yet read by an outside referee. Priority check: no earlier statement found in the
sources searched (arXiv and web search for Dobrushin + axial/maximal-tree gauge + Wilson theory;
Shen-Zhu-Zhu 2204.12737; the Lean project mrdouglasny/LGT, which has beta < 1/(32 n (d-1));
Dimock 1410.2798 on Balaban's axial gauge). Seiler, LNP 159 (1982), and Osterwalder-Seiler (1978)
were not available to the authors; their strong-coupling results are cluster expansions without
an explicit Dobrushin range.

Normalisation as in the paper in `../prd/` (Shen-Zhu-Zhu): weight exp(N beta sum_p ReTr Q_p), beta_W = 4 beta
for SU(2). Ingredients from that paper: the sharp quarter W_1(nu_b, nu_b') <= |b - b'|/4 (chordal
metric on S^3), the comparison lemma (Dobrushin with W_1, rows and columns), Prop. gapdecay
(transfer-matrix gap from time decay in a reflection-positive, translation and hypercubic invariant
state).

## The link set F

Read each site x of Z^4 modulo 2 and set sigma(x) by the table (index x0 + 2x1 + 4x2 + 8x3)
`[0, 0, 2, 3, 3, 1, 2, 1, 1, 2, 1, 3, 3, 2, 0, 0]`. Let F be the set of forward links
(x, sigma(x)). Then

1. every plaquette of Z^4 contains exactly one link of F (96 residue cases; Lean
   `one_fixed_per_plaquette`, by `decide`; independent count on a 6^4 box, `checks/verify_gauge.py`);
2. F is a forest on Z^4: every site has exactly one outgoing F-link and it raises the height
   h(x) = x0 + x1 + x2 + x3 by one, so F has no cycle (a cycle in a graph with out-degree <= 1 is a
   directed cycle, impossible as h increases). For every finite box Lambda (free boundary, not a
   torus), F n E(Lambda) is a finite forest: the outgoing links of the top sites leave the box.
   On a torus F does contain cycles; the argument never uses a torus.
3. every link not in F has exactly 6 x 2 = 12 partners not in F in its six plaquettes (Lean
   `twelve_neighbours`), against 18 without gauge fixing and 14 in temporal gauge.

## The argument

Specification. gamma^F is the local specification on the links of Z^4 \ F: a priori Haar, weight
exp(N beta sum_p ReTr Q_p^F), where Q_p^F is the plaquette holonomy with the links of F replaced
by 1. No infinite-volume gauge transformation is performed: gamma^F is a specification in its own
right.

Dobrushin. The conditional law of a link e of gamma^F is von Mises-Fisher nu_b with b = 4 beta x
(sum of its six staples), each staple a unit quaternion (product of its two partners, F-links being
1). Two distinct links share at most one plaquette, and moving a partner e' moves one staple by an
isometric image of the displacement. So C_{ee'} <= |beta| (sharp quarter), C_{ee'} = 0 unless e, e'
share a plaquette, and every row and column has 12 nonzero entries: alpha <= 12 |beta| = 3 beta_W.
For beta_W < 1/3, gamma^F has exactly one Gibbs measure mu_F, and Lipschitz local functions have
truncated correlations decaying at every rate lambda < -log(3 beta_W) in the distance between link
midpoints.

Finite boxes. For a box Lambda with free boundary conditions, let W_Lambda be the Wilson measure
on the links of Lambda. Gauge-transform along the finite forest F n E(Lambda) (root each tree, as in
Creutz's maximal-tree gauge): the links of F become 1, the other links become g_x U_e g_y^{-1}, which
for fixed values on F is a Haar-preserving bijection. Since the action is gauge invariant, the
image of W_Lambda on the active links is mu^F_Lambda (Haar times exp of the plaquettes inside
Lambda with F = 1), independent of the F-links. Hence for a gauge-invariant O supported in Lambda,
W_Lambda(O) = mu^F_Lambda(O^F), O^F = O with F-links set to 1, which is local and Lipschitz in each
active link (Wilson loops, plaquettes). No leftover holonomy appears: F n E(Lambda) is a forest.

Infinite volume. mu^F_Lambda and mu_F differ in their single-link conditional laws only at links
within distance 1 of the boundary or outside Lambda; Dobrushin's comparison theorem (Georgii, Thm
8.20) bounds |mu^F_Lambda(f) - mu_F(f)| by sum_{e,e'} delta_e(f) D_{ee'} b_{e'}, D = sum C^n, and for
f supported at distance R from the boundary this is O(R^3 e^{-lambda R}). So
omega(O) := lim W_Lambda(O) = mu_F(O^F) exists for every gauge-invariant local O and every sequence
of boxes exhausting Z^4, centred anywhere.

Properties of omega. Translation invariance: omega(tau_v O) = lim W_{Lambda - v}(O), and the boxes
Lambda - v also exhaust Z^4, so the limit is omega(O) (the pattern F is not shifted; only the box
is). Hypercubic invariance: the same with R(Lambda). Reflection positivity: W_Lambda is reflection
positive for boxes symmetric in the plane (Osterwalder-Seiler; free boundary), and weak limits keep
it. Decay: Cov_omega(O, O') = Cov_{mu_F}(O^F, O'^F) <= C e^{-lambda dist}.

Gap. Prop. gapdecay of the paper in `../prd/` then gives a transfer-matrix gap m >= -log(3 beta_W) for beta_W < 1/3.
(The paper's Theorem F used torus limits and uniqueness of the Wilson Gibbs measure; here the
state is the free-boundary limit, which has the properties gapdecay needs.)

## Ceiling

For every forest F' of links (any gauge fixing by a forest), sup_e (row sum) >= 12: in a box of side
n, |F'| <= #sites - 1, so at least 3 n^4 - O(n^3) links are active; the active-link counts a_p of the
6 n^4 - O(n^3) plaquettes have mean 3 - o(1), and sum_e row(e) = sum_p a_p (a_p - 1) >= (number of
active links) x 6 (abar - 1), so the mean row sum is >= 12 - o(1); rows are integers. Staples of a
link involve disjoint link sets, so they can cancel (b = 0) and C_{ee'} = |beta| is attained (Prop.
ceiling of the paper in `../prd/`). Spectral radius >= mean degree, so weighted Dobrushin does not help either.
beta_W = 1/3 is the ceiling of single-link Dobrushin after any forest gauge fixing.

## What it is not

A strong-coupling result at fixed spacing. The window between beta_W = 1/3 and Balaban's
weak-coupling control remains open, and Hypothesis R (the Clay problem on the lattice) is untouched.
