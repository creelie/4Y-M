# Gauge-fixed curvature in four dimensions: |beta| < 1/29 for every SU(N)

Status: proved by hand; not yet read by an outside referee. The local inequality is checked exactly
(sympy, every position of the frozen link), the global identity numerically on tori, the local
inequality adversarially on U(2), U(3), U(4).
Priority check (19 Oct 2026, arXiv and web search for gauge-fixed or axial-gauge Bakry-Emery bounds
for the Wilson action; Shen-Zhu-Zhu 2204.12737, Cao-Nissim-Sheffield 2505.16585 and 2509.04688,
Nissim 2510.22788, 2605.16162, Shen-Smith-Zhu): every curvature argument found uses the ungauged
measure, with the constant 8(d-1) or, in the paper in `../prd/`, the sharp 4d. No gauge-fixed bound found.

Normalisation of the paper in `../prd/`: weight exp(N beta sum_p ReTr Q_p), inner product <X,Y> = -Tr(XY),
Ric_SU(N) = N/2, Ric_SO(N) = (N-2)/4, conventional coupling beta_W = N^2 beta. Plaquette p at x in
the plane mu < nu has links e1 = (x,mu), e2 = (x+mu,nu), e3 = (x+nu,mu), e4 = (x,nu), transported
vectors Y_1..Y_4 (eq. transport), Hessian density f_p (eq. fp) and corner form
c_p = 2<Y1,Y4> - 2<Y1,Y2> + 2<Y2,W Y3 W*> - 2<Y3,Y4> (eq. corner), W = Q_p.

## 1. The line field

sigma(x) is read from x mod 2 by the table [0,0,2,3,3,1,2,1,1,2,1,3,3,2,0,0] (index
x0 + 2x1 + 4x2 + 8x3) and F = {(x, sigma(x))}. Lean `YangMills.StrongCoupling.GaugeFix`
(`YangMills/StrongCoupling/GaugeFix.lean`, `decide` only, axioms `propext` at most):

1. `one_fixed_per_plaquette`: every plaquette contains exactly one link of F;
2. `straight_lines`: sigma(x + e_sigma(x)) = sigma(x), so F is a disjoint union of straight lines,
   one through every site, and at every site the incoming and outgoing F-links point the same way.
   In a box with free boundary F is a finite union of segments, a forest;
3. `two_opposite`: a link e not in F lies in 6 plaquettes; in 2 of them the F-link is the link
   opposite (parallel) to e, in the other 4 it shares an endpoint with e.

## 2. A sum of squares at every vertex

**Lemma 1.** On Z^4, on a torus of side at least 3, or on a box with free boundary (links and
plaquettes of the box only), for every field Q and every X,

    sum_y [ |D_y X|^2 + sum_mu |X_{y,mu} + Q*_{y-mu,mu} X_{y-mu,mu} Q_{y-mu,mu}|^2 ] = 4|X|^2 + sum_p c_p.

Proof. Write o_mu, i_mu for the outgoing and the transported incoming vector in direction mu at y,
so D_y X = sum_mu (o_mu - i_mu). The squares contribute 2|o|^2 + 2|i|^2 at every end, 4|X|^2 in
total. The collinear cross term -2<o_mu, i_mu> of |D_y X|^2 cancels the +2<o_mu, i_mu> of
|o_mu + i_mu|^2. The perpendicular cross terms come only from |D_y X|^2 and are the corner terms
(Lemma divergence of the paper in `../prd/`). In a box a perpendicular pair of ends at y belongs to a plaquette of the
box, because a box is a product of intervals; so the bijection between perpendicular pairs and
corners holds there too. QED

This is Lemma divergence with the collinear term kappa written as a square, so nothing is lost by
Cauchy-Schwarz.

## 3. The local inequality with one frozen link

**Lemma 2.** Let k in {1,2,3,4} and let k-1, k+1 (mod 4) be the positions adjacent to k and k+2
the opposite one. For every N, every W in U(N) and all Y_1..Y_4 in u(N) with Y_k = 0,

    +-f_p + (9/10) c_p <= (21/10) (|Y_{k-1}|^2 + |Y_{k+1}|^2) + (5/4) |Y_{k+2}|^2 .

Since Y_k = 0, c_p consists of the two corners away from the frozen link.

Proof. As in Lemma local of the paper in `../prd/`: W -> -W changes the sign of f_p and fixes c_p; conjugation lets
W = diag(w_a); and the form F = RHS - f_p - (9/10)c_p splits into a block for each diagonal
position a, in y_i = -i (Y_i)_aa, and a block for each unordered pair a != b, in z_i = (Y_i)_ab
(using (Y_i)_ba = -conj z_i), with |Y_i|^2 contributing y_i^2 and 2|z_i|^2.

*Diagonal block* (w_a = e^{i phi}): sum alpha_i y_i^2 + cos(phi) (y1 + y2 - y3 - y4)^2
- (9/5)(y1 y4 - y1 y2 + y2 y3 - y3 y4) with y_k = 0, affine in cos(phi). At cos(phi) = +-1 the
3 x 3 matrices have leading principal minors

| frozen | cos phi = 1 | cos phi = -1 |
|---|---|---|
| e1, e4 | 31/10, 673/200, 8421/2000 | 11/10, 53/200, 21/2000 |
| e2 | 31/10, 861/100, 8421/2000 | 11/10, 21/100, 21/2000 |
| e3 | 9/4, 673/200, 8421/2000 | 1/4, 53/200, 21/2000 |

all positive, so the block is positive definite for every phi.

*Pair block* (w_a = e^{i(chi + theta)}, w_b = e^{i(chi - theta)}): the terms of f_p carry
w_a + conj(w_b) = 2 cos(chi) e^{i theta} or its conjugate, and c_p carries w_a conj(w_b) =
e^{2 i theta}, so the Hermitian matrix M(m, theta), m = cos chi, is affine in m with
M(-1, theta) = M(1, theta + pi). It suffices to take m = 1. For the frozen link e2, in the
variables (z1, z3, z4), with c = cos theta,

    M = [[ 21/5 + 2c,          -2 e^{i theta},       -9/5 - 2 e^{i theta} ],
         [ -2 e^{-i theta},     21/5 + 2c,            9/5 + 2 e^{-i theta} ],
         [ -9/5 - 2 e^{-i theta}, 9/5 + 2 e^{i theta}, 5/2 + 2c ]],

and with u = cos^2(theta/2) in [0,1] its leading principal minors are

    (20u + 11)/5,   (20u + 1)(20u + 21)/25,   (16000u^3 - 9200u^2 + 1600u + 21)/250 .

The determinant is 400u(40u^2 - 23u + 4) + 21 over 250, and 40u^2 - 23u + 4 has discriminant
529 - 640 < 0, so it is positive. The other positions give the same determinant and the leading
minors (20u + 11)/5 or (8u + 1)/2 and (800u^2 - 180u + 53)/50 (discriminant 32400 - 169600 < 0).
Every block is positive definite. QED

Checks: `checks/curvature/exact4.py` (all principal minors of all positions, both signs, m = +-1, exactly;
output `exact4.txt`), `checks/curvature/minors_table.py` (the table above), `checks/curvature/local_un.py`
(Nelder-Mead over W in U(2), SU(3), U(3), U(4): the least eigenvalue of the form is 0.014922 > 0,
the value of the exact blocks).

## 4. The gauge-fixed Hessian bound

**Theorem 1.** Let Lambda be Z^4 or a box of Z^4 with free boundary. For every field Q and every
square-summable X with X_e = 0 for e in F,

    +-Hess S_0(v,v) + (9/10) sum_y [ |D_y X|^2 + sum_mu |X_{y,mu} + Q* X_{y-mu,mu} Q|^2 ] <= (29/2) |X|^2 .

In particular |Hess S_0(v,v)| <= (29/2)|v|^2 on the gauge-fixed directions, against the sharp
constant 16 of Theorem A without gauge fixing.

Proof. Every plaquette has exactly one frozen link and X vanishes on it. Summing Lemma 2 over the
plaquettes, each active link receives 21/10 from at most four plaquettes and 5/4 from at most two
(exactly four and two on Z^4, fact 3), so sum_p (+-f_p + (9/10)c_p) <= (109/10)|X|^2. By Lemma 1,
(9/10) sum_p c_p = (9/10) sum_y [...] - (18/5)|X|^2. Adding, the constant is
109/10 + 36/10 = 29/2. QED

The identity behind it, sum_p S_p + (9/10) sum_y [...] = (29/2)|X|^2 -+ Hess S_0 with
S_p >= 0 the form of Lemma 2, is checked on the torus of side 4 for random SU(2), SU(3), U(3)
fields (`checks/curvature/verify_global.py`: differences below 3e-10 on totals near 10^5, least
eigenvalue of every S_p positive, link counts 4 and 2).

## 5. Consequences

**Corollary 2 (curvature regime).** Let G = SU(N) and |beta| < 1/29, or G = SO(N), N >= 3, and
|beta| < (N-2)/(58N). For every box Lambda with free boundary the gauge-fixed measure mu^F_Lambda
(Haar on the links of Lambda outside F, density proportional to exp(N beta sum_{p in Lambda}
ReTr Q_p) with Q_F = 1) satisfies the Bakry-Emery condition with K = N/2 - (29/2) N |beta|
(resp. (N-2)/4 - (29/2) N|beta|) > 0, uniformly in the box and the field; hence log-Sobolev with
constant 2/K and Poincare with constant 1/K.

Proof. {Q_F = 1} is a product factor, so its Ricci curvature is Ric_G on every link and the Hessian
of the restricted action is the restriction of Hess S_0 to X_F = 0. Theorem 1 and
[BGL14, Cor. 5.7.2]. QED

**Theorem 3 (uniqueness, decay, mass gap).** In the ranges of Corollary 2, d = 4:

1. for every gauge-invariant local O, the free-boundary Wilson states W_Lambda(O) converge as Lambda
   increases to Z^4 along any sequence of boxes, to a limit omega(O) that does not depend on the
   sequence;
2. omega is invariant under translations and the hypercubic group;
3. |Cov_omega(O, O')| <= C e^{-c dist(O,O')} with c > 0 depending only on N, K (C on the
   supports and norms of O, O');
4. if beta >= 0, omega is reflection positive and its transfer matrix has a positive spectral gap.

Proof. (a) Gauge fixing (as in `forest.md`): F n E(Lambda) is a forest, so W_Lambda(O) =
mu^F_Lambda(O^F), O^F = O with the F-links set to 1, a smooth local function of active links.
(b) Decay in boxes: mu^F_Lambda is the Gibbs measure of a nearest-neighbour interaction on a
product of copies of G; the commutator bound [SZZ23, Lemma 4.10] holds with the same constants
(fixing the F-links only removes terms) and the Poincare constant is 1/K, so the proof of
[SZZ23, Cor. 4.11] gives Cov(f,g) <= c_1 e^{-c d(f,g)}(|||f||| |||g||| + |f|_2 |g|_2) uniformly in
the box. (c) Limit: for boxes Lambda in Lambda' let S_s = S_Lambda + s(S_Lambda' - S_Lambda) on the
active links of Lambda', s in [0,1]. Hess S_s is a convex combination of the two Hessians, each
bounded by Theorem 1 (links outside Lambda carry no term of S_Lambda), so mu_s has the same K and
the same commutator constants. Since d/ds mu_s(f) = sum_{p in Lambda' \ Lambda} Cov_{mu_s}(f,
N beta ReTr Q^F_p) and mu_0(f) = mu^F_Lambda(f), we get |mu^F_Lambda'(f) - mu^F_Lambda(f)| <=
C_f R^3 e^{-cR}, R the distance from the support of f to the complement of Lambda. So the limit
exists along every increasing sequence, and two boxes are compared through their intersection,
which is a box. (d) Invariance: W_Lambda(tau_v O) = W_{Lambda - v}(O) and the limit does not depend
on the boxes (the pattern F is a device for computing W_Lambda, never shifted); the same for the
hypercubic group. (e) Decay passes to the limit through (a), (b). (f) For beta >= 0, W_Lambda is
reflection positive for boxes symmetric under the reflection (Osterwalder-Seiler), and so is the
limit along such boxes, which is omega by (c). Prop. gapdecay of the paper in `../prd/`, applied to smooth cylinder
F, gives a gap m >= c > 0. QED

The rate c is of order K over the commutator constants, as in [SZZ23, Rem. 4.12]; it is not
computed here.

What is new: for SU(N), N >= 3, the curvature range |beta| < 1/32 of Theorem B (1/48 in SZZ23) goes
to 1/29; for SU(3), beta_W < 9/29 = 0.310 against 9/32 = 0.281. For SO(N) the range goes from
(N-2)/(64N) to (N-2)/(58N). Theorem C says no Bakry-Emery argument for the ungauged measure passes
1/32; the gauge fixing changes the measure, not the metric, which is why it can. For SU(2) the
range is beta_W < 4/29 = 0.138, below the sharp quarter 2/9 and the forest Dobrushin bound 1/3, so
SU(2) gains nothing from this part.

## 6. How far this gauge can go

Let h_F be the supremum over fields of +-Hess S_0 on the gauge-fixed directions. Uniform
Bakry-Emery for mu^F reaches |beta| < 1/(2 h_F) and no further.

* Lower bound, exact: at the identity field -Hess S_0 restricted to X_F = 0 has the eigenvalue
  9 + sqrt(17) = 13.1231 on the torus of side 2 (factor lam^2 - 18 lam + 64 of the characteristic
  polynomial, `checks/curvature/ident_exact.py`), hence h_F >= 9 + sqrt(17) on Z^4 by periodic extension;
  for SU(2) the field Q = -1 on the links of the frustrated pattern gives the same with the sign +.
* Upper bound, proved: h_F <= 29/2 (Theorem 1).
* Numerical value: h_F = 13.97401 (gradient search over SU(2) fields of period 2, 4 and 6, both
  signs; best 13.9740139454 on the torus of side 4, 13.9740139 on side 6). The extremal field
  singles out one lattice direction: in planes containing it the plaquette half-traces are -0.996
  where the frozen link points along it and -0.439 otherwise, in the other planes -0.920; the top
  eigenvector puts 4 percent of its weight on that direction (`checks/curvature/analyse_best.py`,
  field `checks/curvature/extremal_field_L4.pt`).
* The maximum is attained by an abelian field. With every link in one U(1) subgroup of SU(2) and
  one angle for each residue class and direction (64 angles, extended periodically to the torus of
  side 4), a gradient search from 8 random starts reaches 13.9740139454 (`checks/curvature/abelian.py`,
  `checks/curvature/abelian.txt`, angles in `checks/curvature/abelian_extremal_L4.pt`), the same value as the
  search over all SU(2) fields. The plaquette angles of the maximiser are not rational multiples
  of pi in any visible way, so the exact value of h_F is not identified.
* Certificates of the type used here (any plaquette-local form plus any vertex form with
  perpendicular entries, coefficients free on the period-2 cell) cannot go below 14.3839
  (semidefinite program `checks/curvature/cert.py` with 72 angles: 14.38387; the three-parameter family
  `checks/curvature/sym.py` reaches the same at a = 2.0544, b = 1.2471, t = 0.9180). 29/2 is a rational
  point next to that optimum.

So with this gauge the method stops between |beta| = 1/29 and 1/(2(9 + sqrt 17)) = 1/26.25, and
numerically at 1/27.95. Reaching 1/27.95 needs certificates that couple neighbouring plaquettes.
There are exactly 8 period-2 choices of one frozen forward link at every site with one frozen link
in every plaquette (`checks/curvature/enum_sig.py`), all of the same kind as F. A pattern with one frozen
link at every site and one in every plaquette needs 2(d-1) = d(d-1)/2 (each frozen link lies in
2(d-1) plaquettes, each site carries d(d-1)/2 of them), so it exists only in d = 4; in d = 2, 3 a
pattern with one frozen link per plaquette has density d/4 per site and is not a line field.
