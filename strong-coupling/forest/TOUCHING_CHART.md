# Touching constraints: a local graph and an explicit witness loss

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

Continues `RELATIVE_CHART.md`. **Proved** means a deduction from the stated
paper lemmas, not a Lean verification. The relative decay theorem, (I1), D7
and (D) remain open.

## 1. A local graph retaining incoming averages (proved)

Work with SU(2), in the paper's fixed tree gauge. For a finite window W let
K consist of its outgoing pairs and the incoming pairs (y,mu) with y outside
W and y+e_mu in W. Put

    L = A_W union {all eight private links of each incoming pair},
    Dpath = {d_k : k in K},     Z = L minus Dpath.

Freeze all links off L and all averages c. By the private-link property in
`main.tex` before Lemma `privatesubmersion`, the only averages that involve
L are K. In particular enlarging L this way introduces no further pairs.
This describes an enlarged conditional window, rather than the conditional
law obtained by freezing every link off A_W.

Assume this fibre contains a background U0 for which the hypotheses of
Lemma `smallframes` hold for every pair in K, with curvature bound
kappa <= 10^(-8). This is an assumption, not a selection theorem from the
exterior. For 0 < epsilon <= 10^(-4), set

    D = product over e in Z of B(U0_e, epsilon),
    B_k = closed B(U0_dk, 44 epsilon).

There is a unique smooth graph Gamma over D with Gamma_e(z)=z_e for e in Z,
Gamma_dk(z) in B_k, all averages c_k fixed, and all other links frozen.
Uniqueness is only among configurations in these path-link balls.
Each solved link depends only on the non-path links entering its own pair.
Thus this graph has finite range in physical non-path coordinates.

*Proof.* Repeat the contraction in the proof of Proposition `calmfibre`,
now for every k in K. Set a=(7/8)Lambda, where
Lambda=alpha/sin(alpha), alpha=2 arcsin(1/(2 sqrt(2))). The path products
at U0 satisfy d(v_x,I)<=8 pi kappa. On D times B_k they satisfy

    d(v_x,I) <= 8 pi kappa + 46 epsilon < 0.01.

The cut-off is identically one there. With N_k the mean of the fourteen
nonzero logarithms divided by sixteen, solve

    ell = c_k exp(-N_k(ell,z)).

The map on the right has contraction constant a<0.9561, and moves the
centre by at most a(44+2)epsilon <44 epsilon. The closed balls are complete;
the contraction has a unique fixed point. No other path link enters this
pair, so these equations solve simultaneously without iteration between
pairs. The strict bounds give smoothness by the implicit function theorem
and preserve the finite support of each equation. This proves the claim.

For z,z' in D, let t=max_e d(z_e,z'_e). Each path uses at most two
non-path links. The logarithm bound used in the same contraction gives

    d(psi_k(z),psi_k(z')) <= a d(psi_k(z),psi_k(z')) + 2a t,
    d(psi_k(z),psi_k(z')) <= [2a/(1-a)] t < 44 t.

These constants do not depend on the number or shape of blocks. They do
not establish a Hessian estimate for the action in these coordinates.

## 2. The density in this graph (proved)

Let nu_c be the unnormalised Haar disintegration with respect to product
Haar measure of the averages. On this graph, relative to product Haar
measure dz on Z,

    d nu_c(Gamma(z)) = product over k in K of |det J_k(z)|^(-1) dz,

where J_k is the right-trivialised derivative of V^sharp_k with respect to
the path link d_k alone. The determinant is in dimension three. As in the
proof of `calmfibre`, J_k is an isometry times (1-T_k), with ||T_k||<=a.
Consequently

    (1-a)^3 <= |det J_k| <= (1+a)^3.

*Proof.* In coordinates (z,g), the derivative with respect to g is block
diagonal because path links are private. Change variables from g to the
averages, using Haar volume in both group variables. The graph area and
coarea factors cancel exactly as in the calculation following
`eq:fibrevolume`; only |det partial_g C| remains. The singular-value bounds
follow from ||T_k||<=a. Multiplying proves the formula.

**Sketch (remaining analytic work).** In exponential coordinates on D the
effective potential is

    V_eff = lambda Q(Gamma(z)) + sum_k log |det J_k(z)|
            - sum_e log j_e(z_e),

where j_e is the Haar Jacobian of the exponential coordinate. Its Hessian
needs bounds on second derivatives of Gamma and of log|det J_k|, including
terms involving first derivatives of Q. The determinant bounds alone give
no such derivative bounds. A calm background selected measurably from the
exterior, with quantitative control, is still missing.

The original physical good event restricts solved links too. Its pullback
is D intersect {d(psi_k(z),R_dk(c))<=r for the labelled good path links},
with the corresponding restrictions on the other links. This is not the
product D. Restricting the law to Gamma(D) changes the event; the graph
construction does not prove convexity or covariance decay for the actual
label-conditioned law. Enlarging the free set may also free links assigned
to the bad collar, whose label restrictions must be retained.

## 3. Input (b): explicit finite-intersection loss (proved)

Use the bad blocks b_y of `RELATIVE_CHART.md`, Corollary 1, at
0<r<=50*10^(-8)rho_G. On an even fine torus let zeta(r/50) denote the
single-event chessboard bound called z in `rare/RARE.md`, and put

    M = 7^4 = 2401,     Dsep = 7^4 = 2401,
    A = 6 M zeta(r/50).

If A<1, for every finite set S of coarse blocks,

    P(intersection over y in S of b_y) <= A^ceil(|S|/Dsep)
                                       <= [A^(1/Dsep)]^|S|.

*Proof.* A plaquette with a vertex in a block of {y}^+ has its fine base
in 2y+[-3,3]^4. Hence b_y implies a bad fine site in this box, at threshold
r/50. There are at most M possible bases. Select coarse blocks greedily,
deleting all blocks at periodic l-infinity distance <=3 from each selected
block. Each deletion removes at most 7^4 blocks. This selects at least
ceil(|S|/Dsep) blocks. Their fine-base boxes are disjoint: a common fine
base would imply coarse coordinate distance <=3 in every direction.
For each selected block take the union over its at most M witness bases.
Every choice consists of distinct fine sites. The finite-intersection
bound `rare/RARE.md`(i) is (6 zeta)^n for these sites, irrespective of the
chosen plaquette orientation. Union bounding gives A^n, and A<1 gives the
stated bound. Periodic boxes can wrap; their cardinalities only decrease,
and the periodic-distance argument still applies.

This is annealed torus rarity, including averaging over c. It gives no
uniform conditional probability after fixing c or an exterior. The factor
Dsep belongs in the exponent and cannot be folded into an additive entropy
constant. No new island threshold is asserted. D8 is unchanged.

**Computation:** `d5/touching_checks.py` (standard library only), with output
in `d5/touching_checks.out`, checks the contraction arithmetic, enumerates
the witness offsets, and checks periodic separation for cycle sides 1..64.
These checks support the bookkeeping; they are not a Hessian certificate.

