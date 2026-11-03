# Coarse-curvature labels and a sharper finite-intersection bound

All statements below labelled **proved** are mathematical deductions from
`smallframes`, `RELATIVE_CHART.md`, and the stated finite-intersection
chessboard lemma in `rare/RARE.md`. They are not formal verifications of
those inputs. These estimates are annealed on the SU(2) torus.

## 1. A local witness for rough coarse curvature (proved)

Let c=V^sharp(U), in tree gauge. If every fine plaquette needed by
`smallframes` for the four sides of a coarse plaquette P has chordal
defect at most kappa<=10^(-8), then

    delta(c_P,I) <= (4+32 pi) kappa < 105 kappa.

*Proof.* Write each coarse link c_k=d_k exp(M_k), with d_k the path link
of that pair. `smallframes` and `logbound` give ||M_k||<=8 pi kappa,
so delta(c_k,d_k)<=8 pi kappa. The four path links around P are the
holonomy of the fine square of side two: the other four links of its
straight sides are tree links and equal I. Its defect is at most 4 kappa
by inserting its four fine plaquettes with their parallel transports.
The bi-invariant chordal distance and the triangle inequality give the
bound on c_P. No commutativity is used.

Define a coarse block to be augmented-bad if either

- some active link based in B(y) is farther than r from R_e(c), or
- some coarse plaquette based at y has delta(c_P,I)>eta.

For 0<r<=50*10^(-8) and 0<eta<=105*10^(-8), put

    rho_bad = min(r/50, eta/105).

An augmented-bad block y forces a bad fine site at threshold rho_bad
in Q_y=2y+[-5,5]^4. The relative-link witness is `RELATIVE_CHART.md`.
For the coarse witness, the four pair bases are y,y+e_mu,y+e_nu. All
fine plaquette bases needed by their one-block collars lie in the stated
box. Taking the contrapositive of the preceding inequality proves the
claim. A bad fine site means one of the six fine plaquettes based there
has defect greater than rho_bad.

In particular, coarse-background calmness is now an explicit part of
the environment. The counterexample `ROUGH_BACKGROUND.md` is assigned
bad coarse labels rather than being silently included in the calm chart.

## 2. A coverage-counting lemma (proved)

Suppose fine bad sites satisfy P(T all bad)<=p^|T| for every finite T.
For a family of coarse labels indexed by y, assume label y being bad
forces a bad fine site in a specified set Q_y. Suppose

    |Q_y|<=M,
    each fine site belongs to at most H of the sets Q_y.

For nonempty S of n labels and k=ceil(n/H),

    P(S all bad) <= binomial(M n,k) p^k <= (e M H p)^k.

If A=e M H p<1, this gives P(S all bad)<=A^(n/H).

*Proof.* When all n labels are bad, the union of their witness sets
contains at least k distinct bad sites: any one site can witness at most
H labels. The union has at most Mn sites. Take a union bound over its
k-subsets and apply the assumed finite-intersection bound to each subset.
Finally binomial(Mn,k)<=(eMn/k)^k and n/k<=H. This proof requires neither
independence nor a simultaneous measurable choice of witnesses. The empty
set has probability one and is handled separately.

## 3. Explicit consequences (proved)

Let z(rho) be the **single-event** chessboard quantity in `rare/RARE.md`,
so p=6 z(rho), using its exact Haar-mass definition rather than an
asymptotic expansion. The bounds from section 2 give:

| labels | M | H | A |
|---|---:|---:|---|
| relative links only, radius r | 2401=7^4 | 256=4^4 | 6e*2401*256*z(r/50) |
| relative links plus coarse curvature | 14641=11^4 | 1296=6^4 | 6e*14641*1296*z(rho_bad) |

In either row, P(S all bad)<=A^ceil(|S|/H)<=A^(|S|/H) when A<1.
For the first row the witness box is 2y+[-3,3]^4, from
`TOUCHING_CHART.md`. Along any coordinate a fine site belongs to at most
four such boxes. For the second row there are at most six. Products give
the listed overlap counts, also on short periodic tori (wrapping can
reduce the number of distinct blocks, not increase it past these bounds).

The first row improves the exponent denominator 2401 from greedy
separation to 256, at the price of a larger factor inside A. The improved
bound need not be smaller near the boundary of its applicability; it is
stronger in the small-z regime. The second row pays for an additional
hypothesis that the original relative-link split failed to control.

**Sketch / open.** This supplies a rarity input for an augmented split,
not a mixture comparison theorem. The graph law in `CANONICAL_CHART.md`
is not automatically the law conditional on these labels. Hard-domain
geometry, constrained bad-side response, and the environment covariance
remain to be proved. No numerical threshold for (I1) is established, and
D8 still requires multiscale control.
