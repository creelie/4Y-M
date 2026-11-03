# A locality norm for one constrained Gaussian elimination

**Proved (finite matrix lemma).** Let a positive definite quadratic form
be expressed in retained and eliminated coordinates as

    H = [ A  B ; B*  D ],    S=A-B D^(-1) B*.

S is its constrained minimum form on the retained variables. Give both
index sets a common location metric. For each 1-Lipschitz real function
phi and fixed alpha>=0, conjugate blocks by W=diag(exp(alpha phi)) on
their respective input/output sets. Write ||M||_alpha for the supremum
of these operator norms. For rectangular blocks use the same phi at
input and output. These norms are submultiplicative on compatible
blocks. Because the family includes -phi, ||B*||_alpha=||B||_alpha.

Suppose

    ||A||_alpha<=a,  ||B||_alpha<=b,  ||D^(-1)||_alpha<=g.

Then ||S||_alpha<=a+b^2 g. In particular its individual entries obey

    |S_ij| <= (a+b^2 g) exp[-alpha dist(i,j)].

*Proof.* Conjugation preserves block products, so the norm bound follows
from the Schur formula and submultiplicativity. Choose phi(x)=dist(x,j)
for the entry estimate; then phi(i)-phi(j)=dist(i,j).

## A quantitative perturbation ledger

**Proved (finite matrix deduction).** Let H_0 have blocks A_0,B_0,D_0 and
Schur complement S_0. Assume both B norms are at most b, both inverse-D
norms at most g, and each block difference A-A_0, B-B_0, D-D_0 has
weighted norm at most epsilon. The transpose differences have the same
norm. The resolvent identity gives

    ||S-S_0||_alpha <= epsilon (1+b g)^2.

Indeed split B D^(-1) B* - B_0 D_0^(-1) B_0* into

    (B-B_0) D^(-1) B*
      + B_0 D^(-1) (B*-B_0*)
      + B_0 [D^(-1)-D_0^(-1)] B_0*,

and use D^(-1)-D_0^(-1)=-D^(-1)(D-D_0)D_0^(-1).
The four resulting bounds are epsilon, epsilon b g, epsilon b g and
epsilon b^2 g^2. These hypotheses control locality as well as size;
ordinary unweighted form comparability cannot replace them.

## Meaning and remaining RG inputs

**Sketch / open.** For a general surjective coarse map Q, a right inverse
E and a kernel basis F give x=Ev+Fz and blocks of [E F]* H [E F]. Apply
the lemma only after verifying locality norms of these coordinate maps
and the eliminated inverse. Arbitrary choices of E,F can be nonlocal.

At repeated scales this formula identifies the one-step amplification
factor (1+b_k g_k)^2. Bounds must be expressed in consistently rescaled
coordinates and metrics, with accumulated perturbations controlled;
summability of bare input errors alone is not a bound on their amplified
sum. No Yang--Mills effective action, all-scale bounds or non-Gaussian
large-field estimates are constructed here. This refines the finite D8
interface of `FRONTIER_LEMMAS.md`; D8 remains open.
