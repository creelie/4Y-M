# Finite-box interacting-measure insertions: an exact Wick prototype

**Proved (finite-volume coefficient formula).** On a simply connected finite
open SU(2) box in a full maximal-tree gauge, the order-g^2 connected magnetic
correlation can be computed by the Gaussian contractions below, including
action, observable and Haar corrections.

**Computation (exact rational on the smallest box; floating on larger boxes).**
The implementation evaluates this
formula on open boxes and obtains nonzero corrections to normalised
correlations. Independent CPU derivatives reproduce the B200 results, and
a change of maximal tree preserves the complete small-box result while
changing its sector decomposition. This does **not** close B2: the observable
is unflowed, the box is finite with open boundaries, and neither the
infinite-volume interacting vacuum nor cutoff-uniform bounds are established.

## 1. Exact finite-volume setup

A full maximal tree fixes all vertex gauge freedom except the root; its
links are I. Gauge-invariant integrals become product Haar integrals on
the remaining n edges, with no Faddeev--Popov factor. Put

    U_e=exp(g A_e),     beta_W=4/g^2,
    S=Q(U)/g^2.

The three colour components use an HS-orthonormal su(2) basis. On an open
contractible box, flatness and the tree condition force U=I, so the action
has a unique minimum. The reduced quadratic form H=D^T D is positive
definite. Finite-dimensional Laplace expansion around this minimum gives
asymptotic coefficients; the complement of a fixed coordinate neighbourhood
is exponentially suppressed as g tends to zero. This argument is at fixed
box size. It does not give bounds uniform in box size or a torus expansion
around an isolated minimum.

Expand Q= g^2 Q2+g^3 Q3+g^4 Q4+..., with

    Q2=sum_colour A_colour^T H A_colour,
    C=(2H)^(-1).

The Haar coordinate density is

    log j(g A_e)=-g^2 ||A_e||^2/6+O(g^4).

Thus V3=Q3 and V4=Q4+sum_e ||A_e||^2/6 in the effective density. The
observable is the sum of spatial plaquette Q at one time slice, divided
by g^2, with its own F2,F3,F4 retained.

## 2. Gaussian source formula

For two slice sources s,t let M_s,M_t be their quadratic plaquette forms,
K be the antisymmetric edge cubic tensor, and K_s,K_t its slice parts:

    Q3=6 sum_ijk K_ijk A_i^0 A_j^1 A_k^2,
    C(s,t)=[2H-2s M_s-2t M_t]^(-1),
    B(s,t)=-K+s K_s+t K_t.

For oriented plaquette edges let v_i=sign_i A_i/sqrt(2). Quaternion
multiplication gives the cubic vertex

    Q3_p=4 sum_(i<j<k) det(v_i,v_j,v_k).

This fixes the coefficient sqrt(2)/6 in the antisymmetric K. Independent
quaternion-polynomial multiplication checks this coefficient in the CPU
reference.

Write q4_p(C) for the Gaussian expectation of the quartic plaquette
coefficient. It is a quadratic polynomial in the four-edge covariance.
For d=3 colour components and B_ij=sign_i sign_j C_ij/2, its ingredients
are the elementary Wick identities

    E||v_i||^4=15 B_ii^2,
    E(||v_i||^2 ||v_j||^2)=9 B_ii B_jj+6 B_ij^2,
    E[(v_i dot v_j)||v_k||^2]=9 B_ij B_kk+6 B_ik B_jk,
    E[(v_i dot v_j)(v_k dot v_l)]
      =9 B_ij B_kl+3(B_ik B_jl+B_il B_jk).

They are applied to the quartic scalar part of the product of the four
quaternion exponentials. The cubic square has no same-vertex colour
tadpoles, and

    (1/2) E[Q3^2]=18 <K,C^tensor3 K>.

Consequently the order-g^2 term in the logarithm of the source partition
function is, up to source-independent constants,

    W(s,t)=sum_p[-1+s 1_(p in slice s)+t 1_(p in slice t)] q4_p(C(s,t))
           -(1/2)Tr C(s,t)+18 <B(s,t),C(s,t)^tensor3 B(s,t)>.

The connected coefficients are

    C0=6 Tr(M_s C M_t C),     C2=partial_s partial_t W(0,0).

This follows either by Gaussian generating functions or the cumulant
identity in `FRONTIER_LEMMAS.md`. All sectors, including mixed observable/
action cubic insertions and the Haar term, are present. Omitting them
would change the coefficient.

## 3. An independent derivative with bounded memory

**Proved algebra.** The explicit covariance jets are

    C_s=2 C M_s C,    C_t=2 C M_t C,
    C_st=4 C M_s C M_t C+4 C M_t C M_s C.

For F(B,C)=<B,C^tensor3 B>, with B=-K, B_s=K_s and B_t=K_t,

    F_st=2 <B_s,(C,C,C) B_t>
         +6 <B_s,(C_t,C,C) B>+6 <B_t,(C_s,C,C) B>
         +3 <B,(C_st,C,C) B>+6 <B,(C_s,C_t,C) B>.

The factors three and six follow by permuting tensor legs: both cubic
tensors are antisymmetric, so their permutation signs cancel. The quartic
expectation is quadratic in C and can be differentiated by its polarised
bilinear form. These formulas require only a small number of n-cubed
tensors, rather than the large collection of automatic-Hessian intermediates.

`jet.py` uses NumPy and these explicit derivatives. It reproduces all GPU
coefficients with worst scaled discrepancy below 1.4*10^(-13). The largest
case took about four seconds per target in this CPU implementation. For
these tested boxes, this improved CPU method is preferable to further GPU
runs; a new larger GPU proposal would need a fresh cost comparison.

## 4. Results and interpretation (computation)

For R_g(t)=C_g(t,t_ref)/C_g(t_ref,t_ref), write R=R0+g^2 R2+O(g^4).
The table shows the adjacent earlier slice; all separations are in lattice
units and all boxes have open boundaries.

| temporal x spatial sides | independent edges | R0 | complete R2 |
|---|---:|---:|---:|
| 3 x 2^3 | 29 | 0.03135055622 | 0.009603869231 |
| 5 x 3^3 | 244 | 0.02386137106 | 0.008764521430 |
| 7 x 3^3 | 352 | 0.02354246977 | 0.008599442566 |

The larger cases are deterministic floating Wick computations, not sampling
error bars or rigorous enclosures. The smallest case now has an exact
standard-library integer/Fraction computation, independently constructing
two maximal trees. Both give

    C0(1,1)=10877/2646,       C0(0,1)=341/2646,
    C2(1,1)=596780647/1960368480,
    C2(0,1)=2402567063/49009212000,
    R0(0,1)=341/10877,
    R2(0,1)=26843325164/2795053172625 > 0.

Thus the finite-volume reduction plus exact coefficient computation yields
this explicit nonzero order-g^2 term for the actual normalised finite-box
Wilson observable. The exact result is a trusted integer computation,
not a Lean verification of the Wick reduction. It establishes neither a
mass estimate nor a continuum extrapolation. The computational artifacts are maintained separately.
The irrational cubic factors cancel exactly: K=(sqrt(2)/6)T, T integer,
and 18 times their squared factor is one. Common covariance denominators
then let all cubic sums use integer arithmetic. The two-tree check took
less than one second of one CPU.

A second maximal tree on the smallest box preserves the complete R2 to
roundoff but changes the **measure-only** R2 from approximately 0.6716171
to 0.01378355. The complete R2 remains 0.009603869231. Thus the sector split
is gauge-coordinate-dependent. The computation supports retaining the
whole gauge-invariant correction rather than treating one sector as a
physical vacuum effect or assuming it is entirely a multiplicative Z-factor.
It does not contradict a specifically proved cancellation in the old
free-vacuum calculation; it shows why that cancellation cannot simply be
assumed for the full interacting correction.

## 5. Validation and next steps

The finite coefficients were checked through independent quaternion
polynomial multiplication, Gaussian moment diagnostics, independent
mixed-derivative formulas, and alternative maximal trees. The complete
smallest-box coefficient and its two-tree equality also passed exact
integer/Fraction arithmetic. Operational tooling, logs, source archives
and raw computational data are maintained outside this repository.

**Sketch / next:** reorganise the complete correction into gauge-invariant
connected contributions that can be bounded before taking large volume;
then incorporate the spatial-flow observable and address infrared and
renormalisation limits. Strong sector cancellations make a bound on each
large gauge-dependent term separately unattractive. This is a prototype
for that analysis, not a replacement for the missing B2 theorem or D8.

## 6. Nonlinear-flow follow-up (2026-11-03)

**Proved finite reductions / computation:** `NONLINEAR_FLOW_CHECK.md`
gives full-link spatial-flow Taylor/Duhamel formulas through the required
order, including the exponential-coordinate metric correction. A solvable
one-variable SU(2) model shows that omitting nonlinear flow changes its
normalised order-g^2 coefficient by a factor of three. Independent Haar
quadrature checks that calibration. A separate finite spatial-action
gradient check on a 3x2x2x2 box validates the local vertices and detects
the error from freezing tree links during flow. The actual flowed Wick
covariance and the uniform B2 estimate remain open.

**Proved conditional B/E interface:** `POSITIVE_TIME_TIGHTNESS.md` turns a
lower bound on C_a(s)/C_a(s/2) at fixed positive physical s into a uniform
energy bound for the time-smoothed normalised state. The full interacting
correlation ratio, not merely a finite perturbative coefficient, must
supply that lower bound.
