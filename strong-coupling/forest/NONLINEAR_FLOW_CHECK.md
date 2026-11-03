# Nonlinear flow contributes at the same order as the interacting correction

## 1. An exactly solvable compact-group calibration

**Proved (single-variable SU(2) model).** This is a flow of one group
variable, not the four-dimensional lattice spatial Wilson flow. It tests
which perturbative terms a valid extension of `B_WICK.md` must retain.

Write U=cos(theta)I+i sin(theta)n.sigma, with Hilbert--Schmidt metric
2 dtheta^2 along a radial line, and set Q(U)=4(1-cos(theta)). Use Haar
probability weighted by exp[-Q/g^2]. The negative Riemannian gradient
flow of Q satisfies

    theta'=-2 sin(theta),
    tan(theta_tau/2)=c_tau tan(theta/2),    c_tau=exp(-2tau).

For theta=g|A|/sqrt(2), expanding the exact flow gives

    theta_tau=c_tau theta+c_tau(1-c_tau^2)theta^3/12+O(theta^5),
    F_tau := Q(U_tau)/g^2
      = c_tau^2 X + g^2 d_tau X^2 + O(g^4),
    X=|A|^2,    d_tau=c_tau^2/12-c_tau^4/8.

The small-field Gaussian density is exp(-|A|^2) dA in three dimensions.
The first relative density correction, including Haar, is

    w_2=X^2/24-X/6.

All coefficient statements here are finite-dimensional asymptotic
expansions as g tends to zero. The compact integral has a unique action
minimum at I; its complement outside a fixed neighbourhood contributes
exponentially small terms. Standard local Taylor expansion and Gaussian
domination justify extracting these coefficients at fixed flow times.

## 2. Complete connected coefficient and ratio

**Proved (Gaussian moment calculation).** X has Gamma shape 3/2, scale 1.
Its moments give

    Var(X)=3/2,    Cov(X,X^2)=15/2,
    kappa(X,X,w_2)=7/16.

Consequently Cov_g(F_tau,F_upsilon)=C0+g^2 C2+O(g^4), where

    C0=(3/2)c_tau^2 c_upsilon^2,
    C2=c_tau^2 c_upsilon^2
         [27/16-(15/16)(c_tau^2+c_upsilon^2)].

For the ratio to the equal-reference covariance,

    R=Cov_g(F_tau,F_r)/Cov_g(F_r,F_r)=R0+g^2 R2+O(g^4),
    R0=c_tau^2/c_r^2,
    R2=(5/8)R0(c_r^2-c_tau^2).

If instead one replaces the exact flow by A_tau=c_tau A before
evaluating Q, its quartic coefficient is -c_tau^4/24. Repeating the
same complete measure calculation gives only

    R2_linear_flow=(5/24)R0(c_r^2-c_tau^2).

Thus omitting the nonlinear flow changes this coefficient by a factor
of three when it is nonzero. Normalisation does not remove the omitted
term. At equal flow times both ratio corrections vanish, as they must.

## 3. Scope and next B calculation

**Sketch / next.** A lattice spatial-flow expansion must propagate the
quadratic and cubic field corrections to the required order and include
their contributions to the flowed observable along with the Wilson,
Haar and connected-normalisation terms. Flowing tree-gauge variables
while silently holding tree links fixed changes the flow metric; one
must use the actual full-link flow or derive the compensating gauge
description. This one-variable calibration does not supply those lattice
vertices, an infinite-volume estimate or B2.

Flow time tau here is distinct from the physical Euclidean separation
in `POSITIVE_TIME_TIGHTNESS.md`. The above ratio is a perturbative flow
diagnostic, not a transfer-matrix mass measurement.

**Computation / validation.** Exact rational Gaussian moments give the
displayed 7/16 and 15/2 factors. Independent direct Haar quadrature at
five small couplings and three fixed flow times extrapolates to R2 within
2.1e-10; doubling quadrature from 256 to 512 nodes changes the ratios by
less than 5e-16. At reference flow time zero, R2 at tau=0.1,0.25,0.5 is
approximately 0.1381194262, 0.1453400987 and 0.0731372777. These floating
comparisons are not interval certificates or remainder estimates. The
analytic expansion establishes the coefficient; the quadrature checks
its implementation. Standalone mathematical source is in
`wick/flow_calibration.py`; raw outputs remain outside Git.

## 4. Full-link finite-lattice flow jets

**Proved (local finite-dimensional Taylor reduction).** Use orthonormal
Lie-algebra coordinates U_e=exp(g A_e) on every link, the product
bi-invariant metric, and negative gradient flow of a spatial action Q_s
with expansion Q_s(gA)=g^2 Q_2(A)+g^3 Q_3(A)+g^4 Q_4(A)+O(g^5).
This specifies the flow normalisation; another convention rescales tau.
Write ad_A linkwise. The pulled-back metric divided by g^2 is

    G(A)=I+(g^2/12) ad_A^2+O(g^4).

To see the coefficient, the left exponential Jacobian is
I-(g/2)ad_A+(g^2/6)ad_A^2+O(g^3); multiply it by its transpose and use
ad_A*=-ad_A. Its inverse metric is I-(g^2/12)ad_A^2+O(g^4).
The A-coordinate flow is therefore

    A' = L A + g B(A,A) + g^2 T(A,A,A) + O(g^3),
    L A=-grad Q_2(A),
    B(A,A)=-grad Q_3(A),
    T(A,A,A)=-grad Q_4(A)+(1/12)ad_A^2 grad Q_2(A).

Here B and T are the symmetric polarisations of their homogeneous
polynomials. The Q_k are spatial-action coefficients, not automatically
the full four-dimensional measure-action coefficients of `B_WICK.md`.
For A(0)=A independent of g, put A_tau=A_0+g A_1+g^2 A_2+O(g^3).
Variation of constants gives

    A_0(tau)=exp(tau L) A,
    A_1(tau)=integral_0^tau exp((tau-s)L) B(A_0(s),A_0(s)) ds,
    A_2(tau)=integral_0^tau exp((tau-s)L)
       [2 B(A_0(s),A_1(s))+T(A_0(s),A_0(s),A_0(s))] ds.

For an observable F_2+g F_3+g^2 F_4+O(g^3), its flowed coefficients are

    Phi_2=F_2(A_0),
    Phi_3=DF_2(A_0) A_1+F_3(A_0),
    Phi_4=DF_2(A_0) A_2+(1/2)D^2F_2(A_0)[A_1,A_1]
          +DF_3(A_0) A_1+F_4(A_0).

These follow by equating powers of g in the finite smooth ODE and Taylor
expanding the observable. Phi_2, Phi_3, Phi_4 have degrees two, three and
four in the initial field, so the existing connected Wick/cumulant
reduction applies after substituting their actual coefficients.

**Sketch / implementation target.** Embed the initial tree-gauge Gaussian
field into all links, initially zero on tree edges, and evolve all links
in these equations. Do not freeze tree edges during the flow. Construct
the actual spatial Q_k and observable tensors, check the infinitesimal
flow against direct group evolution, then apply the complete interacting
measure correction. This reduction does not perform those contractions
or bound the large-volume, long-flow-time or continuum remainders.

**Computation (actual finite spatial action).** On a 3x2x2x2 open box
with all 52 links and 18 spatial plaquettes, the vector-field reduction
was checked at a deterministic nonabelian field initially zero on its
23 tree edges. Differentiating a separate quaternion action evaluation
(trigonometric series through degree 16, well within floating precision
at these fields) and using the exact radial/transverse exponential metric
gave full truncation errors approximately

    g=0.2: 3.59e-4,  g=0.1: 4.48e-5,
    g=0.05: 5.59e-6, g=0.025: 6.99e-7.

The factor of about eight per halving agrees with an O(g^3) remainder.
Omitting the metric term instead leaves errors 4.00e-3, 9.98e-4,
2.50e-4 and 6.24e-5, consistent with O(g^2). The initially zero tree
edges have nonzero velocity norm about 2.67. This tests the local
full-link flow vertices; the ODE jets and flowed Wick contractions have
not yet been computed. Runtime was about five seconds after imports on
CPUs, so this check warrants no GPU allocation. These are floating
diagnostics, not rigorous enclosures.
