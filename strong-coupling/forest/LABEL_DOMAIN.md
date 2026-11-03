# The physical label domain is nonconvex even over the flat coarse field

**Proved (local analytic counterexample).** In the physical non-path
exponential coordinates of `CANONICAL_CHART.md`, the pullback of the event
that all physical links are in B(I,r) is nonconvex for every sufficiently
small r>0, already on the fibre of the identity coarse field. This rules
out direct use of the convex-product argument for that particular label
domain. It does not rule out a different estimate on nonconvex domains,
a different conditioning scheme, or a different coordinate system.

## 1. Reduce to one pair

Take c=I everywhere, freeze every other pair and every inside-block link
to I, and change just the eight outgoing private links of one pair k.
The direction can be chosen mu=0. The seven non-path outgoing links have
logarithms z_1,...,z_7. The path link has logarithm X, solved by

    H(X,z) = X+(1/8) sum_(j=1)^7 log(exp(-X) exp(z_j)) = 0.

Each outgoing link occurs in two paths; hence the factor 1/8. The other
two paths are the path link itself and give zero in the sum. Near zero
the cut-off equals one, H_X=I/8, and this equation is exactly the full
nonlinear coarse constraint. Private links enter no other average.

Choose HS-orthonormal A,B in su(2), with [A,[A,B]]=-2B and
[B,[B,A]]=-2A. For small positive a and real h put

    z_1=-(a/7)A+hB,      z_2=-(a/7)A-hB,
    z_3=...=z_7=-(a/7)A.

The two noncommuting terms have cancelling B and [A,B] components when
X is along A (this also follows directly from quaternion multiplication).
The implicit solution therefore has X=x(a,h)A, an analytic, even function
of h, with x(a,0)=a.

## 2. The boundary bends in the wrong direction

Let L(x,b,h) be the A-component of log(exp(-xA)exp(bA+hB)).
The BCH expansion through degree three gives

    partial_h^2 L(x,b,0) = x/3 + O((|x|+|b|)^2).

To check the sign and normalisation explicitly, the only axial h^2 term
at degree three is

    (1/12)[hB,[hB,-xA]] = (x h^2/6) A.

The degree-two commutator is transverse. Higher terms are analytic near
zero, so differentiating the Taylor remainder twice gives the stated
order. At h=0, the axial equation is exactly (x-a)/8=0. Differentiating
twice in h, with b=-a/7, therefore gives

    partial_h^2 x(a,0) = -2 partial_h^2 L(a,-a/7,0)
                       = -2a/3 + O(a^2) < 0

for all sufficiently small a>0. This strict negative sign is an analytic
asymptotic result; no floating-point calculation is needed for it.

Fix sufficiently small r>0. Let c_r=partial_h^2 x(r,0)<0. For small h set

    a(h)=r-(c_r/4)h^2 > r.

The two parameter points z(a(h),+h) and z(a(h),-h) have solved path radius

    x(a(h),+/-h) = r+(c_r/4)h^2+o(h^2) < r.

Their independent coordinate midpoint is z(a(h),0), whose solved path
radius is a(h)>r. All seven independent link radii are strictly below r
at the endpoints and at the midpoint, since they tend to r/7 as h tends
to zero; all other links are I. Thus the endpoints lie strictly inside
the physical all-good label event and their coordinate midpoint lies
outside it. This proves nonconvexity. Any independent-coordinate product
large enough to contain these three points inherits this obstruction.
The cut-off is inactive throughout for sufficiently small r,h.

## 3. Implication for the current route

**Proved (scope).** Adding coarse-curvature labels fixes the rough-background
obstruction but does not make this pullback convex: c=I in this example.
The strongly convex potential on the product graph chart of
`CANONICAL_CHART.md` is still useful, but restricting its support to the
physical label event introduces a nonconvex boundary. Product containment
and a near-isometry do not remove that boundary.

**Sketch (next analytic target).** One can seek boundary estimates for the
nonconvex pullback, or redesign labels in a fixed set of independent
coordinates and prove a local curvature witness for those labels. Both
require new work. Even a good-side estimate would leave the constrained
bad-side update and environment comparison to establish; D8 is unchanged.

**Computation (diagnostic only).** `d5/domain_probe.py` solves the quaternion
constraint numerically at r=0.1,0.03,0.01,0.003. In each case it finds two
good endpoints and a bad midpoint, with equation residual below 10^(-12).
These examples support the sign calculation but do not certify a numerical
value of the asymptotic validity radius. Output: `d5/domain_probe.out`.
