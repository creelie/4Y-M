# Connected normalisers inside the good chart

**Proved (deduction from `CANONICAL_CHART.md`).** Continuous boundary-source
changes inside a calm independent-coordinate product chart have an
exponentially small connected normaliser difference at large separation.
This is a step toward the missing environment comparison, not a bound for
changing a hard good/bad label.

## 1. A finite-dimensional identity

Split independent coordinates of a good graph chart into interior z and
exterior parameters eta,zeta. Let D be the fixed interior product of balls
and U(z;eta,zeta) the effective potential of the larger chart. Define

    Z(eta,zeta)=integral_D exp[-U(z;eta,zeta)] dz.

For scalar source directions u,v, differentiation under the finite smooth
integral gives

    partial_u partial_v log Z
      = -E[partial_u partial_v U]
        +Cov(partial_u U,partial_v U).

If U is a sum of local terms and the source supports are farther apart
than their interaction range, the first term is zero. This is an identity
for a fixed domain. It does not omit a boundary-transport term for moving
domains: moving domains are excluded by hypothesis.

## 2. Apply the explicit graph estimates

Use the larger canonical graph with its strengthened calmness condition
kappa_R+352r<=10^(-10), lambda>=2*10^11. All independent coordinates,
including the source parameters, lie in D_r. The coefficient estimates in
`CANONICAL_CHART.md` give

    Hess U >= lambda/104,
    ||Hess U|| <= 10^6 lambda,
    Hess U has coarse-coordinate range at most 3.

Conditioning on any set of exterior independent coordinates leaves a
product domain; the principal Hessian retains these bounds. The same
weighted covariance proof therefore applies, with alpha=10^(-14) and
amplitude 208/lambda. This also follows by repeating that proof on the
remaining coordinate labels; no probability comparison is being assumed.

Let the two exterior source supports have coarse distance d>6. Their
mixed direct Hessian vanishes. The interior gradients of their scores
partial_u U and partial_v U have support within distance three of their
respective source sets, and norms at most 10^6 lambda ||u|| and
10^6 lambda ||v||. Weighted covariance gives

    |partial_u partial_v log Z|
      <= 2.08*10^14 lambda ||u|| ||v|| exp[-10^(-14)(d-6)].

All constants are independent of the number of interior coordinates.
The source norm is the Euclidean norm of its independent-coordinate
variation. If there are no interior coordinates the covariance is zero
and the same assertion holds.

Integrate this estimate along two straight source segments remaining in
the same good product domain. Their connected finite difference obeys

    |log Z(eta1,zeta1)-log Z(eta1,zeta0)
       -log Z(eta0,zeta1)+log Z(eta0,zeta0)|
      <= 2.08*10^14 lambda ||eta1-eta0|| ||zeta1-zeta0||
         *exp[-10^(-14)(d-6)].

The proof is simply the preceding derivative identity, finite support of
the scores, the established weighted covariance estimate, and two
applications of the fundamental theorem of calculus.

## 3. Meaning for the constrained law

The coordinates belong to a specified larger calm graph, so all prescribed
averages c are fixed while the sources vary; solved path links are varied
with them. This avoids varying an arbitrary incompatible physical exterior.
The independent-link matching of `INDEPENDENT_LABELS.md` identifies the
appropriate enlarged good-window conditional laws. Exterior-only local
terms in the larger potential do not affect the connected difference
between separated supports. The statement is conditional on remaining
inside this chart and its fixed product support throughout the variation.

**Sketch / open:** a hard label change can change the integration domain,
enter a noncalm island, or add several path branches. None is covered by
this derivative identity. The large-field partition weights and their
connections to the good collar still need estimates. The hidden-phase
example in `ENVIRONMENT_OBSTRUCTION.md` shows why rarity alone does not
supply them. No theorem (I1), D7, D5 or D8 is asserted here.
