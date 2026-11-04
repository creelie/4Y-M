# Boundary-preserving transport for constrained fibres

## A sufficient criterion

**Proved (conditional finite-dimensional criterion).** Let M be an oriented
smooth manifold with corners, equipped with a positive C¹ density m, and let
f:M -> R be C². On the region f^{-1}(J), where J is an open interval, suppose
there is a C¹ vector field V satisfying df(V)=1. Assume its flow exists for
the parameter displacements considered, preserves M and every hard label
restriction, and gives diffeomorphisms between the corresponding level sets.
Suppose the coarea disintegration m=ds nu_s is finite, positive, and varies
regularly enough to justify differentiating integrals below; for example this
holds on a compact flow box with C¹ transported density. Write

    S = div_m V,       Z(s) = nu_s(1),       mu_s = nu_s/Z(s).

For a C¹ observable F, with derivatives taken at fixed F,

    d/ds mu_s(F) = mu_s(VF) + Cov_{mu_s}(F,S),
    d/ds log Z(s) = mu_s(S).                                  (1)

Indeed, the flow transports nu_s to a measure whose density relative to
nu_{s+h} has logarithmic derivative -S: this follows from the Lie derivative
L_V m=S m and df(V)=1. Differentiating the integral transported in the
opposite direction gives d nu_s(F)/ds=nu_s(VF+FS); subtracting the
normalizer derivative gives (1). The formula also follows by writing the
flow box density as a(s,y) ds dy and noting partial_s log a=S.

If |V| <= v, osc_{f^{-1}(s)} S <= k, and every fibre has diameter <= D in
an ambient metric for which |VF| <= v for every smooth 1-Lipschitz F, then

    W_1(mu_s,mu_t) <= (v + D k/4)|s-t|.                       (2)

The covariance bound is |Cov(X,Y)| <= osc(X) osc(Y)/4. Integrate (1) and
take the Lipschitz dual supremum; standard smooth approximation is required
when applying this last step on a particular manifold. The estimate is
independent of the absolute size of Z, but the hypotheses on the normalized
score and on the flow must hold uniformly. If |S| <= K, then k <= 2K and
|log Z(t)-log Z(s)| <= K|t-s|. The latter requires an absolute score bound;
the oscillation bound alone does not control Z.

**Proved (boundary feasibility test).** At a point of a smooth active
boundary face h=0, any differentiable flow preserving that face has
dh(V)=0. At a corner it must satisfy this for every active face. Thus for a
vector-valued constraint f:M -> R^r, the restriction of df to the common
tangent space T of the active faces must be onto to lift every parameter
direction. If A=df|_T is onto, the minimum-norm pointwise lift of u is

    V_u = A* (A A*)^{-1} u,       |V_u| <= |u| / sigma_min(A).   (3)

This is elementary finite-dimensional linear algebra. Uniform positive
sigma_min gives a pointwise bound, not by itself a globally C¹, face-preserving
field: the active tangent spaces change at corners, and compatible lifts must
be constructed across strata. Nor does it bound div_m V or its oscillation.

## What the criterion rules out

**Proved (application to the scalar obstruction).** In the example
f(q,t)=q²+t on [0,1]² in `TRANSVERSE_FOLD.md`, the active faces at (0,0)
have common tangent space {0}, so no vector there can simultaneously preserve
both faces and satisfy df(V)=1. Even on the face t=0 with q>0, tangency
forces V_q=1/(2q), so any such lift blows up as q decreases to zero. This
matches the square-root response and excludes a uniform application of (2).
The original example is a warning about a proposed proof route, not a
counterexample for the SU(2) model.

## Gauge-fibre application remains open

**Sketch (testable route).** For each actual hard good/bad label and each
active stratum, form the derivative of the full coarse constraint restricted
to private-link directions tangent to *all* active label faces. Check its
smallest singular value uniformly in the exterior, including coincident
boundaries and corners. Even success leaves the construction of compatible
measurable lifts, bounds on their derivatives and Wilson-density scores,
and a comparison of observables across changing environments. In particular,
the known unconstrained private-link submersion does not imply (3) with a
uniform lower bound after restriction to active tangents.

**Heuristic (alternative).** Smooth or graded labels may avoid an exactly
preserved hard boundary, but derivatives of label weights then enter the
score. Their normalized expectations can still grow when a sector is rare;
the scalar obstruction shows that this trade must be estimated rather than
assumed harmless. Neither this criterion nor the alternative proves (I1),
uniform island control, D8, or continuum reconstruction.
