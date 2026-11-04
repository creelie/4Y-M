# A uniform radial hazard bound on SU(2)

**Proved (a weaker bound than the conjectured Gaussian-scale estimate).**
Identify SU(2) with the unit three-sphere with its standard angular
distance. Let U have density proportional to exp(kappa a·u) relative to
surface measure, with kappa>=0 and a any unit vector. For any centre g,
let R=d(U,g), let f_g be its density, and let 0<r<=pi/2. Then

    f_g(r) <= e sec²(1) (kappa+1) P(R>=r).                       (1)

This bound is uniform in both a and g and in the normalizing constant. The
distance convention in `CALM_LINK.md` must be converted if its group metric
uses a different fixed normalization; this changes only fixed factors.

*Proof.* Parametrize u=cos(s)g+sin(s)omega, where omega lies on the unit
two-sphere tangent to g. Up to a normalization independent of s,

    f_g(s)=sin²(s) integral exp[kappa a·(cos(s)g+sin(s)omega)] d omega.

Put h=1/(kappa+1), so 0<h<=1. For 0<=t<=h, the unit-speed great-circle
parametrization gives |a·u(r+t,omega)-a·u(r,omega)|<=t. Also, since
r<=pi/2, sin(r+t)>=sin(r) cos(t)>=sin(r) cos(1), and r+h<pi.
Therefore

    f_g(r+t) >= cos²(1) exp(-kappa h) f_g(r)
              >= e^(-1) cos²(1) f_g(r).

Integrate over 0<=t<=h and use P(R>=r)>=integral_r^{r+h} f_g(s) ds.
This proves (1), including kappa=0. No Gaussian approximation, small
radius, or assumption on the relative positions of a and g is used.

## Relation to the shell problem

**Proved conditional consequence.** If, after freezing all other variables,
a shell radius is the angular distance from a fixed centre of one SU(2)
link whose conditional density has the above form with kappa<=c lambda,
then its radial density at any fixed r<=pi/2 is bounded by
C(1+c lambda) times its conditional tail, with C=e sec²(1). Integrating
over the frozen variables preserves that inequality. For lambda>=1 this
has the O(lambda) order used in `CALM_LINK.md`'s fixed-radius flux ledger.

**Open.** The sharper proposed O(sqrt(kappa)+kappa r) bound is not proved
by (1). More importantly, applying the estimate to the actual maximum of
path distances requires checking which link is varied, that the reference
centre is fixed under that conditioning, and how ties, cut-locus points,
constraint disintegration and moving kink surfaces enter the flux. The
subsequent bad-event displacement H and environment control are still
open. Thus this calculation does not establish Conjecture B, (I1) or D.

**Computation.** `python3 d5/su2_hazard_check.py` integrates the explicit
radial density by composite Simpson quadrature for four concentrations,
five relative centres and three radii. The largest observed ratio
f(r)/[(kappa+1)P(R>=r)] is approximately 1.163, versus the proved bound
9.312. This floating-point grid is a diagnostic, not a proof or a uniform
optimization.
