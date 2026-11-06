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

**Proved (maximum of dependent radii).** Let R_1,...,R_N be random
angular distances, and suppose each R_i, after conditioning on an
appropriate sigma-field, obeys the hypotheses of (1), with concentration
at most K and a centre measurable with respect to that sigma-field. No
independence between the R_i is assumed. For M=max_i R_i, its law on
(0,pi/2) is absolutely continuous and its density satisfies, for almost
every r there,

    f_M(r) <= N e sec²(1) (K+1) P(M>=r).                         (2)

Indeed, (1) survives conditional averaging and bounds each marginal density
f_i(r) by e sec²(1)(K+1)P(R_i>=r). For every interval I,
P(M in I)<=sum_i P(R_i in I), since M is attained by one of the R_i.
Consequently f_M<=sum_i f_i almost everywhere, and each tail on the right
is at most P(M>=r). The same domination on null sets proves absolute
continuity. For chordal distance delta=2 sqrt(2) sin(R/2), apply the change
of variables; on R<=pi/2 its derivative is bounded away from zero, so only
a fixed factor changes.

**Proved (fixed-radius one-sided shell bound).** Under the same hypotheses,
for every r in (0,pi/2),

    limsup_{h downarrow 0} P(r<=M<=r+h)/h
       <= N e sec²(1) (K+1) P(M>=r).                            (3)

The marginal densities f_i are continuous: their conditional spherical
densities are continuous in r and uniformly bounded when kappa<=K, so
dominated convergence applies. The interval union bound therefore gives
the limsup at the specified r bounded by sum_i f_i(r), to which the
pointwise conditional hazard bound applies. This statement controls an
unweighted shell probability; a flux with an additional moving-boundary
weight still requires a separate comparison.

**Proved conditional application to one block pair.** In the block tree
gauge, write w_a, a in {0,1}^3, for the eight outgoing links in the pair's
direction, with w_0 the straight path P_0. Each of the sixteen two-link
paths contains exactly one w_a. The two paths at transverse offset zero
equal P_0: their remaining links are tree links in the two blocks. For
a!=0, after all links except w_a are frozen, each of the other fourteen
path distances d(P_x,P_0) is the distance from w_a to a fixed SU(2)
element, by bi-invariance. Therefore (2) and (3) apply with N=14 whenever
the conditional law of each free w_a is a von Mises--Fisher law with
concentration at most K. The fourteen variables need not be independent;
each pair sharing a private link is included separately.

For the ordinary Wilson law conditioned on a frozen exterior that leaves
these outgoing links free, and then on every other active link, the
single-link conditional is of that form: the terms containing w_a are
linear in its SU(2) quaternion coordinates. A finite bound K=O(beta_W)
follows from the finite number of adjacent plaquettes, uniformly in their
staples. This deduction is for the unconstrained law before conditioning on
the coarse averages. It does not assert that a fibre-conditional link has
the same law. The frozen exterior must indeed leave each chosen w_a free.

**Open.** The sharper proposed O(sqrt(kappa)+kappa r) bound is not proved
by (1). The block-pair deduction above does not itself identify the
constraint-disintegrated shell flux with the unweighted shell probability,
or establish the needed averaged transport speed on a moving fibre. The
inequality controls the distributional density of the
maximum, including ties; an identification of that density with the
constraint-disintegrated moving-surface flux remains to be proved. The
subsequent bad-event displacement H and environment control are still
open. Thus this calculation does not establish Conjecture B, (I1) or D.

**Computation.** `python3 d5/su2_hazard_check.py` integrates the explicit
radial density by composite Simpson quadrature for four concentrations,
five relative centres and three radii. The largest observed ratio
f(r)/[(kappa+1)P(R>=r)] is approximately 1.163, versus the proved bound
9.312. This floating-point grid is a diagnostic, not a proof or a uniform
optimization.
