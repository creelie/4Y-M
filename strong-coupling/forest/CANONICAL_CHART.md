# A quantitative graph chart around a calm reconstruction

**Matching-law follow-up (proved deduction):** `INDEPENDENT_LABELS.md`
labels only non-path physical links and coarse curvature. A global branch
bound d(d_k,c_k)<=14r makes path restrictions redundant. For good windows
with the stated collar, the actual conditional law equals this graph law.
The original all-physical-link label domain remains nonconvex
(`LABEL_DOMAIN.md`). Bad-side and environment estimates remain open.

**Proved (deduction from the paper's `smallfibre` and elementary analytic
bounds).** This note constructs a measurable background for an enlarged
window, and proves strong convexity on an explicitly defined product graph
chart. It does **not** identify this chart with a physical label-conditioned
law. In view of `ROUGH_BACKGROUND.md`, a coarse-curvature assumption is
necessary for this particular Hessian route.

## 1. Setting and constants

Use SU(2), Hilbert--Schmidt norm, and the tree gauge. Let W be finite.
As in `TOUCHING_CHART.md`, K is the outgoing and incoming pair set,
L contains A_W and all eight private links of each incoming pair,
Dpath={d_k:k in K}, and Z=L minus Dpath. Let B be the set of coarse bases
of pairs in K. Thus L is contained in A_B. Averages based in B but not in
K do not involve L and remain fixed automatically.

Fix **all** coarse averages c and write R=R(c). The frozen exterior must
already have these averages for every pair outside K, whose support is
disjoint from L. This compatibility is automatic for exteriors obtained
from the full conditional law; it is necessary for a nonempty full fibre. For every non-tree link
that occurs in a plaquette with a vertex in a block of B^+, or in a pair
of K, write the fixed exterior, when the link is off L, as

    U_e = R_e exp(z_e),      ||z_e|| <= t.

This is a hypothesis on the exterior near the window, not on distant
links. The following extension has the same proof and constants: frozen
non-path links must satisfy the displayed bound t, while frozen path
links may be within 44t of R. Path links do not enter any other pair's
constraint, so the fixed-point and determinant estimates only use the
non-path bound. The curvature estimate below uses the maximum 44t for
all physical links. This extension is used in `INDEPENDENT_LABELS.md`. Let D_t be the product ||z_e||<t over the free non-path links e in Z.
Assume

    0 < t <= rho/2,       rho=10^(-4),
    kappa_R + 176 t <= 10^(-10),

where every plaquette of R with a vertex in B^+ has chordal defect at most
kappa_R. A plaquette of R is either I or a coarse plaquette of c: checking
the four cases of the two offsets o_mu,o_nu proves this directly. Hence
this condition on R can equivalently be checked on a specified finite
set of coarse plaquettes. It is not implied by relative-link goodness.

The constrained measure below is Haar-disintegrated given c and the
exterior off L, then restricted to the particular graph Gamma(D_t).
No assertion about other components of the fibre is made.

## 2. Explicit inhomogeneous solution (proved)

For a pair k write its path link as ell=c_k exp(X_k). The two offset-zero
paths equal ell. Each of the other fourteen paths has the form
P_x(z), independent of every path variable, and

    A_x(z)=c_k^(-1)P_x(z)

is a product of at most two conjugated exponentials of non-path link
coordinates. The reference products are I. Thus, with the cut-off one,

    N_k(X,z) = (1/16) sum over fourteen paths log(exp(-X) A_x(z)),
    X_k = F_k(X_k,z) := -N_k(X_k,z).

The sums count repeated path functions with their multiplicities.
This equation implies V^sharp_k=c_k exactly.

For clarity the following estimates hold on the **complexified**
traceless matrix space, with its Hilbert--Schmidt norm. Put

    R0=44 rho,   s=exp(46 rho)-1.

On ||X||<=R0 and every ||z_e||<rho, each matrix inside the logarithm is
I+M with ||M||<=s. The series for exp and log give

    ||F_k|| <= (7/8)[-log(1-s)] < R0,
    ||partial_X F_k|| <= (7/8)(1+s)/(1-s) < 0.9.

Here ||D log_(I+M)-I||<=s/(1-s), and
||D_X[exp(-X)A]+I||<=s. These also give the sharper comparison

    ||partial_X N_k + (7/8)I|| <= (7/4)s/(1-s).

The fixed point is therefore unique in the closed R0 ball and is
holomorphic in the non-path coordinates. The iterates starting at X=0
are explicit continuous functions of c and the frozen exterior on the
admissible set. Their uniformly convergent limit is measurable. For real
su(2) inputs it stays in su(2). Each pair solves independently, since path
links are private. This proves existence without assuming a completion
of the conditional fibre in advance.

Apply the same inequalities with rho replaced by t (monotonicity of the
scalar series suffices). Then ||X_k||<=44t. The real solution stays well
inside the region where the cut-off is identically one. The map Gamma
uses these solved path links and the prescribed non-path links. It is
injective, smooth and has finite range. In particular Gamma(0), with the
frozen exterior retained, is the promised measurable background.
Every relevant link is within 44t of R, so every plaquette required by
`smallframes` has defect

    kappa <= kappa_R+176t <= 10^(-10).

This last estimate includes the frozen links in those plaquettes.

## 3. Derivatives and Haar disintegration (proved)

The complex bound ||X_k||<=R0 on the rho product gives Cauchy bounds at
real points where all non-path coordinates have norm at most rho/2.
For a direction v and m_k=max ||v_e|| over the non-path coordinates of k,

    ||D X_k[v]|| <= 88 m_k,
    ||D^2 X_k[v,v]|| <= 3,520,000 m_k^2.

The second bound uses 2 R0/(rho/2)^2. Norms here are Hilbert--Schmidt norms;
one-variable Cauchy along the complex line z+u v/m_k gives the bounds
also for vector-valued X. Each non-path coordinate belongs to at most two
pair supports, so sum_k m_k^2<=2 sum_e ||v_e||^2.

Let H_k=I+partial_X N_k, evaluated at the solution, as an operator on the
three-dimensional Lie algebra. On the whole complex rho product,

    H_k=(1/8)(I+E_k),    ||E_k||<=14s/(1-s)<0.07.

Define the holomorphic logarithm

    ell_k=Tr log(I+E_k).

It satisfies |ell_k|<=3[-log(1-0.07)]<1. Its real value differs from
log det H_k by the constant 3 log(8). Cauchy, with the same radius rho/2,
therefore gives

    |D^2 ell_k[v,v]| <= (8/rho^2) m_k^2,
    |D^2 sum_k log det H_k[v,v]| <= 1.6*10^9 ||v||_2^2.

The sign of the constant difference is immaterial to derivatives:
log det H_k = -3 log(8)+ell_k.

**Density identity.** Relative to Euclidean measure in the free z variables,
up to a constant independent of z, the graph-restricted constrained law is

    exp(-lambda Q(Gamma(z))) * product_e j(z_e) / product_k det H_k(z),
    j(z) = [sin(||z||/sqrt(2))/(||z||/sqrt(2))]^2.

Indeed the path-link Jacobian J_k in `TOUCHING_CHART.md` has determinant
exactly det H_k on the solution. To see the cancellation, write
A_X=integral_0^1 Ad(exp(sX)) ds and Psi_X=Ad(exp(-X)) A_X.
Differentiating exp(X)exp(N) at N=-X gives A_X(dX+dN), while the right
trivialised path derivative is Psi_X dX. Thus J_k is, up to fixed
adjoint isometries, A_X H_k Psi_X^(-1). Since det A_X=det Psi_X,
det J_k=det H_k. Positivity follows from ||E_k||<1.
Coarea cancels the graph area exactly as in `eq:fibrevolume`.

The effective potential is consequently

    V_eff = lambda Q(Gamma(z)) + sum_k log det H_k(z) - sum_e log j(z_e).

The last term has nonnegative Hessian. For example, with u=||z||/sqrt(2),
its radial second derivative is csc(u)^2-1/u^2>0; its tangential derivative
is positive since 1/u-cot(u)>0. The limiting values at zero are positive.
These assertions apply on the small real balls used here.

## 4. Strong convexity, uniformly in the window (proved)

Along a straight coordinate line z+hv let xi be the physical right
velocity of Gamma, and let a_e be its covariant acceleration for the
bi-invariant group metric. The differential of exp has norm at most one
on su(2), and differs from the identity, after trivialisation, by at most
exp(2||z||)-1<0.01. Hence

    0.98 ||v||_2^2 <= T(xi) <= (1+2*88^2) ||v||_2^2.

Differentiating the integral formula for the trivialised differential
of exp bounds a path-link acceleration by

    ||a_dk|| <= ||X_k''|| + 2||X_k'||^2 < 4*10^6 m_k^2.

For an independent link ||a_e||<=2||v_e||^2, so

    sum_e ||a_e|| <= (2+8*10^6)||v||_2^2.

Each link meets six plaquettes. If their defects are at most kappa,
its action gradient has norm at most 12 kappa. All constraints based in
B are fixed, since those outside K do not involve L. Extend xi by zero
to A_B and apply `smallfibre`(ii) at Gamma(z):

    Hess_ambient Q[xi,xi] >= (1/26-5.4*10^4 kappa) T(xi).

This uses the ambient geodesic Hessian on a constrained tangent vector;
it is not an identification with the Hessian in z coordinates. Adding
the acceleration term yields the actual coordinate estimate

    D^2(Q o Gamma)[v,v]
      >= {0.98(1/26-5.4*10^4 kappa)
           -12 kappa(2+8*10^6)} ||v||_2^2
      > (1/52)||v||_2^2.

The final inequality holds for kappa<=10^(-10). Together with section 3,

    Hess V_eff >= (lambda/52 - 1.6*10^9) I,
    Hess V_eff >= (lambda/104) I    if lambda>=2*10^11.

The paper's SU(2) convention is beta_W=4 lambda. This is a deterministic
chart estimate, not a working threshold for (I1) or an island expansion.
The product D_t is convex, so the convex-domain Brascamp--Lieb inequality
used in the proof of `chartdecay` gives

    Var(f) <= (104/lambda) integral ||grad_z f||^2 dmu_chart.

It also gives the corresponding spectral gap for this **restricted chart
law**. No variance or covariance bound for a larger containing event is
inferred from this restriction.

**Computation:** `d5/chart_bounds.py` checks the scalar constants above
using rational Taylor bounds with explicit remainders. Its output is
`d5/chart_bounds.out`. It is bookkeeping for this proof, not a formal
verification of the analytic or measure-theoretic steps.

## 5. Spatial decay for this specified chart law (proved)

For this section impose the slightly stronger condition

    kappa_R+352t <= 10^(-10).

Then the previous Hessian bounds hold throughout D_(2t), giving room
for the product-domain approximation used in `chartdecay`. Distances
below are periodic coarse l-infinity distances between the base blocks
of physical links or non-path coordinates.

First obtain an absolute upper Hessian bound. The calculation in the
proof of `smallfibre`(ii) bounds the error in the ambient geodesic Hessian
from 2Q_U in absolute value by 84 kappa T, and Q_U<=24T. The velocity
and acceleration bounds above therefore give

    ||Hess(Q o Gamma)|| <= 8*10^5.

For completeness the independent Haar term also has an absolute bound:
its analytic continuation is -2 log S(q), where
S(q)=sum_j (-1)^j(q/2)^j/(2j+1)! and q=sum_a z_a^2. On a complex rho ball,
|S(q)-1|<1/4 and |-2 log S(q)|<1, by their series. Cauchy gives an
absolute Hessian bound 8/rho^2=8*10^8 per independent coordinate.
Together with the logdet bound,

    ||Hess V_eff|| <= 8*10^5 lambda + 2.4*10^9 < 10^6 lambda.

All terms have range at most three in the coarse base-block metric:
a solved link uses non-path links from the two blocks of its pair;
physical plaquettes involve blocks at mutual distance at most one.
Thus a pair of independent coordinates entering the same local term
has distance at most three. At most 45 independent links are based in
any block, so each Hessian block row has at most

    N0=45*7^4=108045

possibly nonzero blocks. A block is a three-dimensional Lie-algebra
coordinate. Its operator norm is at most the norm of the full Hessian.
For any coarse 1-Lipschitz function phi of these coordinate labels, and
alpha=10^(-14), Schur's test gives

    ||e^(alpha phi) Hess V_eff e^(-alpha phi)-Hess V_eff||
       <= N0*10^6 lambda*(exp(3 alpha)-1)
       < lambda/208.

Apply the weighted Helffer--Sjostrand argument in the proof of
`chartdecay`, now with m=lambda/104. Its scalar generator acts diagonally
on the coordinate label and commutes with these weights. The conjugated
vector operator is accretive by at least lambda/208, as is its adjoint.
For product balls use the same smooth radial truncation and convex
penalty approximation as that proof, inside D_(2t). The truncation
Jacobian and additional diagonal Hessian commute with the weights; the
Schur bound on the off-diagonal part persists. Passing to the restricted
law on D_t gives, for smooth F,H,

    |Cov_chart(F,H)| <= (208/lambda)
       ||e^(-alpha phi) grad F||_L2 ||e^(alpha phi) grad H||_L2.

This is a deduction using the same finite-dimensional covariance
representation and product-domain approximation as the paper, not a
formal verification of that representation.

If physical observables have supports X,Y and gradients at most one per
physical link, their pullbacks under Gamma have gradient norm at most
sqrt(15489)*sqrt(|X|) and sqrt(15489)*sqrt(|Y|). Each support expands by at
most one coarse block. Choosing phi to be the distance to the expanded
support of H yields the explicit consequence

    |Cov_chart(F,H)| <= (3.3*10^6/lambda) sqrt(|X||Y|)
       * exp(-10^(-14) max(dist_coarse(X,Y)-2,0)).

The rate is very crude but independent of the window and its shape.
It concerns only Gamma(D_t), with the above exterior and coarse-curvature
hypotheses. It is not a comparison between different exteriors, nor a
covariance estimate conditional on the original physical good/bad labels.

## 6. What moved and what did not

- **Proved:** measurable inhomogeneous background selection for this enlarged
  conditional window when the exterior is close to a calm R(c).
- **Proved:** a finite-range product graph chart, its exact density, and a
  volume-independent strong-convexity estimate on the stated domain, and
  weighted covariance decay under the strengthened collar condition.
- **Sketch / open:** physical label conditioning also restricts solved path
  links. Pulling these restrictions back need not give D_t. The graph law
  above is a deliberately specified restriction of the fibre law.
- **Sketch / open:** good/bad labels for coarse curvature, constrained bad-side
  coefficients, and the environment covariance must still fit together.
  The specified chart estimate does not prove spatial decay for the actual
  mixture, (I1), D7, D5, (D), or multiscale composition D8.

**Computation (diagnostic, not proof):** `d5/chart_probe.py` solves a
nonabelian pair with a noncentral c and checks by finite differences that
the physical path-link determinant equals det H_k. The discrepancy is
below 10^(-7) at each tested scale. The largest test scale exceeds the
certified analytic radius; it is only a check of the determinant identity.
