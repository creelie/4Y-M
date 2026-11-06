# A transverse coordinate removes the explicit scalar path fold locally

**Proved (local scalar example).** The fold in `BAD_BRANCHES.md` is a
failure of the path-only coordinate, not a singularity of the scalar
two-variable constraint. Integrating one suitable non-path variable
regularises a smooth local contribution to its conditional measure.
This is not a uniform bad-side SU(2) transport estimate.

Let the exceptional non-path angle be t and the other six angles be 0.49.
On the rectangle

    0.04 <= q <= 0.18,    -0.54 <= t <= -0.46,

the unique maximal relative angle is q-t, and the cut-off stays in its
smooth transition branch. Thus the same average as in `BAD_BRANCHES.md` is

    f(q,t)=q+[2-4 sqrt(2) sin((q-t)/2)](t+2.94-7q)/8.

Here (q-t)/2 lies in [1/4,9/25], and S=t+2.94-7q >= 57/50.
For the upper transition endpoint, sqrt(2)<99/70 and the alternating
Taylor upper bound for sine give

    2 sqrt(2) sin(9/25)
      < 1946158731/1953125000 < 1.

The lower endpoint exceeds 1/2 by sqrt(2)>7/5 and
sin(1/4)>1/4-(1/4)^3/6. Hence 0<chi<1 throughout the rectangle.
Differentiation gives

    f_t = [chi+2 sqrt(2) cos((q-t)/2) S]/8
        > 466431/1250000 > 1/3.

The last rational bound uses cos u >= 1-u^2/2. In particular, t can
serve as the solved variable even where f_q=0.

For any smooth nonnegative density w(q,t) supported in the rectangle,
with a smooth zero extension, change variables (q,t) -> (q,gamma=f(q,t)). For
fixed q there is at most one t=T(q,gamma). The pushforward density and
the conditional observable numerator are

    Z(gamma)=integral w(q,T(q,gamma))/f_t(q,T(q,gamma)) dq,
    N_F(gamma)=integral F(q,T(q,gamma))
                         w(q,T(q,gamma))/f_t(q,T(q,gamma)) dq.

Interpret the integrands as zero outside the image of the support. The
smooth zero extension makes these extensions smooth. Thus Z and N_F are smooth
for smooth F, and N_F/Z is smooth wherever Z>0. No 1/|f_q| factor appears.

**Computation.** At t=-1/2 the path derivative vanishes near
q=0.074119614463, gamma=0.169782627016. With the product bump density

    w=bump((q-0.11)/0.07) bump((t+0.5)/0.04),
    bump(x)=exp[-1/(1-x^2)] for |x|<1, otherwise zero,

128-, 256- and 512-node quadrature give Z approximately 0.016032740833.
The 256/512 difference is below 3e-13. Central differences for the
conditional mean of q approach a finite slope near 1.04053 as the
gamma increment decreases. An independent two-dimensional Gaussian
delta approximation approaches the same density; halving its width
reduces the leading error by about four. These are floating diagnostics,
not rigorous numerical enclosures. Standalone source is in
`d5/frontier_diagnostics.py`; raw data remain private.

**Sketch / next.** Cover a full bad SU(2) fibre using the eight-private-link
submersion, choose nondegenerate coordinate minors and a partition of
unity, then bound the integrated scores and chart changes. Hard label
boundaries introduce flux terms that the compact smooth bump avoids.
The scalar example neither controls those terms nor the full branch sum,
and gives no volume- or exterior-uniform Wasserstein coefficient.

## Hard-boundary obstruction even with a transverse coordinate

**Proved (scalar counterexample).** A uniformly nonzero transverse
derivative by itself does not provide a uniform response bound after a
hard restriction. Let (q,t) range over [0,1]^2 with Lebesgue density and
condition on f(q,t)=q^2+t=gamma. Here partial_t f=1 everywhere, so f is a
submersion and the coarea integral can be taken in t without a singular
Jacobian. For 0<gamma<1 the unnormalised density is

    Z(gamma)=integral_0^1 1_{0<=gamma-q^2<=1} dq=sqrt(gamma).

The conditional expectation of q is sqrt(gamma)/2. Hence both Z' and
the derivative of this bounded observable diverge as gamma decreases to
zero, despite complete transverse solvability. For gamma in (0,1), the
conditional law is uniform on [0,sqrt(gamma)], so its Wasserstein-1
distance to the law at gamma=0, interpreted by the weak limit delta_0,
is sqrt(gamma)/2. Therefore even uniform Lipschitz transport fails at this
hard boundary. The mechanism is the moving endpoint q=sqrt(gamma), where
the level surface meets t=0 tangentially; it is not a bulk coarea fold.

**Scope / next.** This toy example does not establish a failure for the
actual SU(2) fibre. It shows that a full-constraint submersion and a
bounded coarea factor alone cannot settle the hard-label transport task.
An application needs quantitative transversality to the active boundary,
a cancellation or a bound after averaging over the boundary, with a
normalizer controlled on the relevant parameter range. In particular,
refining a numerical branch solver near a fold cannot certify such a
uniform estimate without controlling these boundary terms.
