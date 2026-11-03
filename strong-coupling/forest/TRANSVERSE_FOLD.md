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
