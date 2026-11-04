# Averaging a shell weight over coarse fibres

## Normalizers cancel in the disintegration

**Proved (measure-theoretic lemma).** Let (X,mu) be a probability space,
f:X->Y measurable, rho=f_*mu, and {mu_y} a disintegration of mu over f.
Let h:X->R and q:X->R be measurable, with 0<=q<=Q. For a radius r and
epsilon>0 define

    M_epsilon(y) = epsilon^(-1) integral 1_{r<=h<=r+epsilon} q dmu_y.

Then

    integral M_epsilon(y) rho(dy)
      = epsilon^(-1) integral 1_{r<=h<=r+epsilon} q dmu
      <= Q epsilon^(-1) mu(r<=h<=r+epsilon).                      (1)

Consequently, for any fixed sequence epsilon_n downarrow 0 and
M_*(y)=liminf_n M_{epsilon_n}(y),

    integral M_*(y) rho(dy)
      <= Q liminf_n epsilon_n^(-1) mu(r<=h<=r+epsilon_n).         (2)

This follows from disintegration and Fatou's lemma. If the actual fibre surface measure has
total mass M(y) and one proves M(y)<=M_*(y) almost everywhere, the same
bound applies to M. In a smooth transverse setting, equality with a
coarea-defined one-sided surface mass is a standard additional regularity
condition; it is not automatic for arbitrary measurable h and fibres.
There is no lower bound on fibre normalizers in (1)--(2): they are already
absorbed in rho, the actual distribution of the coarse averages.

## The path maximum and its speed

**Proved (deterministic speed bound).** Let h be the maximum, over a fixed
pair's sixteen paths, of the Hilbert--Schmidt chordal distances from the
straight path. Regard h as a function of all active links with the product
Euclidean metric formed from their Hilbert--Schmidt tangent norms. Each
two-link path and the straight path contain at most four link occurrences
in total. By bi-invariance and the triangle inequality, changing the
configuration by link increments of lengths d_e changes each path
distance by at most the sum of those four lengths. Each length is at most
the product distance, giving Lipschitz constant at most 4 even if an
occurrence is repeated.
Taking a maximum preserves this bound. Thus along any absolutely
continuous path with product speed <=V,

    |d h/dt| <= 4 V                                                (3)

at almost every differentiability time. In particular, if the transported
velocity is supported in these links and has product norm <=2.20, its
shell-normal speed is at most 8.80 in this metric. This constant is a
coarse bound; the normalization of the velocity norm must be checked
against the norm used in the cited transport result.

**Proved conditional combination.** Suppose the assumptions of the
block-pair result in `SU2_RADIAL_HAZARD.md` hold for the unconstrained
exterior law mu, with concentration K, and the chordal shell radius r
corresponds to angular radius in (0,pi/2). Then (3) of that note and the
chordal change of variables give

    limsup_{epsilon downarrow 0}
      mu(r<=h<=r+epsilon)/epsilon
       <= 14 e sec²(1) (K+1) mu(h>=r).                            (4)

Here the angular-to-chordal derivative is >=1 on that range. If q is a
measurable shell-speed weight bounded by 8.80, (2) yields

    integral M_* d rho
      <= 123.2 e sec²(1) (K+1) mu(h>=r).                         (5)

For an O(beta_W) bound on K this is the O(beta_W) averaged shell-mass
estimate targeted by `CALM_LINK.md`, subject to the stated exterior,
normalization and metric assumptions. It is a bound for one-sided shell
approximations. It does not by itself identify the moving constraint's
signed flux, establish that its total variation is controlled by M_*, or
bound the displacement of a fibre law conditioned on the shell. Those
are separate open steps needed for (I1) and D.

## A local-observable fallback

**Proved (limited consequence).** Let nu be a finite signed measure on a
fibre and mu_y a probability measure. For bounded F,

    |integral (F-mu_y F) d nu| <= osc(F) ||nu||_TV.                (6)

If F depends on k SU(2) links and is 1-Lipschitz for the sum of their
Hilbert--Schmidt chordal distances, osc(F)<=2 sqrt(2) k. Thus, if the
physical signed shell flux is shown to obey the mass bound in (5), its
contribution to any such local F is bounded uniformly in the ambient
volume by 2 sqrt(2) k times the right side of (5). No Wasserstein bound
on the entire fibre is needed for this limited claim.

This does not provide decay with separation from the shell. A volume-free
local bound alone cannot be summed over arbitrarily many distant shells
and does not imply the cube coefficient or environment mixing required for
D. The claim is conditional on the missing signed-flux identification.
