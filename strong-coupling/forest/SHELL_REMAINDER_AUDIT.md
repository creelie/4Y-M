# Audit: shell mass is not the whole large-field derivative

**Proved (scope correction).** Section 4 of the source
[`CALM_DATA.md`](https://github.com/romanov360/OpenGluon/blob/claude/cube-coefficient-calm-data/yang-mills/cube-coefficient/CALM_DATA.md)
defines the remainder Xi as the derivative of the *whole* large-path
contribution, centred by the original fibre expectation. It does not give
the representation as a sum only of kink-surface measures asserted in the
earlier sketch in `CALM_LINK.md`. The derivation of that representation
must therefore be regarded as open.

**Proved (elementary obstruction).** On [0,1], choose a smooth nonzero
real function phi supported in (1/2,3/4) with integral zero, and let
mu_t have density 1+t phi for sufficiently small |t|. Let A=[1/2,1],
and set F=phi. Then E_0 F=0, A and its boundary are fixed, and

    d/dt [integral_A F dmu_t - E_0(F) mu_t(A)] at t=0
       = integral phi² > 0.

There is no moving-boundary flux. The entire derivative is a bulk density
term supported in A. Thus a bound on surface flux alone cannot control a
general large-region remainder. A correct argument must account for the
bulk derivative of the action and density, any observable motion, the
derivative of the cutoff, and the normalization; cancellations must be
shown rather than presumed.

**Sketch (next check).** Starting from the exact source definition, write
one common change of variables for the large-path contribution on each
regular fibre. Identify the bulk score and transport terms and their
support, and separately identify any surface terms. An O(beta) local
pointwise bound on the bulk score times the shell mass would suffice for
bounded local F, but it does not automatically give separation decay or
the full Wasserstein bound. `AVERAGED_SHELL_MASS.md` controls an averaged
one-sided surface approximation under explicit assumptions; it does not
prove the missing bulk or fibre differentiation estimates. D remains open.
