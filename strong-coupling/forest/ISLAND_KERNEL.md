# Integrating hard island weights through a common collar kernel

## 1. Positive-kernel integration lemma

**Proved (finite-measure lemma).** Let K(x,y)>0 be measurable. Let a_0,a_1
and b_0,b_1 be nonnegative finite measures, and suppose all four numbers

    Z_ij = integral K(x,y) da_i(x) db_j(y)

are finite and positive. Suppose, for all relevant x_0,x_1,y_0,y_1,

    exp(-delta) <= K(x_1,y_1) K(x_0,y_0)
                    / [K(x_1,y_0) K(x_0,y_1)] <= exp(delta),
    delta >= 0.

Then

    |log Z_11 + log Z_00 - log Z_10 - log Z_01| <= delta.

*Proof.* Multiply the upper pointwise inequality by the product measure
da_1(x_1) da_0(x_0) db_1(y_1) db_0(y_0) and integrate. Tonelli factors
the two sides as Z_11 Z_00 and exp(delta) Z_10 Z_01 respectively.
The lower inequality follows in the same way. Positivity permits taking
logarithms. No interpolation between the a_i or b_j is required.

The measures may have disjoint hard supports, contain finite or countable
branch sums, or have large and uneven weights. Their masses do not enter
the estimate. They must be separate x-only and y-only measures sharing
one kernel. Finiteness and positivity of the four Z_ij are essential.

## 2. How the calm-chart estimate can supply the kernel bound

**Proved conditional deduction.** Suppose K(x,y) is the partition function
of one fixed calm collar domain, and all x and y range in convex source
sets X,Y within the same larger good chart. Suppose the mixed derivative
bound of `CONNECTED_NORMALIZERS.md` holds uniformly there:

    |partial_u partial_v log K| <= L(d) ||u|| ||v||,
    L(d)=2.08*10^14 lambda exp[-10^(-14)(d-6)],   d>6.

Integrating over the two source segments gives the lemma with

    delta = L(d) diam(X) diam(Y).

Thus arbitrary hard changes of separate island weights preserve this
connected-normaliser bound, **provided the common collar kernel and the
uniform calm-source description have actually been established**. This
is stronger than interpolating the island weights through good fields:
only the collar sources are interpolated, while island interiors are
integrated without such a restriction.

## 3. Binary-label consequence

**Proved (finite probability deduction).** Normalise the four Z_ij to a
positive binary joint law. Writing R=Z_11 Z_00/(Z_10 Z_01), the two
conditional odds for the second label differ by the factor R. The
largest difference between logistic probabilities whose log odds differ
by s is tanh(|s|/4), attained at log odds -s/2 and s/2. Consequently

    |P(B=1 | A=1)-P(B=1 | A=0)| <= tanh(delta/4),
    |Cov(A,B)| <= (1/4) tanh(delta/4).

This concerns a pair of labels with all other specified data frozen.
Integrating a shared environment afterwards can introduce an additional
covariance; the hidden-phase example in `ENVIRONMENT_OBSTRUCTION.md`
still applies. This lemma is not the environment bootstrap.

## 4. Exact gap in the Yang--Mills application

**Sketch / open.** Enlarge each bad island to absorb all rough constraint
supports and branch-dependent Jacobians. Cut through good layers, retain
the interface coordinates and all constraints, and seek a disintegration

    Z_ij = integral K_collar(x,y) da_i(x) db_j(y).

One must prove that (i) the same fixed product collar domain works in all
four sectors, (ii) all its source coordinates remain in the larger calm
chart, (iii) island measures include every branch and depend only on their
own interface, and (iv) the collar estimate is uniform in the frozen rest.
Neither factorisation nor these gauge-geometry properties is established
here. If the label change also changes the collar kernel or a common
constraint couples both island measures, section 1 cannot be applied.
The five-block good collar does not automatically discharge these checks.
Even after them, island entropy, uniform-exterior comparison and D8 remain.

## 5. CPU diagnostic

**Computation (floating, not a gauge-theory certificate).** A bounded scalar
chain with potential z^2+0.1z^4, nearest-neighbour coupling 1/2, and a
fixed product domain [-1,1]^n was integrated by 48- and 96-node Gauss
quadrature. Separate endpoint weights had disjoint hard domains and
two-peak mixtures. For n=1,...,9 the integrated log odds stayed below
the computed pointwise kernel cross-ratio range. At 96 nodes, for n=7,8,9
the kernel ranges were approximately 9.11e-7, 1.13e-7 and 1.41e-8;
the integrated log odds were 3.87e-7, 4.82e-8 and 6.01e-9.
An unbounded Gaussian chain supplied a separate exact-formula control.
Hard-domain quadrature errors are distinct from the smooth kernel check.
The lemma in section 1 is proved by integration, not by these numbers.
Standalone source is in `d5/frontier_diagnostics.py`; operational data
and raw outputs are retained outside Git.
