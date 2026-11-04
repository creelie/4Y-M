# Counting weak-coupling steps: a D8 subproblem

## A conditional discrete estimate

**Proved (scalar recursion).** Fix b,C,U>0 with U>b+C/U and C/U<b.
Suppose u_0>U and, while u_k>U,

    u_{k+1}=u_k-b+r_k,              |r_k| <= C/u_k.             (1)

Let n be the first index with u_n<=U, and put delta=b-C/U>0. Then n is
finite, U-b-C/U < u_n <= U, and

    |n-(u_0-u_n)/b|
      <= (C/b) [1/U + log(u_0/U)/delta].                       (2)

*Proof.* Each step decreases u by at least delta and at most b+C/U, so the
first exit is finite and the exit interval follows. The telescoping identity
gives b n=u_0-u_n+sum_{k<n} r_k. For the error, monotonicity and
u_k-u_{k+1}>=delta imply, for k<n-1,

    1/u_k <= delta^{-1} integral_{u_{k+1}}^{u_k} dx/x.

The remaining term is at most 1/U; sum these inequalities and use u_{n-1}>U.
This proves (2).

For a fixed block factor L>1, (2) implies the explicit two-sided bound

    |log(L^n)-(log L)(u_0-u_n)/b|
      <= (C log L/b) [1/U + log(u_0/U)/delta].                 (3)

Thus the number of steps grows linearly in u_0 with a logarithmic correction;
the corresponding length is exponential in u_0 up to polynomial factors.
The constant inside those factors depends on the hypotheses and threshold.
The terminal effective coupling belongs to a fixed compact interval, but
this does not identify its effective action or prove mixing there.

## An explicit D8 decomposition

**Sketch (proof obligations, not established for Yang--Mills).** To apply
such counting to the actual blocking construction one needs, in order:

1. A gauge-covariant scale map on a domain of local interactions, sources
   and large-field polymers, with an exact relation to the original measure.
2. A nonperturbative running-coupling coordinate and a controlled remainder
   of the form (1), with constants uniform in volume and cutoff; extracting
   a coefficient from a formal expansion is not this estimate.
3. A stable induction bounding all other relevant and irrelevant directions,
   large-field weights, and source derivatives up to the exit scale. The
   dependence of their errors on u_k must be propagated by the actual RG
   derivative, not simply summed without amplification.
4. A crossover theorem covering the entire terminal set of effective
   actions, including all allowed defects and boundary environments. A
   finite interval in the scalar coupling does not make this set compact
   in the required interaction norm, nor imply an infrared mixing estimate.
5. A physical-units matching theorem: the lattice spacing, block factor,
   exit scale, observable renormalizations and infrared decay rate must
   yield a strictly positive cutoff-independent physical mass.

**Proved warning (accumulation).** Even if every one-step error is bounded
by C/u_k, the bound obtained from (1) is only O(log u_0), not uniformly
summable. For the exact recursion u_{k+1}=u_k-b, the sum of 1/u_k over
pre-exit steps grows at least a positive multiple of log u_0. Hence an
argument needing a uniformly bounded accumulated error must obtain extra
decay, cancellations, counterterm extraction, or a contraction-weighted
sum; one-loop-size remainders cannot simply be declared negligible.

**Proved (summable higher-order errors).** Under the same recursion, for
every p>1,

    sum_{k<n} u_k^(-p)
      <= U^(-p) + U^(1-p)/[delta (p-1)].                      (4)

For k<n-1, monotonicity gives
u_k^(-p) <= delta^(-1) integral_{u_{k+1}}^{u_k} x^(-p) dx;
bound the last term by U^(-p) and integrate to infinity. In particular,
if u=1/g² and a suitably normalized per-step defect is bounded by
A g^4=A/u², its unamplified sum is at most
A(U^(-2)+1/(delta U)), uniformly in the initial u_0. This makes extracting
the lower-order marginal term before summing the remainder a concrete
target. It does not control propagation factors from later scales: the
bound applies directly only if those factors are uniformly bounded, or
after an appropriate weighted estimate is proved.

The recursion is a diagnostic model. Neither (1) nor the existence of an
appropriate Yang--Mills coordinate is proved here; D8 remains open.

**Computation.** `python3 d5/crossover_count_check.py` checks the exit and
accumulation bounds for three initial values and remainders -C/u, 0, C/u.
For u_0=100, 1000, 10000 with b=1, C=1/2 and U=10, the zero-remainder
inverse-coupling sums are approximately 2.258, 4.557 and 6.859. This is a
floating-point diagnostic of the scalar statement, not evidence for a
Yang--Mills recursion.
