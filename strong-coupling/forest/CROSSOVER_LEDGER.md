# A finite multiscale error ledger

**Proved (abstract finite-step estimate).** Let (X_k, ||.||_k),
0 <= k <= N, be normed spaces, F_k and Fbar_k maps X_k -> X_{k+1},
and suppose trajectories u_{k+1}=F_k(u_k), v_{k+1}=Fbar_k(v_k) remain
in domains on which

    ||F_k(x)-F_k(y)||_{k+1} <= L_k ||x-y||_k,
    ||F_k(v_k)-Fbar_k(v_k)||_{k+1} <= epsilon_k,

with L_k,epsilon_k >= 0. Then, writing e_k=||u_k-v_k||_k,

    e_N <= (prod_{j=0}^{N-1} L_j)e_0
           + sum_{k=0}^{N-1} epsilon_k prod_{j=k+1}^{N-1} L_j.  (1)

*Proof.* The triangle inequality gives e_{k+1} <= L_k e_k+epsilon_k.
Iteration gives (1). Empty products equal one. The norms may change with
scale; (1) requires the stated operator bounds in those actual norms.

If L_k <= 1+ell_k, ell_k >= 0 and sum_k ell_k <= A independently of N,
then all products are <= exp(A), giving

    e_N <= exp(A)(e_0+sum_k epsilon_k).                             (2)

A finite crossover interval J contributes a multiplier prod_{j in J} L_j.
The multiplier must itself be uniformly bounded as the cutoff varies;
calling an interval 'finite' does not give such a bound if its length or
coefficients depend on the cutoff. More generally, (1) is the appropriate
weighted summability condition; sum epsilon_k < infinity alone is not.
For L_k=2, epsilon_0=1 and other errors zero, e_N=2^(N-1) is attainable.

## Relevant and contracting directions

**Proved (abstract split criterion).** Let nonnegative errors r_k and s_k
satisfy

    r_{k+1} <= (1+a_k) r_k + b_k s_k + u_k,
    s_{k+1} <= theta s_k + c_k r_k + v_k,  0 <= theta < 1.

Suppose sum a_k <= A, sum b_k <= B, sup c_k <= C,
U=sum u_k, V=sum v_k are finite on the interval and

    K := exp(A) B C/(1-theta) < 1.

Then, with R=max_k r_k and S=max_k s_k,

    R <= exp(A)[r_0 + B(s_0+V/(1-theta)) + U]/(1-K),
    S <= s_0 + (C R+V)/(1-theta).

*Proof.* Iteration of the second inequality yields the stated bound on S.
Applying (2) to the first gives R <= exp(A)(r_0+B S+U).
Substitute the S bound and solve for R. These conservative estimates
permit arbitrary finite lengths under uniform constants. They illustrate a
possible bootstrap; the hypotheses and the splitting are not verified for
a Yang--Mills RG map. Couplings often have marginal or relevant directions,
so blindly imposing contraction on the entire effective action is unjustified.

## Gauge-theory application

**Open / next.** Define normed spaces of gauge-invariant local and polymer
interactions, a blocking map, the relevant/marginal extraction and source
observables, and prove the actual constants in (1), or the split bounds,
including bad-domain changes and scale-dependent rescaling. Establish the
uniform crossover multiplier and match physical units to the infrared
mixing estimate. `WEIGHTED_SCHUR.md` is one Gaussian one-step input, not
a nonlinear RG map. No all-scale Yang--Mills control follows yet.
