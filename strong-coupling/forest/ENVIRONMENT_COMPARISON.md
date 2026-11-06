# An environment comparison criterion from conditional odds

This note is a **proved finite binary-model criterion**, with an **open**
Yang--Mills application. It connects the two-island cross-ratio estimate
in `ISLAND_KERNEL.md` to a bound on the shared label environment. Its
hypothesis is uniform over all other labels; torus rarity alone does not
supply it.

## 1. Criterion and proof

Let V be a finite nonempty set, let mu be a strictly positive probability
law on {0,1}^V, and let

    c_ij = sup |mu(B_i=1 | xi) - mu(B_i=1 | xi')|,  i != j,
    c_ii = 0,

where xi and xi' differ only at j. Write C=(c_ij), and suppose
q=max_i sum_j c_ij < 1. For f define its coordinate oscillation
osc_i(f)=sup_{xi=xi' off i}|f(xi)-f(xi')|. Then

    |Cov_mu(f,g)| <= (1/4) sum_ij osc_i(f) [(I-C^T)^(-1)]_ij osc_j(g).

The inverse is the convergent nonnegative Neumann series. In particular,
for single labels x,y this gives |Cov(B_x,B_y)| <= (1/4)[(I-C)^(-1)]_xy.

*Proof.* Let P_i resample B_i from its conditional distribution, and let
P be the uniform average of the P_i. The local conditional covariance
identity and the one-coordinate range inequality imply

    |mu[f (I-P_i)h]| <= (1/4) osc_i(f) osc_i(h).

For j != i, couple the two i-site conditionals optimally when the exterior
changes only at j. This gives

    osc_j(P_i h) <= osc_j(h)+c_ij osc_i(h),
    osc_i(P_i h) = 0.

Thus osc(P h) <= M osc(h), where
M=(1-1/n)I+C^T/n and n=|V|. The infinity norm of C is q<1, hence the
spectral radius of M is <1; P^t h converges uniformly to mu(h), since its
coordinate oscillations vanish and mu is invariant. Telescoping yields

    Cov(f,g)=sum_{t>=0} mu[f (I-P)P^t g].

Apply the local estimate and sum M^t/n = (I-C^T)^(-1).
All sums are finite-dimensional and converge absolutely by the same
spectral-radius bound. This also fixes the orientation of C.

## 2. From hard-island kernels to C

**Proved conditional deduction.** Freeze all labels except i,j. If the
four positive conditional weights have a log cross-ratio bounded in
absolute value by delta_ij uniformly in the frozen labels, then

    c_ij <= tanh(delta_ij/4).

This is the binary logistic estimate of `ISLAND_KERNEL.md`. Thus, if the
common-collar factorisation and its cross-ratio estimate hold for every
frozen configuration and the resulting row sums are below one, section 1
controls covariance of the *fully integrated label law*, including the
shared environment. No integration over a random environment is treated
as if it preserved the original two-label odds.

For a graph distance, if sum_j c_ij exp(m d(i,j)) <= q_m < 1 uniformly,
then the path expansion and triangle inequality give

    [(I-C)^(-1)]_xy <= exp(-m d(x,y))/(1-q_m).

Consequently label covariances decay at least at that rate. More generally,
the first displayed estimate controls observables whose conditional means
have summable, spatially local label oscillations.

## 3. What remains to apply this

**Open.** The required bound must hold for *every* conditioning of all
other labels, including adjacent and large bad clusters. The existing
calm-window result does not establish this, and the row contraction is not
a consequence of chessboard rarity. A large cluster can invalidate a
common calm collar. Likewise, obtaining local oscillations of the
conditional expectations of physical observables requires branch-complete
transport or another comparison estimate. These are substantial inputs,
not bookkeeping. A graded cluster expansion may handle sparse exceptional
regions without requiring this uniform single-site condition, but its
geometric and mixing hypotheses would themselves need proof for this
constraint system. See `RESEARCH_DIRECTIONS_2026_11_04.md`.

## 4. Exact diagnostic

**Computation.** `python3 d5/environment_check.py`, run from this directory,
enumerates an asymmetric three-label law using rational arithmetic, checks
the covariance bound and verifies the different orientations of C_01 and
C_10. Its row sums are 2/15, 8/45 and 8/153. It also reproduces the
hidden-phase covariance 1/256 at q=1/4. This finite example checks arithmetic
and indexing; it does not verify the gauge hypotheses or prove the lemma.
