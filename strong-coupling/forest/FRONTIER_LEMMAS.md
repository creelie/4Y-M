# Small analytic advances across A, B, C, D8 and E

This note collects finite analytical interfaces. **Proved** means the
finite mathematical deduction stated here, not a proof of the Yang--Mills
continuum conclusions. Computational artifacts are maintained separately.

## A. An exact attainable block response at zero exterior (proved)

For the three-link quaternion block with all exterior fields zero, write
kappa=beta_W=4 beta and w=y_2 y_3. The joint law of (y_1,w) is proportional
to exp(kappa Re(y_1 w)) times product Haar measure. The marginal of y_1
is uniform, and

    E[conj(w) | y_1] = A(kappa) y_1,
    A(kappa) = I_2(kappa)/I_1(kappa).

For a unit quaternion direction h, the test function

    f = <y_1,h> + <conj(y_2 y_3),h>

is 1-Lipschitz for the sum of the three link chordal distances: multiplication
and conjugation preserve those distances, and each link is charged at most
its own displacement. Its response to a field in direction h at y_1 is

    Cov(f,<y_1,h>) = [1+A(kappa)]/4.

This is an **attainable lower bound**, not the supremum over all Lipschitz
functions or exterior fields. Differentiating the S^3 partition function,
or its Bessel series, gives the stated A. The expansion is

    A(kappa)=kappa/4-kappa^3/96+O(kappa^5).

Thus the frequently used (1+beta)/4 is the first-order approximation to
this test response, not an exact finite-coupling identity. Rational positive
Bessel-series tails give approximately 0.271825746653093 at beta_W=0.351
and 0.272071902546725 at 0.355. Neither number establishes the upper bound
0.281690140845... required for the proposed improved endpoint. The full
supremum and the previously identified shared Lipschitz budget remain open.
This is a useful calibration, not an extension of the proved coupling range.

## B. Keep the interacting-measure insertions (proved identity)

Let a finite reference measure have sufficient dominated moments to
expand expectations, and let

    dmu_g proportional to exp(-g V3-g^2 V4) dmu_0,
    F_g=F2+g F3+g^2 F4,   G_g=G2+g G3+g^2 G4.

Assume parity makes F2,G2,V4 even and F3,G3,V3 odd. With cumulants under
mu_0, the connected covariance is C_g=C0+g^2 C2+o(g^2), where

    C0 = cumulant(F2,G2),
    C2 = cumulant(F4,G2)+cumulant(F2,G4)+cumulant(F3,G3)
         -cumulant(F3,G2,V3)-cumulant(F2,G3,V3)
         -cumulant(F2,G2,V4)+(1/2)cumulant(F2,G2,V3,V3).

*Proof.* Expand log E_0 exp(s F_g+t G_g-g V3-g^2 V4), take the mixed s,t
derivative at zero, and collect powers of g. Terms of odd total parity
vanish. Normalisation by the partition function is already included in
logarithmic cumulants. Domination justifies the derivatives.

For a ratio C_g(t)/C_g(t0), its order-g^2 coefficient is

    [C2(t)-R0(t) C2(t0)]/C0(t0).

Only a contribution proportional to C0 at **both** separations cancels.
The last two action-insertion terms need not have this form.

**Proved finite matrix example.** In a number-basis oscillator truncated
above level six, let H_lambda=H0+lambda x^4 and use the centred x^2
observable. At lambda=0 only the level-two excitation contributes to its
connected ground-state correlation. First-order perturbation gives

    E0'=3/4,    E2'=39/4,    (E2-E0)'=9.

The normalised correlation therefore has expansion

    R_lambda(t)=exp(-2t)[1-9 lambda t+O(lambda^2)].

New excited-state amplitudes are O(lambda), so their squared weights start
at O(lambda^2); the level-two amplitude cancels in the ratio. These statements
follow from finite-dimensional analytic perturbation and four ladder steps.
They illustrate an energy change that survives normalisation, not a theorem
about the Yang--Mills vacuum. The actual B2 infinite-volume, flowed-observable
calculation still needs its own vertices, renormalisation and limit control.

## C. What coverage counting cannot improve by itself (proved)

In the abstract witness lemma with overlap at most H, let H coarse labels
share one Bernoulli witness of probability p, and let each label be bad
exactly when that witness is bad. Then

    P(all H labels bad)=p.

Any universal domination q^|S| for this abstract class must satisfy
q>=p^(1/H). Thus the exponent loss 1/H cannot be removed using the witness
implication and overlap bound alone. Better gauge-field estimates need
additional geometry or joint information. This is optimality for the
abstract incidence model, **not** a statement that the current gauge
collar constants are optimal or that one bad plaquette causes every label.

## D8. Form stability is inherited; locality is an additional input (proved)

Let H0,H be positive definite finite matrices with

    (1-epsilon) H0 <= H <= (1+epsilon) H0,    0<=epsilon<1.

For a surjective linear coarse map Q define the effective quadratic forms

    S(v)=inf_{Qx=v} x^T H x,    S0(v)=inf_{Qx=v} x^T H0 x.

Taking infima proves (1-epsilon)S0<=S<=(1+epsilon)S0. This needs no
locality assumption. Repeating such bounds over scales costs products
of (1+epsilon_k) and (1-epsilon_k); a scale-independent fixed epsilon is
insufficient to keep both products uniformly bounded away from zero and
infinity. Summable errors, bounded uniformly below one, are a sufficient
condition. This is a perturbation interface, not an RG construction.

Form stability alone does not preserve exponential locality. On a path
of n sites let H0=I, u=(1,...,1)/sqrt(n), H=I+epsilon u u^T. Then

    H^(-1)=I-[epsilon/(1+epsilon)] u u^T.

Its endpoint off-diagonal entry has magnitude epsilon/[(1+epsilon)n],
not bounded by a uniform constant times exp(-m n) for any fixed m>0.
The perturbation is small in operator norm but spatially nonlocal.

A usable additional hypothesis is weighted operator control: if for
weights W=exp(alpha phi),

    ||W H0^(-1) W^(-1)||<=G0,
    ||W (H-H0) W^(-1)||<=epsilon/G0,

then the Neumann identity gives ||W H^(-1) W^(-1)||<=G0/(1-epsilon).
This identifies a norm an effective-action construction must control;
plain energy comparison does not supply it. No such all-scale Yang--Mills
estimate is claimed here.

## E. A concrete nontriviality input connected to B (proved conditional lemma)

Suppose a centred observable has a finite positive spectral measure sigma_a with

    total mass v_a >= v>0,
    integral E dsigma_a(E) <= M<infinity.

For t>=0, Jensen and monotonicity of x exp(-Mt/x) give

    C_a(t)=integral exp(-tE) dsigma_a(E)
          >= v_a exp[-Mt/v_a] >= v exp[-Mt/v] > 0.

Thus if these correlations converge, their limit is nonzero. This criterion
connects B's control of spectral leakage to E's nontriviality requirement.
It assumes the positive spectral representation and uniform physical-energy
moment bound. It proves neither existence of a subsequential field limit
nor the full OS axioms. A normalisation with fixed variance but uncontrolled
energy moments can still lose all its mass to arbitrarily high energies.

A uniform physical gap alone does not imply nontriviality: take a symmetric
two-state continuous-time Markov process with decay rate m>0, sample it at
time step a, and observe a times its +/-1 state. The transfer eigenvalues
are 1 and exp(-ma); the gap in physical units is m, while

    C_a(t)=a^2 exp(-mt) -> 0.

All positive-order moments of this bounded vanishing observable tend to
zero. This is a reflection-positive time-process example, not a 4D
Yang--Mills construction or a claim about its correct renormalisation.
It isolates why a nontrivial spectral-weight statement must accompany
any proposed continuum-gap argument.
