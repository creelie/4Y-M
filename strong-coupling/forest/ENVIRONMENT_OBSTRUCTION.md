# Finite-intersection rarity does not control the environment covariance

**Proved (abstract counterexample).** Arbitrarily strong finite-intersection
rarity bounds do not imply spatial decay of the label environment, even
when the covariance within every fully label-conditioned law is zero.
This is a limitation of the proposed inputs alone, not a counterexample
to mixing of the Yang--Mills specification.

Fix q in (0,1). Sample a fair hidden bit S. Given S, let the labels B_x
be independent Bernoulli variables with parameter q/2 if S=0 and q if S=1.
For every finite set T of distinct sites,

    P(B_x=1 for all x in T)
      = [(q/2)^|T|+q^|T|]/2 <= q^|T|.

Thus the same finite-intersection domination used for the rarity input
holds at every order. Yet for distinct x,y at any distance,

    Cov(B_x,B_y)=Var(P(B_x=1 | S))=q^2/16.

The correlation has a positive distance-independent floor. If the physical
spins themselves are B_x, conditioning on the entire label set makes the
law deterministic, so every quenched covariance is zero. All the original
covariance then sits in Cov_omega(mu^omega F,mu^omega H).

The hidden-bit law is not asserted to be a finite-range lattice gauge law.
Its role is to show why the original Gibbs structure, connected partition
functions, or some further mixing input must be used in an environment
argument. Neither the chessboard moment bound nor a generic probability
lemma can supply that missing information automatically.

## An exact diagnostic identity (proved)

For a binary label B, condition on the remaining labels xi. Write
p=P(B=1 | xi), assume 0<p<1, and define
Delta_F=mu(F | B=1,xi)-mu(F | B=0,xi). Then

    Cov(F,1_B | xi)=p(1-p) Delta_F,
    Cov(E[F | B,xi],E[H | B,xi] | xi)=p(1-p) Delta_F Delta_H.

These follow by expanding the two possible label values. Small p by itself
only gives a small floor if Delta_F and Delta_H are bounded; it does not
make their product decay with the separation of the physical observables.
Normalising a covariance by a rare probability can also remove the apparent
small factor. A valid bootstrap must control the label-change expectations,
not just count rare labels.

For two binary labels with strictly positive joint conditional weights
Z_ij (unnormalised conditional partition functions), the connected odds
ratio is

    R = Z_11 Z_00/(Z_10 Z_01),
    log R = log Z_11+log Z_00-log Z_10-log Z_01.

The covariance numerator is p_11 p_00-p_10 p_01. Bounding nonnegative
intersection probabilities alone does not force log R to be small.
Controlling such connected normaliser differences through the good region
is a concrete next target. `CONNECTED_NORMALIZERS.md` proves one version
for continuous source changes inside the good product chart, not for hard
label flips or bad islands.

**Sketch / open:** prove connected normaliser bounds for actual label or
island changes, with the branch contributions and hard domains retained.
The example does not prove this impossible; it specifies an indispensable
extra ingredient. (I1), the environment bootstrap and D8 remain open.
