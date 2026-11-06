# A diagonal and reflection-positivity ledger for continuum limits

## What bounded countable data give

**Proved (countable diagonal extraction).** Let A be a countable complex
*-algebra over Q+iQ with an antilinear involution Theta preserving A. Let
ell_n:A->C be normalized linear functionals, positive in the sense

    ell_n(F* F)>=0,             ell_n((Theta F)* F)>=0,  F in A_+,

where A_+ is a specified countable positive-time subalgebra, and suppose
the first inequality in fact holds for every F in A. If for every F in A
the numbers ell_n(F) are bounded, then a common subsequence converges for
every F in A to a normalized linear functional ell with the same two
positivity properties. If one fixed H in A has

    ell_n(H* H)-|ell_n(H)|² >= c > 0,

then the limiting variance is at least c.

*Proof.* Enumerate A. Successively extract convergent subsequences of the
bounded complex sequences, and take a diagonal subsequence. Algebraic
identities and each nonnegative real inequality pass to the limit. The
variance statement follows from convergence of the three entries. For
finite families F_i in A_+, the limiting reflection matrices
ell((Theta F_i)* F_j) are positive semidefinite over Q+iQ by the same
argument and hence over C by continuity of a finite quadratic form.

**Scope.** This proves only a positive functional on a countable algebra.
Boundedness for each lattice approximation, without bounds uniform in n,
does not suffice. A variance witness gives a nonconstant random variable
for ordinary positivity; it need not give a nonzero OS Hilbert-space vector,
whose squared norm is ell((Theta H)* H), not ell(H* H). A nonzero OS
witness needs a positive lower bound for that reflection norm, after
subtracting its vacuum component where appropriate.

**Proved (explicit reflection-positive example).** Let X_+ and X_- be
independent symmetric signs, let reflection exchange them, and let A_+
consist of functions of X_+ alone. Then for every F in A_+,

    E[conjugate(Theta F) F] = |E F|² >= 0.

The observable H=X_+ has ordinary centered variance one but reflection
norm E[X_- X_+]=0. Thus it represents zero in the OS quotient despite
being nonconstant in the original probability space. This finite example
isolates the distinction; it is not a lattice Yang--Mills model.

## What must still be supplied

**Open (continuum Yang--Mills).** Identify a common algebra of
renormalized gauge-invariant observables with stable positive-time support;
prove the required uniform moments and a positive centered reflection norm.
Extend the limit continuously to the appropriate test-function space,
prove full Euclidean covariance, symmetry, clustering, regularity and the
growth hypotheses of an applicable reconstruction theorem, and identify
the resulting theory with the desired Yang--Mills dynamics. The countable
algebra, diagonal extraction and positivity alone do not give tightness of
a random distribution or any of those analytic properties. Nor do they
guarantee uniqueness of the subsequential limit.

The distinction between ordinary and reflection norms matters for flowed
fields: a smearing kernel crossing the reflection plane can invalidate
positive-time support even when the unflowed lattice measure is reflection
positive. Compare `CONTINUUM_SYMMETRY_TEST.md` and
`CONTINUUM_GAP_INTERFACE.md`.
