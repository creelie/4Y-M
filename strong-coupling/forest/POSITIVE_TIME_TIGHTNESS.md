# A positive-physical-time ratio gives a normalised energy bound

## 1. A B/E interface without an unsmeared energy moment

**Proved (conditional spectral lemma).** Let sigma_a be positive spectral
measures on [0,infinity), and let

    C_a(t)=integral exp(-tE) dsigma_a(E).

Fix a physical time s>0, independent of the cutoff a. Suppose

    0<C_a(s/2)<infinity,    C_a(s)/C_a(s/2) >= c>0

uniformly in a. Define the probability measure

    dnu_a(E)=exp(-sE) dsigma_a(E)/C_a(s).

Then

    integral E dnu_a <= 2/(e s c) =: M_s,
    nu_a([R,infinity)) <= M_s/R,       R>0,
    C_a(s+t)/C_a(s) >= exp(-M_s t),    t>=0.

*Proof.* Since u exp(-u)<=1/e for u>=0,

    E exp(-sE) <= [2/(e s)] exp(-sE/2).

Integrate and divide by C_a(s). Markov's inequality gives tightness at
large energies; Jensen's inequality for the probability nu_a gives the
last bound. In particular, an unbounded unsmeared first moment does not
prevent this positive-time criterion from applying.

When the spectral representation comes from a centred observable O_a and
a nonnegative Hamiltonian, nu_a is the spectral law of

    exp(-s H_a/2) O_a Omega_a / sqrt(C_a(s)).

This vector is normalised. It is a time-smoothed state, not itself a
proof that a local continuum field exists. The estimates prevent its
spectral mass from escaping to infinite physical energy and make any
convergent normalised correlations nonzero.

## 2. What B must supply, and what this does not prove

**Sketch / open.** B can target a cutoff-uniform lower bound on this
specific fixed-physical-time ratio for a renormalised gauge-invariant
observable, with a positive spectral representation. Overall positive
multiplicative normalisation cancels from the ratio. The full interacting
and flowed correction, its remainder and its limits still need control;
the finite-box unflowed perturbative calculation does not establish c.

If s instead tends to zero in lattice units converted to physical units,
M_s can diverge. A lattice-step ratio is not enough. For sigma_a a unit
atom at E=1/a, the half-time ratio is exp[-s/(2a)] and violates the
uniform lower bound at fixed s. Normalising the variance alone therefore
does not repair energy escape.

**Scope.** No convergence of the Schwinger distributions, local-field
normalisation, OS axioms, vacuum uniqueness or Yang--Mills identification
is proved. Time smoothing is a spectral diagnostic and a conditional
nontriviality input. It does not replace the missing local continuum
construction. D8 and E remain open.

**Computation.** Finite positive mixtures at energies 1,3,1/a, for three
cutoffs, satisfy the moment and correlation inequalities; a pure escaping
atom gives the stated vanishing ratio. These floating checks are
calibrations only; the displayed inequalities follow analytically above.
