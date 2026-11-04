# A necessary continuum symmetry test

**Proved (counterexample to symmetry-by-subsequence).** Define a centered
stationary real Gaussian field on R^4 with covariance

    K(x) = (1/4) sum_{j=1}^4 cos(x_j).

This is positive definite: it is the Fourier transform of equal positive
point masses at the eight frequencies +/-e_j. It is invariant under
translations, coordinate permutations and sign changes, hence the full
hypercubic symmetry group, but it is not invariant under arbitrary rotations.
For example K(t,0,0,0)=(3+cos t)/4, while rotating the vector through 45
degrees in the first two coordinates gives K(t/sqrt(2),t/sqrt(2),0,0)
=(1+cos(t/sqrt(2)))/2; their fourth Taylor coefficients differ.

The constant sequence of its restrictions to increasingly fine lattices
therefore has a tight continuum limit with exact lattice rotational
symmetry, but that limit is not Euclidean rotation invariant. This is an
obstruction to deducing full rotations from discrete symmetry and
subsequential convergence alone. It is not asserted to satisfy reflection
positivity or to describe Yang--Mills.

**Sketch (constructive target).** In addition to compactness, prove an
asymptotic rotation-defect estimate for a separating class of renormalized
local observables: for each fixed rotation and finite collection of test
functions, the difference between the rotated and unrotated lattice
correlations tends to zero. Separately establish reflection positivity for
observables with the required positive-time support, regularity/growth,
clustering and a nonzero renormalized state. Gradient-flow smearing across
the reflection plane does not automatically satisfy the support hypothesis.
No reconstruction or nontrivial continuum theory is proved here.
