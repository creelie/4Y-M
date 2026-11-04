# Finite large-field diagnostics and certificates

The scripts in this directory support the accompanying mathematical notes.
The small tracked `.out` files are mathematical certificate/diagnostic
summaries. Raw arrays, job logs and submission tooling are not included.

**Standard-library checks:** `bad_branch.py`, `chart_bounds.py`,
`coverage_checks.py`, `rough_certificate.py`, `touching_checks.py` and
`scales.py` run without numerical packages. Their current outputs were
reproduced byte-for-byte during upstream integration. The rough-field
certificate checks exact integer identities; scalar inequalities in
`chart_bounds.py` use rational arithmetic. These checks do not substitute
for the analytic arguments or extend their hypotheses.

**Numerical diagnostics:** the other chart/domain scripts use NumPy and,
where imported, SciPy; threshold/ledger scripts use mpmath. The stored
threshold computations predate later domain and environment obstructions;
they are conditional ledger calculations, not certified working island
thresholds. See the latest section of `../ROADMAP.md` for current status.

Three additional standalone scientific programs accompany the newer work:

- `branch_controls.py --device cpu --reference --batch 64 --batches 8`
  uses NumPy/JAX to check constraint Jacobians against finite differences.
  A path-only singularity does not establish divergent integrated transport.
- `frontier_diagnostics.py` uses NumPy/SciPy for the scalar transverse fold,
  separate hard island weights through a bounded scalar collar, and finite
  spectral-mixture checks. These are the diagnostics of
  `../TRANSVERSE_FOLD.md`, `../ISLAND_KERNEL.md` and
  `../POSITIVE_TIME_TIGHTNESS.md`.
- `frontier_calibration.py` uses the standard library for the Bessel,
  oscillator, overlap, locality and normalisation calibrations of
  `../FRONTIER_LEMMAS.md`.

These three programs each write `result.json` in their current working
directory; use separate scratch directories to keep generated data out of
the source checkout. They contain no scheduler or cluster transport code.
Numerical diagnostics alone close neither the bad-side transport problem
nor the environment recursion or multiscale estimate.

**Island thresholds (mpmath, standard library for the enumerations):**
`island_threshold.py` computes the block-scale threshold of the island
theorem (the largest root of the chessboard exponent against the entropy
of the collar graph); `composed_threshold.py` the thresholds at every
scale of the composed cut-off average with the box sizes of the fine
witnesses, and the scale of the defect propagation at the same couplings;
`witness_geometry.py` enumerates the witness boxes and their overlaps;
`abelian_growth.py` iterates the cut-off average on an abelian field of
constant curvature and reproduces the factor four per step;
`grid_check.py` compares the two scales on a grid of couplings. The
tracked `.out` files are their outputs.
