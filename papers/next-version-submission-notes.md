# Next-version submission notes

Title: Lattice Mass Gaps and Continuum Criteria for Four-Dimensional Yang-Mills Theory

Authors in manuscript order: Deep Bhattacharjee; Ushashi Bhattacharya; James Montgomery Lockwood; Tristan Romanov

Submission route: Hexagon, v7 to v8. The platform version is independent of GitHub release numbering.

Existing record: [v7](https://hexagonmath.org/2610.00341v7). The v8 identifier is not claimed to exist before submission is accepted by the platform.

Research release: [v4.41.1](https://github.com/creelie/4Y-M/releases/tag/v4.41.1). Version DOI: [10.5281/zenodo.24045649](https://doi.org/10.5281/zenodo.24045649).

## Abstract to paste

We prove strong-coupling lattice mass gaps for compact simple gauge groups and a conditional transfer from a physical-scale cube estimate. An integrated response criterion controls spatial covariance row sums without a uniform relaxation gap. The required weak-coupling Wilson bounds and nontrivial four-dimensional continuum construction remain open.

## Changes in this revision

Adds a space-time integrated response estimate for weighted force covariance row sums. Earlier manuscript material is retained. The title is unchanged. The abstract is shortened. Current submission files replace superseded compiled copies on the main branch.

## Scope that must remain visible

The required Wilson response bound, the uniform cube estimate along continuum tuning, and a nontrivial continuum construction remain unproved. A finite-volume response implication is not a proof of those hypotheses.

## Files

- `ym4d.pdf` is the current compiled manuscript.
- `ym4dtex.zip` has `main.tex` at its root, with required figures and ancillary checks.
- `aiprovenance.md` records the origin and checking scope of this revision.
- `release-details.md` identifies the current release and DOI.

Compile with pdfLaTeX at least three times. Upload the TeX ZIP when the platform compiles sources. Do not upload superseded PDFs or use a new-version suffix as an already assigned identifier.

## Checks and declarations

The current source compiles locally. Finite ancillary checks do not prove the remaining existence or uniform estimates. This revision has not supplied a new Lean formalization of its analytic or geometric arguments. The submitting authors must confirm the human-understanding, authorship, licence and platform declarations for the current additions.


## Bibliography and layout in v4.41.1

Removes the two unpublished bibliography entries and their citation markers at the author’s request. Consolidates a duplicate record for the same cited preprint version. The surrounding text, mathematical statements and contributor credit are preserved. The compiled manuscript has 392 pages. See `bibliography-audit.md` for the source of every retained reference.
