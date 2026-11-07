# Lattice Mass Gaps and Continuum Criteria for Four-Dimensional Yang-Mills Theory

Current research release [v4.41.0](https://github.com/creelie/4Y-M/releases/tag/v4.41.0). Version DOI: [10.5281/zenodo.24044466](https://doi.org/10.5281/zenodo.24044466).

We prove strong-coupling lattice mass gaps for compact simple gauge groups and a conditional transfer from a physical-scale cube estimate. An integrated response criterion controls spatial covariance row sums without a uniform relaxation gap. The required weak-coupling Wilson bounds and nontrivial four-dimensional continuum construction remain open.

## Current manuscript

- [PDF](papers/ym4d.pdf)
- [Compilable TeX ZIP](papers/ym4dtex.zip)
- [Next-version submission notes](papers/next-version-submission-notes.md)
- [AI provenance](papers/aiprovenance.md)
- [Release details](papers/release-details.md)
- [Hexagon v7 to v8 form text](papers/hexagon-submission-form.md)

## Scope

The required Wilson response bound, the uniform cube estimate along continuum tuning, and a nontrivial continuum construction remain unproved. A finite-volume response implication is not a proof of those hypotheses.

## Source and checking

Current source is in `strong-coupling/forest/`. Compile `main.tex` with pdfLaTeX at least three times. Ancillary scripts in `strong-coupling/forest/anc/` document their own finite scope. Historical verification records and independent research material are retained outside `papers/`. Earlier manuscript versions are available through release tags, rather than duplicate active PDFs and ZIPs.

Adds a space-time integrated response estimate for weighted force covariance row sums.
