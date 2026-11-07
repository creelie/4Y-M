# Lattice Mass Gaps and Continuum Criteria for Four-Dimensional Yang-Mills Theory

The complete two-column REVTeX manuscript incorporates the slow Gaussian heat-bath mode, slow marked-force response, an explicit operator counterexample and a corrected gap-transfer lemma at equal physical times. The weak-coupling Wilson bounds and nontrivial continuum construction remain open.

[Complete PDF](revision-v4.40.0/ym4d.pdf) · [TeX ZIP](revision-v4.40.0/ym4dtex.zip) · [Submission package](revision-v4.40.0/ym4d.zip) · [Release](https://github.com/creelie/4Y-M/releases/tag/v4.40.0) · [Zenodo DOI](https://doi.org/10.5281/zenodo.24028596)

<!-- integrated-round26 -->

# Lattice Yang–Mills Mass Gaps and the Continuum Cube Criterion

The Gaussian block action has a sharp spatial decay rate of order b^(-1/2). A constrained heat-bath estimate gives a quantitative sufficient condition for the Wilson covariance bound. James Montgomery Lockwood is credited for the internal-fibre formulation that motivated this approach. The weak-coupling Wilson estimates and the nontrivial continuum construction remain open.

[PDF](revision-v4.39.0/YangMills_Hexagon.pdf) · [TeX ZIP](revision-v4.39.0/hexagontex.zip) · [Release](https://github.com/creelie/4Y-M/releases/tag/v4.39.0) · [Version DOI](https://doi.org/10.5281/zenodo.24025857)

<!-- current-2026-11-07-revision -->

## Earlier revision and repository guide

## Earlier revision and repository guide

# One Cube Estimate Implies the Yang–Mills Mass Gap

Deep Bhattacharjee, Ushashi Bhattacharya, James Montgomery Lockwood and
Tristan Romanov.

The paper is `strong-coupling/forest/main.tex` (REVTeX, two columns; PDF at
`preprints/preprintshexagon/YangMills_Hexagon.pdf`), with a one-column
version in `strong-coupling/forest/cmp/main.tex`. Both build with
`latexmk -pdf main.tex`. An earlier version is the Hexagon preprint
[2610.00341](https://hexagonmath.org/2610.00341v6).

It proves a unique state and an explicit transfer-matrix gap for lattice
Yang–Mills theory in a forest gauge (SU(2) to β_W = 0.348, SU(3) to 0.52,
every compact simple group at small coupling), and that a Dobrushin–Shlosman
estimate on one cube of fixed physical size gives a mass gap in every
continuum limit. Towards that estimate it labels the blocks of Bałaban's
averages and controls the islands of bad blocks: their weights factorise
through kernels on the good layers between them, and past an explicit
coupling they are rarer than their entropy. The cube estimate along the
continuum tuning and a nontrivial continuum limit are not proved.

- `strong-coupling/forest/`: the paper, its figures, scripts and notes.
- `strong-coupling/gauge-fixed/`: the certified constants (`blocks/`, `su3/`).
- `YangMills/`: the Lean 4 / Mathlib files (`lake exe cache get && lake build`).
- `preprints/`: the Hexagon PDF (`preprintshexagon/`), the MDPI PDF
  (`preprintsmdpi/`) and their TeX sources (`texsources/`).
- `YangMills_masterzipfinal2.zip`: the repository at the current release with
  the submission notes of the paper (`hexagon-package/`).

Archived on Zenodo under the concept DOI
[10.5281/zenodo.23444338](https://doi.org/10.5281/zenodo.23444338), which
resolves to the latest release. License: MIT.
