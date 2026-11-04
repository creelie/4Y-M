# Yang–Mills Mass Gap at Strong Coupling and One Cube

Deep Bhattacharjee, Ushashi Bhattacharya, James Montgomery Lockwood and
Tristan Romanov.

The paper is `strong-coupling/forest/main.tex` (REVTeX, two columns; PDF at
`preprintYangMillsForestGauge/YangMillsForestGauge.pdf`), with a one-column
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
- `YangMills_masterzipfinal2.zip`: the repository at the current release with
  the Hexagon package of the paper.

Archived on Zenodo under the concept DOI
[10.5281/zenodo.23444338](https://doi.org/10.5281/zenodo.23444338), which
resolves to the latest release. License: MIT.
