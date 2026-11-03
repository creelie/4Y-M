# Research continuation: status, source and reproducibility

This directory includes the public research continuation from
`romanov360/4Y-M:claude/yang-mills-mass-gap-wf5h8x`. Start with the last
section of `ROADMAP.md`, then `D_STATUS.md` and `FRONTIER_PLAN.md`.
Older entries are retained as a research history; later corrections
supersede earlier proposed thresholds and heuristic diagnoses.

The complete upstream contribution is
[draft PR #160](https://github.com/creelie/4Y-M/pull/160).

[PR_SUMMARY.md](../../../PR_SUMMARY.md) combines the eight-program overview with
the item-by-item A1--A4/B/C/D/E/I progress and its relation to the Clay
problem.

## Results and open fronts

| Front | Contribution | Remaining obligation |
|---|---|---|
| A/G1 | four kernel-checked Lean modules for finite comparison, cost minimality, the measure bridge and compact-space continuous-cost minima | measurable conditional coupled updates, the specification-dependent estimates and the requisite limit arguments |
| A2/A3 | response diagnostics, exact zero-field Bessel test and block-shape computations | the first-order ceiling estimates and Monte Carlo responses do not prove a global optimum or a new threshold |
| B | lattice free-vacuum analysis, exact finite-box interacting coefficients, nonlinear flow formulas and finite local flow checks | complete flowed interacting coefficients, uniform remainders, infrared and continuum limits |
| C/D/I | chart construction, independent-label law matching, curvature witnesses, branch/domain obstructions and conditional collar-kernel integration lemma | bad-side transport, hard-boundary contributions, actual island/collar factorisation and environment comparison |
| D8 | finite weighted locality and Schur-complement interfaces | all-scale nonlinear large-field estimates with controlled accumulated errors |
| E | conditional positive-time spectral tightness and nontriviality criteria | Yang--Mills correlation-ratio bounds, local-field construction, convergence and reconstruction |

`A4_LEAN.md` and `COUPLING_INTERFACE.md` delimit the formal results.
`B_WICK.md`, `NONLINEAR_FLOW_CHECK.md`, `INDEPENDENT_LABELS.md`,
`BAD_BRANCHES.md`, `ENVIRONMENT_OBSTRUCTION.md`, `ISLAND_KERNEL.md`,
`WEIGHTED_SCHUR.md` and `POSITIVE_TIME_TIGHTNESS.md` describe the latest
analytical interfaces and their hypotheses. No new (I1), D8 or continuum
existence theorem is asserted by this contribution.

## Reproducible scientific programs

- `d5/README.md` indexes the exact/rational checks and numerical chart,
  branch, collar and spectral diagnostics. Six standard-library scripts
  reproduce their small mathematical output summaries byte-for-byte.
- `wick/README.md` documents the finite-box exact and floating Wick
  implementations, the independent covariance jets, nonlinear compact-group
  calibration and actual finite spatial-action derivative check. These
  source programs are included; generated `result.json` files are not.
- `rare/` contains finite rarity calculations; `slice/` and `slice/sprime/`
  contain the slice-ratio and free-vacuum calculations. Some historical
  analyses require separately supplied raw arrays or the external
  `big-bang-monotile/yang-mills-continuum` checkout. They are not represented
  as standalone reruns from this repository alone. The external module
  for `leak_split.py` must be provided on `PYTHONPATH`.
- `../gauge-fixed/blocks/sharp/`, `feedback/` and `shapes/` retain the
  fork's certificate checks and numerical experiments alongside upstream's
  canonical published programs and certificates. Treat Monte Carlo and
  first-order shape ceilings according to their stated limitations.

Run generated-output programs from separate scratch directories. The
small tracked mathematical certificate summaries are distinct from raw
arrays, job logs, submission tooling and machine-specific configuration;
the latter remain excluded from the repository.

## Upstream reconciliation

The full contribution is based on upstream `main` at `fdd2dc5` and includes
the research branch through `e8c3554`. Both current paper sources, release
metadata, PDFs and the release archive are preserved from upstream.
Shared SU(2)/SU(3) endpoint improvements were already upstream and are not
claimed again. Upstream's newer fixed-set, coarse-defect and covariant-tilt
discussion remains intact; the continuation's conditional good-chart
estimates do not discharge that paper's open uniform-exterior/all-scale
input. Other new notes remain research-level material for review.

Canonical SU(3) certificate logs remain upstream's versions. Alternate
fork summaries are retained as `consts_log_fork.txt` and
`cover_log_fork.txt`, with their individual run summaries alongside them.
The fork's one-line alpha output is already contained verbatim in
upstream's two-line alpha output, so no duplicate is needed.

**Validation scope.** The four Lean modules were checked against the same
pinned v4.33.0 dependencies as upstream; the complete default target was
not rebuilt. Integration checks parse all contributed Python sources,
reproduce selected existing certificates, check source coverage and
relative links, and compare protected upstream files byte-for-byte.
They do not independently reprove every analytical note or rerun the
expensive historical Monte Carlo and interval sweeps.

**Integration checks completed.** All 123 paths changed by the research
branch since the common ancestor are present: 104 retain the fork's exact
bytes and 19 are reconciled for upstream content, current scope corrections
or source locations, including one trailing-blank-line cleanup. No fork
path is missing. The alternate SU(3) logs
retain the fork's original bytes under their new names. The eight newly
included scientific programs ran successfully from their public layouts;
both floating Wick methods matched the exact rational small-box result,
and the independent core Lean block certificate compiled. Six existing
standard-library diagnostics reproduced their committed summaries exactly.
The program validation used two CPUs for about 18 seconds and no GPU.
All Python sources parse and the contributed prose links resolve.
