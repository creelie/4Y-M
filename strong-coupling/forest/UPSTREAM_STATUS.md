# Upstream comparison and first PR candidate (2026-11-03)

**Superseding full-contribution preparation.** At the user's request the
full public research continuation is now prepared on
`romanov360/4Y-M:upstream/all-research`, based on upstream `fdd2dc5` and
the research branch through `e8c3554`. The six-file proposal below is a
historical first candidate, not the requested final scope.
`RESEARCH_INDEX.md` describes the complete contribution, its scientific
source programs, validations and conflict resolutions. Current paper
sources and release artifacts remain upstream's versions. Private
submission tooling, machine configuration, raw data and audit records
remain excluded. GitHub authentication became available after preparation,
and the full contribution is now open as
[draft PR #160](https://github.com/creelie/4Y-M/pull/160).

[Review the complete contribution](https://github.com/creelie/4Y-M/pull/160).

**Assessment.** Compared research head `6e0608b` with
`creelie/4Y-M:main` at `fdd2dc5`. Their common ancestor is `5ca904f`;
upstream has 16 commits and the research branch 43 commits beyond it.
These counts describe that snapshot, before this comparison note.

Upstream is at v4.31.0 and already includes the SU(2) 0.351 and SU(3)
0.77 endpoint improvements, newer paper sections on coarse defects and
covariant forest tilts, and release metadata/PDF updates. Those should
not be replaced by the older paper and release state on the research
branch. The full branch merge has seven conflicting paths: README,
both paper sources, the SU(3) block README and three certificate/output
logs. The broad research diff also includes exploratory work requiring
its own assessment against the newer paper.

## Focused formalisation candidate

**Prepared / pushed:** `romanov360/4Y-M:upstream/g1-coupling`, commit
`e7a52e4`, starts directly from upstream `fdd2dc5` and adds one commit.
It has six changed files, 723 added lines and a clean test merge:

- four Lean modules: `DobrushinComparison`, `CouplingMinimality`,
  `WassersteinBridge` and `CouplingMinimizer`;
- their imports in `YangMills.lean`;
- a self-contained `strong-coupling/forest/A4_LEAN.md` describing the
  checked results, validation and undischarged G1 inputs.

**Validation.** All proof declarations match the sources already compiled
against Lean/Mathlib v4.33.0. The only bridge difference is a comment
removing references to fork-only notes. Upstream's toolchain, lakefile
and dependency manifest are byte-identical to those checked. The complete
default target was not rebuilt. The PR excludes the exploratory notes,
paper edits, numerical outputs, private tooling and operational data;
its history is upstream plus one focused commit.

**Status:** the candidate branch is public and ready for upstream review,
but no PR has been opened. `gh pr create --draft` could not authenticate:
GitHub CLI has no API login on this node. Git SSH push succeeded. The
prepared title is *Formalize finite Dobrushin comparison and compact
coupling minima*. The commit body includes the proposed description.

[Open the upstream comparison / PR form](https://github.com/creelie/4Y-M/compare/main...romanov360:upstream/g1-coupling?expand=1).

## Remaining research upstreaming

**Sketch / next.** Review the chart, bad-branch and island notes against
upstream's added fixed-set, coarse-defect and forest-tilt sections before
proposing a separate analytical contribution. Preserve the distinction
between proved conditional finite lemmas, numerical diagnostics and the
open uniform-exterior / all-scale estimates. The finite Wick/flow and
spectral-interface work can form a further independent contribution once
its intended scope and reproducibility package are agreed. No new
threshold, (I1), D8 or continuum theorem is claimed by the first PR.
