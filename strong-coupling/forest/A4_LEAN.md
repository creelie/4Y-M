# Direction A4: porting the block lemmas onto mrdouglasny/LGT (scoping)

`DIRECTIONS.md` (A4) proposes to port Lemmas `tilt`, `onesite`, `blockK` and Proposition `blockstate` of `main.tex` into
Lean on top of [mrdouglasny/LGT](https://github.com/mrdouglasny/LGT). The goal is a formally verified SU(2) lattice gap
up to β_W ≈ 0.35.

This note records what LGT provides, what our repository already has in Lean, and what is missing. The LGT checkout
read for this note is the public default branch, toolchain v4.30.0.

**Conclusion.**
- The port is feasible, but it is not a port. LGT's Dobrushin layer is built on **total variation**, while our
  block route uses **Kantorovich (W₁) coefficients**. The W₁ Dobrushin theorem is the critical-path item, and
  neither LGT nor its dependency markov-semigroups has it. Their own plan lists it as "Phase A: Wasserstein (LATER)".
- Everything specific to SU(2) blocks is either already in Lean in this repository or is elementary measure theory:
  - the finite combinatorics;
  - the arithmetic of the constants;
  - the algebraic halves of the sharp quarter.
- The exception is the analytic enclosures of g = log Z₀ (Lemma `vmfconst`), which are certified in Arb, not in
  Lean.

## 1. What LGT proves

The main statements are `ym_mass_gap_exponential_decay` and `ym_mass_gap_rate_exists` in
`LGT/MassGap/StrongCoupling.lean`. Both have zero sorries and no project axioms, according to its README.
- **Setting:** U(n) Wilson action S_p = β(n − Re Tr U_p) on the torus (ℤ/Nℤ)^d.
- **Condition:** β < 1/(32n(d−1)).
- **Method:** single-link total-variation Dobrushin. The influence is 1 − e^{−4nβ}, with 8(d−1) neighbours.

For SU(2) (n = 2) in d = 4, the convention gives β_LGT = β_W/2. Their range is therefore **β_W < 1/96 ≈ 0.0104**.

| result | β_W range (SU(2), d = 4) | metric | status |
|---|---|---|---|
| LGT | < 0.0104 | TV, single link | Lean, complete |
| `Thresholds.lean` (this repo) | < 2/9 = 0.222 | W₁, single link, sharp quarter | Lean arithmetic only |
| Theorem `blocks` (`main.tex`) | ≤ 0.351 | W₁, plaquette blocks | paper + Arb + Lean arithmetic |

A port would raise the formal range by a factor of 34.

## 2. What this repository already has in Lean

All of the following files have zero sorries.

| file | content | role in the block proof |
|---|---|---|
| `YangMills/StrongCoupling/QxCoefficients.lean` | transverse half of the sharp quarter (coefficient positivity) | Theorem `quarter`, used in `onesite`(i) |
| `YangMills/StrongCoupling/ZonalBound.lean` | zonal half of the sharp quarter (Riccati subsolution) | same |
| `YangMills/StrongCoupling/GaugeFix.lean`, `BlockGauge.lean` | finite gauge-fixing patterns and counts | the forest and the block geometry (Lemma `plaqblocks`) |
| `YangMills/StrongCoupling/Thresholds.lean` | single-link W₁ Dobrushin arithmetic | template |
| `strong-coupling/gauge-fixed/blocks/sharp/lean/BlockCheck.lean` | kernel-checked chain w → s → C → K → α ≤ 0.99997 from rational enclosures | eq. (`blockconst`) and α(0.351) |

## 3. The gaps, in dependency order

| # | item | size | where it lives |
|---|---|---|---|
| G1 | **W₁ Dobrushin theorem** for a specification on a compact metric spin space with block updates: uniqueness and covariance decay from Σ_B C_{Bj} ≤ α < 1 (Lemma `blockcomparison`, Proposition `blockstate`) | large, roughly 1–2k lines: Kantorovich duality on compact metric spaces, coupling of conditional laws, Neumann series | new; markov-semigroups' "Phase A" is the natural home |
| G2 | **Forest gauge-fixing identity** (Lemma `gaugefix`) for the 4D forest | medium | LGT's `GaugeFixing.lean` proves the axial (tree) version; the forest version is the same argument segment by segment |
| G3 | **Lemma `tilt`** (covariance under a bounded tilt, Cauchy–Schwarz) | small, about 100 lines | pure Mathlib measure theory |
| G4 | **Poincaré on S³ for vMF laws** with constant 1/3 (Lemma `q-lem:vmf`) | medium: needs Bakry–Émery or a direct spectral bound on S³ | markov-semigroups has `Diffusion/BakryEmery.lean` (abstract), but no sphere instance |
| G5 | **Sharp quarter as a measure statement.** The analytic glue: covariance = derivative of ν_{c+th}(f) and its reduction to the zonal and transverse coefficients | medium | the algebraic halves exist (§2) |
| G6 | **Lemma `onesite`** (i) and (ii), and **Lemma `blockK`** (the two-site coupling argument) | medium | elementary given G3–G5 and a W₁ coupling API from G1 |
| G7 | **Analytic enclosures of g = log Z₀** (Lemma `vmfconst`): Hessian eigenvalues, ‖D³g‖ and \|∇g\| on balls, from the Bessel series with a tail bound | large if done in Lean (interval arithmetic for series); small if stated as checked hypotheses | Arb now; LGT has no interval-arithmetic layer |

The critical path is G1. Without it the TV framework applies, and TV coefficients cannot reach the block range: the
single-link TV influence at β_W = 0.35 (β_LGT = 0.175) is 1 − e^{−4nβ} = 1 − e^{−1.4} ≈ 0.75 per neighbour, against
the W₁ coefficient β_W/4 ≈ 0.09.

## 4. Suggested staging

1. **Statement-level port.** Restate Theorem `blocks` in LGT's types, taking G1 and G7 as explicit hypotheses
   (G7 as numeric hypotheses already discharged by `BlockCheck.lean`). Prove G2, G3, G5 and G6 against them. This
   isolates the two big items and gives a reviewer the full logical chain.
2. **G1** as a contribution to markov-semigroups' Phase A, generic in the spin space and the block structure.
3. **G4** as a sphere instance of the Bakry–Émery machinery: S³ has Ricci curvature 2 and the vMF potential is
   λ-convex on the relevant ball.
4. **G7** last, or never: the Arb certificate with outward rounding is a reasonable trust boundary, and `BlockCheck.lean`
   already removes the arithmetic from it.

## 5. Practical notes

- Toolchains: LGT pins Lean v4.30.0 and Mathlib v4.30.0 (its roadmap says v4.29.0, while its lakefile says v4.30.0).
  This repository pins v4.33.0, and big-bang-monotile's `SliceRatio` pins v4.19.0. A port should live in a fork of LGT
  at its toolchain. Our existing files use only basic Mathlib and core, and should port with renaming only.
- Nothing was built for this note. Building LGT needs the Mathlib cache for v4.30.0, about 5 GB, which this session
  did not download.
- The discrepancy noted in `DIRECTIONS.md` (α′(0.351) ≤ 0.99997 against ≤ 0.99995) should be reconciled before the
  port, since the formal statement will quote one of them. `BlockCheck.lean` certifies ≤ 0.99997; `main.tex` quotes
  ≤ 0.99995 from the direct Arb evaluation. Both are below 1, so either suffices for the statement at β_W = 0.351.

**Status of A4:** scoped; not started in Lean. The critical-path item is a W₁ (Kantorovich) Dobrushin theorem, which
the Lean ecosystem does not yet have.
