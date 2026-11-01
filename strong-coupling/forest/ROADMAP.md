# Roadmap after the A/B/C directions (November 2026)

This file records:
- where every direction of `DIRECTIONS.md` ended up;
- how the directions relate to the two remaining walls, (D) the cube estimate and (E) existence;
- which path we bet on, and why.

The detail behind each line is in the linked files. The split of (D) into pieces D1–D8 is in `D_STATUS.md`.

## 1. Track A: strong coupling, fixed β

| item | outcome | where |
|---|---|---|
| A1 SU(3) blocks | **done.** β_W ≤ 0.77 (was 0.76); α₃(0.77) ≤ 0.99978 and m ≥ 1.1·10⁻⁴ there. Certified cover at r = 77/50 with 2825 boxes; constants on 826 boxes | `../gauge-fixed/blocks/su3/README.md` |
| A2 SU(2) 0.351 → 0.355 | **diagnosed, not proved.** The true block response is (1+β)/4 in every outside field. The slack is the 1/(1−C) factor (4.1%). A second-order expansion would reach about 0.355, but it is not written | `../gauge-fixed/blocks/feedback/README.md` |
| A3 larger blocks | **done, negative.** Built-up blocks reach a ceiling of about 0.369, against 0.366: not worth proving | `../gauge-fixed/blocks/shapes/README.md` |
| A4 Lean port onto LGT | **scoped.** The critical path is a W₁ (Kantorovich) Dobrushin theorem, which the Lean ecosystem lacks. LGT itself proves only β_W < 0.0104 (total variation) | `A4_LEAN.md` |

Current proved ranges:
- SU(2): β_W ≤ 0.351.
- SU(3): β_W ≤ 0.77.

**Relation to (D): none.** Track A is about fixed coupling, while (D) is posed along β_W → ∞. Pushing the range by
0.01 leaves the gap to the scaling region near β_W ≈ 2.2 essentially unchanged. A4 raises trust in proved results; it
does not add content.

## 2. Track B: the flowed-slice route

| item | outcome | where |
|---|---|---|
| B1 3D next order | **done** | `slice/THREE_D.md` |
| B2 lattice S′ at O(g²) | **done in the free vacuum**, uniformly in a. The vertex bound is proved (Lemma 3). The O(g²) vacuum correction and the constant c are open | `slice/sprime/LATTICE_SPRIME.md` |
| B3 stochastic quantization (CCHS) | **done** (literature). Spacetime-flowed slices are smoothed observables; the 3D measure the route needs does not exist | `slice/CCHS.md` |
| B4 3D unflowed Monte Carlo | not pursued (low priority) | |

**Relation to (D): a parallel route, not a prerequisite.** If it worked, it would go around the cube estimate.
Everything it has produced so far is perturbative. Turning it into a statement about the measure needs the same
nonperturbative control of large fields and of all scales that D7 and D8 ask for. It also needs a 3D measure that
nobody has constructed (B3).

## 3. Track C: the cube estimate

| item | outcome | where |
|---|---|---|
| C1 sharper cube constants | not pursued. Sharper local constants change nothing without the composition step | |
| C2 rare bad fibres | **half done.** "Rare" is proved on the torus by chessboard estimates, with explicit thresholds. "Harmless" cannot come from disagreement percolation at weak coupling | `rare/RARE.md` |
| C3 tree bounds, continuum | only reachable through the 3D route; not started | |

**Relation to (D): part of it.**
- C2's "harmless" half is D7.
- C3 is D8 and beyond.

## 4. (D), the cube estimate

The estimate: D_c(n, α) for block variables, with a supremum over boundary conditions, uniformly along
β_W ∼ b log(1/a). Pieces, from `D_STATUS.md`:

- **Proved:**
  - D1: the Gaussian version at every scale (n ≥ 10⁸²⁷/α);
  - D2: one composition step;
  - D3: decay of the constrained law on a calm fibre, at rate 2.5·10⁻⁴;
  - D6: large fields are rare.
- **Open:**
  - D4: the mass off the chart on fibres with no calm completion;
  - D5: fibre decay ⇒ cube coefficient, with the supremum over boundary conditions;
  - D7: large fields are harmless;
  - D8: the multiscale composition, the coarsest averages and the crossover.

New since v4.29.0 (`D_STATUS.md` §3) concerns D7. The single-link W₁ coefficients are bounded uniformly in the field,
by β_W per entry and 4.5 β_W per row (Lean). So a bad island of n sites amplifies influence by at most
(4.5 β_W)^{4n}, and its rarity beats this once β_W ≥ 1.1·10⁷. If three caveats are met, D7 reduces to D5:
1. an annealed comparison theorem;
2. rarity under an arbitrary exterior;
3. the transport of the island bound to chart coordinates.

## 5. (E), existence

(E) asks for a continuum limit that satisfies the Osterwalder–Schrader axioms and has a gap in physical units
(`OPEN.md` §3–4). The paper proves (D) ⇒ a gap ≥ log(1/α)/(2ℓ) for every continuum limit. What is left after (D) is
the existence of a non-trivial limit with the axioms, which is the Clay problem itself. It is not attacked; it waits
on (D).

## 6. The bet

**Do (D) itself, starting with D5. Keep B as a cheap diagnostic, not a route.**

- **Why not B.** B postpones the wall rather than avoiding it. Its upgrade from perturbation theory to the measure is
  D7 + D8 under another name, and it also needs a 3D measure that does not exist. B is still useful as a fast source
  of yes/no answers on which estimates should be true. It should run when a numerical check would settle a question
  raised by D.
- **Why D5.** It is the next step from a proved theorem (`chartdecay`). It is a finite argument with known constants,
  and the only new work is controlling the chart's boundary terms. Each proved piece of D narrows the problem
  permanently. Through §4, D5 may bring D7 with it.
- **The real risk is D8.** Composition across all scales is Bałaban's territory: that program ran for years without
  reaching a gap, and the proved Gaussian constant (10⁸²⁷) shows the distance. The realistic goal is to make D5 and D7
  theorems and isolate D8 cleanly. Closing (D) is not the expected outcome.

## 7. Order of work

1. **D5 on calm exteriors:** the cube coefficient from Theorem `chartdecay`(ii), with the boundary terms of the chart.
   **Done for links in no average (`D5.md`):** C_G ≤ 3·10¹³ plus an off-chart error. Global calmness is
   volume-limited, so D5 must be local, and then it merges with D7 through the missing local off-chart bound (L).
2. **Annealed Dobrushin comparison:** random coefficients dominated in finite intersections by a Bernoulli field
   (caveat 1). A self-contained probability lemma; with §4 it turns D7 into D5.
   **Done, negative (`ANNEALED.md`):** one-point rarity gives only an O(q) floor (approximate uniqueness). A bootstrap
   lemma reduces the rest to quenched bounds (I1) for laws conditioned on the large-field set. Those are Bałaban
   territory, not a probability lemma.
3. **Local rarity under an arbitrary exterior** (caveat 2), in the constrained-law form (L) of `D5.md` §4.
   **Partly done (`ISLANDS.md` §1):** Proposition L (sketch) gives the off-chart bound in window form with
   volume-independent constants, by adding the private links of the straddling pairs to the free set. Rarity of calm
   windows under an arbitrary exterior in finite-intersection form is still open (only the density bound `localcalm`).
4. **Island bound in chart coordinates** (caveat 3). The section S has row sums ≤ 14 and range one.
   **Done, crude (`ISLANDS.md` §2):** with a one-block collar, an island of n blocks costs at most
   (4.5β_W C_G)^{324n}. With item 3 the thresholds balance at β_W ≈ 10¹⁶⁴, independently of the volume, but only
   given (I1), and only for a gap in lattice units.
5. **D4 and D8:** multiscale analysis, no short route visible.
   **Scoped (`D4_D8.md`):**
   - D4 in window form is the large-field part of (I1).
   - The base case of D8 lands at an effective coupling of about 10–700 (heuristic, one loop). That is far below every
     weak-coupling threshold here (≥ 10³⁸), so the crossover is unavoidable.
6. **Side items, if cheap:** the A2 second-order proof (0.351 → 0.355), and the Lean W₁ Dobrushin theorem (A4).
   Neither bears on (D).
   **Outcome:**
   - A2 is parked. It needs a new Arb certificate (`../gauge-fixed/blocks/feedback/README.md`).
   - The arithmetic half of the Lean W₁ Dobrushin theorem is proved
     (`YangMills/StrongCoupling/DobrushinComparison.lean`, no `sorry`). The coupling half is still open.

### After §7: the calm-data link (`CALM_LINK.md`)

- The work here is merged with the OpenGluon calm-data note. The constant for the first kind of link is C_* ≤ 7.4·10¹². The second kind
  has C_b ≤ 9.2·10¹³ plus a shell term.
- The shell term, averaged over v, reduces to the localised off-chart covariance (sketch, §2).
- A window bootstrap gives that covariance with decay rate α_G/2, down to a floor, at β_W ≳ 10⁶⁴ (§3).
- What is left of the calm-data cube coefficient is one statement, (R1): one-point calm rarity under an arbitrary
  exterior. All three routes to it fail (§4), and it is the large-field half of a renormalisation-group expansion.
  **D is blocked there**, and at D8. Work moves to other paths.
