# The two walls, (D) and (E): where they stand (November 2026)

A summary from a reader (Lockwood) of release v4.29.0:

> (D), the one-cube estimate: the bounds on a single block have to hold uniformly as the lattice spacing goes to zero.
> They are proved only at strong coupling (up to 0.345 for SU(2)), not along the weak-coupling path to the continuum.
> The large-field cases are only shown to be typical at physical scale, not ruled out. (E), existence: building the
> continuum theory as a limit and showing it satisfies the axioms with a gap in physical units. Everything else in
> v4.29.0 is proved. Closing (D) would reduce the problem to (E), and (E) is essentially the Clay problem itself.

This note checks the summary against the paper and splits (D) into the pieces the paper leaves open. It records one
observation made after the release (§3) and ranks the next targets (§4). Nothing here closes (D).

## 1. Corrections to the summary

- **Strong-coupling range.**
  - SU(2): Theorem `blocks` gives β_W ≤ 0.351, with α(0.345) ≤ 0.9805 and α(0.351) ≤ 0.99995. The value 0.345 is
    the earlier endpoint, from before the sharpening in `gauge-fixed/blocks/sharp/`.
  - SU(3): Theorem `su3blocks` gave β_W ≤ 0.76 in v4.29.0. It now gives β_W ≤ 0.77, with α₃(0.77) ≤ 0.99978
    (direction A1 of `DIRECTIONS.md`).
- **"Everything else is proved."** True of the paper's theorems, with two trust boundaries:
  - the numerical certificates are in Arb ball arithmetic, and only the final arithmetic is rechecked by the Lean
    kernel (`A4_LEAN.md`);
  - nothing has been refereed.
- **The large fields.** The summary is accurate. Large fields are shown to be rare and small; they are not shown to be
  harmless. In more detail:
  - the bad plaquettes form clusters of microscopic physical size (Theorem `defects`, `rare/RARE.md`);
  - they are present along the tuning (Proposition `manydefects`);
  - fibres without a calm completion are the typical case on cubes of fixed physical size (Proposition `calmrare`);
  - no estimate shows that the bad regions do not spoil the cube coefficient.
- **(E).** The summary is accurate. Section "The continuum step at one physical scale" proves (D) ⇒ a gap ≥
  log(1/α)/(2ℓ) for every continuum limit. The existence of a limit satisfying the axioms is `OPEN.md` §3, the Clay
  problem itself.

## 2. (D) split into its pieces

(D) is posed for block variables (`OPEN.md` §4, subsection "Block-constrained cubes"): the constrained cube condition
D_c(n, α), with a supremum over boundary conditions, uniformly along β_W ∼ b log(1/a).

| piece | content | status | where |
|---|---|---|---|
| D1 | Gaussian (quadratic) version at every scale | **proved**: gap 537/5000 (Proposition `blockgapsharp`); D_c(n, α) for n ≥ 10⁸²⁷/α at every scale k (Theorem `constrainedscales`) | `main.tex` |
| D2 | one composition step for the averages | **proved**: Corollary `compose`; Q_k ≤ S_k ≤ 256 Q_k (Theorem `effectiveupper`); locality (Theorem `effectivelocal`) | `main.tex` |
| D3 | decay of the constrained law on one fibre with a calm background | **proved**: rate α_G = 2.5·10⁻⁴, amplitude O(1/β) (Theorem `chartdecay`) | `main.tex` |
| D4 | mass off the chart on fibres with no calm completion | **open**; bounded only in mean over exteriors (Theorem `offchart`) | Remark `chartdecayopen`(a) |
| D5 | fibre decay ⇒ cube coefficient, with the supremum over boundary conditions | **open**. Calm case done in `D5.md` (C_G ≤ 3·10¹³, links in no average); the global form is volume-limited, and the local form needs (L) | Remark `chartdecayopen`(b), `D5.md` |
| D6 | large fields: rare | **proved** on the torus, with explicit constants | Theorem `defects`, `rare/RARE.md` |
| D7 | large fields: harmless | **open** | `rare/RARE.md` §4 |
| D8 | the cube condition for the averages ⇒ the law one scale up; the coarsest averages; the crossover | **open** | `OPEN.md` §4 |

D4, D5, D7 and D8 are the open content. D5 and D7 are the two halves of Lockwood's sentence: "proved only at strong
coupling", and "large fields not ruled out".

## 3. An observation: islands are cheap compared with their rarity

The Bassalygo–Dobrushin scheme for random bad sites [BD86] needs three things:
- (G) contraction on good regions;
- (B) bad islands amplify influence by at most C^{|island|};
- (P) islands have tails beating C and the entropy.

`rare/RARE.md` gives (P). Item (B) is cheaper than it looked:
- The single-link Kantorovich coefficients are bounded **uniformly in the field**. With the sharp quarter (modulus 1/4), each
  Dobrushin entry is ≤ 4|β| · 1/4 = |β| = β_W/4, and each row ≤ 6(d−1)|β| = 4.5 β_W in d = 4 (`Thresholds.lean`; both halves of
  the quarter are in Lean).
- So an island of n sites, which has at most 4n links, amplifies W₁ influence by at most (4.5 β_W)^{4n}, whatever the
  field inside it.
- The rarity lemma at ρ = β_W^{−1/4} with ∗-adjacency (D = 80) gives cluster tails qⁿ with
  log q ≤ −β_W^{1/2}/16 + 6 log β_W + 24.7 + log(6e·80).
- Hence q·(4.5 β_W)⁴ < ½, and the expected island weight per site is summable, once **β_W ≥ 1.1·10⁷**. This is the
  same scale as the rarity thresholds of `RARE.md`.

What this does and does not give:
- It does **not** close D7. It reduces "harmless" to (G), in a metric in which the trivial single-link bound for
  islands can be inserted, under three caveats:
  1. **The bad set depends on the spins.** In [BD86] the potential is quenched. Here the chessboard domination
     controls the expected weight of paths through islands (finite intersections), which is what an *annealed*
     comparison argument needs. An annealed Dobrushin comparison theorem with configuration-dependent coefficients is
     not standard and is not proved here.
  2. **Torus only.** The rarity lemma uses reflection positivity. The supremum over exteriors in (D) is not covered
     (`RARE.md` §4).
  3. **The metric of (G).** At weak coupling the single-link coefficients do not contract even in calm regions
     (`OPEN.md` §4: bounds on single lines cannot give (D)). The contraction of D3 is in the chart coordinates of one
     fibre. To use the island bound there, it has to be transported from links to chart coordinates. The section S of
     Theorem `chartdecay` has row sums ≤ 14 and range one, which suggests a loss of a constant factor per site, but
     this is not checked.
- If the three caveats are met, D7 follows from D3 + D5. **The large fields are then not a separate wall: D5 is the
  wall.**

## 4. Next targets, ranked by what they would decide

1. **D5 on calm exteriors.** Prove the constrained cube coefficient for boundary conditions that admit a calm
   completion, from Theorem `chartdecay`(ii). The covariance bound there is of Helffer–Sjöstrand type with a weight
   e^{αφ}. A W₁ coefficient for the cube needs the response of the cube law to a boundary change, which is a
   covariance with the boundary interaction. This is the most concrete open step: a finite argument with the
   paper's constants, and no new idea is needed beyond controlling the boundary terms of the chart.
2. **Annealed comparison (caveat 1).** A Dobrushin comparison theorem whose coefficients are random and dominated in
   finite intersections by a Bernoulli field. This is a self-contained probability lemma, and with §3 it would turn
   D7 into D5.
3. **Local rarity (caveat 2).** Extend the rarity lemma from the torus to boxes under an arbitrary exterior. The
   energy–entropy count of Proposition `anyexterior` is local up to the boundary layer. The obstacle is that scattered
   bad sets pay a boundary layer of size comparable to their own volume (`RARE.md` §4).
4. **D4 and D8** need the multiscale analysis. They are Bałaban's territory, and no short route is visible.

**(E)** is not attacked. With (D) it is reduced to the existence of a limit with the axioms (`OPEN.md` §3), which is
the Clay problem.

## Reference

[BD86] L. A. Bassalygo, R. L. Dobrushin, *Uniqueness of a Gibbs field with a random potential — an elementary
approach*, Theory Probab. Appl. 31 (1986) 572–589.
