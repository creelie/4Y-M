# One-scale large-field polymer expansion: what it must prove (Nov 2026)

Labels: **proved**, **sketch**, **heuristic**, **open**. Nothing here closes D.

## 0. Why the architecture changes
Every attempt to get (R1) under a fixed exterior failed: localcalm, chessboard, Peierls (`CALM_LINK.md` §4),
Dobrushin–Pechersky, and BKMP (`FRONTIER.md`). The calm-collar variant fails as well, because rarity has to hold given the
coupling's past, and that is the hidden-bit trap. In Bałaban/Dimock, large-field regions are **summed** on the torus, never conditioned
on. So the target is a polymer representation, not a cube condition.
Reference: Dimock arXiv:1212.5562 (φ⁴₃ UV). It covers the scalar case only, where the large-field cost comes from λφ⁴. In
the gauge case the cost is the action itself, e^{−βK²/C} per large plaquette.

## 1. Statement needed (open)
On the torus Λ at one block scale, write Z = Σ_Ω Z(Ω). Here Ω is the closure of the set of non-calm blocks plus a collar of width w,
and Z(Ω) integrates over configurations that are non-calm exactly on Ω and calm on Λ∖Ω.
(P1) **Peierls by summation:** Z(Ω)/Z_calm(Ω) ≤ e^{−c β K² |Ω|/w^d}. Here Z_calm(Ω) is the same integral with Ω also calm,
     and the bound is uniform in the configuration on Λ∖Ω **only through its calm collar**.
(P2) **Analytic calm region:** log Z(calm on X | collar data) has a convergent polymer expansion with
     exponentially decaying activities (the chart / `chartdecay` machinery, already in the paper at sketch level).
(P1)+(P2) together give a standard polymer gas, and a cluster expansion then gives decay of correlations on the torus. That
replaces the cube coefficient and D4/D7 at one scale. D8 (composition across scales) stays open.

## 2. The calm-completion step inside (P1) (heuristic)
(P1) needs a calm configuration on Ω compatible with the calm collar, with action ≤ C|Ω|·(small).
- Holonomy: the collar is calm, but loops around Ω can have O(1) holonomy. Spread over area R², that is plaquettes of size ~θ/R²,
  so it stays calm for large Ω. For Ω of a few blocks, θ ≤ π gives plaquettes of size ~π/b², which is calm only if K ≳ π/b².
  **Potential obstruction for small polymers**, where the threshold K ≈ 4.8·10⁻¹⁸ is far below π/b².
- Topology: collar boundary data can carry instanton winding. A fill of size R has plaquettes ~1/R², so the fill is again
  impossible to make calm for small Ω.
- In Bałaban, small polymers with non-trivial collar holonomy simply *are* large-field, because Ω is enlarged until the
  collar is flat-ish. So the polymer definition has to be "closure under fillability", and its entropy has to be counted.
  **This counting is the concrete open lemma.** It has to beat e^{−cβK²} per block. Fillability-enlargement can multiply
  |Ω| only by a constant if holonomy ≤ C·(collar plaquettes)·(perimeter), which is a non-abelian Stokes bound that holds on
  calm collars (sketch).

## 3. Next step
Prove the non-abelian Stokes/fillability lemma at one scale. Statement: a w-collar with all plaquettes ≤ κ around a region
Ω of diameter R has a calm fill with plaquettes ≤ C κ (R/w)^{?}, in the absence of instanton winding. Then check whether
the exponent in (R/w) is compatible with the threshold K. This is testable numerically on small SU(2) lattices
(minimise the action of the fill given random calm collars).

## 4. Adversarial correction to §2 (heuristic)
- **Holonomy obstruction: spurious for ball-like Ω in d=4.** The collar ≅ S³×[0,w] is simply connected. Every loop in
  it bounds a disk inside the collar, so non-abelian Stokes gives holonomy ≤ C·κ·(area in the collar). A loop does not need
  to pass through Ω to be bounded. The calm data on the collar is therefore close to a pure gauge g: S³→SU(2).
- **Winding is the real obstruction.** g has a degree k (π₃(SU(2)) = ℤ). By Lüscher admissibility, if all plaquettes are
  small, the lattice charge is an integer, so a fill must carry charge k. Its action is ≥ 8π²|k|/g², i.e. ≳ β_W·const per
  unit of winding, and the fill can be calm only if its size is ≳ K^{−1/2}.
- So non-trivial winding is suppressed by the action itself, e^{−cβ_W}: it is just one more large-field cost. Non-ball Ω
  (handles, i.e. non-simply-connected collars) bring the holonomy problem back on non-contractible loops. That needs the
  enlargement rule "fill the holes of Ω", and its entropy is a constant factor per block.
- Revised open lemma: **calm fill for k = 0 with plaquettes ≤ Cκ, C independent of R.** The danger is C growing with R
  (elliptic estimates for the fill are fine for Laplace-type problems, but the fill is non-linear). This can be tested
  numerically.
