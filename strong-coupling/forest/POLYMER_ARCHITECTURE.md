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

## 5. Literature lead for the fill lemma: quantitative null-homotopy (sketch, not checked in detail)

Chambers–Manin–Weinberger, *Quantitative nullhomotopy and rational homotopy type*, GAFA 2018 (arXiv:1611.03513):
a null-homotopic L-Lipschitz map S^m → S^n has a null-homotopy of thickness C(m,n)(L+1) and width C(m,n)(L+1)².
Companion paper: Chambers–Dotterrer–Manin–Weinberger, *Quantitative null-cobordism*, JAMS 2018.

How this maps onto the fill. The collar's gauge function is g: S³_R → SU(2) ≅ S³ with degree k = 0. Rescaled to the unit
sphere it has Lipschitz constant L ≈ κR, where κ is the per-link deviation.
- **L ≲ 1 (κR ≲ 1).** g stays near a constant up to O(L). A geodesic radial contraction then fills with link
  deviation O(κ), so C is independent of R. This is elementary.
- **L ≫ 1 (κR ≫ 1).** Put the CMW homotopy in the radial direction. Spatial gradient ~ L/R ~ κ, which is fine. But the
  width is ~ L², which forces a radial Lipschitz of order L²/R ~ κ²R. That is **not** uniform in R.
  Either a better (linear-width) bound is needed for maps S³ → S³, or the polymer counting has to pay for
  κR ≫ 1 collars separately. Gromov's linear conjecture is open in general. For S³ → S³, which is degree theory only,
  linear null-homotopies may be easier: **open, to check**.
- The numerical test (job 2023a16d…, 3D cube) cannot see the L² issue, because π₂(SU(2)) = 0 and its collars have L = O(1).

## 6. Correction to §5: the parameter is κR² (curvature), not the gauge Lipschitz constant (sketch / open)

**Scalability lead.** Berdnikov–Manin, *Scalable spaces* (Invent. Math. 2022, arXiv:1912.00590), and Manin, *Plato's cave*
(arXiv:1801.00335), give O(L)-Lipschitz null-homotopies X → Y for scalable Y. Confirmed from the abstract: the scalable class includes spheres and Lie groups. Precise theorem form (linear null-homotopy
constant, domain hypotheses) still to be read from the paper. If S³ = SU(2) is scalable, the L² width in §5 disappears for the *flat* part of a collar.

**Adversarial point: the gauge function is the wrong object.**
- **Exactly pure-gauge collar.** If U_xy = g_x g_y⁻¹ on the collar, then any extension of g to interior sites gives
  plaquettes exactly 1. Lipschitz bounds are irrelevant, and no winding question arises, because k counts curvature
  and the collar's curvature is zero.
- **Calm collar.** Plaquettes ≤ κ only. In a radial gauge from a pole of S³_R, the collar connection has |A| ≲ κ·dist ≤ κR.
  The cone fill A_t = tA has F = tF_A + (t² − t)A∧A, which gives **|F| ≲ κ + (κR)²**.
  - So a naive fill is calm iff κR ≲ 1, i.e. κR² ≲ R. This regime is not uniform.
  - Uhlenbeck-type gauge needs ∫_{S³_R}|F|^{3/2} ~ κ^{3/2}R³ small, which also fails for large R.
- **Consistency check (no contradiction).** For any fill, ∫_B tr F∧F ≡ CS(collar) mod ℤ. Both sides are O(κ²R⁴), and a
  calm fill has ∫|F|² ≤ κ²R⁴. So there is no dimensional obstruction to a calm k = 0 fill.
- **Revised open lemma (multiscale).** For a collar on S³_R with |F| ≤ κ and k = 0, there is a fill with |F| ≤ Cκ and C
  independent of R. A plausible proof is a dyadic decomposition: on each scale-ℓ sub-block, κℓ² ≪ 1 holds once ℓ ≪ κ^{-1/2},
  which allows local Uhlenbeck gauges. The transition functions then form a Lipschitz map into SU(2), which needs scalability to
  be null-homotopic with linear constants. **Open.**
- **Relevance check for (P1).** Calm regions in the polymer expansion have κ ~ g·polylog, so κR² ≫ 1 holds for large
  regions. The hard regime is generic, not an edge case.

**Numerics.** Job 2023a16d (κ ≤ 0.01, R ≤ 10, so κR² ≤ 1) only tests the easy regime. The follow-up needs κ = 0.1–0.3 at R = 10
(κR² = 10–30).

## 7. Calm-fill numerics, 3D SU(2) cooling (heuristic)

Jobs 2023a16d (κ = 0.001, 0.01) and d2257246 (κ = 0.1, 0.3), with R = 3…10. Collar links are a smooth gauge transform composed
with iid noise of size κ. The table gives interior max plaquette / boundary max plaquette:

| κ | R=3 | 4 | 6 | 8 | 10 |
|---|---|---|---|---|---|
| 0.001 | 0.57 | 0.53 | 0.43 | 0.44 | 0.36 |
| 0.01 | 0.57 | 0.53 | 0.43 | 0.44 | 0.36 |
| 0.1 | 0.53 | 0.38 | 0.40 | 0.47 | 0.36 |
| 0.3 | 0.57 | 0.42 | 0.44 | 0.53 | 0.39 |

No growth with R: C ≤ 1 in every run.

**Adversarial caveat.** With iid noise, the curvature integrated over the collar is a random walk of size ~κR, not
κR², so these collars are nearly the easy case. Coherent abelian fields are trivially fillable by the same field.
The real test is coherent *noncommuting* curvature. That is the hedgehog collar: each face carries field κ with
generator n̂·σ. It runs as job 0e0d65f7 (κ = 0.05 and 0.2, so κR² reaches 20).

**Hedgehog v1 (job 0e0d65f7): invalid as a test.** The generator jumps where faces meet, so the edge plaquettes are of order
κR. bd_max therefore grew linearly in R, from 0.08 to 0.36 at κ = 0.05, and the collar was not calm. Measured against the
intended κ, the interior grows like ≈ 0.3κR. That is ambiguous: it may only reflect the non-calm edges.

**The real threat, identified.** The test case is a reducible connection with direction field n̂ of degree N, covariantly constant, and
uniform curvature F = κ n̂·σ. Total flux is 4πN, with N ≈ κR²/4π. This collar is calm. Filling it in the ball means unwinding
charge N through non-abelian configurations, and **I know no bound on that cost**: the "monopole charge" is not
topological in SU(2) because π₁ = 0, but non-topological does not mean cheap.
- **3D check: job af45d8e1** (`fill_mono.py`). Valid only if flux_err ≈ 0. Its output includes the depth profile.
- **4D analogue.** On S³ the corresponding object is a calm collar carrying a large, non-quantised Chern–Simons
  density. Same open question.
