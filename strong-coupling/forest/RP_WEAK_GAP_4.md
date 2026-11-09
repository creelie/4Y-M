# Direction #4: weak mass gap from reflection positivity (one focused session)

Status: **outcome is a no-go for the correlation-level version, plus a reduction of the chessboard
version to the existing open problem (D).** No new lemma. #4 is closed; see the end for what remains.

Setting: SU(2) Wilson action on a periodic torus Λ = (Z/L_tZ) × (Z/L_sZ)³, reflection positive
(Osterwalder–Seiler) with positive self-adjoint transfer matrix T, 0 ≤ T ≤ ‖T‖ = 1 after normalization.
The target ("weak gap"): truncated correlations of local gauge-invariant observables decay as e^{−m|t|},
uniformly in L_s, in the periodic state, with no supremum over boundary conditions.

## 1. No-go: finitely many periodic-state correlators give only *upper* bounds on the gap

Let O be a gauge-invariant observable on a time slice with ⟨O⟩ = 0, and let a_n = ⟨O, T^n O⟩ (L_t → ∞).
By the spectral theorem a_n = ∫ λ^n dμ_O(λ), with μ_O ≥ 0 on [0, 1). Cauchy–Schwarz gives
a_n² ≤ a_{n−1} a_{n+1}, so

  m_eff(n) := log(a_n / a_{n+1}) is non-increasing in n, and m_eff(n) ↓ m_O := −log sup supp μ_O.

So any effective mass measured at finite separation is an **upper** bound on the true decay rate of O.
A criterion of the form "κ(r) < 1 at some finite r in the periodic state ⇒ gap" cannot follow from RP and
spectral theory alone. A finite stretch of a log-convex sequence is compatible with any m_O ≥ 0, including 0
(take μ_O with a small atom near λ = 1). Consistency check: 4D U(1) is RP and gapless at weak coupling,
yet at any fixed r its periodic-state correlators can be made to look gapped by taking β just above the transition.

Consequence for `BALL_CERTIFICATE.md` (weak-form variant) and for #3 (`strong-coupling/forest/experiments/linresp/`):
L* = 2 is a statement about a finite stretch of a log-convex sequence. It bounds the gap from above, not below.
The "chessboard estimate transporting κ(r) < 1 to decay" asked for there does not exist in that form:
chessboard estimates bound probabilities of **block events** (via tiled partition functions), not covariances.

## 2. What RP does give: rarity of bad blocks, which reduces to (D)

The one volume-uniform input RP provides is the chessboard estimate: for a block event A reflected over
the torus, ⟨Π_{b∈B} A_b⟩ ≤ z(A)^{|B|}, with z(A) a tiled-partition-function ratio. Combining this with a
coupling of two copies (whose marginals are both the periodic state, so the chessboard bounds apply in each
copy) gives a disagreement-path bound of the form

  P(disagreement path of k blocks) ≤ μ^k (z^{1/2} + r_good^{1/2})^k   (Hölder, p = q = 2),

where r_good is the contraction of the block coupling **with both boundaries good**. This is exactly
direction C2 (`rare/RARE.md`):
- (P), bad-block rarity, is already proved there by chessboard estimates, with explicit thresholds;
- r_good < 1 is (G) + (B), and RARE.md §3 shows it fails in TV at weak coupling (conditional laws are at
  TV distance O(1) even for good neighbours), so it needs a Wasserstein/chart-metric contraction — the open
  composition problem (D).

So the RP route replaces sup over all ω by sup over good ω, which is a genuine weakening of certificate (C),
but its remaining input is the same open problem (D), and it must be supplied at the block scale where the
running coupling is O(1) (the window). Sanity check: for 4D U(1) at weak coupling, r_good ≥ 1 at every
block scale (massless photon), so the reduction does not prove anything false.

## 3. Verdict and what remains

- Correlation-level RP criterion: **no-go** (§1). Drop.
- Chessboard + good-boundary coupling: **= C2 = (D)** (§2). Nothing new beyond restating (C) with good
  boundaries. Record this as an optional weakening of `BALL_CERTIFICATE.md` (C): sup over ω with all
  boundary plaquettes good, plus a chessboard rarity term.
- Next per `../CLAY_NEXT_DIRECTIONS.md`: #6 (audit of Tomboulis 2007), a literature task plus small numerics.
