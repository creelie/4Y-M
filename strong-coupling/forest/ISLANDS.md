# Items 3 and 4: local rarity and the island bound in chart coordinates (November 2026)

This is items 3 and 4 of `ROADMAP.md` §7. It has two results and a ledger of thresholds.

1. **Proposition L (§1).** This is the local off-chart bound (L) of `D5.md` §4, in window form. It follows from the
   proof of Proposition `calmfibre`, once the free links of a window are enlarged by the private links of the pairs
   that straddle its boundary. The constants depend on the window and **not on the volume**. Status: the steps of
   `calmfibre` have been checked to transfer, but the proof is not written in the paper's detail.
2. **Island transport (§2).** The single-link island bound of `D_STATUS.md` §3 passes to chart coordinates at the
   cost of a collar of one block around each island and a factor C_G (from `D5.md` §2) at each exit. Status: crude
   bound, proved given D5-calm.
3. **Thresholds (§3).** With 1 and 2, the calm windows and the islands balance at β_W ≈ 10¹⁶⁴, independently of
   the volume. This is a ledger, not a theorem. It still needs the annealed step (I1) of `ANNEALED.md`, which is not
   available.

## 1. Proposition L: the local calm-fibre bound

**Setting.** W is a set of blocks (the cube), with constrained law μ_{η,v}. W₁ ⊂ W is a window of m blocks. The
**straddling pairs** are the pairs k = (y, μ) with y ∉ W₁ and y + e_μ ∈ W₁. Their averages involve the links w_a
inside B(y + e_μ) ⊂ W₁, by the formulas P_{x′_a} = ℓ_a w_a that precede Lemma `privatesubmersion`.

Let L₁ be A_{W₁} together with the eight private links ℓ₀, …, ℓ₇ of every straddling pair. There are at most 4m
straddling pairs, so |L₁| ≤ 49m + 32m.

**Claim.** Condition μ_{η,v} on every active link outside L₁. Then:
1. The conditional law is the constrained law on the set of fields on L₁ for which:
   - the averages of the pairs of W₁ are fixed;
   - the averages of the straddling pairs are fixed;
   - everything else is fixed.
2. Each of these pairs has its path link d_k = ℓ₀ in L₁, and that link is private to it.
3. Suppose this conditional fibre contains a field with δ(U_p, I) ≤ K/(4√n_P) on every plaquette p with a vertex in a
   block of W₁⁺, where n_P ≤ 576m counts the plaquettes that contain a link of L₁. Then the conditional law gives
   mass at most

     ε_loc(m, K, λ) = (3·10⁸ n_P^{3/2}/K³)^{73m} e^{−3λK²/4}

   to the fields with such a plaquette above K.

**Proof (transfer of `calmfibre`).**

(1) and (2). A private link enters only the average of its pair (paragraph after Lemma `privatesubmersion`). The
pairs whose averages involve L₁ are those of W₁ and the straddling pairs. The private links of distinct pairs are
distinct.

(3) The proof of `calmfibre` uses four facts, and each holds here with the pair set K₁ = (pairs of W₁) ∪ (straddling
pairs):
- **(a) Lemma `smallframes`(i) applies at every pair of K₁.** Its hypothesis concerns plaquettes in B(y) ∪ B(y + e_μ),
  and for a straddling pair y ∈ W₁⁺.
- **(b) The contraction Ψ solves each constraint for its path link.** Its constant (7/8)Λ, and the radius bounds
  ρ = K/(188√n_P), R = 44ρ, are local to the pair.
- **(c) The volume bounds of `eq:fibrevolume` and h(B_ρ).** These go link by link. The exponent is the number of free
  links that are not path links: |L₁| − |K₁| ≤ 45m + 28m = 73m. The Jacobian factor 1.957^{−3|K₁|} has
  |K₁| ≤ 8m, and is absorbed exactly as (M⁴1.957¹²)^{1/45} < 2.93 is in the paper.
- **(d) Every plaquette has at most one path link, and n_P counts the plaquettes through free links.** Each of the
  32m straddling private links lies in six plaquettes, so n_P ≤ 384m + 192m.

**What has not been checked.** The integration step (`eq:fibrevolume`) is stated in the paper for A_W. It holds for a
general free set with private path links by the same coarea computation, but this has not been rewritten line by
line. ∎ (sketch)

**What it changes.** The constant depends on m and K only. Run inside the window, the argument of `D5.md` §2 has the
error 320λ|V₁|ε_loc instead of 320λ|V|ε_W. So the volume obstruction of `D5.md` §3 disappears **for cubes of fixed
lattice size**. It does not disappear for cubes of fixed physical size, which is what (D) asks for at the coarsest
scale; that is D8. The hypothesis of (3) is still a property of the conditioning configuration. It holds for most
window positions only at very large β_W (§3).

## 2. The island bound in chart coordinates

Caveat 3 of `D_STATUS.md` §3: islands are bounded link by link, while the contraction of D3 lives in the chart
coordinates of a calm region. The interface is handled as follows.

- **Collar.** The hypothesis (H) of `D5.md` §2 asks that every plaquette with a vertex in a block of W⁺ be small. A
  calm window that touches an island fails it. So treat the ∗-neighbours of an island's blocks (a collar of one
  block) as part of the island. An island of n blocks then becomes at most 81n blocks.
- **Inside the island and the collar.** Use the single-link coefficients, uniformly bounded in the field: each row
  ≤ 4.5β_W (`Thresholds.lean`). Counting links per site as in `D_STATUS.md` §3 (at most 4 per site, one site per
  block at the block scale), the island and its collar contain at most 4·81n = 324n links, each with row ≤ 4.5β_W.
  Crude count: (4.5β_W)^{324n}.
- **Exits.** A step from the island into a calm window is the response of the window law to one of its exterior links,
  which is at most C_G ≤ 3·10¹³ by `D5.md` §2 (case (a); case (b) needs its uncomputed constant). Bound the number of
  exit steps by the island-and-collar sites.

So an island of n blocks amplifies influence by at most (4.5β_W · C_G)^{324n}. This is a crude bound with the
collar included.

**Status.** Given D5-calm, this is proved as a bound on the coefficient matrix of a mixed block/link Dobrushin scheme:
- calm windows are rows bounded through §2 of `D5.md`;
- island links are rows bounded by the single-link coefficients.

It is caveat 3, answered with a constant per site. Two weaknesses:
- the constant is huge;
- case (b) of `D5.md` §2 enters at the straddling exits.

## 3. Ledger of thresholds (`d5/islands.py`, `d5/islands.out`)

These are order-of-magnitude numbers for SU(2), with K = 4.8·10⁻¹⁸ and C_G = 3·10¹³.

| condition | window side n₀ (α = 0.9) | β_W needed |
|---|---|---|
| D_c(n₀, α) from C_G (`D5.md` §2) | 5.7·10¹⁴ | — |
| off-chart error 320λ\|V₁\|ε_loc small | 5.7·10¹⁴ | ≥ 3.8·10⁹⁷ |
| windows calm in the sense of `calmexterior` (ρ_c = K/(8n_P) = 1.6·10⁻⁷⁹), with islands, collars and exits beaten by their rarity: q(ρ_c)(4.5β_W C_G)^{324} < ½ | 5.7·10¹⁴ | ≥ 8.3·10¹⁶³ |

All three thresholds are **independent of the volume**. So along the tuning β_W ∼ b log(1/a), they are passed for
small a. What they would give, if (I1) of `ANNEALED.md` were available, is the cube condition D_c(n₀, α) at the
lattice scale, with n₀ ≈ 6·10¹⁴. That is a gap in lattice units, which tends to zero in physical units. The physical
gap needs the multiscale composition D8. The chessboard rarity used in the last row needs the torus. Under an
arbitrary exterior, only the density bound of Corollary `localcalm` is available, and it does not give the
finite-intersection domination that the path expansion uses.

## 4. Summary for (D)

| piece | status after items 1–4 |
|---|---|
| (G) calm contraction | `D5.md` §2: C_G ≤ 3·10¹³, proved for links in no average; sketch otherwise |
| (L) local off-chart | Proposition L: window form, volume-independent; sketch-level transfer of `calmfibre` |
| (B) island cost in chart coordinates | §2: (4.5β_W C_G)^{324} per block, with a collar |
| (P) rarity | torus: proved (`rare/RARE.md`); arbitrary exterior: density only (`localcalm`) |
| (I1) annealed / conditioned laws | open (`ANNEALED.md`); Bałaban-type large-field expansion |
| D8 multiscale | open |

The walls that remain for (D) are therefore:
- (I1), together with rarity under an arbitrary exterior in finite-intersection form;
- D8.

Everything else on the list has at least a crude, volume-independent bound.
