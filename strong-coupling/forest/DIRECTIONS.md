## Status (November 2026)

| item | status | where |
|---|---|---|
| A1 SU(3) blocks 0.77 | interval cover at r = 77/50 running (level 3) | `../gauge-fixed/checks/su3/deg2/` |
| A2 SU(2) true feedback | diagnosed, not proved. The true response is (1+β)/4 in all fields. The slack is the 1/(1−C) factor (4.1%), not a doubled C, so the premise below is corrected | `../gauge-fixed/blocks/feedback/README.md` |
| A3 block shapes | done (negative: built-up blocks reach a ceiling of about 0.369, against 0.366) | `../gauge-fixed/blocks/shapes/README.md` |
| A4 Lean bridge | scoped; the critical path is a W₁ Dobrushin theorem, which LGT (TV only, β_W < 0.0104) lacks | `A4_LEAN.md` |
| B1 3D next order | done | `slice/THREE_D.md` |
| B2 lattice S′ at O(g²) | done in the free vacuum, uniformly in a. The vertex bound is proved (Lemma 3); the O(g²) vacuum correction is open | `slice/sprime/LATTICE_SPRIME.md` |
| B3 CCHS | done (literature): spacetime-flowed slices are smoothened observables; the 3D measure is missing | `slice/CCHS.md` |
| B4 | not pursued (low priority, as below) | |
| C1 cube constants | not pursued. Without a structural idea for the composition step, sharper local constants do not change what is proved | |
| C2 rare bad fibres | rarity proved on the torus (chessboard); harmlessness is not reachable by disagreement percolation | `rare/RARE.md` |
| C3 | only via 3D, as below | |
| (D) overview | the open pieces of the cube estimate, and the observation that islands cost (4.5β_W)^{4n} against their rarity | `D_STATUS.md` |

## Track A: certain gains at strong coupling (days each, fully checkable)

**A1. SU(3) blocks, 0.76 → about 0.77.**
- Theorem su3blocks (`main.tex:379`, constants in Sec. su3blocks, line 2152) has the same structure as the SU(2) lemma we sharpened.
- Port both changes: interpolating exponents, and the Poincaré inequality for the block law.
- Certify with the same arb, mpmath and Lean pipeline as `blocks/sharp/`.
- This is the most predictable win on the list.

**A2. SU(2), 0.351 → about 0.355.**
- Our C ≈ 0.038 is about twice the true feedback. *(Corrected in `../gauge-fixed/blocks/feedback/`: C itself is close to
  the true conditional feedback; the slack is that the true response does not grow by 1/(1−C) at all.)*
- Bound the actual y₁–y₃ correlation instead of the two-step estimate.
- This is idea R2 in `forest/IDEAS.md`.

**A3. Compute the ceiling before proving larger blocks.**
- Before writing any proof for two-plaquette or cube blocks, measure each shape's first-order sensitivity K by Monte Carlo, the way we found (1+β)/4 for plaquettes. That gives each shape's ceiling β\*.
- This takes CPU-minutes per shape. It tells us whether a bigger block can reach, say, 0.45, or only 0.37.
- Then we only write proofs for shapes worth the effort.

**A4. Formal bridge.**
- Port the block lemmas (Lemmas 13–14 and Proposition 17) into Lean on top of mrdouglasny/LGT.
- That gives a formally verified SU(2) gap up to about β ≈ 0.35, which is new as a *formal* result.
- It's slower, but it's the item an outside reviewer can check without trusting us.

The two copies of the SU(2) certificate disagree slightly:

| | α′(0.351) bound | Chain fails at |
|---|---|---|
| Earlier write-up (big-bang-monotile) | ≤ 0.99997 | 0.352 |
| What I pushed to 4Y-M | ≤ 0.99995 | 0.3511 |

The difference is probably in the series truncation or the number of radius cells. One of us should reconcile it before anything is quoted.

## Track B: the slice ratio (I2), the real frontier

What we now know:
- **Free 4D:** the ratio at r > 0 converges for the paper's magnetic slice field, flowed or not, even where the one-step energy diverges.
- **2D SU(N):** proved.
- **4D at order g²:** the one-step energy has a log of size 8/11·2b₀, but the r > 0 ratio only sees a component with small, ultraviolet-finite norm.
- **Proposition S′** (bound the soft part, control the hard part's norm) is the right target, and it is proved as an abstract operator inequality in Lean.

**B1. Do the next-order leakage computation in three dimensions. This is my top new idea.**
- In 3D, g² has dimensions of mass, so the leakage should be suppressed by a power (about g²a) rather than a log. The one-step bound would then hold outright at next order.
- If that checks out, 3D Yang–Mills is the natural proving ground for I2:
  - it is superrenormalizable, with finitely many divergent diagrams;
  - it already has a repo, `yang-mills-3d-forest-gauge`;
  - its continuum limit is itself a major open result.
- **Test:** rerun the `LEAKAGE.md` free-vacuum calculation (`leak_time.py`) with d = 3 and check that the coefficient goes to 0 like a as the spacing shrinks.
- **Cost:** a few CPU-hours, no Modal.

**B2. Prove the soft-part bound S′ at order g² on the lattice, not just in the continuum.**
- This upgrades the continuum argument in LEAKAGE §5 into a lattice statement, with explicit constants from the lattice propagator.
- It uses the same formulas as `slice/composite_free.py` plus one vertex.
- It is a concrete perturbative theorem: it doesn't solve I2, but it shows the obstruction is purely nonperturbative.

**B3. Connect to stochastic quantization.**
- Chandra–Chevyrev–Hairer–Shen build the 3D Yang–Mills(–Higgs) gauge-orbit state space using the Yang–Mills heat flow. That is the same object as the paper's flowed observables.
- Shen–Zhu–Zhu (arXiv:2204.12737) give Langevin mass-gap results for the lattice theory at strong coupling.
- **Concrete question:** are the paper's spatially flowed slice fields continuous functions on the CCHS orbit space? If so, I2-type statements for short Euclidean times could be imported in 3D.
- This step is a literature check, then possibly a lemma.

**B4. The 3D unflowed Monte Carlo test (I11).**
- It came back negative: at our statistics the signal is pure noise from two lattice steps on.
- Only worth redoing with a multilevel estimator and about 100× the statistics.
- Lower priority than B1, which answers a sharper question for less.

## Track C: the other two inputs

**C1. Cube condition (I1).** The proved size thresholds (10¹⁵ and up) are far from the numerical estimate (about 64).
- Certified local eigenvalue computations would shrink the constants.
- But the composition step (Remark 101) is a genuine gap. Work on it only if someone has a structural idea for it.

**C2. Rare bad fibres.** Restated as "bad fibres are rare and harmless", it looks provable by disagreement percolation from the paper's Theorems 122 and 131. That's a self-contained combinatorial lemma and good to hand off. *(Outcome in `rare/RARE.md`: "rare" holds by chessboard estimates; disagreement percolation fails at weak coupling.)*

**C3. Tree bounds (I3) and continuum existence.** These are Bałaban territory. I'd only attack them through B1: if 3D works, its tree bounds are the next target there.
