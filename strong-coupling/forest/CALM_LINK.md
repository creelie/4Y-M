# Linking the calm-data cube coefficient to the large fields (November 2026)

This note does three things:
- §1 merges the D5 work here with the OpenGluon note on calm data;
- §2 and §3 push the two open items of that note as far as the tools here allow;
- §4 states the wall that stops both pushes.

**The OpenGluon note.** romanov360/OpenGluon, branch `claude/cube-coefficient-calm-data` (commit cd76a3c),
`yang-mills/cube-coefficient/CALM_DATA.md`. It is cited below as **[CD]**. That repository does not use the labels D1–D8;
its "cube coefficient", "non-calm fibres" and "harmless large fields" are D5, D4 and D7 here.

Status labels:
- §1 is a comparison;
- §2 is a reduction (sketch), with one lemma at the level of a Gaussian computation;
- §3 is a proved abstract step, with a sketch-level input and an order-of-magnitude ledger (`d5/bootstrap.py`, `d5/bootstrap.out`);
- §4 is a diagnosis.

## 1. Merge with [CD]

| item | here | [CD] | combined status |
|---|---|---|---|
| links coupled to the cube by a plaquette, in no average | `D5.md` §2: C_G ≤ 2.95·10¹³ (c_Q = 4) | Theorem 1: C_* ≤ 7.4·10¹², with an explicit count of the 18 links sharing a plaquette and S(α_G) ≤ 2.46·10¹⁶ | **proved, C_* ≤ 7.4·10¹²** (given `chartdecay`(ii) and `explicitdecay`). The factor 4 is the sharper sharing count of [CD] |
| links entering one cut-off average | `D5.md` §2 case (b): sketch, no constant | Theorem 2: transport along the eight private links, \|Y\|₂ ≤ 2.20, leading constant C_b ≤ 9.2·10¹³ proved; the shell term Ξ_j open (Conjecture B) | **reduced to Ξ_j.** §2 below reduces averaged Ξ_j to the localised statement of §3 |
| split of exterior links | implicit | Lemma 1, checked by enumeration for n = 4, 6 | proved |
| cube condition on calm data | n ≥ 17C_G/α, i.e. 5.7·10¹⁴ at α = 0.9 | Corollary 2: n ≥ 2.3·10¹⁵(1 + C_Ξ/C_b)/α, conditional on Conjecture B and λε_W\|A_W\| ≤ 1 | conditional; not uniform in the volume |
| off-chart term | `D5.md` §3: volume-limited. `ISLANDS.md` §1: Proposition L, window form ε_loc (sketch) | Conjecture L: localised off-chart covariance, weighted by e^{−αd} and π(not calm at x) | Proposition L is the **mass** half of Conjecture L. §3 below derives the covariance form from it, up to a floor |
| averaging over the averages v | not used | Lemma `mixture` needs only ∫W₁(κ(·\|v), κ′(·\|v)) ρ(dv) | used in §2 and in `D4_D8.md` §1 |

The two notes agree on every constant they share, up to the sharing count. Neither has a bound on the non-calm fibres.

## 2. The shell term Ξ_j, averaged over v (reduction, sketch)

**Setting ([CD] §4).** Fix a link j of the second kind, with its block y, the eight private links w₀, …, w₇ of y, and
m_y, the largest distance of a path of y from straight. The transport of Theorem 2 is cut off by a piecewise-linear χ(m_y)
that equals 1 for m_y ≤ ρ_G/8 and 0 for m_y ≥ ρ_G. The remainder Ξ_j(F) is the part of d/dt E_t F that comes from the
shell {ρ_G/8 ≤ m_y ≤ ρ_G}: the kinks of χ, and the flux of the moving constraint across the kink surfaces
{m_y = r_k}.

**Step 1 (structure of Ξ_j).** On a fibre at fixed v, Ξ_j(F) = Cov_{π_v}(F, g_v), where g_v is the shell term. It is
supported on the shell and is a sum of surface measures on {m_y = r_k}, with speed |Y|₂ ≤ 2.20. Split
g_v = g_v⁺ − g_v⁻ into positive parts with masses M_v^±, and let M_v = M_v⁺ + M_v⁻. Since Cov(F, g^±) =
M_v^±(E_{π_v^±}F − E_{π_v}F),

  |Ξ_j(F)| ≤ M_v⁺ |E_{π_v^+}F − E_{π_v}F| + M_v⁻ |E_{π_v^−}F − E_{π_v}F|,  π_v^± = g_v^± π_v / M_v^±.

So **|Ξ_j(F)| ≤ M_v · H_v**, where H_v is the largest W₁ displacement of π_v under reweighting by a probability
density supported on the shell of block y. H_v is the response of the fibre law to **conditioning on a large-field
event at one block**, which is the "harmless large fields" statement at one block.

**Step 2 (the flux mass, averaged).** By Lemma `mixture` (averaged form, [CD] §5), only ∫ M_v H_v ρ(dv) is needed.
Bound it by sup_v H_v · ∫ M_v ρ(dv). Disintegrating along the map to the averages, ∫ M_v ρ(dv) is the flux of the
transport across the kink surfaces under π_E, the law given the exterior **without** fixing v. Under π_E the private links
are not determined by constraints. Each w_a has, given all other links, a von Mises–Fisher-type law on SU(2) with
concentration κ_a ≤ cλ, and each path of y depends on exactly one private link (P_{x′_a} = ℓ_a w_a).

**Lemma (radial hazard rate; Gaussian-level computation, not written for SU(2)).** Let U have density ∝ e^{κ Re tr(u₀*U)/2}
on SU(2), κ ≥ 1, and let g₀ ∈ SU(2) and 0 < r ≤ π/2. The density f of d(U, g₀) satisfies

  f(r) ≤ C(κ^{1/2} + κr) P(d(U, g₀) ≥ r),

with C absolute and independent of u₀ and g₀.

*Reason.* In normal coordinates, U is close to a Gaussian of width κ^{−1/2} around u₀, and d(U, g₀) is a non-central
χ variable with 3 degrees of freedom. Its hazard rate is at most C/σ before the centre and C(r − a)/σ² after it. ∎ (sketch)

The surfaces {m_y = r_k} lie in ∪_a {d_a = r_k}, and P(m_y ≥ r) ≥ max_a P(d_a ≥ r). Apply the lemma to each private link
given all the others, and average:

  ∫ M_v ρ(dv) ≤ 8 · 2.20 · C(√(cλ) + cλρ_G) Σ_k π_E(m_y ≥ r_k) ≤ C′ λ π_E(m_y ≥ ρ_G/8).

**Conclusion of §2 (sketch).** With H := sup_v H_v,

  ∫ |Ξ_j| ρ(dv) ≤ C′ λ H · π_E(m_y ≥ ρ_G/8).

This is Conjecture B of [CD] in the averaged form that Lemma `mixture` needs, **provided H = O(1)**. The averaging
over v is essential: at fixed v the private links are fixed by the constraints, and the hazard-rate lemma has nothing to
act on.

H = O(1) is not proved. It asks that conditioning on a large field at one block moves the law of a 1-Lipschitz
function by O(1), uniformly in the volume. That is Conjecture L applied to the event at block y, so **Conjecture B
(averaged) follows from Conjecture L plus the hazard-rate lemma**, and the two open items of [CD] become one.

## 3. Conjecture L by a window bootstrap: decay down to a floor

**Setup.** (E, c) is calm. F is 1-Lipschitz with support S_F, and G is a function of links within distance 1 of j,
with per-link gradient ≤ 1. Let Γ(d) be the supremum of |Cov_π(F, G)| over such pairs with d(S_F, j) ≥ d and
|S_F| ≤ s. The cube coefficient involves λ·Cov, so the quantity that matters is λΓ.

Take the window W_R(j), a cube of side 2R + 1 blocks around j, with the free set L₁ of `ISLANDS.md` §1. Condition on
the configuration "out" outside L₁:

  Cov_π(F, G) = E[Cov(F, G | out)] + Cov(F, E[G | out]).

1. **Inside.**
   - If the conditioned window has a calm completion, then Theorem `chartdecay`(ii), with Proposition L for the off-chart
     part, gives |Cov(F, G | out)| ≤ C₀ e^{−α_G d} + c·s·ε_loc(R).
   - Otherwise, use the trivial bound c·s.
   - So E[Cov(F, G | out)] ≤ C₀ e^{−α_G d} + c·s·(ε_loc + p_nc), where p_nc = π(W_R(j) has no calm completion).
2. **Outside.** E[G | out] depends on out only through the N_∂(R) links at the window's boundary. On calm windows, its
   sensitivity to each of them is at most C e^{−α_G R}, by the window cube coefficient (§1, with C = max(C_*, C_b)).
   Each such link is at distance ≥ d − R from S_F. So |Cov(F, E[G | out])| ≤ N_∂(R) C e^{−α_G R} Γ(d − R), plus the
   same non-calm term.

So Γ(d) ≤ C₀ e^{−α_G d} + η + ε Γ(d − R), with ε = N_∂(R) C e^{−α_G R} and η = 2c·s(ε_loc + p_nc). Lemma 3 of
`ANNEALED.md` applies verbatim with C₀e^{−md} replaced by C₀e^{−md} + η: the extra constant passes through the
iteration as Σ_k 2^{−k}η.

**Proposition (window bootstrap; proved given the two inputs of steps 1 and 2).** If N_∂(R) C e^{−α_G R/2} ≤ ½, then

  Γ(d) ≤ (2C₀ + Γ₀ e^{α_G R/2}) e^{−α_G d/2} + 2η.

Its inputs are chart decay on calm windows with a calm completion (`chartdecay`(ii)), Proposition L, and the window
cube coefficient at the window's boundary. That last input is item (a) of §1 plus the reduction of item (b) to Ξ_j, so it
inherits §2.

**Numbers (`d5/bootstrap.out`; SU(2), λ = β_W/4, K = 4.8·10⁻¹⁸, C = 9.2·10¹³).**

| quantity | value |
|---|---|
| bootstrap radius R | 6.5·10⁵ blocks |
| window | m = 2.9·10²⁴ blocks, \|L₁\| ≤ 2.4·10²⁶ |
| λ\|L₁\|ε_loc ≤ 10⁻³ | β_W ≥ 1.1·10⁶⁴ |
| λ\|L₁\|p_nc ≤ 10⁻³ with torus chessboard rarity (`rare/RARE.md`) | β_W ≥ 2.1·10⁹⁵ |
| λ\|L₁\|p_nc with the density bound `localcalm` (any exterior) | ≥ 10¹⁴⁸ at every β_W (grows like log β_W) |

**What this gives.**
- Decay of covariances at rate α_G/2 down to a floor 2λη, with constants independent of the volume.
- The window is far smaller than the cube of the cube condition (m ≈ 3·10²⁴ blocks against 6·10⁵⁷). So the off-chart
  threshold drops from 3.8·10⁹⁷ (`ISLANDS.md` §3) to 1.1·10⁶⁴.
- **Only one-point rarity of a non-calm window is used, not finite-intersection rarity.** The floor needs p_nc, the
  probability of one event at the window of j. This is weaker than what the annealed comparison of `ANNEALED.md` asked
  for.

**What it does not give.** The floor is not of the form of Conjecture L. Conjecture L weights the non-calm events by
e^{−αd(x,j)}, so a non-calm site far from j costs nothing. In the bootstrap, a non-calm window around j costs its full
probability, and the floor does not decay. For the cube condition at the lattice scale, the floor is harmless if λ|L₁|p_nc
times the number of exterior links of the cube is small. That is a volume-independent requirement, met on the torus at
β_W ≳ 10⁹⁵. Under an arbitrary exterior it would need the one-point bound on p_nc, which is §4.

## 4. The wall, stated exactly

After §§2–3, the calm-data cube coefficient with a volume-independent error needs one statement:

> **(R1) One-point calm rarity under an arbitrary exterior.** For every exterior of a cube, and every window W_R(j) in it
> with R = 6.5·10⁵, the law of the cube gives the event "W_R(j) has no calm completion" probability at most
> 10⁻³/(λ|L₁|).

The rest of the chain then holds at β_W ≳ 10⁶⁴ (sketch level): §3 for Conjecture L with a floor, §2 for Conjecture B, and
§1 for the cube condition. There are three ways to try (R1), and each fails.

1. **The density bound `localcalm`** is proved under every exterior. It controls the *fraction* of window positions that
   are non-calm, not every position. Even read as a one-point bound, its best value is
   p ≈ 4096m · 32(6 log β_W + 31)/(β_W ρ²), with ρ = K/(3072m). Then λ|L₁|p ≥ 10¹⁴⁸ for the bootstrap window,
   **independently of β_W** up to log β_W (third row of the table). Density rarity costs a power of β_W, and the factor λ eats it.
2. **Chessboard estimates** give e^{−β_Wρ²/16} rarity of finite intersections. They need reflection positivity, which
   needs the torus (or a reflection-symmetric boundary). An arbitrary exterior breaks this.
3. **A localised Peierls argument on calm data.** The aim is to compare a non-calm field on the window with a calm one
   in the same fibre, at action cost proportional to the large plaquettes removed. Interpolating between them costs:
   - Fibre fields can differ by up to 92.3κ per link (`offchart`(i), main.tex line 8560).
   - Block coercivity loses D_v² ≤ 224Q (main.tex line 6536).
   - So the interpolation can raise the action by about 224·cap(S)·K² on the region S it touches. The gain is of order K²
     per removed large plaquette.

   This is a constant-factor loss, but it is multiplied by the capacity of the region. It cannot be won by taking β_W
   large, because both sides scale with λ. *(Order-of-magnitude computation; constants not optimised.)*

**What would cross the wall.** In Bałaban's large-field expansion the thresholds increase from scale to scale (and
within a scale, in the multi-threshold versions). Each large-field site pays e^{−p(g)} under **every** exterior. This
uniformity is exactly what (R1) needs and what route 3 lacks, and Dimock (arXiv:1212.5562) extracts it in a
self-contained form for a model problem. So (R1) is not a lemma within reach of the tools here. It is the large-field
half of a renormalisation-group expansion, the same object as (I1) of `ANNEALED.md` and D4 of `D4_D8.md`.

**Lead (b) of the OpenGluon survey: Dobrushin–Pechersky (assessed, heuristic, November 2026).** The criterion in the
form of Conache–Kondratiev–Kozitsky–Pasurek (arXiv:1501.00673, Theorem 2.6) asks for two things:
- a contraction (2.9) with Σκ < 1, only on boundaries whose neighbours have h ≤ K;
- a Lyapunov bound under **every** boundary, π_ℓ^x(h) ≤ 1 + Σ_{ℓ′} c_{ℓℓ′} h(x_{ℓ′}), with
  c̄ = sup_ℓ Σ_{ℓ′} c_{ℓℓ′} < 1/Δ^χ, where Δ is the degree of the block graph and χ its chromatic number.

The first is the calm contraction we have (`D5.md`, C_* above). The second is a quantitative form of (R1): the badness h
of a block (a function of its plaquette sizes) has to be **attenuated** by a factor Δ^χ ≥ 80^{16} from the boundary
inward. At the linearised level it is not:
- a constant abelian-embedded field strength B is a classical solution of the lattice equations;
- with boundary links taken from it, the interior minimiser is the same constant field, and the fluctuations around it
  are O(β_W^{−1/2});
- so h(interior) ≈ h(boundary), i.e. c̄ ≈ 1.

The criterion therefore needs a **nonperturbative** attenuation mechanism. For SU(2) in a constant field the candidate is the
Nielsen–Olesen unstable mode: it fits once the region is wider than about B^{−1/2}, and with it the uniform field is
not the minimiser. That mechanism is the same multiscale, large-field question as before.

The same example bears on (R1) as stated: (R1) has to fail for exteriors carrying a uniform field B ≫ K on any window of side
≲ B^{−1/2}, unless the nonabelian instability rescues it. So the honest form of (R1) is not "under every exterior" but
**averaged over exteriors under the Gibbs measure**, or under a class of exteriors that excludes such backgrounds. Lemma
`mixture` and §2 already use only averaged statements. Whether the cube coefficient itself, with its supremum over
boundary conditions, survives uniform-field boundaries is part of D4. *(Heuristic: classical solution plus linearised
fluctuations; not a theorem.)*

**Conclusion for D.** Under (R1), the three remaining pieces of the calm-data cube coefficient reduce to sketch-level
statements with explicit, volume-independent thresholds:
- the shell term (§2);
- the localised off-chart covariance (§3);
- the non-calm windows, through the floor of §3.

(R1) itself is not available under an arbitrary exterior, and D8 is untouched. **D is blocked at the large-field
expansion.** The next steps on D are research-level, not refinements: a multi-threshold large-field expansion in the
style of Dimock, worked out for the block variables here. This is where the work on D stops. Other paths are in
`ROADMAP.md` §7.
