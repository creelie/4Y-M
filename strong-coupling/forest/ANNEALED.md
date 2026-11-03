# Item 2: an annealed Dobrushin comparison (November 2026)

This is item 2 of `ROADMAP.md` §7, caveat 1 of `D_STATUS.md` §3. The question: can Dobrushin's comparison
(Lemma `comparison` of `main.tex`) be run with coefficients that are **good off a random bad set B(σ)** and only
**trivially bounded on it**, when B depends on the configuration and is rare only in the chessboard sense
P(S ⊂ B) ≤ q^{|S|}?

The outcome is **negative for the cheap versions**:
- §2 proves that every argument which uses only the one-point rarity of B has a floor of order q. It gives closeness
  of two Gibbs measures to O(q) per site, but no decay.
- §3 explains why iterating loses the rarity at every step.
- §4 computes that smoothing the bad set with an auxiliary variable trades Lipschitz cost against rarity, and loses.

One route survives (§5): a mixture over the exact large-field set, plus a bootstrap. Its abstract lemma is proved
there. Its inputs are a quenched coefficient bound for laws conditioned on the large-field set, which is not available,
and an expansion of the large-field weights, which is Bałaban's. So caveat 1 does not reduce to a soft probability
lemma. §6 records what this changes in `D_STATUS.md` §3.

## 1. Setting

Index the sites by e (links, or blocks V_x). For two configurations σ, σ′, take the coupled single-site (or
single-block) updates of Lemma `comparison`, with W₁ responses

  W₁(γ_e(·|σ), γ_e(·|σ′)) ≤ Σ_{e′} c_{ee′}(σ, σ′) d(σ_{e′}, σ′_{e′}).

The random coefficients satisfy:
- **(G)** c_{ee′} ≤ a_{ee′} when e is good in both σ and σ′, with sup_e Σ_{e′} a_{ee′} e^{λd(e,e′)} ≤ θ < 1;
- **(B)** c_{ee′} ≤ A_{ee′} always, with Σ_{e′} A_{ee′} ≤ Ā. For links, Ā = 4.5β_W (`D_STATUS.md` §3);
- **(P)** "e is bad" is the event b_e(σ), which depends on σ near e, and μ(b_e for all e ∈ S) ≤ q^{|S|}.

Let D = diam G.

## 2. Proposition 1: rarity in one-point form gives an O(q) floor (proved)

**Claim.** Assume (G), (B), and only the one-point bound μ_i(b_e) ≤ q. Then the minimising coupling of Lemma
`comparison` for μ₁, μ₂ (both Gibbs for γ) satisfies

  sup_e ρ(e) ≤ 2ĀDq/(1−θ),  ρ(e) = E d(σ_e, σ′_e).

The bound does not decay with the distance from the place where μ₁ and μ₂ are allowed to differ.

**Proof.** By minimality, as in Lemma `comparison`,

  ρ(e) ≤ E Σ_{e′} c_{ee′} d_{e′} ≤ Σ_{e′} a_{ee′} ρ(e′) + E[(1_{b_e(σ)} + 1_{b_e(σ′)}) Σ_{e′} A_{ee′} d_{e′}],

and the last term is at most 2ĀDq. Take the supremum over e. ∎

The floor is not an artefact of using D as the bound for d. Hölder's inequality gives
E[1_b d] ≤ (qD)^{1/p} ρ^{1/p′}, which leads to the sublinear recursion ρ ≤ θρ + Ā(qD)^{1/p}ρ^{1/p′}. That recursion has
the positive fixed point ρ* ≍ Ā^p qD/(1−θ)^p. Weighted versions behave the same way: e^{−md} becomes e^{−md/p′} at
each step, so the rate degrades geometrically. **Any bound that sees B only through its one-point marginals has an
O(q) floor.** Removing the floor needs the joint law of the bad set and the disagreement.

What Proposition 1 does give is approximate uniqueness. Every two infinite-volume Gibbs measures of γ are within
2ĀDq/(1−θ) per site in W₁. Through the coupling, local covariances satisfy |Cov(F, H)| ≤ ‖H‖_∞ L_F (e^{−λd}·C +
C′q|S_F|). That is a correlation floor of order q, not a gap.

## 3. Why iterating loses rarity (structural)

The quenched proof [BD86] expands ρ along paths e₀ → e₁ → … with the product Π c_{e_k e_{k+1}}, and pays A only at
the bad sites of the path. In the annealed case the expansion needs

  E[1_{b_{e₀}} d_{e₁}] ≤ E[1_{b_{e₀}} Σ c_{e₁e₂} d_{e₂}],

and the right side requires minimality at e₁ with the event b_{e₀} held fixed. But b_{e₀} depends on σ near e₀, which
includes e₁ (e₁ is a neighbour of e₀ wherever c_{e₀e₁} ≠ 0), and the update at e₁ resamples exactly those spins.

Dropping from the bad set the sites near the update keeps the inequality valid. But then the bad site added at step
k is dropped at step k + 1, because consecutive path sites are neighbours. The chessboard factor q^{|S|} never
accumulates. Hence the floor of §2 at every order.

In the TV setting, disagreement percolation [vdBM94] escapes this by exploring the disagreement cluster from the
boundary, so that the explored spins are fixed and the badness of explored sites is quenched for the rest. It needs
a coupling that agrees exactly with probability bounded below at good sites. At weak coupling no such coupling exists
(`rare/RARE.md`: "harmless" cannot come from disagreement percolation).

## 4. Proposition 2: an auxiliary soft split loses (computation)

The natural fix makes the bad set quenched by adding an auxiliary label. Fix a small-field profile b(x) ∈ [0, 1]:
- b = 1 for δ(U_p, I) ≤ K;
- b = 0 for δ(U_p, I) ≥ 2K;
- b is L-Lipschitz, with L ≥ 1/K.

Write 1 = θb + (1 − θb), and sample ω_p = "good" with probability θb(U_p). Given ω, σ is Gibbs with the extra factors
θb on good plaquettes and 1 − θb on bad ones. With δ = 1 − θ:
- **rarity:** P(p bad) ≥ δ, so q ≥ δ;
- **cost on ω:** the factor 1 − θb ∈ [δ, 1] has log-derivative ≤ L/δ. It adds up to L/δ·D to each single-link
  coefficient, so the per-site island cost grows by a factor (Ā + 6LD/δ)⁴.

A site of an island then costs at least (L/δ)⁴·δ = L⁴/δ³ ≥ L⁴ ≥ K⁻⁴ ≫ 1, at every choice of δ ≤ 1. **The (P)-versus-(B)
balance of `D_STATUS.md` §3, which needed per-site cost × q < ½, cannot hold.**

The hard split (δ → 0, b an indicator) has the same problem in another form. Conditioning on "p large" or "p small"
restricts a link to an intersection of balls and ball complements that move with the neighbours. Its W₁ response is
the boundary density over the mass of the allowed set, and that is unbounded when the set is a thin lens. So the
coefficient is not bounded uniformly in the field, and (B) fails. With a randomised threshold τ_p uniform in [K, 2K],
the bad set is ω = {δ_p > τ_p}, and the factor on ω is a linear ramp g. Then E|g′/g| ≍ (1/K)·μ(shell)/μ(δ_p > τ_p).
When λK² ≫ 1 the single-link density decreases across the shell, so this ratio is at least of order one, and the
cost is at least 1/K per link. Four links per site then cost ≳ K⁻⁴ against a rarity factor that is at most one, so the
cost again wins.

These are order-of-magnitude computations. The constants in the exponents are not optimised, and the conclusion holds
for every choice.

## 5. The surviving route: mixture over the large-field set, plus a bootstrap

Write μ = Σ_ω μ(B = ω) μ^ω, where μ^ω = μ(· | B(σ) = ω) and B is the exact bad set. Then

  Cov_μ(F, H) = E_ω[Cov_{μ^ω}(F, H)] + Cov_ω(μ^ω F, μ^ω H).

- **The first term is quenched.** For fixed ω, the [BD86] path expansion applies if μ^ω has coefficients (G) off ω
  and (B) on ω. The annealed average then uses only P(S ⊂ B) ≤ q^{|S|}, which is (P). This is where chessboard
  estimates enter correctly.
- **The second term is a covariance of the environment.** μ^ω F depends on ω near S_F up to exponentially small
  corrections, by the quenched decay. So it is a covariance of (smoothed) local functions of the bad set under μ.
  Those are again local observables of μ, so the argument is circular. The circle closes as a bootstrap.

**Lemma 3 (bootstrap; proved).** Let Γ : ℕ → [0, Γ₀] satisfy

  Γ(d) ≤ C₀ e^{−m d} + ε Γ(d − r)  for d ≥ r.

If ε e^{m′r} ≤ ½ for some m′ < m, then Γ(d) ≤ (2C₀ + Γ₀ e^{m′r}) e^{−m′d} for every d.

*Proof.* Let M(d) = sup_{d′≤d} Γ(d′) e^{m′d′}. For d ≥ r,

  Γ(d) e^{m′d} ≤ C₀ e^{−(m−m′)d} + ε e^{m′r} Γ(d − r) e^{m′(d−r)} ≤ C₀ + ½ M(d),

while Γ(d) e^{m′d} ≤ Γ₀ e^{m′r} for d < r. Hence M(d) ≤ max(Γ₀ e^{m′r}, C₀ + ½M(d)), which gives
M ≤ 2C₀ + Γ₀ e^{m′r}. ∎

Here Γ(d) would be the supremum of |Cov_μ| over local observables at distance d with unit Lipschitz constants and
fixed supports. ε would come from replacing the bad-set indicators by Lipschitz approximations (cost 1/K per
indicator) and from the rarity of bad sets near the supports: ε ≍ q·|N|²/K², with N the support neighbourhood.
Since q is e^{−c√β_W} at ρ = β_W^{−1/4} (`rare/RARE.md`), the smallness ε e^{m′r} ≤ ½ is not the problem.

**The missing inputs.**
1. **(I1) Quenched coefficients for μ^ω.** These are laws conditioned on the exact large-field set: the thin-lens
   problem of §4 at the boundary of ω, and the off-chart problem of `D5.md` §3 inside ω^c, now with ω^c of arbitrary
   shape. The small-field part of (I1) is (L) of `D5.md` §4. The large-field part is not bounded by any estimate
   available here.
2. **(I2) Locality of μ^ω F in ω.** This follows from the quenched expansion once (I1) holds.

In Bałaban's work the large-field regions are not handled by conditioning. They are handled by an expansion in which
the large-field set is a summation variable, its rarity is a weight e^{−cβK²} per site, and the coupling to the
small-field region goes through analyticity. **So (I1) is a piece of Bałaban's large-field expansion, not a
probability lemma.**

## 6. What changes in `D_STATUS.md` §3

The island observation stands. The island cost (4.5β_W)^{4n} against the rarity q^n holds once β_W ≥ 1.1·10⁷. But
caveat 1 is harder than "a self-contained probability lemma":
- with rarity used only through probabilities, the comparison has an O(q) floor (Proposition 1);
- the honest way to use rarity is quenched, through the mixture of §5, and that needs bounds on the conditioned laws
  (I1).

So:
- **proved (new):**
  - approximate uniqueness at O(q) per site (Proposition 1);
  - the bootstrap reduction (Lemma 3);
- **open:** (I1). Its small-field part is (L), and its large-field part belongs with D4.

The large fields remain a separate wall. "D5 is the wall" in `D_STATUS.md` §3 is too optimistic, and that file is
updated to say so.

## References

- [BD86] L. A. Bassalygo, R. L. Dobrushin, Theory Probab. Appl. 31 (1986) 572–589.
- [vdBM94] J. van den Berg, C. Maes, *Disagreement percolation in the study of Markov fields*, Ann. Probab. 22 (1994)
  749–763.

**Follow-up (proved scope statement, November 2026).**
`ENVIRONMENT_OBSTRUCTION.md` gives an explicit rare-label hidden-phase
example satisfying every finite-intersection bound while retaining a
distance-independent environment covariance, even with zero quenched
covariance. `CONNECTED_NORMALIZERS.md` proves the connected normaliser
estimate for continuous changes inside the good product chart; extending
it to hard label flips and bad islands remains open.
