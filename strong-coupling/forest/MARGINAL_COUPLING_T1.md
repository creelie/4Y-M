# T-1: near/far decomposition and exact SU(2) marginalization (d = 4)

Status: **draft.** §1–§3 are elementary and exact. §4 is a crude bound that shows
where the difficulty lies. §5 is a plan. Nothing here is interval-certified yet.
dump10 deliverable T-1 (Tristan, first 7 days); feeds T-2/T-3 and
direction #2 of `../CLAY_NEXT_DIRECTIONS.md`.

Conventions: SU(2) ≅ unit quaternions S³, normalized Haar measure dU, Wilson weight
exp(β ⟨U_p⟩) per plaquette with ⟨U⟩ = ½ Tr U = Re U (so β = β_W/… as in
`gauge-fixed/`; **check the normalization against `blocks/sharp.py` before using
any number below**). Box V (interior links), boundary ∂V = links outside V that
share a plaquette with V, exterior configuration ω.

## 1. Single-link identity (exact)

For a quaternion S ∈ R⁴,

  ∫_{SU(2)} exp(β ⟨U, S⟩) dU = Z_β(|S|),  Z_β(s) := 2 I₁(βs)/(βs),

where ⟨U,S⟩ = Re(U S̄) and |·| is the quaternion norm. Proof: by invariance, reduce to
S = s·1; then ∫ e^{βs cos θ} (2/π) sin²θ dθ = 2I₁(βs)/(βs). If link e lies in
plaquettes p₁…p₆, its Boltzmann factor is exp(β⟨U_e, S_e⟩), where S_e = Σ_k A_{p_k}
is the sum of staples. Its conditional law is von Mises–Fisher vMF(β S_e), and
integrating it out leaves the weight Z_β(|S_e|).

## 2. The near set of j

Let j = (x, μ) be a changed boundary link. Its near set is

  P_j = {(x + σν, μ) : ν ≠ μ, σ = ±1}   (6 links in d = 4).

(a) Each e ∈ P_j shares exactly one plaquette with j, namely p(j,e) = the (μ,ν)
plaquette spanned by j and e.
(b) P_j is plaquette-independent. Two parallel μ-links share a plaquette only if
their base points differ by a unit vector. Here x+σν and x+σ′ν′ differ by
σν − σ′ν′, which is never a unit vector unless (σ,ν) = (σ′,ν′).
(c) Every plaquette containing j contains exactly one element of P_j. The six
plaquettes through j are p(μ, ν, ±) for ν ≠ μ, and each contains one parallel
partner (x ± ν, μ).

## 3. Marginalization identity (exact)

Let N = P_j ∩ V (near links in the box) and F = V \ N (far links). By (b), the
links of N are conditionally independent given F and ω, so

  ρ^F_ω(U_F) ∝ Π_{e∈N} Z_β(|S_e(U_F, ω)|) · Π_{p ∩ N = ∅} exp(β⟨U_p⟩).

By (c), j enters ρ^F_ω only through:
- one staple of each S_e, e ∈ N, changing by a quaternion of norm ≤ d(Q_j, Q′_j);
- the plaquettes p(j,e) with e ∈ P_j \ V (partner fixed by ω). In these, j enters
  directly through the Wilson factor. Call their number m_j ≤ 6.

The conditional law of N given F is Π_{e∈N} vMF(β S_e), also explicit.

Lipschitz constants. d/ds log Z_β(s) = β I₂(βs)/I₁(βs), which increases in s and is
< β. Since |S_e| ≤ 6:

  |Δ log Z_β(|S_e|)| ≤ β R(β) d_j,  R(β) := I₂(6β)/I₁(6β) < 1.

So the j-dependence of ρ^F has Lipschitz constant at most

  Λ_j = β (|N| R(β) + m_j)

per unit of d_j, against 6β before marginalization. The gain comes from interior
near links (|N| large, m_j small), which happens when j sits on a face of the box
rather than on an edge or corner.

## 4. A crude coupling bound, and where it fails

With TV(ρ^F_ω, ρ^F_{ω^j}) ≤ tanh(Λ_j d_j / 2), maximally couple F. Then couple N
through the vMF conditionals, which are equal unless S_e changed. The expected
disagreement is

  E d_V ≤ TV · 2|F| + Σ_{e∈N} W₁(vMF(βS_e), vMF(βS′_e)).

Summed over j ∈ ∂V, the first term is ~ |∂V| · Λ · |V|. That is not < |V| at any
useful β, so **a maximal TV coupling of the whole far marginal pays the volume
factor**. This is presumably why the naive coupling gave r = 1.48 at β_W = 0.30.
The marginalization must therefore be combined with a W₁ (Dobrushin-type) treatment
of the far marginal, not a TV one.

## 5. Plan for T-2 (all-boundary certification)

1. Reduced specification. On F, the measure ρ^F_ω has local interactions: Wilson
   plaquettes plus the Z-terms. Each Z-term couples the up to 18 links of the
   staples of one e ∈ N. Compute single-link Dobrushin influences C_{f,g} for this
   reduced specification, using vMF-type conditionals with Bessel-ratio
   Lipschitz bounds. The new boundary influence of j is ≤ β R(β) per Z-term.
   Target: Σ_j C_{·,j} + Σ_g C_{·,g} < 1 in the GDS sense for the box.
2. Supremum over boundaries. Every bound above depends on ω only through |S_e| ≤ 6
   and through the monotone ratios I₂/I₁. The worst case is therefore |S_e| = 6
   (aligned staples). This is an explicit, uniform bound, with no interval covering
   needed at this level of crudeness. A refined version would cover |S_e| ∈ [0,6].
3. Certification. Bessel values via python-flint/Arb (available in the cluster
   venv), then rational bounds checked in Lean as in `blocks/lean/`.
4. Benchmark: d = 4, L = 2, β_W = 0.30. Then try to pass 0.351 (frozen benchmark,
   `BALL_CERTIFICATE.md`).

Risk: the Z-terms enlarge the interaction range (18-link cliques). The Dobrushin sum
for f ∈ F may get worse even though the boundary term improves. Step 1 is exactly
the computation that decides this, and it should be done before any certification
effort.

## 6. Step 5.1 result (job 1679093, 2026-11-06): raw form fails

Raw potential-Lipschitz sums (β = β_W, no vMF W₁ factor), `research/experiments/reduced_dobrushin/`:

| L | β_W | R(6β) | max row unred → red | max ∂-column unred → red |
|---|---|---|---|---|
| 2 | 0.30 | 0.40 | 5.40 → 6.54 | 0.90 → 1.08 |
| 2 | 0.351 | 0.45 | 6.32 → 7.95 | 1.05 → 1.42 |
| 3 | 0.30 | 0.40 | 5.40 → 6.54 | 1.80 → 1.80 |
| 3 | 0.351 | 0.45 | 6.32 → 7.95 | 2.11 → 2.37 |

So the risk named in §5 is real. The 18-link Z-cliques add more interior partners than
the R(6β) < 1 factor saves, and both the row sums and the boundary columns get worse.
Marginalizing the whole near set P_j for every j at once is not a win in raw form.
Redesign options:
(a) Keep the Z-terms but bound them by the *second* derivative of log Z_β, since the
    cross-influence of two staple links through one Z-term is O(β²), not βR.
(b) Marginalize only the near links of the single j being tested, one j at a time
    (that is what T-1 actually needs).
(c) Include the vMF W₁ conversion factor before comparing.
Option (a) is the cheapest next computation.

### 6a. Variant (a), Hessian bound (job 1679518)

Numerically sup_{|S|≤6} ‖Hess_S log Z_β‖ = β²/4 exactly (the small-|S| limit), i.e. 0.0225 at
β_W = 0.30 against βR = 0.12. With that weight for Z-term cross influences:

| L | β_W | max row unred / red | max ∂-col unred / red | max_j col ratio red/unred |
|---|---|---|---|---|
| 2 | 0.30 | 5.40 / 5.40 | 0.90 / 0.90 | 1.00 |
| 2 | 0.351 | 6.32 / 6.32 | 1.05 / 1.05 | 1.00 |
| 3 | 0.30 | 5.40 / 5.40 | 1.80 / 1.80 | 1.13 |
| 3 | 0.351 | 6.32 / 6.32 | 2.11 / 2.11 | 1.32 |

So the reduced specification is no longer worse in the maxima, but it is not better either.
The maxima are attained by links and boundary links that the marginalization does not touch.
At L = 3, some boundary links still get a larger column (up to ×1.32). Under the decision
rule this is a **fail**: no gain in the quantities that enter the Dobrushin sum.
Next: variant (b), marginalizing only the near set of the single tested j and reporting
per-j column ratios (min/median/max), not just the global maxima.

### 6b. Variant (b), one j at a time (job 1680007): park #2

| L | β_W | total ∂-influence unred → red | per-j ratio min / q25 / median / q75 / max | j improved |
|---|---|---|---|---|
| 2 | 0.30 | 86.4 → 77.0 (−11%) | 0.675 / 0.675 / 0.675 / 1.00 / 1.00 | 96/160 |
| 2 | 0.351 | 101.1 → 94.0 (−7%) | 0.79 / 0.79 / 0.79 / 1.00 / 1.00 | 96/160 |
| 3 | 0.30 | 388.8 → 366.1 (−6%) | 0.675 / 0.675 / 0.90 / 1.00 / 1.125 | 384/648 |
| 3 | 0.351 | 454.9 → 449.6 (−1%) | 0.79 / 0.79 / 1.00 / 1.05 / 1.32 | 192/648 |

Face-type j gain a factor 0.675 (β_W = 0.30), but edge and corner j gain nothing and some
get worse. The total gain is 11% at the benchmark, under the 20% threshold, and it shrinks
with β_W and with L. Decision: **near-set marginalization gives no usable gain. #2 is parked**
and #4 (weak gap in the RP state) is next, as the #3 result also suggested.
