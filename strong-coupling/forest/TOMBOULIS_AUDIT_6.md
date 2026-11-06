# Direction #6: audit of Tomboulis 2007 (arXiv:0707.2179)

Status: **localization done from the literature, plus one exact 2D computation.** No numerics run yet.
Sources: Tomboulis 0707.2179; Ito–Seiler 0711.4930; Tomboulis reply 0712.2620; Ito–Seiler 0803.3019.

## 1. The argument in one paragraph

Periodic box with sides a·b^{n_i}. MK-type decimations give upper and lower bounds on Z in terms of a
decimated box Λ⁽¹⁾ with character coefficients c_j(1); for d ≤ 4, c_j(n) → 0 (Ito). Tomboulis
interpolates between the bounds with α ∈ (0,1], t > 0, h(α,t) = exp[−t(1−α)/α], c̃_j = α c_j(1):

  Z = F₀(1)^{h(α(t),t)|Λ⁽¹⁾|} Z₁(c̃(α(t))),   Z⁺ = F₀(1)^{h(α⁺(t),t)|Λ⁽¹⁾|} Z₁⁺(c̃(α⁺(t))),

with Z⁺ = ½(Z + Z⁻), Z⁻ the 't Hooft-twisted (vortex) partition function. If some t* has
**α(t*) = α⁺(t*)**, then Z⁺/Z = Z₁⁺/Z₁ exactly. Iterating to the strong-coupling regime gives a positive
't Hooft string tension, hence confinement at all β.

## 2. The missing link, stated as one inequality

Existence of t* (Ito–Seiler Problems 1–3) is obtained by continuing an implicit-function branch in λ ∈ [0,1].
The branch continues provided Ψ_t ≠ 0 along it, which reduces to

  (★) A(α) > A⁺(α),   A = ∂_α log Z₁ / (|Λ⁽¹⁾| log F₀(1)),  A⁺ the same with Z₁⁺.

Since Z₁⁺/Z₁ = ½(1 + ρ) with ρ := Z₁⁻/Z₁, (★) is equivalent to

  (★′) ∂_α log(1 + ρ(α)) < 0:  **the vortex ratio is decreasing under the uniform rescaling c_j ↦ α c_j**,

on the decimated lattice, along the path, at whatever couplings the path visits. Tomboulis proves (★) only
where the high-temperature expansion converges, and argues that it is invoked only there. Ito–Seiler dispute
this: Problems 2–3 arise at the first decimation step, where the couplings are the original ones. They also
note that implicit-function local solvability does not give global continuation (Kanazawa counterexample
Ψ = e^{−t} − 1 + 2λ). Our reading: the dispute is exactly whether (★′) is needed at weak coupling. It is a
**monotonicity of the vortex free energy in the coupling**, a correlation-inequality statement with no known proof
for SU(2) in d = 4.

## 3. Exact 2D computation: (★′) holds trivially, so 2D cannot discriminate

On a 2D torus with N plaquettes, Z₁ = Σ_j (α c_j)^N (c₀ = 1), and the twist gives
Z₁⁻ = Σ_j σ_j (α c_j)^N with σ_j = (−1)^{2j} for SU(2) (or the center character for SU(N) / U(1)).
With x_j = (α c_j)^N ≥ 0,

  1 + ρ = 2(1 + Σ_{σ_j=+1, j≠0} x_j) / (1 + Σ_{j≠0} x_j),
  ∂_α log(1 + ρ) = (N/α) [ y_+/(1 + y_+) − y/(1 + y) ] ≤ 0,  y_+ ≤ y,

since y ↦ y/(1+y) is increasing. So (★′) holds in 2D for **every** compact group and every coupling.
That includes U(1), for which the 4D conclusion is false (Coulomb phase). The plan's 2D test is therefore
uninformative, and the 3D/4D content of (★′) is where the obstruction lives.

## 4. The U(1) consistency test

The argument's conclusion is false for 4D U(1) at weak coupling (Guth; Fröhlich–Spencer), so some step must
fail for U(1). Ito–Seiler state that t* does not exist for U(1). Two cases:
- (★′) **fails** for 4D U(1) at weak coupling. Then (★′) is the decisive step, and an SU(2) proof of (★′)
  at weak coupling would carry the full non-abelian content of confinement. That is a sharp, publishable
  reformulation, but not easier than the original problem.
- (★′) **holds** for 4D U(1) too. Then the failure is elsewhere: the global continuation (Problem 3) or
  eq. (5.21). In that case (★′) is a red herring, and the audit must move to those steps.

Deciding which case holds is a numerical question: ρ(α) on 4⁴–6⁴ tori with twisted boundary conditions,
for SU(2) and U(1), with α scaling the Wilson coupling. It is a small Monte Carlo job (CPU, ≤ 16 cores).

## 5. Verdict for Clay

#6 localizes to (★′) or to the continuation step. Neither is a repair route with lower difficulty than the
window problem; at best it is a clean obstruction statement. Value: moderate as a write-up, low as a route
to the gap. Optional next step: the §4 Monte Carlo test.
