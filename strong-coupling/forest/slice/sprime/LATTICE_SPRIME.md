# Proposition S′ at order g² on the lattice (direction B2)

`../../DIRECTIONS.md` (B2) asks for the robust slice-ratio statement S′ of
`big-bang-monotile/yang-mills-continuum/LEAKAGE.md` §5 at order g², uniformly in the lattice spacing a, in d = 4.

S′ is proved in Lean (`lean/SliceRatio/Operator.lean`, `slice_robust`): if u = v + w and ‖Sv‖² ≥ c‖v‖², then
‖Sⁿu‖ ≥ c^{n/2}‖v‖ − ‖w‖. Here S is the one-step transfer matrix, T = S², and u = F|Ω⟩ is the flowed slice state.
Two inputs are needed:

- (i) a one-step bound for a soft part v, uniform in a;
- (ii) a bound on ‖w‖/‖v‖ for the hard remainder, uniform in a.

**Result.** Both hold at order g² in the free vacuum, uniformly in a. The vertex bound (V) that the hard-norm
estimate needs is proved below (Lemma 3); it holds for elementary reasons on the lattice. The
logarithm of LEAKAGE.md §1 sits entirely in the energy of w, which S′ never uses.

## Setting

- Lattice units a = 1; three-dimensional slices of side L = 4√τ, so the soft modes are the same at every spacing.
- Free transverse vacuum: ⟨|A(p)|²⟩ = 1/(2p̂), p̂² = Σ 4 sin²(pᵢ/2).
- F = F₂ + F₃ + O(g²), where F₂ is quadratic (leading order) and F₃ cubic, the soft–hard–hard fusion term of LEAKAGE.md §2.
- Soft/hard split at Λ = M/√τ: w_Λ is the part of F₃|Ω⟩ with a leg of momentum |k̂| ≥ Λ; v = F₂|Ω⟩ + (the rest).

## Lemma 1 (soft one-step bound, leading order)

F₂|Ω⟩ is a superposition of pairs (P, −P) of energy 2P̂, with weights w_P = P̂⁴ e^{−4τP̂²}/(4P̂²). By Jensen,

    −log(‖S v‖²/‖v‖²) = −log⟨e^{−H}⟩_v ≤ ⟨H⟩_v = Σ_P w_P 2P̂ / Σ_P w_P =: Ē/√τ.

`hardsum.py` evaluates Ē on the lattice:

| τ/a² | 1 | 4 | 16 | 64 | 256 | 1024 | 4096 |
|---|---|---|---|---|---|---|---|
| Ē | 2.830 | 3.062 | 3.122 | 3.137 | 3.141 | 3.142 | 3.142 |

Ē increases to its continuum value ≈ π. Ē is bounded by its continuum value up to 0.01% at every spacing tested, so
c = e^{−π/√τ} (in physical units, per lattice step a: e^{−πa/√τ}) is a uniform one-step constant. The soft components
of F₃ with at most three legs of momentum < Λ have energy < 3Λ = 3M/√τ in the free theory, so the bound extends to
all of v with Ē replaced by max(Ē, 3M).

## Lemma 2 (hard norm)

Use the vertex bound

**(V)** the bilinear vertex R_lat(k, l) of the lattice flow drift obeys |R_lat| ≤ c_V(|k̂| + |l̂|) (Lemma 3).

The flow-time integral obeys |I(k, l)| ≤ 2e^{−τP̂²}/(k̂² + l̂²). With (V), each hard pair
contributes at most c g² e^{−2τP̂²}P̂⁴ k̂^{−4} times the soft-leg weight, and

    ‖w_Λ‖² / ‖v‖² ≤ c g² · [(1/V) Σ_{|k̂|≥Λ} k̂^{−4}] · [Σ_P P̂⁴e^{−4τP̂²}/(2P̂)] / [Σ_P w_P] = c g² · N(M),

    N(M) = H(M) Ē / M,   H(M) = Λ · (1/V) Σ_{|k̂|≥Λ} k̂^{−4}.

The second bracket is again the weighted mean of 2P̂, i.e. Ē/√τ. The sum Σ k̂^{−4} converges in the UV for
three-dimensional slices: in the continuum it is ∫d³k/(2π)³ k^{−4} over |k| ≥ Λ, which equals 1/(2π²Λ).

| τ/a² | 4 | 16 | 64 | 256 | 1024 | 4096 | continuum |
|---|---|---|---|---|---|---|---|
| H(2) | .0610 | .0589 | .0568 | .0555 | .0548 | .0544 | .0507 |
| H(4) | .0436 | .0543 | .0524 | .0504 | .0491 | .0484 | .0507 |
| N(2) | .093 | .092 | .089 | .087 | .086 | .086 | .080 |
| N(4) | .033 | .042 | .041 | .040 | .039 | .038 | .040 |
| N(8) | 0 | .016 | .022 | .021 | .021 | .021 | .020 |

At τ/a² ≤ 4, H(8) = 0 because Λ exceeds the zone corner √12; there are no hard modes. N(M) is bounded uniformly in
a and decreases like 1/M.

## Lemma 3 (the vertex bound (V))

Let Z(A) be the drift of the lattice flow written in the exponential chart U = e^{A}: the Wilson force
vec(U_μ(x) Σ staples), plus, if one uses it, a lattice DeTurck term D_μω with ω = α₀ ∂̂·A. Its bilinear part is
Z₂(A, A)(P) = Σ_{k+l=P} R_lat(k, l) A(k) A(l), with R_lat a trigonometric polynomial on T³ × T³, because Z is a
translation-invariant, finite-range, real-analytic function of the links. Then:

- **R_lat(0, 0) = 0.** For constant links U_μ = e^{εX_μ}, each term U_μ·(staple) is a plaquette holonomy. The forward
  and backward plaquettes in the plane (μ, ν) are 1 ± ε²[X_μ, X_ν] + O(ε³), so their ε² parts cancel in the sum. The
  DeTurck term vanishes on constant fields because ∂̂·A = 0 there and ω stays 0 to all orders.
- **Hence (V).** R_lat is smooth and vanishes at the origin, so |R_lat(k, l)| ≤ ‖∇R_lat‖_∞ (|k| + |l|) with
  representatives k, l ∈ [−π, π]³. Since 2 sin(t/2) ≥ 2t/π on [0, π], |k| ≤ (π/2)|k̂|. So
  c_V = (π/2)‖∇R_lat‖_∞, a number independent of a and τ. ∎

On the lattice a quadratic piece l̂² − k̂² is harmless: in lattice units |k̂| ≤ √12, so it is bounded by √12(|k̂| + |l̂|).
What would be fatal is a vertex not vanishing at zero momentum, and gauge structure excludes that.

`vertex_lat.py` checks this numerically for the SU(2) Wilson force. It computes the exact bilinear part by a symmetric
second difference of plane waves in different colours (`vertex_lat.out`):

| lattice | momenta | pairs | max \|R\|/(\|k̂\|+\|l̂\|+\|P̂\|) |
|---|---|---|---|
| L = 8 | whole zone | 400 | 0.354 |
| L = 32 | \|nᵢ\| ≤ 2 (continuum regime) | 200 | 0.161 |
| L = 32 | \|nᵢ\| ≤ 6 | 200 | 0.132 |

For constant fields the bilinear part is 0 to machine precision.

The bound on I(k, l) uses damping e^{−sk̂²} of every leg. For the pure Wilson flow this holds only for the transverse
legs: the longitudinal modes are zero modes of the linearised flow. Since F is gauge-invariant, one may flow with the
DeTurck-modified drift instead (Lüscher's α₀ trick), which damps all legs and leaves F unchanged. Lemma 3 covers
that drift too.

## Corollary (S′ at order g², free vacuum)

For every M ≥ 2 and every a,

    ‖Sⁿ F|Ω⟩‖ ≥ e^{−n max(Ē, 3M)/(2√τ)} ‖v‖ − (c g² N(M))^{1/2} ‖v‖ + O(g⁴),

with Ē ≤ π and N(M) ≤ 0.1·2/M from the tables. Fix a physical time t = na and choose M with (c g² N(M))^{1/2} < ½ e^{−t max(π, 3M)/(2√τ)}. This is possible at small
g, since N(M) does not depend on g or a. Then ‖Sⁿu‖/‖u‖ is bounded below by a constant independent of a. The
logarithm of LEAKAGE.md §1 is the energy of w_Λ; it does not appear.

## What this does not cover

1. **O(g²) corrections to the vacuum** (the interacting ground state) are not included, as in LEAKAGE.md §3. They
   change both ‖v‖ and the weights; they are tied to the coupling renormalisation.
2. **The constant c** is not tracked. It contains c_V² (Lemma 3, not evaluated in closed form), the colour factor N(N²−1)/16 and the 4/3 angular factor of
   LEAKAGE.md §2.
3. **`leak_norm.out`.** The measured total order-g² change of Var(F) per g², 0.785 → 1.035 over τ/a² = 1 → 16, drifts.
   Lemma 2 predicts that the hard part is uniform. The drift is therefore in the soft part (lattice artefacts at
   small τ/a², where the soft modes sit near the cutoff) or in F₃ × F₃ soft–soft terms that Lemma 2 does not separate.
   This is not resolved here; a split of that measurement at Λ = M/√τ would decide it.
4. Order g² only. At realistic couplings g₀² ≈ 1.5 these are 5–15% effects (LEAKAGE.md §5).

**Status of B2:** S′ holds at order g² in the free vacuum, uniformly in a, with the hard norm controlled by
Lemmas 2–3 and the uniform lattice sums. Open: the O(g²) vacuum correction, the explicit constant, and the
interpretation of the `leak_norm` drift.

## Files

| File | Contents |
|---|---|
| `hardsum.py` | Lattice sums Ē, H(M), N(M) at τ/a² = 1, …, 4096 |
| `hardsum.out` | Output |
| `vertex_lat.py` | Bilinear part of the SU(2) Wilson-flow drift vs \|k̂\|+\|l̂\|+\|P̂\| (`python3 vertex_lat.py L NPAIRS NMAX`) |
| `vertex_lat.out` | Output for (8, 400, all), (32, 200, 2), (32, 200, 6) |
