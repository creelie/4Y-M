# Flowed slice fields and the CCHS state space (direction B3)

Question (`../DIRECTIONS.md`, B3): are the paper's spatially flowed slice fields continuous functions on the
gauge-orbit state space of Chandra–Chevyrev–Hairer–Shen? If so, slice-ratio statements could perhaps be imported
from stochastic quantization in d = 3.

**Short answer.**
- For the full three-dimensional field: **yes**. Flowed gauge-invariant observables are exactly CCHS's
  "smoothened observables". They are continuous on the orbit space and separate its points.
- For the slice fields the paper actually uses (spatial flow of one time slice): **not covered by any existing
  theorem**. In d = 3 the slice sits exactly at the regularity CCHS handle, but in two dimensions, and that version has
  not been written. In d = 4 the slice is half a derivative too rough for any deterministic state space.
- **Nothing can be imported yet.** CCHS construct the state space and the Langevin dynamics, not the invariant
  measure. The three-dimensional Yang–Mills measure on that space is still open.

## 1. What CCHS prove

Source: A. Chandra, I. Chevyrev, M. Hairer, H. Shen, *Stochastic quantisation of Yang–Mills–Higgs in 3D*,
Invent. Math. 237 (2024) 541–696, arXiv:2201.03487. References below are to the arXiv version.

- **Theorem 1.2 (state space).** For η ∈ (−½ − κ, −½) there is a complete metric space (S, Σ) of distributional
  connections on T³ with S ↪ C^η(T³).
  - (ii) The DeTurck–Yang–Mills flow F_t, taken modulo gauge as F̃_t, extends continuously to all of S for every
    t > 0.
  - (iii) X ∼ Y :⇔ F̃(X) = F̃(Y) extends gauge equivalence, and the orbit space O = S/∼ is separable and
    completely Hausdorff.
- **Remark 1.3 and §2.8 (Definition 2.62, Proposition 2.67).** Let O be a gauge-invariant observable of smooth
  connections, e.g. a Wilson loop f(hol(A, ℓ)). Then O_t = O ∘ F̃_t is defined on all of S and is continuous by (ii).
  - Holonomy observables separate O (Proposition 2.67, after Sengupta).
  - For many groups, Wilson loops alone do so (Remark 2.65).
- **The metric (Remark 1.5).** Σ = Θ + |||·|||, where Θ is built to make F̃_t continuous in the initial condition.
  The space is defined so that "flow for a positive time, then observe" is continuous.
- **Relation to Cao–Chatterjee (Remark 2.14).** Their space 𝒳 (S. Cao, S. Chatterjee, arXiv:2111.10652) consists of
  solutions of the same flow modulo gauge on (0, ∞). The CCHS orbit space embeds into it.
- **Construction.** CCHS build the Langevin dynamics on O, with local well-posedness and gauge covariance.
  Chevyrev–Shen, arXiv:2503.03060, settle uniqueness of the gauge-covariant renormalisation. The invariant measure,
  i.e. the 3D Yang–Mills measure itself, is not constructed.

## 2. Matching the paper's slice field

The paper's slice field is a gauge-invariant function, e.g. the magnetic energy density, of the spatial links of
one time slice x₀ = t. It is taken after a **spatial** gradient flow of that slice for time τ.

Regularity comes from Gaussian power counting, my estimate. A d-dimensional field with propagator 1/k² restricted to
a hyperplane has spectrum ∫dk₀/(k₀² + p²) = π/|p| in the D = d − 1 slice dimensions. That is Hölder–Besov regularity
(1 − D)/2 −.

| | full field | time slice | deterministic flow well-posed for |
|---|---|---|---|
| d = 3 | C^{−1/2−} on T³ (CCHS's η) | C^{−1/2−} on T² | ν > −½ |
| d = 4 | C^{−1−} on T⁴ | C^{−1−} on T³ | ν > −½ |

- **d = 3, full field, spacetime flow.** Any Wilson loop of the flowed field at t > 0 is a CCHS smoothened
  observable. It is continuous on O, so its law under any measure on O is well defined. Correlations between flowed
  observables in two different time slices fall in this class *if the flow is the full 3D flow*.
- **d = 3, spatial flow of a slice.** This is not of the CCHS form. The flow acts on a two-dimensional restriction,
  and restriction is not continuous on S (S ↪ C^{−1/2−} has no trace on a plane).
  - The restricted field has exactly the regularity C^{−1/2−} at which CCHS's construction works, but on T².
  - The known two-dimensional state spaces need regularity 0−: CCHS's 2D paper (Publ. IHES 2022, arXiv:2006.04987) and
    Chevyrev's 2D space (Comm. Math. Phys. 2019, arXiv:1808.09196). So they do not apply.
  - What is missing is a T² copy of (S, Σ) at η ≈ −½, together with a statement that the time-slice trace of the
    3D measure lives in it. The first part looks like a routine adaptation, since CCHS's estimates are heat-kernel
    estimates and the threshold ν > −½ does not depend on dimension. The second part needs the measure.
- **d = 4.** The slice is in C^{−1−}(T³), half a derivative below the threshold of the deterministic flow. The flowed
  slice field exists only as a renormalised random variable on the lattice, not as a continuous function on a fixed
  space of distributions. CCHS-type state spaces do not reach this regime; even the full 4D theory is beyond
  stochastic quantisation's subcritical methods, since d = 4 is critical.

## 3. Consequences for the programme

1. **No import is possible now.** Every statement about slice ratios is a statement about the measure, and the 3D
   measure on O is open. B3 gives a target formulation, not a shortcut.
2. **The natural 3D version of (I2) is in spacetime-flowed form.** Replace the paper's spatial flow by the full 3D
   flow at time τ. Then the slice correlators D(r) are correlations of CCHS smoothened observables. A future
   construction of the 3D measure on O, e.g. as the invariant measure of the CCHS dynamics, would make them defined
   and continuous automatically.
   - Check before switching: the spacetime flow smears the slice over a time width √τ. The transfer-matrix
     argument of Proposition "slice ratio" uses that the slice field depends on one time slice only (reflection
     positivity across x₀ = t ± ½). A flowed observable is no longer strictly localised, so (ii) of the proposition
     would need the version for observables with Gaussian tails in time. This is a real but probably manageable change.
3. **Lattice side.** Shen–Zhu–Zhu (arXiv:2204.12737) prove a Langevin mass gap and uniqueness for lattice
   Yang–Mills at strong coupling by Bakry–Émery. That is the same regime as our Dobrushin bounds, and gives nothing at
   weak coupling. Their method is the lattice analogue of the CCHS dynamics. It is the most direct route from our
   finite-β statements to the stochastic-quantisation literature.

**Status of B3:** literature check done. The clean statement is "spacetime-flowed slice observables are smoothened
observables in the sense of CCHS §2.8". Spatially flowed slices in d = 3 need a T² state space at regularity −½, and
in d = 4 no existing space applies.
