# Three dimensions as the proving ground for the slice ratio (direction B1)

Question: is the obstruction found in d = 4 at next order (a logarithm in the first-step energy of the flowed slice
state, `LEAKAGE.md` in big-bang-monotile/yang-mills-continuum) absent in d = 3? If so, three-dimensional SU(N)
Yang–Mills, which is superrenormalizable, is the natural place to try to prove the slice-ratio input (I2) first.

Answer from the data: **yes at next order, and the nonperturbative Monte Carlo agrees.** Nothing here is a proof.

## 1. Next order in d = 3: no logarithm, a finite correction

Setup (`leak.py`, `analyse_leak.py`): free transverse vacuum, exact nonlinear spatial Wilson flow, exact f-sum
rule for the mean energy of the flowed slice state, ⟨E⟩(g₀) = E0 + g₀²E2 + O(g₀⁴), SU(2), g₀² the lattice
coupling (β = 4/g₀²). Raw table: `data3d/leak_table.out`.

| √τ/a | 1 | 1.5 | 2 | 3 | 4 | 6 |
|---|---|---|---|---|---|---|
| E2/E0 | 0.321(35) | 0.800(64) | 1.240(82) | 1.85(21) | 2.36(24) | 4.03(71) |
| (E2/E0)/(√τ/a) | 0.32 | 0.53 | 0.62 | 0.62 | 0.59 | 0.67 |

For √τ/a ≥ 2 a linear fit E2/E0 = c·√τ/a + d gives c = 0.61(10), d = 0.02(24), χ² = 0.43/2; a logarithmic fit
is worse (χ² = 1.87/2). In d = 3 the lattice coupling is g₀² = g²a with g² the dimensionful continuum coupling, so the
physical correction is

    g₀² E2/E0 = 0.61(10) · g²√τ + O(a),

finite as a → 0. In d = 4 the same quantity grows like (N/3π²)·g²·ln(√τ/a). The reason is power counting: the
hard-pair leakage spectrum is ∝ 1/q³ in both dimensions (`LEAKAGE.md` §2), and it is integrated over a slice of
dimension d − 1: d³q/q³ is logarithmic, d²q/q³ converges in the ultraviolet and is dominated by q ≈ 1/√τ.

So in d = 3 the one-step bound (Proposition S of `SLICE.md`) survives next order without any appeal to
asymptotic freedom, and the split S′ (soft part plus small hard remainder) is not even needed at this order.

Not covered: order g₀⁴. Power counting gives at most (g²√τ)² times a possible ln(√τ/a) (two-loop logarithms exist in
three-dimensional gauge theory, in the vacuum energy), which would again be multiplied by a vanishing power only if it
attaches to g₀⁴ = g⁴a². This should be checked before relying on it.

## 2. Nonperturbative Monte Carlo, SU(2), d = 3

First-step energy E_k = −log(D(a)/D(0))/a of the spatially flowed slice field, in units of √σ, over a factor 5.3
in lattice spacing (`data3d/d3_summary.out`, `data3d/jax3d_summary.out`; L√σ ≥ 3.1):

| a√σ | 0.298 | 0.224 | 0.158 | 0.112 | 0.084 | 0.056 |
|---|---|---|---|---|---|---|
| τσ = 0.03 | 7.00(31) | 7.39(30) | 7.24(29) | 7.49(13) | 7.79(21) | 7.67(29) |
| τσ = 0.06 | 6.28(21) | 6.21(22) | 5.92(21) | 6.14(10) | 6.12(14) | 5.86(23) |
| τσ = 0.1 | 6.04(18) | 5.66(19) | 5.43(19) | 5.60(11) | 5.39(12) | 5.09(20) |
| unflowed | 10.3(7) | 12.8(1.0) | 18.4(1.6) | 28.9(2.2) | 41.1(2.3) | 56.4(4.5) |
| unflowed × a√σ | 3.07 | 2.85 | 2.91 | 3.23 | 3.45 | 3.15 |

- **With the flow E_k is bounded.** Linear extrapolations in a give E_k(a → 0)/√σ = 7.85(20), 5.99(14), 5.15(13)
  at τσ = 0.03, 0.06, 0.1, with χ² ≤ 4/4. A logarithmic growth is excluded at τσ = 0.06 and 0.1 (fitted slopes
  −0.14(15) and −0.44(13) per e-fold, i.e. if anything decreasing); at τσ = 0.03, where the flow radius is below one
  lattice spacing on the coarsest lattices, the slope is 0.43(21).
- **Without the flow E_k ≈ 3.1/a**, as in the free magnetic field (`composite_free.out`: 3.92/a), so the flow is
  what makes the first step bounded, and the interacting theory behaves like the free one.
- The slice ratio at r > 0 without the flow is not resolvable at this statistics (the correlator is noise from
  two steps on; see `../IDEAS.md`).

## 3. What this suggests

- The d = 3 version of (I2) is the right next target: there is no next-order obstruction and the Monte Carlo
  shows the first-step bound directly. The companion repo `yang-mills-3d-forest-gauge` provides the forest-gauge
  framework in d = 3.
- The proof would need: (a) the first-step bound of the flowed slice field at weak coupling in d = 3 to all orders in
  a convergent expansion (superrenormalizability leaves finitely many divergent
  diagrams; whether they all drop out of this ratio is the first thing to check); (b) a link between the
  first-step bound and the slice ratio, which is Proposition "slice ratio" (ii) of the paper (log-convexity).
- A useful intermediate theorem: for d = 3 Yang–Mills on a fixed lattice with β = 4/(g²a), the first-step energy of
  the flowed slice field is bounded uniformly in a at every fixed order of perturbation theory. Order g² is §1.

Reproduce: the scripts are in big-bang-monotile, `yang-mills-continuum/` (commit 0ef7a65):
`python3 analyse_leak.py data/leak/*.npy` and `python3 analyse.py`, `python3 analyse_jax.py` for the Monte Carlo summaries.
