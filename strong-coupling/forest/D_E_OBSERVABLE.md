# D. One concrete (E) programme: observable, normalization, and every use of the gap

Status: **direction D deliverable, part 1** (`../CLAY_NEXT_DIRECTIONS.md`). Line references are to `forest/main.tex`, §`ssec:limitE`.
"MIX" means any use of the averaged-transfer hypothesis (Thm `avgtransfer`, Cor `avggap`): the transfer bound ‖𝕋_k|_{Ω^⊥}‖ ≤ e^{−m_* a_k}, with m_* = 16ϖ/(625e²ℓ). The rule of direction D is that no step of (E) may use MIX as a premise unless it also cites an infrared input that is proved independently of the gap.

## 1. The choice

- **Observable.** Φ = Σ_{p ∋ 0} (N − Re Tr U_p), the sum over the 24 plaquettes with a corner at 0, so r = 1. It is the example in §`ssec:limitE`. It is B_4-invariant, local, and has no flow.
- **Packets.** A_k(f) = Σ_x a_k⁴ f(a_k x)(Φ_x − c_k).
- **Normalization.** Z_k = 𝒩_k(ε)^{−1/2}, the two-point normalization of Thm `limitEtwo`. It is preferred to C_k(σ_0)^{−1/2} (Prop `limitEnorm`(ii)) because Thm `limitEtwo` proves the two-point moment bound for it, instead of assuming (M_s).
- **Targets.** The four targets (nonzero limit, uniform moments, rotations, convergence) are mapped to results in the table in §2.

## 2. Every MIX use

| Result | Uses MIX for | Can it run without MIX? |
|---|---|---|
| Lem `limitElattice`(iii) | the decay e^{−m_* a_k t} of reflected correlations | RP part (ii, positivity, monotonicity): **yes**. Decay bound: **no**. |
| Thm `limitEsep`(i)–(iii) (existence of S_n on separated tuples, OS axioms) | compactness along a subsequence uses only (M_s); the transfer bound enters the spectral statements | Existence and RP: **yes**, given (M_s). |
| Thm `limitEsep`(iv)–(v) (σ(H_μ) ⊂ {0} ∪ [m_*,∞), clustering K e^{−m_*λ|y|}) | the gap itself | **No.** This *is* the gap transferred to the continuum. |
| Prop `limitEnorm`(i)–(ii) (C_k(t) ≤ e^{−m_* t}, ‖Ψ‖² ≥ e^{m_*σ_0}, nontriviality) | **nontriviality** via ‖Ψ_0(f_0)‖² ≥ e^{m_*σ_0} | **No.** The nonzero limit is currently derived *from* the gap. |
| Prop `limitEnorm`(iii) | a *lattice-unit* gap m_lat implies triviality | Warning result, not a use. |
| Cor `limitEexact` | convergence of separated moments | Depends on avgtransfer; to be checked line by line (open). |
| Cor `limitEstrong` | gap ≥ ϖ/(2e²q_2) in physical units | **No.** |
| Prop `limitEsing` | none (|x−y|^{−8} incompatible with (M)) | **Yes.** It is a no-go result and needs no input. |
| Lem `limitEbox`, Thm `limitEtwo`(i)–(iv) (two-point moment bound for Z_k = 𝒩_k(ε)^{−1/2}) | γ_k ≤ e^{−m_* a_k} in the box reduction | **Yes (D1 done).** The paper states that (i)–(iii) use only translation/B_4 invariance and RP. Checked: Lem limitEbox uses ‖𝕋_k^{n_0}‖ ≤ 1 and (ii) uses γ_k^{⌊d/a_k⌋} ≤ 1, so γ_k enters only (i)'s decay factor and (iv). |
| Thm `limitEtwo`(v), Prop `limitEsd`(i),(iii) | spectral decay, σ(H^{(1)}) ⊂ [m_*,∞) | **No.** |
| Cor `limitEtwo`(a),(b) existence, RP, nontriviality ⇔ the (NE)-type ratio | the decay estimate and σ(H^{(1)}) ⊂ [m_*,∞) only | **Yes** for existence and for nontriviality. (SD) ⇔ (M_s) for n = 2, and (SD) ⇔ (UV) (Prop limitEsd(ii)) are gap-free statements. |

## 3. Reading

1. **Nontriviality.** In my first pass I flagged this as resting on MIX; that is **corrected** here. Prop `limitEnorm`(ii) does derive nontriviality from the gap, but that is for the C_k(σ_0) normalization. For the chosen Z_k = 𝒩_k(ε)^{−1/2}, Cor `limitEtwo`(b) shows: some S_2 ≠ 0 ⇔ limsup_k ω_k(ϑ(Ā_k(g))A_k(g))/𝒩_k(ε) > 0 for some g, and (NE) gives every limit ≠ 0. This is MIX-free. The D rule is therefore met for n = 2, provided (SD)/(UV) and (NE) are proved as independent infrared/UV inputs. Those two are now the honest open inputs, and the gap is not one of them.
2. **Uniform moments** for n = 2 are MIX-free, modulo (SD) at smaller separations (D1, done). For n ≥ 3, (M_s) is assumed, not proved.
3. **Rotations** come from B_4 plus one rotation of infinite order (Prop `rotations`(iv)). The clustering step that upgrades this uses MIX.
4. **Convergence** is through a subsequence by compactness, which is MIX-free given (M_s). Full-sequence convergence (Cor `limitEexact`) is still to be audited.

## 4. Next items
- **D1.** Done: no γ_k in the bounds.
- **D2.** Replaced by: prove (NE) and (SD)⇔(UV) for the 24-plaquette density. These are the remaining non-gap inputs for n = 2. (UV) is the natural target: the spectral weight of near-plane packets must sit at bounded physical energy.
- **D3.** Audit Cor `limitEexact`.
