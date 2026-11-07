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
| Lem `limitEbox`, Thm `limitEtwo`(i)–(iv) (two-point moment bound for Z_k = 𝒩_k(ε)^{−1/2}) | γ_k ≤ e^{−m_* a_k} in the box reduction | The bound K_L K_{4L/δ}² looks **combinatorial plus RP**. Whether γ_k is used for the bound or only for decay is to be checked (open item D1). |
| Thm `limitEtwo`(v), Cor `limitEtwo`, Prop `limitEsd` | spectral decay, σ(H^{(1)}) ⊂ [m_*,∞) | **No.** |

## 3. Reading

1. **A nonzero limit currently rests on MIX** (Prop `limitEnorm`(ii), and Cor `limitEtwo` "some S_2 not identically 0"). This is the one dependency that violates the D rule outright. An independent infrared input is needed, for example a lower bound on 𝒩_k(ε) relative to 𝒩_k(ε') that comes from RP and the spectral measure alone (Prop `limitEsd` with γ_k ≤ 1 instead of e^{−m_* a_k}).
2. **Uniform moments** for n = 2 may be MIX-free through Thm `limitEtwo`(i)–(iii) (item D1). For n ≥ 3, (M_s) is assumed, not proved.
3. **Rotations** come from B_4 plus one rotation of infinite order (Prop `rotations`(iv)). The clustering step that upgrades this uses MIX.
4. **Convergence** is through a subsequence by compactness, which is MIX-free given (M_s). Full-sequence convergence (Cor `limitEexact`) is still to be audited.

## 4. Next items
- **D1.** Audit the proofs of Lem `limitEbox` and Thm `limitEtwo`(i)–(iii) for γ_k. Decide whether the two-point bound for the normalization survives γ_k ≤ 1.
- **D2.** A MIX-free nontriviality argument for S_2, built from the RP spectral measure alone.
- **D3.** Audit Cor `limitEexact`.
