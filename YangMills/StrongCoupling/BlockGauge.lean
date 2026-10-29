/-!
# The block-rooted axial gauge on blocks of side two

Offsets in a block `2y + {0,1}⁴` are read as `o : Fin 16`, bit `μ` of `o` being `o_μ`.  A link
`(o, μ)` is a tree link if `o_μ = 0` and `o_ν = 0` for every `ν > μ`.  The theorems check the
finite counts of Lemma `lem:blocktree` and of the proof of Theorem `thm:blockgap` of the paper in
`strong-coupling/forest/`, and the twelve classes of Proposition `prop:invisible`.
-/

namespace YangMills.StrongCoupling.BlockGauge

def bit (o : Fin 16) (μ : Nat) : Nat := (o.val >>> μ) % 2

def flip (o : Fin 16) (μ : Nat) : Fin 16 := ⟨(o.val ^^^ (1 <<< μ)) % 16, Nat.mod_lt _ (by decide)⟩

def isTree (o : Fin 16) (μ : Fin 4) : Bool :=
  bit o μ.val == 0 && (List.range 4).all (fun ν => !(decide (μ.val < ν)) || bit o ν == 0)

/-- Fifteen tree links per block. -/
theorem tree_count :
    ((List.finRange 16).flatMap fun o => (List.finRange 4).filter fun μ => isTree o μ).length = 15 := by
  decide

/-- Seventeen active links inside a block and thirty-two between blocks. -/
theorem active_counts :
    ((List.finRange 16).flatMap fun o => (List.finRange 4).filter fun μ =>
        bit o μ.val == 0 && !isTree o μ).length = 17 ∧
    ((List.finRange 16).flatMap fun o => (List.finRange 4).filter fun μ =>
        bit o μ.val == 1).length = 32 := by
  decide

/-- The highest set bit of a nonzero offset. -/
def top (o : Fin 16) : Nat := if bit o 3 = 1 then 3 else if bit o 2 = 1 then 2 else if bit o 1 = 1 then 1 else 0

/-- Every offset other than the corner is joined by a tree link to an offset with a smaller value,
so the tree links connect the block and, being fifteen, form a spanning tree. -/
theorem tree_connects :
    ∀ o : Fin 16, o.val ≠ 0 →
      isTree (flip o (top o)) ⟨top o, by unfold top; split <;> (try split) <;> (try split) <;> decide⟩ = true ∧
      (flip o (top o)).val < o.val := by
  decide

/-- The chain of an active link inside a block: the directions above `μ` in which the offset is one. -/
def ones (o : Fin 16) (μ : Nat) : List Nat := (List.range 4).filter fun ν => μ < ν && bit o ν == 1

/-- The plaquettes `(x', μ, ν)` used by the chain of `(z, μ)`, with `x'` the offset after clearing
the bits of `z` from the top down to `ν`. -/
def clearFrom (o : Fin 16) (ν : Nat) : Fin 16 :=
  ⟨(o.val % (1 <<< ν)) % 16, Nat.mod_lt _ (by decide)⟩

def load (x : Fin 16) (μ ν : Nat) : Nat :=
  ((List.finRange 16).filter fun z =>
      bit z μ == 0 && (ones z μ).length > 0 && bit z ν == 1 && decide (μ < ν) &&
      clearFrom z ν == x).foldl (fun s z => s + (ones z μ).length) 0

/-- Chains have length at most three and every plaquette carries load at most eight. -/
theorem chain_bounds :
    (∀ z : Fin 16, ∀ μ : Fin 4, (ones z μ.val).length ≤ 3) ∧
    (∀ x : Fin 16, ∀ μ ν : Fin 4, load x μ.val ν.val ≤ 8) := by
  decide

def table : List Nat := [0, 0, 2, 3, 3, 1, 2, 1, 1, 2, 1, 3, 3, 2, 0, 0]

def sig (x : Fin 16) : Nat := table.getD x.val 0

def invisible (x : Fin 16) : Bool := x.val == 14 || x.val == 15

/-- Exactly twelve classes of active links of the line field have one end on an invisible line
(residues `0111`, `1111`), and none has both. -/
theorem invisible_classes :
    ((List.finRange 16).flatMap fun x => (List.finRange 4).filter fun μ =>
        sig x != μ.val && (invisible x != invisible (flip x μ.val))).length = 12 ∧
    (∀ x : Fin 16, ∀ μ : Fin 4, sig x ≠ μ.val → ¬(invisible x ∧ invisible (flip x μ.val))) := by
  decide

end YangMills.StrongCoupling.BlockGauge
