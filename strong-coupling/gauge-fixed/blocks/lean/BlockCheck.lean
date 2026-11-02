/-
Kernel check of the chain of constants of the SU(2) plaquette blocks (Lean core, no Mathlib).

BlockData.lean, written by ../sharp.py, holds Arb enclosures rounded outward to 10^-12:
  lamLo <= every eigenvalue of Hess g,  tauHi >= every eigenvalue,  mHi >= |D^3 g|  on |c| <= 24 beta,
  aHi   >= |grad g(c)|                                                           on |c| <= 20 beta,
with g = log Z0, beta = beta_W / 4 and beta_W = 0.351 fixing the balls.

Three facts about real numbers enter as rational inequalities:
  e^w <= 1 / (1 - w)   for 0 <= w < 1     (E below stands for 1 / (1 - W)),
  sqrt(tau) <= q       because q^2 >= tau,
  1 / sqrt 3 <= r      because 3 r^2 >= 1.
Every quantity of the chain is increasing in each input bound, so the rational values bound the
true ones.  The chain is

  W = 8 beta^2 (tau - lambda_*),   S = 1/4 + (E - 1) E q r,
  C = 4 beta S aHi + (320/3) beta^3 mHi E,   alpha = 40 beta (1 + beta) S / (1 - C).
-/
import BlockData

/-- nonnegative rationals, numerator over a positive denominator -/
structure Q where
  n : Nat
  d : Nat

def Q.add (a b : Q) : Q := ⟨a.n * b.d + b.n * a.d, a.d * b.d⟩
def Q.mul (a b : Q) : Q := ⟨a.n * b.n, a.d * b.d⟩
/-- truncated difference, used only where `b ≤ a` has been checked -/
def Q.sub (a b : Q) : Q := ⟨a.n * b.d - b.n * a.d, a.d * b.d⟩
def Q.div (a b : Q) : Q := ⟨a.n * b.d, a.d * b.n⟩
def Q.le (a b : Q) : Bool := a.n * b.d ≤ b.n * a.d
def Q.lt (a b : Q) : Bool := a.n * b.d < b.n * a.d
def Q.nat (k : Nat) : Q := ⟨k, 1⟩

def lam : Q := ⟨lamLo, den⟩
def tau : Q := ⟨tauHi, den⟩
def mM : Q := ⟨mHi, den⟩
def aA : Q := ⟨aHi, den⟩
/-- q >= sqrt(tau): q^2 >= tau is checked below -/
def q : Q := ⟨500000000008, den⟩
/-- r >= 1/sqrt 3: 3 r^2 >= 1 is checked below -/
def r : Q := ⟨577350269190, den⟩
def one : Q := Q.nat 1

def sq_ok : Bool := tau.le (q.mul q)
def r_ok : Bool := one.le ((Q.nat 3).mul (r.mul r))
def lam_le_tau : Bool := lam.le tau

/-- the chain at beta = bW / 4; returns (W, C, alpha) -/
def chain (bW : Q) : Q × Q × Q :=
  let b := bW.div (Q.nat 4)
  let W := (Q.nat 8).mul (b.mul (b.mul (tau.sub lam)))
  let E := one.div (one.sub W)
  let S := (⟨1, 4⟩ : Q).add ((E.sub one).mul (E.mul (q.mul r)))
  let C := ((Q.nat 4).mul (b.mul (S.mul aA))).add
           ((⟨320, 3⟩ : Q).mul (b.mul (b.mul (b.mul (mM.mul E)))))
  let alpha := ((Q.nat 40).mul (b.mul ((one.add b).mul S))).div (one.sub C)
  (W, C, alpha)

/-- W < 1 and C < 1 (so the truncated differences are honest) and alpha <= bound -/
def holds (bW bound : Q) : Bool :=
  let (W, C, alpha) := chain bW
  W.lt one && C.lt one && alpha.le bound

theorem data_ok : sq_ok && r_ok && lam_le_tau = true := by decide
theorem alpha_third : holds ⟨1, 3⟩ ⟨9428, 10000⟩ = true := by decide
theorem alpha_345 : holds ⟨345, 1000⟩ ⟨9805, 10000⟩ = true := by decide
theorem alpha_35 : holds ⟨35, 100⟩ ⟨99671, 100000⟩ = true := by decide
theorem alpha_351 : holds ⟨351, 1000⟩ ⟨99997, 100000⟩ = true := by decide
/-- the same rational chain exceeds one at beta_W = 0.352 -/
theorem alpha_352_fails : holds ⟨352, 1000⟩ one = false := by decide
