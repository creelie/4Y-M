/- Written by sharp.py: Arb enclosures rounded outward, denominator 10^12. -/

def den : Nat := 1000000000000
/-- lower bound for every eigenvalue of Hess g on |c| <= 24 beta -/
def lamLo : Nat := 156381938687
/-- upper bound for every eigenvalue of Hess g on |c| <= 24 beta -/
def tauHi : Nat := 250000000007
/-- upper bound for the norm of D^3 g on |c| <= 24 beta -/
def mHi : Nat := 61517136568
/-- upper bound for |grad g(c)| on |c| <= 20 beta -/
def aHi : Nat := 391442796709
