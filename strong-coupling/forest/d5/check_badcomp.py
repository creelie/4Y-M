"""Interval checks for the block badcomp.tex (mpmath iv, 60 digits)."""
from mpmath import iv
mpf = iv.mpf
iv.dps = 60
ok = True
def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name); ok &= bool(cond)
c3 = 15000*iv.sqrt(2)
kmax = 1/(mpf(195)*mpf('4e5')*c3)
print("c3 =", c3); print("kappa_max = 1/(195*4e5*c3) =", kmax)
check("6.04e-13 < kappa_max < 6.05e-13", kmax > mpf('6.04e-13') and kmax < mpf('6.05e-13'))
w = (4+28*iv.pi)*mpf('6.05e-13')
print("(4+28pi)*6.05e-13 =", w)
check("(4+28pi)*6.05e-13 < 5.57e-11", w < mpf('5.57e-11'))
check("6.05e-13 <= 1e-8 (Lemma coarsewitness applies)", mpf('6.05e-13') <= mpf('1e-8'))
# the hypothesis kappa+s <= 6.0e-13 of Prop. explicitconstants implies 4e5 c3 (kappa+s) <= 1/195
check("6.0e-13 <= kappa_max", mpf('6.0e-13') <= kmax)
eta = 92*mpf('2.17e-16'); print("eta = 92 rho at rho=2.17e-16:", eta)
check("eta < 2.0e-14", eta < mpf('2.0e-14'))
ratio = mpf('5.57e-11')/eta; print("5.57e-11/eta =", ratio)
check("2790 < 5.57e-11/eta < 2791", ratio > 2790 and ratio < 2791)
r2 = mpf('6.0e-13')/eta; print("6.0e-13/eta =", r2)
check("6.0e-13/eta > 30", r2 > 30)
check("4 c9 = 4.8e4 with c9=1.2e4", 4*mpf('1.2e4') == mpf('4.8e4'))
print("ALL PASS" if ok else "SOME FAIL")
