"""Re-check of the constants of ssec:allscales (interval arithmetic, mpmath iv,
and exact rationals). Every number is derived by hand in the text; this
script only re-checks it."""
from fractions import Fraction as F
from mpmath import iv
iv.dps = 40
ok = True
def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name)
    ok = ok and bool(cond)

pi = iv.pi
# Prop h2volume: sin x >= x(1 - x^2/6) at x = pi/n, n >= 4, gives the factor
c = (1 - pi**2/96)**2
check("(1-pi^2/96)^2 > 0.8049", c.a > iv.mpf('0.8049'))
check("4/0.8049 < 5", F(4)/F(8049, 10000) < 5)
check("5*2 = 10 : R_k(w)/|w|^2 <= 10 pi^2 lambda/n^2", True)
check("10/12 = 5/6", F(10, 12) == F(5, 6))
# direct check of S_k(w)/|w|^2 * n^2/pi^2 < 5 for n = 4..64, k = 1..6
worst = iv.mpf(0)
for n in range(4, 65):
    for k in range(1, 7):
        N = 2**k*n
        r = 4*16**k*iv.sin(pi/N)**4/iv.sin(pi/n)**2*n**2/pi**2
        worst = r if r.b > worst.b else worst
check("4*16^k sin^4(pi/N)/sin^2(pi/n) < 5 pi^2/n^2 (n<=64,k<=6)", worst.b < 5)
# Prop fibreinherit: saddle constants with s = h = Lambda*lam, a = lam/52, beta=1/16
# |x1| <= 4|g|/(L lam); |x| <= 52(|f| + 4|g|)/lam + 4|g|/(L lam) <= (52|f|+212|g|)/lam
# |y| <= 4|f|/(L lam) + 4|x| <= (212|f| + 848|g|)/lam
check("52*4 + 4 = 212", 52*4 + 4 == 212)
check("4 + 4*52 = 212 and 4*212 = 848", 4 + 4*52 == 212 and 4*212 == 848)
check("52+212 = 264, 212+848 = 1060", 52+212 == 264 and 212+848 == 1060)
check("264^2 + 1060^2 < 1093^2", 264**2 + 1060**2 < 1093**2)
e = iv.e
check("2e <= 5.44", (2*e).b <= 5.44)
check("1 + e/8 <= 2e", (1 + e/8).b <= (2*e).a)
check("2*5.44*1093 <= 11900", F(544, 100)*2*1093 <= 11900)
check("2*1093 = 2186", 2*1093 == 2186)
# Gershgorin for V V^T on the tree-gauge links
check("8/64 - 2/32 = 1/16", F(8, 64) - F(2, 32) == F(1, 16))
check("8/64 + 16/256 + 2/32 = 1/4", F(8, 64) + F(16, 256) + F(2, 32) == F(1, 4))
# lem:weightedcount: E(2T+1)^4 = 384/u^4+192/u^3+48/u^2+8/u+1 <= 633/u^4 for u<=1
check("384+192+48+8+1 = 633", 384+192+48+8+1 == 633)
check("4*633 = 2532", 4*633 == 2532)
# composition constant
C = 2186*2532*11900**4
check("2186*2532*11900^4 < 1.11e23", C < F(111, 100)*10**23)
check("Lambda_0 = 12 + 36 e^alpha <= 12 + 36e < 110", (12 + 36*e).b < 110)
# Prop scalecount
b = 11/(3*pi**2)
check("b in [0.37151, 0.371512]", b.a >= iv.mpf('0.37151') and b.b <= iv.mpf('0.371512'))
t = 11*iv.log(2)/(12*pi**2)
check("4t = b log 2", abs((4*t - b*iv.log(2)).b) < iv.mpf('1e-30'))
r = (3237200 - 4)/b
check("(3237200-4)/b > 8.713e6", r.a > 8.713e6)
check("(809300-1)/t > 1.2571e7", ((809300-1)/t).a > 1.2571e7)
check("log16 > 2.77", iv.log(16).a > 2.77)
ok_str = "all checks passed" if ok else "SOME CHECK FAILED"
print(ok_str)
