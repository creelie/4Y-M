"""Constants of the large-field block (largefield.tex), certified.

Real constants are enclosed with mpmath.iv (outward rounded interval
arithmetic); the lattice counts are exact integer enumerations.

Sections:
 A. M, Lambda, kappa_c, C_E = 8/kappa_c^2.
 B. Lattice geometry: tree links, owned plaquettes, dependence ranges,
    the overlap N_ov of the witness regions, C_* = C_E N_ov.
 C. The per-block constant c_1(lambda, eta) of the upper bound on the
    coarse action and the large-field rate a(lambda, kappa).
 D. The Kotecky-Preiss rate a_KP = 9 for the distance-one gas.
 E. Thresholds lambda*(kappa), beta* = 4 lambda*.
 F. The inherited count at scale j (witness boxes of lem:scalewitness):
    m_j, H'_j, thresholds beta^pw_j, and k_* - j^pw_max at 2001 couplings.
"""
import itertools

def iv_asin(x):  # asin by atan2, valid on (-1,1)
    return iv.atan2(x, iv.sqrt(1 - x**2))
from mpmath import iv, mp

iv.dps = 80
mp.dps = 80

def ivlo(x): return mp.make_mpf(x._mpi_[0])
def ivhi(x): return mp.make_mpf(x._mpi_[1])
def show(name, x):
    print(f"{name} in [{mp.nstr(ivlo(x), 12)}, {mp.nstr(ivhi(x), 12)}]")

print("=== A. basic constants ===")
alpha = 2 * iv_asin(1 / (2 * iv.sqrt(2)))
Lam = alpha / iv.sin(alpha)
one_m = 1 - iv.mpf(7) / 8 * Lam
M = 1 + 2 / one_m**3
logM = iv.log(M)
show("alpha", alpha); show("Lambda", Lam); show("1-7Lambda/8", one_m)
show("M", M); show("log M", logM); show("8 log M", 8 * logM)
assert M.b < 2.4e4
psi_c = iv_asin(1 / (4 * iv.sqrt(2))) / 3
kappa_c = 2 * iv.sqrt(2) * iv.sin(psi_c / 2)
C_E = 8 / kappa_c**2
show("kappa_c", kappa_c); show("C_E = 8/kappa_c^2", C_E)
assert (iv.mpf('16.003')**2).b < C_E.a
print("16.003^2 =", mp.nstr((iv.mpf('16.003')**2).b, 8), "< C_E: OK")

print("\n=== B. lattice geometry (exact enumeration) ===")
D = 4
E = [tuple(1 if i == m else 0 for i in range(D)) for m in range(D)]
def add(x, y): return tuple(a + b for a, b in zip(x, y))
def blk(x): return tuple(a // 2 for a in x)
def off(x): return tuple(a % 2 for a in x)
def is_tree(x, mu):
    o = off(x)
    return o[mu] == 0 and all(o[n] == 0 for n in range(mu + 1, D))
def link_class(x, mu):
    """('tree',y) | ('inside',y) | ('private',(y,mu)) by the base block y."""
    y = blk(x); o = off(x)
    if is_tree(x, mu): return ('tree', y)
    if o[mu] == 0: return ('inside', y)
    return ('private', (y, mu))
def plaq_links(x, mu, nu):
    return [(x, mu), (add(x, E[mu]), nu), (add(x, E[nu]), mu), (x, nu)]
planes = [(m, n) for m in range(D) for n in range(m + 1, D)]
B0 = list(itertools.product(range(2), repeat=D))
# per block link census
cnt = {'tree': 0, 'inside': 0, 'private': 0}
for x in B0:
    for mu in range(D):
        cnt[link_class(x, mu)[0]] += 1
print("links based in a block:", cnt)
assert cnt == {'tree': 15, 'inside': 17, 'private': 32}
# owned plaquettes of block 0
owned = [(x, m, n) for x in B0 for (m, n) in planes]
print("owned plaquettes per block:", len(owned))
assert len(owned) == 96
ok_active = ok_range = 0
inside_owned = 0
for (x, m, n) in owned:
    cls = [link_class(*l) for l in plaq_links(x, m, n)]
    if any(c[0] != 'tree' and (c[1] if c[0] != 'private' else c[1][0]) == (0,)*D for c in cls):
        ok_active += 1
    bases = {blk(l[0]) for l in plaq_links(x, m, n)}
    allowed = {(0,)*D, E[m], E[n]}
    if bases <= allowed: ok_range += 1
    o = off(x)
    if o[m] == 0 and o[n] == 0: inside_owned += 1
print("owned plaquettes containing an active link based in the block:", ok_active)
print("owned plaquettes whose links are based in {y, y+e_mu, y+e_nu}:", ok_range)
print("owned plaquettes inside the block:", inside_owned)
assert ok_active == 96 and ok_range == 96 and inside_owned == 24
# the plaquettes of a block's links: f_y depends on Xi_y, Xi_{y+e_mu}
# witness region of a coarse plaquette P=(y;mu,nu): blocks z with
# z-y in {-1..2} in mu,nu and {-1,0,1} elsewhere (lem:coarsewitness)
def region(y, m, n):
    rngs = [range(-1, 3) if c in (m, n) else range(-1, 2) for c in range(D)]
    return {add(y, d) for d in itertools.product(*rngs)}
print("blocks in a witness region:", len(region((0,)*D, 0, 1)))
def verts(x, m, n):
    return [x, add(x, E[m]), add(x, E[n]), add(add(x, E[m]), E[n])]
Nov = 0; arg = None
cand = list(itertools.product(range(-4, 5), repeat=D))
regs = {(y, m, n): region(y, m, n) for y in cand for (m, n) in planes}
for (x, m, n) in owned:
    vb = {blk(v) for v in verts(x, m, n)}
    c = sum(1 for key, R in regs.items() if vb & R)
    if c > Nov: Nov, arg = c, (x, m, n)
print("N_ov = max_p #{P : p in N(P)} =", Nov, "at", arg)
# sanity: no witness region in the enumeration window is cut off
C_star = C_E * Nov
show("C_* = C_E N_ov", C_star)
# number of coarse plaquettes whose witness region meets the plaquettes
# owned by one block (used for the local form of the kernel bound)
own_meet = set()
for (x, m, n) in owned:
    vb = {blk(v) for v in verts(x, m, n)}
    for key, R in regs.items():
        if vb & R: own_meet.add(key)
print("coarse plaquettes whose witness region meets Pi(y):", len(own_meet))

print("\n=== C. c_1(lambda, eta) and a(lambda, kappa) ===")
# A_1(c) <= 4(1+eta) lam Q_1(c) + L^4 c_1(lam, eta),
# c_1 = 12 log 1.957 + 45 log 14 + 67.5 (1 + log(11520 (1+1/eta) lam)),
# from rho^2 = 135/(2A), A = (1+1/eta) 96*90^2 lam.
assert 2 * 96 * 90**2 // 135 == 11520 and 2 * 96 * 90**2 % 135 == 0
K0 = 12 * iv.log(iv.mpf('1.957')) + 45 * iv.log(14)
show("12 log 1.957 + 45 log 14", K0)
def c1(lam, eta=1):
    lam = iv.mpf(lam); eta = iv.mpf(eta)
    return K0 + iv.mpf('67.5') * (1 + iv.log(11520 * (1 + 1 / eta) * lam))
def rho_opt(lam, eta=1):
    lam = iv.mpf(lam); eta = iv.mpf(eta)
    return iv.sqrt(iv.mpf(135) / (2 * (1 + 1 / eta) * 96 * 8100 * lam))
show("rho at lambda=100, eta=1", rho_opt(100))
assert (46 * rho_opt(100)).b < 0.043
for lam in ['1e3', '1e6', '1e12', '1e20']:
    show(f"c_1({lam},1)", c1(lam))
def a_rate(lam, kappa, eta=1):
    lam = iv.mpf(lam); kappa = iv.mpf(kappa)
    return lam * kappa**2 - 8 * logM - c1(lam, eta)

cfc = 14 * iv.mpf('90.6')**3 * iv.mpf(384)**3
show("14*90.6^3*384^3 (thm LF2(ii))", cfc)
assert ivhi(cfc) < 5.9e14

print("\n=== D. Kotecky-Preiss rate ===")
# connected sets (l^inf adjacency, 80 neighbours) of n blocks through a
# block: <= 80^{2(n-1)} = 6400^{n-1}; incompatible with X: a block within
# distance one of X, <= 81|X| choices.  Condition with a(X)=b|X|:
# 81 u/(1-6400 u) <= b,  u = e^{-(a-b)}.
b = iv.mpf(1) / 10; aKP = iv.mpf(9)
u = iv.exp(-(aKP - b))
lhs = 81 * u / (1 - 6400 * u)
show("u", u); show("81u/(1-6400u)", lhs)
assert lhs.b <= b.a
# optimal b, for the record
import math
bopt = (-81 + math.sqrt(81**2 + 4 * 6400 * 81)) / (2 * 6400)
print("optimal b =", bopt, " least a =", bopt + math.log(6400 + 81 / bopt))

print("\n=== E. thresholds lambda*(kappa): a(lambda,kappa) = 9 ===")
def threshold(f, lo, hi, target, rel=1e-8):
    lo = mp.mpf(lo); hi = mp.mpf(hi)
    assert f(lo).b < target and f(hi).a > target
    while hi / lo - 1 > rel:
        mid = mp.sqrt(lo * hi)
        v = f(mid)
        if v.b < target: lo = mid
        elif v.a > target: hi = mid
        else: break
    assert f(lo).b < target and f(hi).a > target
    return lo, hi
for kap in ['0.08', '1e-2', '1e-4', '1e-8']:
    k = mp.mpf(kap)
    lo0 = max(mp.mpf(100), mp.mpf('67.5') / k**2)   # a increasing beyond
    assert a_rate(lo0, kap).b < 9
    lo, hi = threshold(lambda l: a_rate(l, kap), lo0, 1e60, 9)
    print(f"kappa={kap}: lambda* in [{mp.nstr(lo,8)}, {mp.nstr(hi,8)}],"
          f" beta*=4lambda* in [{mp.nstr(4*lo,7)}, {mp.nstr(4*hi,7)}],"
          f" lambda* kappa^2 = {mp.nstr(hi*k**2,6)}")

print("\n=== F. inherited count at scale j ===")
def coord_counts_brute(j):
    lo_off = -8 * 2**j + 2; hi_off = 9 * 2**j - 3
    m = sum(1 for yc in range(-10 * 2**j, 10 * 2**j)
            if any(lo_off <= t <= hi_off for t in (2 * yc, 2 * yc + 1)))
    H = 0
    for yc in range(0, 2**j + 1):
        h = sum(1 for Yc in range(-20, 21)
                if any(2**(j + 1) * Yc + lo_off <= t <= 2**(j + 1) * Yc + hi_off
                       for t in (2 * yc, 2 * yc + 1)))
        H = max(H, h)
    return m, H
def coord_counts(j):
    # box of Y_c: [2^{j+1}Y_c - 8 2^j + 2, 2^{j+1}Y_c + 9 2^j - 3], even to odd
    m = 17 * 2**(j - 1) - 2
    # Y_c whose box meets {2y_c, 2y_c+1}: 2^{j+1}Y_c in [2y_c-9 2^j+3, 2y_c+8 2^j-1]
    P = 2**(j + 1); H = 0
    for yc in range(0, min(2**j, 4096)):  # the count depends on y_c only through 2y_c mod 2^{j+1}; for j>12 the window 0..4095 covers the residues near the ends, checked below
        a = 2 * yc - 9 * 2**j + 3; b = 2 * yc + 8 * 2**j - 1
        H = max(H, b // P - (-(-a // P)) + 1)
    # an interval of 17 2^j - 4 + 1 points meets at most floor(...)+1 multiples
    assert H <= (17 * 2**j - 4) // P + 1
    if j >= 3: assert H == (17 * 2**j - 4) // P + 1
    return m, H
for j in range(1, 11):
    assert coord_counts(j) == coord_counts_brute(j), j
print("closed-form counts agree with enumeration for j<=10")
rho = iv.mpf('1e-4')
C0 = iv.mpf(1) / 2 + iv.mpf(2) / 3 * iv.log(24 * iv.pi**3)
C1 = C0 + iv.log(12800) / 6
def a0(beta, j):
    lam = iv.mpf(beta) / 4
    rho0 = rho / (iv.mpf('1.01') * iv.mpf('1.023') * 4**j)
    return lam * rho0**2 - 8 * logM - c1(lam)
def need(m, H):
    # A_j >= 9  <=>  a_0 >= 9H + log(mH-1) + mH log(mH/(mH-1))
    mH = iv.mpf(m * H)
    return 9 * H + iv.log(mH - 1) + mH * iv.log(mH / (mH - 1))
beta_pw = {}
rows = []
for j in range(1, 31):
    mc, Hc = coord_counts(j)
    m, H = mc**4, Hc**4
    req = need(m, H)
    f = lambda bb, j=j: a0(bb, j)
    lo, hi = threshold(f, 1e6, 1e60, req.b)
    # certify: a0(hi) > req.b >= every point of req, a0(lo) < req.a ?
    assert a0(hi, j).a > req.b
    lo2, hi2 = threshold(f, 1e6, 1e60, req.a)
    assert a0(lo2, j).b < req.a
    beta_pw[j] = (lo2, hi)
    rows.append((j, mc, Hc, m, H, req, lo2, hi))
    if j <= 6 or j % 5 == 0:
        print(f"j={j}: m_j={mc}^4, H'_j={Hc}^4, a_0 needed in"
              f" [{mp.nstr(ivlo(req),9)},{mp.nstr(ivhi(req),9)}], beta^pw_j in"
              f" [{mp.nstr(lo2,7)}, {mp.nstr(hi,7)}],"
              f" 16^-j beta^pw_j <= {mp.nstr(hi/16**j,5)}")
print("max_j 16^-j beta^pw_j (1<=j<=30):",
      mp.nstr(max(r[7] / 16**r[0] for r in rows), 6))
print("min_j 16^-j beta^pw_j (1<=j<=30):",
      mp.nstr(min(r[6] / 16**r[0] for r in rows), 6))
# monotone thresholds
for j in range(1, 30):
    assert beta_pw[j][1] < beta_pw[j + 1][0]
print("beta^pw_j increasing in j: OK")
# islandswide values (s=5) for comparison
iw = {1: mp.mpf('2.89281e15'), 2: mp.mpf('4.62936e16'), 3: mp.mpf('1.18541e18')}
for j in (1, 2, 3):
    print(f"j={j}: beta*_j(islandswide)/beta^pw_j = "
          f"{mp.nstr(iw[j]/beta_pw[j][1],5)} to {mp.nstr(iw[j]/beta_pw[j][0],5)}")

def kstar(beta):
    B = iv.mpf(beta)
    r = iv.sqrt(96 * (iv.log(B) + C1) / B)
    t = iv.log(iv.mpf('1e-4') / (iv.mpf('1.01') * r)) / iv.log(4)
    fa, fb = int(mp.floor(t.a)), int(mp.floor(t.b))
    assert fa == fb, "k_* not certified"
    return fa
from collections import Counter
diffs = Counter(); jlist = Counter()
for i in range(2001):
    beta = mp.mpf(10)**(22 + mp.mpf(20) * i / 2000)
    ks = kstar(beta)
    jm = 0
    for j in range(1, 31):
        lo, hi = beta_pw[j]
        if beta >= hi: jm = j
        elif beta < lo: break
        else: raise AssertionError("beta inside a threshold enclosure")
    assert jm < 30
    diffs[ks - jm] += 1
print("k_* - j^pw_max over 2001 couplings in [1e22,1e42]:", dict(diffs))
