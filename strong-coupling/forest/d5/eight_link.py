"""Certified arithmetic for eightlink2.tex.

Every inequality below is checked in interval arithmetic (mpmath.iv) on a
subdivision of the stated range; a failure raises AssertionError.  Part D is
an exact finite enumeration.  Part E computes the island threshold of
thm:islandthreshold and certifies six digits of each root by the sign of
f - log q* at the two ends of a bracket, in interval arithmetic.

Notation (prop:fold, prop:insidefold): torus point g = e^{iq}; the far
non-path link has angle -1/2, the six near links 49/100; the moved inside
link u_b = e^{it} sits on the path of a near link b, so that the fourteen
path angles are -1/2 (twice), 49/100 (eleven times) and 49/100 + t (once).
"""
from mpmath import mp, iv, mpf, log, pi, e, findroot

mp.dps = 40
iv.dps = 40
X = lambda a, b: iv.mpf([a, b])
TH_FAR, TH_NEAR = iv.mpf('-0.5'), iv.mpf('0.49')
TWO_SQRT2 = 2 * iv.sqrt(2)


def pieces(a, b, n):
    a, b = mpf(a), mpf(b)
    h = (b - a) / n
    return [iv.mpf([a + i * h, a + (i + 1) * h]) for i in range(n)]


def chi_of(d):
    # chi = 2 - 2m, m = 2 sqrt2 sin(d/2), valid for m in [1/2, 1]
    return 2 - 2 * TWO_SQRT2 * iv.sin(d / 2)


def f_torus(q, t):
    """Angle of V on the torus for the configuration with u_b = e^{it}.
    The maximiser is the far link for q >= 0 (checked in A2) and the
    largest of the near paths for q <= -0.01 (checked in A1)."""
    if q.b <= mpf('-0.01'):
        d = TH_NEAR - q + (t if t.a >= 0 else 0)  # shifted path wins iff t > 0
    elif q.a >= 0:
        d = q - TH_FAR
    else:
        raise ValueError
    return q + chi_of(d) * (iv.mpf('4.88') + t - 14 * q) / 16


def P(x, k=10):
    """print an interval endpoint or a real with k significant digits"""
    try:
        return f"{float(x):.{k}g}"
    except ValueError:          # a proper interval: print its upper end
        return f"{float(x.b):.{k}g}"


def report(name, lo, hi):
    print(f"{name}: [{P(lo)}, {P(hi)}]")


# ---------------------------------------------------------------- Part A1
print("Part A1: near the point g_L*, q in QL = [-0.0182, -0.0146], |t| <= 1e-3")
QL = pieces('-0.0182', '-0.0146', 400)
T1 = X('-0.001', '0.001')
mn, mx = mpf(10), mpf(-10)
for q in QL:
    # the near paths are strictly farther from g than the far link,
    # for every value of the shift t in [-1e-3, 1e-3]
    m1 = (TH_NEAR - q + X('-0.001', '0')) - (iv.mpf('0.5') + q)
    assert m1.a > mpf('0.018'), q
    mn = min(mn, m1.a)
print("  min over QL, t of (0.49 + min(t,0) - q) - (0.5 + q) >=", P(mn, 6))
lo, hi = mpf(10), mpf(-10)
for q in QL:
    m = TWO_SQRT2 * iv.sin((TH_NEAR - q + T1) / 2)
    assert m.a > 0.5 and m.b < 1, q
    lo, hi = min(lo, m.a), max(hi, m.b)
report("  m on QL x T", lo, hi)
lo = mpf(10)
for q in QL:
    br = iv.sin(TH_NEAR + T1 / 2 - q)
    assert br.a > mpf('0.47'), q
    lo = min(lo, br.a)
print("  bracket sin(0.49 + t/2 - q) >=", P(lo, 6))
lo = mpf(10)
for q in QL:
    dm = iv.sqrt(2) * iv.cos((TH_NEAR - q + T1) / 2)
    assert dm.a > mpf('1.3'), q
    lo = min(lo, dm.a)
print("  d/ds of 2 sqrt2 sin((0.49 + s - q)/2) at s = 0 >=", P(lo, 6))

# the torus equation near q_L, for t of both signs
GAMMA_LO, GAMMA_HI = mpf('0.16978261'), mpf('0.16978263')
for tt in (X('-0.001', '0'), X('0', '0.001')):
    a = f_torus(iv.mpf('-0.0182'), tt)
    b = f_torus(iv.mpf('-0.0146'), tt)
    assert a.b < GAMMA_LO and b.a > GAMMA_HI
    # derivative of f in q: 1 + chi'(q)(4.88+t-14q)/16 - 14 chi/16,
    # chi' = 2 sqrt2 cos(d/2) (d = 0.49 + t_+ - q decreases in q)
    dlo, dhi = mpf(10), mpf(-10)
    for q in QL:
        d = TH_NEAR - q + (tt if tt.a >= 0 else 0)
        chi = chi_of(d)
        dchi = 2 * iv.sqrt(2) * iv.cos(d / 2)
        fp = 1 + dchi * (iv.mpf('4.88') + tt - 14 * q) / 16 - 14 * chi / 16
        assert fp.a > mpf('1.3')
        dlo, dhi = min(dlo, fp.a), max(dhi, fp.b)
    print(f"  t in [{P(tt.a, 3)},{P(tt.b, 3)}]: f(-0.0182) < gamma*, "
          f"f(-0.0146) > gamma*, f' in [{P(dlo, 5)}, {P(dhi, 5)}]")

# ---------------------------------------------------------------- Part A2
print("Part A2: near the point g*, q in [0.05, 0.1], |t| <= 1e-3")
QS = pieces('0.05', '0.1', 400)
mn = mpf(10)
clo, chi_hi = mpf(10), mpf(-10)
mlo, mhi = mpf(10), mpf(-10)
c2 = mpf(-10)
for q in QS:
    marg = (q - TH_FAR) - (TH_NEAR - q + abs(T1))
    assert marg.a > mpf('0.1'), q
    mn = min(mn, marg.a)
    m = TWO_SQRT2 * iv.sin((q - TH_FAR) / 2)
    assert m.a > 0.5 and m.b < 1
    mlo, mhi = min(mlo, m.a), max(mhi, m.b)
    ch = chi_of(q - TH_FAR)
    clo, chi_hi = min(clo, ch.a), max(chi_hi, ch.b)
    d2 = iv.sqrt(2) * iv.sin((q - TH_FAR) / 2)  # chi_0''
    c2 = max(c2, abs(d2).b)
print("  far margin (0.5 + q) - (0.49 + |t| - q) >=", P(mn, 6))
report("  m", mlo, mhi)
report("  chi_0", clo, chi_hi)
print("  |chi_0''| <=", P(c2, 6), "; |t chi_0''/16| <=", P(c2 * mpf('1e-3') / 16, 4))
for q0 in ('0.05', '0.1'):
    v = f_torus(iv.mpf(q0), T1)
    assert v.a > GAMMA_HI + mpf('9e-4')
    print(f"  f({q0}, t) >= {P(v.a, 10)} > gamma* for |t| <= 1e-3")
for q0, sgn in (('0.0741', -1), ('0.0742', 1)):
    q = iv.mpf(q0)
    d = q - TH_FAR
    fp = 1 - 14 * chi_of(d) / 16 + (-2 * iv.sqrt(2) * iv.cos(d / 2)) * (iv.mpf('4.88') - 14 * q) / 16
    assert (fp.b < 0) if sgn < 0 else (fp.a > 0)
    print(f"  f_0'({q0}) in [{P(fp.a, 6)}, {P(fp.b, 6)}]")

# ---------------------------------------------------------------- Part B
print("Part B: the Kantorovich gap")
dist = TWO_SQRT2 * iv.sin(X('0.0897', '0.0914') / 2)
print("  delta(g*, g_L*) >=", P(dist.a, 6))
assert dist.a > mpf('0.1268')
print("  0.1268 - 2*0.0034 =", P(mpf('0.1268') - 2 * mpf('0.0034'), 6), ">= 0.12")
print("  geodesic length of t -> e^{it}: sqrt2 |t|; sqrt2 =", P(iv.sqrt(2).b, 8))

# ---------------------------------------------------------------- Part C
print("Part C: constants of the chain law")
s2 = iv.sqrt(2)
a1 = (1 + 6 * iv.mpf('1.1986')) / s2
a2 = 6 * iv.mpf('0.5993') / s2
assert a1.b <= mpf('5.793') and a2.b <= mpf('2.543')
print("  (1 + 6*1.1986)/sqrt2 <=", P(a1.b, 8), "<= 5.793")
print("  6*0.5993/sqrt2      <=", P(a2.b, 8), "<= 2.543")
for name, v in (("sqrt2*5.793", s2 * iv.mpf('5.793')), ("sqrt2*64.94", s2 * iv.mpf('64.94')),
                ("sqrt2*2.543", s2 * iv.mpf('2.543')), ("sqrt2*32.47", s2 * iv.mpf('32.47'))):
    print(f"  {name} <= {P(v.b, 8)}")
assert (s2 * iv.mpf('5.793')).b < mpf('8.193') and (s2 * iv.mpf('64.94')).b < mpf('91.84')
assert (s2 * iv.mpf('2.543')).b < mpf('3.597') and (s2 * iv.mpf('32.47')).b < mpf('45.92')
for n in (24,):
    A = n * s2 * iv.mpf('5.793'); Bc = n * s2 * iv.mpf('64.94') + iv.mpf('1.1986')
    assert A.b < 197 and Bc.b < 2206
    A = n * s2 * iv.mpf('2.543'); Bc = n * s2 * iv.mpf('32.47') + iv.mpf('0.5993')
    assert A.b < 87 and Bc.b < 1103
print("  n = 24 reproduces 197 beta_W + 2206 and 87 beta_W + 1103 of thm:site")

# ---------------------------------------------------------------- Part D
print("Part D: exact enumeration of the chains on a line of blocks")


def chains(bad, L, p):
    """Chains of a set of bad blocks on the line 0..L-1 (indices j), with p
    active inside links per block.  Links: ('P', j, a) private links of the
    pair j (a = 0 the path link), ('U', j, a) inside links of block j."""
    js = sorted(bad)
    comps, cur = [], []
    for j in js:
        if cur and j - cur[-1] > 2:
            comps.append(cur); cur = []
        cur.append(j)
    if cur:
        comps.append(cur)
    out = []
    for c in comps:
        j0, j1 = c[0], c[-1]
        C = set()
        for j in c:
            C |= {('P', j, a) for a in range(8)}
        for j in range(j0, j1 + 2):
            C |= {('U', j, a) for a in range(1, p + 1)}
        for j in range(j0 - 1, j1 + 2):
            if j not in bad:
                C.add(('P', j, 0))
        out.append((j0, j1, len(c), C))
    return out


def constraint(j, p):
    return ({('P', j, a) for a in range(8)} | {('U', j, a) for a in range(1, p + 1)}
            | {('U', j + 1, a) for a in range(1, p + 1)})


L = 15
checked = 0
viol_s1 = 0
for p in (7, 6, 4):
    for mask in range(1 << (L - 6)):
        bad = {j + 3 for j in range(L - 6) if mask >> j & 1}   # bad blocks in 3..L-4
        if not bad:
            continue
        ch = chains(bad, L, p)
        allC = set()
        for (j0, j1, nb, C) in ch:
            # (i) disjointness
            assert not (allC & C)
            allC |= C
            # (iv) count
            assert len(C) == 8 * nb + p * (j1 - j0 + 2) + (j1 - j0 + 3 - nb)
            if p == 7:
                assert len(C) == 7 * nb + 8 * (j1 - j0) + 17
            # links of C are based within distance 1 of a bad block
            for (kind, j, a) in C:
                assert min(abs(j - b) for b in bad) <= 1
            # constraints: which pairs see C, and what else they see
            for j in range(L - 1):
                K = constraint(j, p)
                if not (K & C):
                    continue
                assert j0 - 1 <= j <= j1 + 1
                ext = K - C
                for (kind, jj, a) in ext:
                    if kind == 'P':      # non-path private link of a good pair
                        assert a != 0 and jj not in bad and jj == j
                    else:                # inside link of block j0-1 or j1+2
                        assert jj in (j0 - 1, j1 + 2)
                if j in bad:
                    assert ('P', j, 0) in C and not ext & {('P', j, a) for a in range(8)}
            for j in range(j0 - 1, j1 + 2):
                assert constraint(j, p) & C
            assert len(C) <= 23 * nb + 1
        # every bad pair's private links lie in exactly one chain
        for b in bad:
            assert all(('P', b, a) in allC for a in range(8))
        # windows: a block w with dist(w, bad) >= s+1; the window frees the
        # private and inside links of w and the path link of the pair w-1
        # entering it, and imposes the constraints of the pairs w and w-1
        for s in (1, 2, 3):
            for w in range(1, L - 1):
                if min(abs(w - b) for b in bad) < s + 1:
                    continue
                free = ({('P', w, a) for a in range(8)} | {('U', w, a) for a in range(1, p + 1)}
                        | {('P', w - 1, 0)})
                cons = constraint(w, p) | constraint(w - 1, p)
                hit = (free | cons) & allC
                if s >= 2:
                    assert not hit
                elif hit:
                    viol_s1 += 1
        checked += 1
print(f"  {checked} bad patterns (p = 7, 6, 4): chains disjoint, counts exact, "
      f"n <= 23 n_b + 1, exterior links of chain constraints as stated, "
      f"windows of collar s >= 2 "
      f"disjoint from chains")
print(f"  collar s = 1: {viol_s1} (pattern, window) incidences, e.g. the end path link "
      f"of a chain entering a window block at distance 2")
assert viol_s1 > 0

# ---------------------------------------------------------------- Part E
print("Part E: the island threshold of thm:islandthreshold")
M, H = mpf(11) ** 4, mpf(6) ** 4
C0 = mpf(1) / 2 + mpf(2) / 3 * log(24 * pi ** 3)
print("  C0 =", P(C0, 12))


def logq(s):
    D = (4 * s + 3) ** 4 - 1
    return -H * (log(2) + 2 * log(D)) - log(e * M * H)


def logq_iv(s):
    D = iv.mpf((4 * s + 3) ** 4 - 1)
    Hi, Mi = iv.mpf(H), iv.mpf(M)
    return -Hi * (iv.log(2) + 2 * iv.log(D)) - iv.log(iv.e * Mi * Hi)


def f_iv(b, rho):
    C0i = iv.mpf(1) / 2 + iv.mpf(2) / 3 * iv.log(24 * iv.pi ** 3)
    return 6 * (iv.log(b) + C0i) - b * rho ** 2 / 16


for s in (1, 2, 5):
    print(f"  s={s}: Delta={(4*s+3)**4-1}, log q* = {P(logq(s), 9)}")
for rho in ('1e-8', '1e-9'):
    for s in (1, 2, 5):
        r = mpf(rho)
        target = logq(s)
        f = lambda b: 6 * (log(b) + C0) - b * r ** 2 / 16 - target
        b = findroot(f, 16 * abs(target) / r ** 2 * mpf('1.01'))
        lo, hi = b * (1 - mpf('1e-7')), b * (1 + mpf('1e-7'))
        ri, ti = iv.mpf(rho), logq_iv(s)
        assert (f_iv(iv.mpf(lo), ri) - ti).a > 0 and (f_iv(iv.mpf(hi), ri) - ti).b < 0
        assert b > 96 / r ** 2
        print(f"  beta*({rho}, {s}) in [{P(lo, 9)}, {P(hi, 9)}]")
print("all checks passed")
