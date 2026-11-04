"""Constants of lem:lambdatwo (the bound on Lambda_2 of prop:layercross, G = SU(2)).

Part 1 checks by enumeration the combinatorial facts used in the proof.
Part 2 evaluates the analytic constants (mpmath, 40 digits; printed values
are rounded up where they are upper bounds).
Part 3 evaluates Delta of cor:twoislands and the width that makes Delta <= 1.
"""
import itertools
from mpmath import mp, mpf, sqrt, sin, exp, log, ceil

mp.dps = 40
D = 4
E = [tuple(1 if i == m else 0 for i in range(D)) for m in range(D)]
ZERO = (0,) * D


def add(a, b): return tuple(x + y for x, y in zip(a, b))
def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def off(x): return tuple(c % 2 for c in x)
def blk(x): return tuple(c // 2 for c in x)


def is_tree(l):
    x, m = l
    o = off(x)
    return o[m] == 0 and all(o[n] == 0 for n in range(m + 1, D))


def is_path(l):
    x, m = l
    return off(x) == E[m]


def kind(l):
    if is_tree(l): return 'tree'
    if is_path(l): return 'path'
    return 'private' if off(l[0])[l[1]] == 1 else 'inside'


def path_link(k):
    y, m = k
    return (add(tuple(2 * c for c in y), E[m]), m)


def args(k):
    """non-path links of the fourteen paths of nonzero offset of the pair k"""
    y, m = k
    base = tuple(2 * c for c in y)
    res = []
    for o in itertools.product((0, 1), repeat=D):
        if o in (ZERO, E[m]): continue
        x = add(base, o)
        for l in ((x, m), (add(x, E[m]), m)):
            assert not is_path(l)
            if l not in res: res.append(l)
    return res


def mult(k, g):
    """number of paths of nonzero offset of k that contain g"""
    y, m = k
    base = tuple(2 * c for c in y)
    n = 0
    for o in itertools.product((0, 1), repeat=D):
        if o in (ZERO, E[m]): continue
        x = add(base, o)
        n += g in ((x, m), (add(x, E[m]), m))
    return n


def plaqs(l):
    x, m = l
    out = []
    for n in range(D):
        if n == m: continue
        mu, nu = min(m, n), max(m, n)
        for b in (x, sub(x, E[n])):
            out.append(frozenset([(b, mu), (add(b, E[mu]), nu), (add(b, E[nu]), mu), (b, nu)]))
    return out


def pairs_of(g):
    x, m = g
    y = blk(x)
    return [k for k in ((y, m), (sub(y, E[m]), m)) if g in args(k)]


def part1():
    print('Part 1: combinatorics (enumeration over all link types)')
    offs = list(itertools.product((0, 1), repeat=D))
    nonpath = [(o, m) for o in offs for m in range(D) if kind((o, m)) in ('private', 'inside')]
    print('  non-path active links per block:', len(nonpath),
          '(private', sum(kind(g) == 'private' for g in nonpath), ', inside',
          sum(kind(g) == 'inside' for g in nonpath), ')')
    for m in range(D):
        a = args((ZERO, m))
        assert len(a) == 21 and sum(kind(l) == 'private' for l in a) == 7
    print('  each pair: 21 arguments, 7 private non-path, 14 inside (tree ones included)')
    for g in nonpath:
        assert len(set(plaqs(g))) == 6
        K = pairs_of(g)
        if kind(g) == 'private':
            assert len(K) == 1 and mult(K[0], g) == 2
        else:
            assert len(K) == 2 and all(mult(k, g) == 1 for k in K)
            assert path_link(K[0])[1] == path_link(K[1])[1]
    print('  private: one pair, two paths; inside: two pairs (parallel path links), one path each')
    maxpl = 0
    for o in itertools.product(range(4), repeat=D):
        for m in range(D):
            for p in plaqs((o, m)):
                maxpl = max(maxpl, sum(is_path(l) for l in p))
    print('  largest number of path links on one plaquette:', maxpl)
    for m in range(D):
        pl = path_link((ZERO, m))
        two = 0
        for p in plaqs(pl):
            npath = sum(is_path(l) for l in p)
            nnon = sum(kind(l) in ('private', 'inside') for l in p)
            assert npath + nnon <= 4 and nnon <= 2
            two += npath == 2
        assert two == 3
    print('  each path link: six plaquettes, three of them with a second path link,'
          ' at most two non-path links on each')
    worst = 0
    for g in nonpath:
        tot = sum(sum(is_path(l) for l in p) for p in plaqs(g))
        worst = max(worst, tot)
    print('  largest number of path links on the six plaquettes of a non-path link:', worst)
    # exact Schur row sum of the Wilson part, all path links solved, gamma_P = 1
    gP, gI = 1.0, 0.5
    gam = {'private': gP, 'inside': gI, 'tree': 0.0}
    best = 0
    for o in itertools.product(range(2), repeat=D):
        for m in range(D):
            g = (o, m)
            if kind(g) not in ('private', 'inside'): continue
            K = pairs_of(g)
            P = set(plaqs(g))
            for k in K: P |= set(plaqs(path_link(k)))
            T = 0.0
            for p in P:
                A = 1.0 if g in p else 0.0
                S = 0.0
                for l in p:
                    if is_path(l):
                        k = (blk(sub(l[0], E[l[1]])), l[1])
                        if k in K: A += gam[kind(g)]
                        S += sum(gam[kind(h)] for h in args(k) if h != g)
                    elif l != g:
                        S += 0.0 if is_tree(l) else 1.0
                T += A * S
            best = max(best, T)
    print('  exact value of the Wilson row sum at gamma_P = 1 (not used in the proof): %.1f' % (best / 2))
    print()


def constants(r, kappa, label):
    r = mpf(r); kappa = mpf(kappa)
    s2 = sqrt(2)
    thZ = 16 * s2 * r
    m1 = (thZ / 2) / sin(thZ / 2)
    eZ = (m1 - 1) + m1 * thZ / 2
    thX = 14 * s2 * r
    eX = (1 - sin(thX / 2) / (thX / 2)) + thX / 2
    eA = 2 * s2 * r
    epsT = eZ + eA + eX
    m = 8 / (1 - 7 * epsT)
    gP = m1 / (1 - 7 * epsT)
    sig = mpf(1) / 2
    MG = -log(2 - exp(sig))
    B2 = MG / ((sig - 16 * r) / 2) ** 2
    B3 = MG / ((sig - 16 * r) / 3) ** 3
    X2 = m * B2 / 16 * (14 * gP ** 2 + 4 * gP + 1)
    C2 = exp(14 * r) * gP ** 2 + X2
    M1 = B2 * (14 * gP + 2) / 16
    M2 = B3 * (14 * gP ** 2 + 4 * gP + 1) / 16 + mpf(14) / 16 * B2 * X2
    L = 3 * m * M2 + 3 * m ** 2 * M1 ** 2
    wil = 9 + 27 * gP + 63 * gP ** 2
    kap = kappa / 2 * 240 * C2
    cW = wil + kap
    cpsi = 40 * L
    print('Part 2: constants,', label)
    for name, val in [('r', r), ('kappa', kappa), ('m_1', m1), ('e_Z', eZ), ('e_X', eX),
                      ('e_A', eA), ('eps_T', epsT), ('m = 8/(1-7eps_T)', m), ('gamma_P', gP),
                      ('M_G', MG), ('B_2', B2), ('B_3', B3), ('X_2', X2), ('C_2', C2),
                      ('M_1', M1), ('M_2', M2), ('L', L),
                      ('9+27g+63g^2', wil), ('kappa term', kap), ('c_W', cW), ('c_psi', cpsi)]:
        print('  %-18s %s' % (name, mp.nstr(val, 12)))
    print()
    return cW, cpsi


def nblocks(d):
    return (2 * d + 1) ** 4 - (2 * d - 1) ** 4


def part3(cW, cpsi, beta, rho, s, label):
    beta = mpf(beta); rho = mpf(rho)
    r = 43 * rho
    eta = 92 * rho
    aG = mpf('2.5e-4')
    nx = 45 * (nblocks(s + 1) + nblocks(s + 2)) + 17 * nblocks(s + 3)
    Lam2 = beta * cW + cpsi
    pre = 1920 / beta * Lam2 ** 2 * r ** 2 * nx
    wmin = 2 * s + 8 + log(pre) / (2 * aG)
    print('Part 3: Delta of cor:twoislands,', label)
    print('  beta_W =', mp.nstr(beta, 6), ' s =', s, ' rho =', mp.nstr(rho, 6),
          ' r = 43 rho =', mp.nstr(r, 6), 
          ' eta+30r =', mp.nstr(eta + 30 * r, 6))
    print('  single-block islands: N_{s+1},N_{s+2},N_{s+3} =',
          nblocks(s + 1), nblocks(s + 2), nblocks(s + 3), ' n_x = n_y =', nx)
    print('  Lambda_2 <=', mp.nstr(Lam2, 8), ' beta_W r^2 =', mp.nstr(beta * r ** 2, 8))
    print('  Delta <= %s * exp(-2 alpha_G (w-2s-8))' % mp.nstr(pre, 6))
    print('  ln(prefactor) =', mp.nstr(log(pre), 8))
    print('  Delta <= 1 for w >=', mp.nstr(wmin, 8), ' i.e. w >=', int(ceil(wmin)))
    print('  each factor e in beta_W/beta* adds', mp.nstr(1 / (2 * aG), 6), 'layers;',
          'islands of a and b blocks add 1000 ln(ab) layers at most')
    print()


if __name__ == '__main__':
    part1()
    constants('8e-5', '1e-8', 'general lemma: r <= 8e-5, kappa <= 1e-8')
    rcal = mpf('6.0e-13') / 62
    cW, cpsi = constants(rcal, '6.0e-13', 'hypothesis (a) of thm:layerbound: r <= 6.0e-13/62, kappa <= 6.0e-13')
    part3(cW, cpsi, '1.13e37', '2.17e-16', 1, 'at the threshold of cor:twoislands, s = 1')
    part3(cW, cpsi, '1.33e37', '2.17e-16', 5, 'at the threshold of cor:twoislands, s = 5')
