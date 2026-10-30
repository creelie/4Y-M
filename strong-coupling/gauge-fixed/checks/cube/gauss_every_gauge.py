"""Checks for prop:gaussevery: in the quadratic approximation the cube
condition D(n, alpha), alpha < 1, fails for every 2Z^4-invariant gauge given
by a forest F of links whose trees are infinite and cover the sites.

For the direction field of Table tab:sigma and for the axial field sigma = 0:
(A) the 2Z^4-periodic flat fields form a space of dimension 19, and those
    vanishing on F a space of dimension at least 3;
(B) for n = 2, 4 the cube law is nondegenerate, its mean given the exterior
    of a flat field psi vanishing on F is psi itself, and the class sums obey
    sum_c a_c K_c >= N_n sum_c a_c with a_c = |psi| on the class c.
"""
import itertools, sys
import numpy as np

TAB = """0000 0
1000 0
0100 2
1100 3
0010 3
1010 1
0110 2
1110 1
0001 1
1001 2
0101 1
1101 3
0011 3
1011 2
0111 0
1111 0"""
SIGMA_TABLE = {tuple(int(ch) for ch in r.split()[0]): int(r.split()[1]) for r in TAB.splitlines()}
E = [tuple(1 if i == m else 0 for i in range(4)) for m in range(4)]
add = lambda x, y: tuple(a + b for a, b in zip(x, y))
mod2 = lambda x: tuple(a % 2 for a in x)

passed = failed = 0
def check(name, ok):
    global passed, failed
    print(("PASS " if ok else "FAIL ") + name)
    if ok: passed += 1
    else: failed += 1

def plaquettes_at(x):
    for m in range(4):
        for v in range(m + 1, 4):
            # links with signs: (x,m)+, (x+e_m,v)+, (x+e_v,m)-, (x,v)-
            yield [((x, m), 1), ((add(x, E[m]), v), 1), ((add(x, E[v]), m), -1), ((x, v), -1)]

def run(label, sigma):
    inF = lambda l: sigma(mod2(l[0])) == l[1]
    # (A) periodic flat fields on the 2-torus cell: psi(x,mu) = c_mu + h(x+e_mu) - h(x)
    sites = list(itertools.product(range(2), repeat=4))
    links = [(x, m) for x in sites for m in range(4)]
    idx = {s: i for i, s in enumerate(sites)}
    M = np.zeros((64, 20))
    for r, (x, m) in enumerate(links):
        M[r, 16 + m] += 1
        M[r, idx[mod2(add(x, E[m]))]] += 1
        M[r, idx[x]] -= 1
    check(f"{label}: periodic flat fields have dimension 19", np.linalg.matrix_rank(M) == 19)
    Frows = [r for r, l in enumerate(links) if inF(l)]
    check(f"{label}: F has 16 classes", len(Frows) == 16)
    # flat fields vanishing on F: kernel of M[Frows] restricted to the image of M
    _, sv, Vt = np.linalg.svd(M[Frows])
    null = Vt[np.sum(sv > 1e-9):]
    sols = M @ null.T
    dim = np.linalg.matrix_rank(sols)
    check(f"{label}: flat periodic fields vanishing on F have dimension {dim} >= 3", dim >= 3)
    rng = np.random.default_rng(1)
    psi_cell = sols @ rng.standard_normal(sols.shape[1])
    psi_cell[np.abs(psi_cell) < 1e-12] = 0.0
    psi = lambda l: psi_cell[links.index((mod2(l[0]), l[1]))]
    # flatness on a sample of plaquettes
    flat = max(abs(sum(s * psi(l) for l, s in p)) for x in sites for p in plaquettes_at(x))
    check(f"{label}: psi is flat and vanishes on F", flat < 1e-10 and all(abs(psi_cell[r]) < 1e-12 for r in Frows))
    for n in (2, 4):
        V = [(x, m) for x in itertools.product(range(n), repeat=4) for m in range(4) if not inF((x, m))]
        vid = {l: i for i, l in enumerate(V)}
        Vset = set(V)
        # plaquettes meeting V: base points in [-1, n]^4
        P = []
        for x in itertools.product(range(-1, n + 1), repeat=4):
            for p in plaquettes_at(x):
                if any(l in Vset for l, _ in p):
                    P.append(p)
        ext = sorted({l for p in P for l, _ in p if l not in Vset and not inF(l)})
        eid = {l: i for i, l in enumerate(ext)}
        Av = np.zeros((len(P), len(V))); Ae = np.zeros((len(P), len(ext)))
        for r, p in enumerate(P):
            for l, s in p:
                if l in vid: Av[r, vid[l]] += s
                elif l in eid: Ae[r, eid[l]] += s
        H = Av.T @ Av
        lam = np.linalg.eigvalsh(H)[0]
        check(f"{label}, n={n}: cube form nondegenerate (least eigenvalue {lam:.4f})", lam > 1e-6)
        U = -np.linalg.solve(H, Av.T @ Ae)          # columns u_j
        pe = np.array([psi(l) for l in ext]); pv = np.array([psi(l) for l in V])
        check(f"{label}, n={n}: mean given psi outside equals psi", np.max(np.abs(U @ pe - pv)) < 1e-9)
        Cj = np.abs(U).sum(axis=0)                  # C_{V,j} = ||u_j||_1
        K = {}
        for l, c in zip(ext, Cj):
            key = (mod2(l[0]), l[1]); K[key] = K.get(key, 0.0) + c
        Nn = (n // 2) ** 4
        lhs = sum(abs(psi((k[0], k[1]))) * K.get(k, 0.0) for k in K)
        tot = sum(abs(psi(l)) for l in links if not inF(l))
        check(f"{label}, n={n}: sum_c a_c K_c = {lhs:.4f} >= N_n sum_c a_c = {Nn * tot:.4f}", lhs >= Nn * tot - 1e-9)
        cls = [(x, m) for x in sites for m in range(4) if not inF((x, m))]
        kmax = max(K.get(k, 0.0) for k in cls if abs(psi(k)) > 1e-12)
        check(f"{label}, n={n}: max_c K_c / N_n = {kmax / Nn:.4f} >= 1", kmax >= Nn - 1e-9)

run("forest gauge", lambda r: SIGMA_TABLE[r])
run("axial gauge", lambda r: 0)
print(f"{passed} checks passed, {failed} failed")
sys.exit(1 if failed else 0)
