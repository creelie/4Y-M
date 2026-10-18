"""Per-pair form of (2.14) of the OpenAI preprint, by quadrature (numerical evidence).
Correct form:  D^4W(u_i a, conj u_i b, conj u_j a, u_j b) = R(S_i,S_j) + R(S_ij, conj S_ij),
with S_i = S(f_i, conj f_i).  The preprint prints -R(S_i,S_j) + R(S_ij,S_ij); since S_j is purely
imaginary, R(S_i,S_j) = -R(S_i, conj S_j), and the summed identity (2.15) is the same.
Run: python3 -I check_pair.py   (from this directory)"""
import sys, os, itertools, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import check_identity as ci
rng = np.random.default_rng(5)
GEN = ci.GEN
pts, wts = ci.sphere_nodes(20, 24)
nrm = lambda v: np.asarray(v) / np.linalg.norm(v)
pinned = {2: nrm([0.3, -0.5, 0.81]), 3: nrm([-0.7, 0.1, 0.7])}
edges = [(2, 0), (0, 1), (1, 3), (2, 1), (0, 3)]
b = rng.uniform(0.2, 1.0, len(edges))
I, Jg = [m.ravel() for m in np.meshgrid(np.arange(len(wts)), np.arange(len(wts)), indexing='ij')]
W = wts[I] * wts[Jg]
S = {0: pts[I], 1: pts[Jg], 2: np.broadcast_to(pinned[2], (len(W), 3)), 3: np.broadcast_to(pinned[3], (len(W), 3))}
en = sum(be * np.einsum('ni,ni->n', S[x], S[y]) for be, (x, y) in zip(b, edges))
g = W * np.exp(en); g /= g.sum()
def J(gs):
    M = np.zeros((3, 3))
    for p in itertools.permutations(range(len(gs))):
        P = np.eye(3)
        for q in p: P = P @ GEN[gs[q]]
        M += P
    M /= math.factorial(len(gs))
    return np.stack([be * np.einsum('ni,ij,nj->n', S[x], M, S[y]) for be, (x, y) in zip(b, edges)], -1)
def Hk(args):
    pr = np.ones(len(edges), complex)
    for v, _ in args: pr = pr * v
    return J([gg for _, gg in args]) @ pr
E = lambda z: np.sum(g * z)
R = lambda h, k: -E(Hk([(h, 'c'), (k, 'c')]) + Hk([(h, 'c')]) * Hk([(k, 'c')]))
def parts(s):
    if not s: yield []; return
    f, r = s[0], s[1:]
    for p in parts(r):
        for i in range(len(p)): yield p[:i] + [[f] + p[i]] + p[i + 1:]
        yield [[f]] + p
def D4W(args):
    t = 0
    for pa in parts([0, 1, 2, 3]):
        term = np.ones(len(g), complex)
        for bl in pa: term = term * Hk([args[r] for r in bl])
        t = t + term
    return E(t)
d = lambda f: np.array([f[y] - f[x] for x, y in edges])
Sf = lambda f, h: np.array([0.5 * (f[x] * h[y] - h[x] * f[y]) for x, y in edges])
fi = np.array([1 + 2j, -0.5 + 0.3j, 0, 0]); fj = np.array([-0.4 + 1j, 0.8 - 0.6j, 0, 0])
ui, uj = d(fi), d(fj); Si = Sf(fi, fi.conj()); Sj = Sf(fj, fj.conj()); Sij = Sf(fi, fj)
print("D4W                             =", D4W([(ui, 'a'), (ui.conj(), 'b'), (uj.conj(), 'a'), (uj, 'b')]))
print("R(Si,Sj) + R(Sij,conj Sij)      =", R(Si, Sj) + R(Sij, Sij.conj()))
print("-R(Si,Sj) + R(Sij,conj Sij)     =", -R(Si, Sj) + R(Sij, Sij.conj()), "(as printed; differs)")
