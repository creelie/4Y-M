"""Exact check of the dual certificate for the staggered coarse field gamma_B.

Fine torus of side L (L divisible by 4), blocks B(y)=2y+{0,1}^4, tree links
as in Definition def:blockgauge, block average V of that definition.
gamma_B: (gamma)_{y,0}=s(y), (gamma)_{y,1}=-s(y), s(y)=(-1)^{y0+y1}, others 0.
phi*: gamma_{y,mu} on the eight links (x,mu), x in B(y), o_mu=1; zero elsewhere.
Checks, in integer arithmetic:
  (a) V(phi*) = gamma_B;
  (b) d^T d phi* = V^T(32 gamma_B) on every link (tree links included);
  (c) Q(phi*) = 64|T|, Q_1(gamma_B) = 16|T|, |gamma_B|^2 = 2|T|.
V is scaled by 16 to stay in integers: 16 V = sum over the 16 paths.
"""
import itertools, sys
import numpy as np

def run(L):
    n = L // 2
    sh = (L, L, L, L)
    X = np.indices(sh)                      # X[nu] = coordinate nu
    O = X % 2                               # offsets
    Y = X // 2                              # block coordinates
    s = np.where(((Y[0] + Y[1]) % 2) == 0, 1, -1)
    # coarse gamma on (n,n,n,n,4)
    Yc = np.indices((n, n, n, n))
    sc = np.where(((Yc[0] + Yc[1]) % 2) == 0, 1, -1)
    gam = np.zeros((4, n, n, n, n), dtype=np.int64)
    gam[0] = sc; gam[1] = -sc
    # phi* on fine links: phi[mu][x]
    phi = np.zeros((4,) + sh, dtype=np.int64)
    for mu in range(4):
        g = gam[mu][Y[0], Y[1], Y[2], Y[3]]
        phi[mu] = np.where(O[mu] == 1, g, 0)
    # tree links: o_mu=0 and o_nu=0 for nu>mu
    tree = np.zeros((4,) + sh, dtype=bool)
    for mu in range(4):
        t = (O[mu] == 0)
        for nu in range(mu + 1, 4):
            t &= (O[nu] == 0)
        tree[mu] = t
    assert not np.any(phi[tree]), "phi* not in tree gauge"
    roll = lambda a, nu, k: np.roll(a, -k, axis=nu)   # a(x + k e_nu)
    # (a) 16 V_{y,mu}(phi) = sum_{x in B(y)} (phi_{x,mu} + phi_{x+e_mu,mu})
    Vs = np.zeros((4, n, n, n, n), dtype=np.int64)
    for mu in range(4):
        path = phi[mu] + roll(phi[mu], mu, 1)
        Vs[mu] = path.reshape(n, 2, n, 2, n, 2, n, 2).sum(axis=(1, 3, 5, 7))
    ok_a = np.array_equal(Vs, 16 * gam)
    # curvature of fine plaquettes: (dphi)_p at x, plane (mu,nu), mu<nu
    Q = 0
    dTd = np.zeros_like(phi)
    for mu, nu in itertools.combinations(range(4), 2):
        w = phi[mu] + roll(phi[nu], mu, 1) - roll(phi[mu], nu, 1) - phi[nu]
        Q += int((w * w).sum())
        # adjoint: link (x,mu) +w(x); (x+e_mu,nu) +w(x); (x+e_nu,mu) -w(x); (x,nu) -w(x)
        dTd[mu] += w
        dTd[nu] += np.roll(w, 1, axis=mu)
        dTd[mu] -= np.roll(w, 1, axis=nu)
        dTd[nu] -= w
    # V^T alpha with alpha = 32 gamma: V_{y,mu} = (1/16) sum_x (phi_{x,mu}+phi_{x+e_mu,mu})
    # => (V^T alpha)_{(z,mu)} = (1/16)[alpha_{y(z),mu} + alpha_{y(z-e_mu),mu}]
    VT = np.zeros_like(phi)
    for mu in range(4):
        a = 32 * gam[mu][Y[0], Y[1], Y[2], Y[3]]          # alpha at block of z
        a_prev = np.roll(a, 1, axis=mu)                   # alpha at block of z-e_mu
        tot = a + a_prev
        assert np.all(tot % 16 == 0)
        VT[mu] = tot // 16
    ok_b = np.array_equal(dTd, VT)
    T = n ** 4
    Q1 = 0
    for mu, nu in itertools.combinations(range(4), 2):
        rc = lambda a, ax: np.roll(a, -1, axis=ax)
        w = gam[mu] + rc(gam[nu], mu) - rc(gam[mu], nu) - gam[nu]
        Q1 += int((w * w).sum())
    g2 = int((gam * gam).sum())
    print(f"L={L}: |T|={T}  V(phi*)=gamma_B: {ok_a}  d^T d phi* = V^T(32 gamma_B) on all links: {ok_b}")
    print(f"   Q(phi*)={Q} = 64|T|: {Q == 64*T}   Q_1(gamma_B)={Q1} = 16|T|: {Q1 == 16*T}   |gamma_B|^2={g2} = 2|T|: {g2 == 2*T}")
    return ok_a and ok_b and Q == 64*T and Q1 == 16*T and g2 == 2*T

allok = all(run(L) for L in (4, 8, 16, 32))
print("ALL CHECKS PASSED" if allok else "FAILURE")
