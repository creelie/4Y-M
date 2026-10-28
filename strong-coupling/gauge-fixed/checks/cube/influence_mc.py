"""Monte Carlo probe of the cube condition D(n, alpha) for SU(2) in the forest gauge.

The coefficient C_{V,j} of the paper is a supremum over all exteriors; no
simulation can bound it. This script estimates the same quantity for typical
exteriors: an exterior is drawn by heat-bath sampling of the gauge-fixed
measure in a box with trivial outer boundary, the link j next to the cube is
rotated, and two heat-bath chains on the cube, one for each exterior, are run
with common random numbers. The time average of d_V between the two chains is
an upper bound for the Kantorovich distance of the two cube measures (up to
Monte Carlo error), so the ratio to delta(Q_j, Q'_j) estimates the influence
of j from above for that exterior. Summing over the links of a class gives an
estimate of K_c / N_n, the quantity that D(n, alpha) bounds by alpha.

Run: python3 influence_mc.py --beta 2.0 --n 4 [--samples 64 --sweeps 600]
Prints one line per run; nothing in the paper rests on it.
"""
import argparse
import math
import numpy as np

TABLE = """0000 0
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
SIGMA = np.zeros((2, 2, 2, 2), int)  # residue modulo two -> direction of the line
for _row in TABLE.splitlines():
    _r, _s = _row.split()
    SIGMA[tuple(int(c) for c in _r)] = int(_s)

# ---------------------------------------------------------------- quaternions
def qmul(a, b):
    a0, a1, a2, a3 = a[..., 0], a[..., 1], a[..., 2], a[..., 3]
    b0, b1, b2, b3 = b[..., 0], b[..., 1], b[..., 2], b[..., 3]
    return np.stack([a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
                     a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
                     a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1,
                     a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0], axis=-1)

def qconj(a):
    return a * np.array([1.0, -1.0, -1.0, -1.0])

# --------------------------------------------------------------------- lattice
class Box:
    """Links U[..., mu, x0, x1, x2, x3, :] on a box of side L (leading batch axis).

    Links with base point on the outer layer (index 0 or L-1) are the
    identity and never updated; links of the forest F are the identity.
    """

    def __init__(self, L, batch, rng):
        self.L, self.B, self.rng = L, batch, rng
        self.U = np.zeros((batch, 4, L, L, L, L, 4))
        self.U[..., 0] = 1.0
        idx = np.indices((L, L, L, L))
        res = idx % 2
        sig = SIGMA[res[0], res[1], res[2], res[3]]
        inner = np.ones((L,) * 4, bool)
        for d in range(4):
            inner &= (idx[d] >= 1) & (idx[d] <= L - 2)
        self.active = np.stack([(sig != mu) & inner for mu in range(4)])
        par = idx.sum(0) % 2
        self.colour = [[self.active[mu] & (par == p) for p in (0, 1)]
                       for mu in range(4)]

    def staple(self, mu):
        U, S = self.U, np.zeros_like(self.U[:, 0])
        sh = lambda a, d, s: np.roll(a, -s, axis=1 + d)  # a(x + s e_d)
        for nu in range(4):
            if nu == mu:
                continue
            Umu, Unu = U[:, mu], U[:, nu]
            fwd = qmul(qmul(sh(Unu, mu, 1), qconj(sh(Umu, nu, 1))), qconj(Unu))
            Unu_m = sh(Unu, nu, -1)
            bwd = qmul(qmul(qconj(sh(Unu_m, mu, 1)), qconj(sh(Umu, nu, -1))), Unu_m)
            S += fwd + bwd
        return S

    def heatbath(self, beta, rng, group=1, mask_extra=None):
        """One sweep. Rows of the batch are grouped in consecutive blocks of
        size `group`; rows of one block share their random numbers, which
        couples them (the chains of a pair). mask_extra restricts the update."""
        for mu in range(4):
            for p in (0, 1):
                m = self.colour[mu][p]
                if mask_extra is not None:
                    m = m & mask_extra[mu]
                if not m.any():
                    continue
                S = self.staple(mu)[:, m]          # (B, k, 4)
                k = np.linalg.norm(S, axis=-1)
                # density of the new link u: exp(beta <u, conj(S)>)
                mdir = qconj(S / k[..., None])
                self.U[:, mu][:, m] = sample_coupled(beta * k, mdir, group, rng)

def kp_x0(c, rng):
    """Kennedy-Pendleton: x0 with density sqrt(1-x0^2) exp(c x0), c > 0."""
    x0 = np.empty(c.shape)
    todo = np.ones(c.shape, bool)
    while todo.any():
        cc = c[todo]
        r = rng.random((4, cc.size))
        lam2 = -(np.log(r[0]) + np.cos(2 * np.pi * r[1]) ** 2 * np.log(r[2])) / (2 * cc)
        ok = r[3] ** 2 <= 1 - lam2
        idx = np.flatnonzero(todo)
        x0[idx[ok]] = 1 - 2 * lam2[ok]
        todo[idx[ok]] = False
    return x0

def vmf(c, mbar, rng):
    """Unit quaternions with density proportional to exp(c <u, mbar>)."""
    shape = c.shape
    x0 = np.empty(shape)
    small = c < 0.3
    g = rng.standard_normal(shape + (4,))
    g /= np.linalg.norm(g, axis=-1, keepdims=True)
    x0[small] = g[small][:, 0]
    if (~small).any():
        x0[~small] = kp_x0(c[~small], rng)
    v = g[..., 1:]
    v = v / np.linalg.norm(v, axis=-1, keepdims=True)
    rr = np.sqrt(np.maximum(0.0, 1 - x0 ** 2))
    x = np.concatenate([x0[..., None], rr[..., None] * v], axis=-1)
    x[small] = g[small]
    return qmul(x, mbar)

def sample_coupled(a, m, group, rng, max_iter=3000):
    """Rows in blocks of `group` sample u with density exp(a <u, m>) by
    rejection from a common proposal exp(c <u, mbar>): identical rows give
    identical samples, nearby rows agree with high probability."""
    B, k = a.shape
    G = B // group
    a3 = a.reshape(G, group, k)
    m3 = m.reshape(G, group, k, 4)
    msum = m3.sum(axis=1)
    nrm = np.linalg.norm(msum, axis=-1)
    mbar = np.where(nrm[..., None] > 1e-9, msum / np.maximum(nrm, 1e-300)[..., None], m3[:, 0])
    cosang = np.einsum("gjkc,gkc->gjk", m3, mbar)
    c = np.maximum(0.0, (a3 * np.maximum(cosang, 0.0)).min(axis=1))   # (G, k)
    drift = a3[..., None] * m3 - c[:, None, :, None] * mbar[:, None]
    M = np.linalg.norm(drift, axis=-1)                                  # (G, j, k)
    out = np.empty((G, group, k, 4))
    done = np.zeros((G, group, k), bool)
    for _ in range(max_iter):
        open_gk = ~done.all(axis=1)
        if not open_gk.any():
            break
        gi, ki = np.nonzero(open_gk)
        u = vmf(c[gi, ki], mbar[gi, ki], rng)                         # (q, 4)
        w = rng.random(gi.size)
        expo = np.einsum("qjc,qc->qj", drift[gi, :, ki], u) - M[gi, :, ki]
        acc = (np.log(w)[:, None] <= expo) & ~done[gi, :, ki]
        qi, ji = np.nonzero(acc)
        out[gi[qi], ji, ki[qi]] = u[qi]
        done[gi[qi], ji, ki[qi]] = True
    if not done.all():   # exact independent fallback keeps the marginals exact
        gi, ji, ki = np.nonzero(~done)
        out[gi, ji, ki] = vmf(a3[gi, ji, ki], m3[gi, ji, ki], rng)
    return out.reshape(B, k, 4)

# ------------------------------------------------------------------- the probe
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta", type=float, required=True, help="beta_W")
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--margin", type=int, default=3)
    ap.add_argument("--therm", type=int, default=400)
    ap.add_argument("--sweeps", type=int, default=600)
    ap.add_argument("--samples", type=int, default=32)
    ap.add_argument("--angle", type=float, default=math.pi / 2)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    n, w = a.n, a.margin
    L = n + 2 * w + 2
    if (w + 1) % 2:
        raise SystemExit("margin+1 must be even so that the cube starts on an even site")
    c0 = w + 1

    # thermalise the exterior
    box = Box(L, 1, rng)
    for _ in range(a.therm):
        box.heatbath(a.beta, rng)
    ext = box.U[0].copy()

    # the pair of chains lives on the cube with two layers around it
    off = c0 - 2
    P = n + 4
    ext = ext[:, off:off + P, off:off + P, off:off + P, off:off + P]
    idx = np.indices((P,) * 4)
    incube = np.ones((P,) * 4, bool)
    for d in range(4):
        incube &= (idx[d] >= 2) & (idx[d] < 2 + n)
    act = np.stack([SIGMA[tuple(idx % 2)] != mu for mu in range(4)])
    V = act & incube[None]

    # boundary links j = (y, mu): active, not in V, sharing a plaquette with V
    def at(a, shifts):  # a(y + sum s e_d)
        for d, s in shifts:
            a = np.roll(a, -s, axis=d)
        return a
    near = np.zeros_like(V)
    for mu in range(4):
        for nu in range(4):
            if nu == mu:
                continue
            near[mu] |= at(V[nu], [(mu, 1)]) | V[nu] | at(V[nu], [(mu, 1), (nu, -1)]) \
                | at(V[nu], [(nu, -1)]) | at(V[mu], [(nu, 1)]) | at(V[mu], [(nu, -1)])
    J = act & ~V & near
    jlist = np.argwhere(J)
    pick = jlist[rng.choice(len(jlist), size=min(a.samples, len(jlist)), replace=False)]

    B = 2 * len(pick)
    pair = Box(P, B, rng)
    pair.U[:] = ext[None]
    rot = np.array([math.cos(a.angle / 2), math.sin(a.angle / 2), 0.0, 0.0])
    dj = 2 * math.sin(a.angle / 2)
    for b, (mu, *x) in enumerate(pick):
        pair.U[2 * b + 1, mu, x[0], x[1], x[2], x[3]] = qmul(
            pair.U[2 * b + 1, mu, x[0], x[1], x[2], x[3]], rot)
    dist = []
    for t in range(a.sweeps):
        pair.heatbath(a.beta, rng, group=2, mask_extra=V)
        if t >= a.sweeps // 2:
            d = np.linalg.norm(pair.U[0::2] - pair.U[1::2], axis=-1)  # (P,4,L,L,L,L)
            dist.append((d * V[None]).sum(axis=(1, 2, 3, 4, 5)))
    dist = np.array(dist)
    half = len(dist) // 2
    drift = (dist[half:].mean() - dist[:half].mean()) / dj
    infl = dist.mean(axis=0) / dj              # per sampled j
    counts = len(jlist)
    Nn = (n // 2) ** 4
    # sum over all boundary links of all classes divided by 48 N_n: the class
    # average of K_c / N_n
    est = infl.mean() * counts / (48 * Nn)
    print(f"beta_W={a.beta} n={n} |V|={int(V.sum())} boundary={counts} "
          f"mean influence={infl.mean():.3f} (sd {infl.std():.3f}, max {infl.max():.3f}) "
          f"avg K_c/N_n={est:.3f} (late-minus-early {drift:+.3f})", flush=True)

if __name__ == "__main__":
    main()
