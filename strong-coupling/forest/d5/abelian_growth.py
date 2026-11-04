"""
Check of the obstruction mechanism: for an abelian field the composed cut-off
average multiplies a constant curvature by exactly 4 per step (Lemma torusstep
/ Theorem torusscales), so a field with no fine defect above theta has every
scale-j plaquette equal to 4^j theta.

Model: U(1), two directions (the two directions of the plane of the curvature;
in four dimensions the fields are taken independent of the two transverse
coordinates, the transverse tree links are trivial, and the sixteen paths of a
pair reduce to these four), periodic torus of side L=2^m, blocks of side 2,
tree gauge (tree links = 0), cut-off average
   c_{y,mu} = p_0 + (chi/4) sum_x wrap(p_x - p_0),   chi=1 for small angles,
with p_x the angle of the path of offset x. Coarse plaquette angles are gauge
invariant; the coarse field is put back in the tree gauge and the step repeated.

Checks:
 (A) constant curvature theta=2 pi n/L^2: every scale-j plaquette angle is 4^j theta
     (mod 2 pi), j up to log2(L)-2.
 (B) random small field: coarse plaquette angle = mean over the four offsets of the
     fine square holonomies (the statement of Lemma torusstep), one step, to 1e-12.
"""
import numpy as np

def wrap(a):
    return (a + np.pi) % (2*np.pi) - np.pi

def plaq(th):
    # th[mu, x0, x1]; plaquette at x in plane (0,1): th0(x)+th1(x+e0)-th0(x+e1)-th1(x)
    t0, t1 = th
    return wrap(t0 + np.roll(t1, -1, axis=0) - np.roll(t0, -1, axis=1) - t1)

def tree_gauge(th):
    """Gauge transform so that the tree links of every block vanish:
    tree links of block y: (2y, dir 0) and (2y+o0 e0, dir 1), o0 in {0,1}.
    th'_{x,mu} = th_{x,mu} + g(x) - g(x+e_mu)."""
    L = th.shape[1]
    g = np.zeros((L, L))
    t0, t1 = th
    for y0 in range(0, L, 2):
        for y1 in range(0, L, 2):
            g[y0, y1] = 0.0
            # link (2y,0): need t0[y0,y1] + g[y0,y1] - g[y0+1,y1] = 0
            g[y0+1, y1] = t0[y0, y1] + g[y0, y1]
            # links (2y+o0 e0, 1): t1[x] + g[x] - g[x+e1] = 0
            for o0 in (0, 1):
                g[y0+o0, y1+1] = t1[y0+o0, y1] + g[y0+o0, y1]
    new = np.empty_like(th)
    new[0] = wrap(t0 + g - np.roll(g, -1, axis=0))
    new[1] = wrap(t1 + g - np.roll(g, -1, axis=1))
    return new

def cutoff_average(th):
    """One step: returns coarse field c[mu, y0, y1] on the torus of side L/2."""
    L = th.shape[1]
    th = tree_gauge(th)
    t0, t1 = th
    assert np.allclose(t0[0::2, 0::2], 0, atol=1e-9)          # tree links dir 0
    assert np.allclose(t1[:, 0::2], 0, atol=1e-9)          # tree links dir 1
    c = np.zeros((2, L//2, L//2))
    for y0 in range(L//2):
        for y1 in range(L//2):
            # pair (y,0): paths through x = 2y+o, two links in direction 0
            paths = []
            for o0 in (0, 1):
                for o1 in (0, 1):
                    x0, x1 = 2*y0+o0, 2*y1+o1
                    paths.append(t0[x0 % L, x1 % L] + t0[(x0+1) % L, x1 % L])
            p0 = paths[0]
            dev = [wrap(p - p0) for p in paths]
            chi = 1.0 if max(abs(d) for d in dev) < 0.5 else 0.0  # cut-off (inactive here)
            c[0, y0, y1] = wrap(p0 + chi*np.mean(dev))
            paths = []
            for o0 in (0, 1):
                for o1 in (0, 1):
                    x0, x1 = 2*y0+o0, 2*y1+o1
                    paths.append(t1[x0 % L, x1 % L] + t1[x0 % L, (x1+1) % L])
            p0 = paths[0]
            dev = [wrap(p - p0) for p in paths]
            chi = 1.0 if max(abs(d) for d in dev) < 0.5 else 0.0
            c[1, y0, y1] = wrap(p0 + chi*np.mean(dev))
    return c

def square_means(th):
    """Mean over the four offsets of the holonomy of the fine square of side two at 2y+o."""
    P = plaq(th)
    L = th.shape[1]
    S = np.zeros((L//2, L//2))
    for y0 in range(L//2):
        for y1 in range(L//2):
            tot = 0.0
            for o0 in (0, 1):
                for o1 in (0, 1):
                    x0, x1 = 2*y0+o0, 2*y1+o1
                    tot += P[x0 % L, x1 % L] + P[(x0+1) % L, x1 % L] + P[x0 % L, (x1+1) % L] + P[(x0+1) % L, (x1+1) % L]
            S[y0, y1] = tot/4
    return S

if __name__ == '__main__':
    rng = np.random.default_rng(1)
    # (A) constant curvature
    L, n = 256, 1
    theta = 2*np.pi*n/L**2
    th = np.zeros((2, L, L))
    th[1] = theta*np.arange(L)[:, None]*np.ones((1, L))   # th1(x) = theta x0
    th[0][L-1, :] = -theta*L*np.arange(L)                 # transition function on the seam column
    th = wrap(th)
    print(f"(A) L={L}, constant plaquette angle theta={theta:.6e}")
    cur = th
    j = 0
    while cur.shape[1] >= 4 and 8*4**j*theta < 0.5:   # stop before the cut-off can activate
        P = plaq(cur)
        print(f"   scale j={j}: side {cur.shape[1]:3d}, plaquette angle min={P.min():.6e} max={P.max():.6e}, 4^j theta={wrap(4**j*theta):.6e}")
        assert np.allclose(P, wrap(4**j*theta), atol=1e-10)
        cur = cutoff_average(cur)
        j += 1
    # (B) random small field, one step: coarse plaquette = mean of square holonomies
    L = 32
    th = rng.normal(scale=1e-3, size=(2, L, L))
    c = cutoff_average(th)
    Pc = plaq(c)
    S = square_means(th)
    err = np.abs(wrap(Pc - S)).max()
    print(f"(B) L={L}, random field with link angles ~1e-3: max |coarse plaquette - mean of four squares| = {err:.3e}")
    assert err < 1e-12
    print("    growth of the maximal plaquette angle in one step (random field):",
          f"{np.abs(plaq(th)).max():.3e} -> {np.abs(Pc).max():.3e}")
    print("all checks passed")
