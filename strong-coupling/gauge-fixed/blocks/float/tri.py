"""Floating-point bound for the Kantorovich derivative of a plaquette-triple block of SU(2).

Block law on (S^3)^3: exp(<x1,b1>+<x2,b2>+<x3,b3>+4 beta Re(x1 x2 x3)).  Coordinates (x1, w, x2) with
w = conj(x2 x3); the coupling is 4 beta <x1, w>, <x3,b3> = <w, conj(x2 b3)> = <x2, conj(w) conj(b3)>.
For phi = <x1,h> and f 1-Lipschitz for the sum of chordal distances,
Cov(f, phi) <= 1/2 [ E Q_1(T1)^(1/2) + E Q_w(T_w)^(1/2) + 2 Q_2(T_2)^(1/2) ],
T1 = u - Pi_1 u, T_w = Pi_1 u - Pi_{1w} u, T_2 = Pi_{1w} u, u = (phi - E phi) d mu / d sigma^3.
Q(F) = sum_{k>=1} |F_k|^2 / k (harmonic components on S^3)."""
import numpy as np
from scipy.special import iv, roots_legendre

def qmul(p, q):
    a1, b1, c1, d1 = np.moveaxis(p, -1, 0); a2, b2, c2, d2 = np.moveaxis(q, -1, 0)
    return np.stack([a1*a2-b1*b2-c1*c2-d1*d2, a1*b2+b1*a2+c1*d2-d1*c2,
                     a1*c2-b1*d2+c1*a2+d1*b2, a1*d2+b1*c2-c1*b2+d1*a2], -1)
def qconj(p):
    return p * np.array([1, -1, -1, -1])

def Z0(k):
    k = np.asarray(k, float); out = np.ones_like(k); m = k > 1e-8
    out[m] = 2*iv(1, k[m])/k[m]; out[~m] = 1 + k[~m]**2/8; return out
def Z0p(k):  # d/dk Z0 = 2 I_2(k)/k
    k = np.asarray(k, float); out = np.zeros_like(k); m = k > 1e-8
    out[m] = 2*iv(2, k[m])/k[m]; out[~m] = k[~m]/4; return out

class Grid:
    """Product quadrature on S^3 exact for polynomials of degree <= D, and an orthonormal basis of
    the harmonics of degree 1..K (Hopf coordinates, Jacobi polynomials)."""
    def __init__(self, D, K):
        ns = D//4 + 2; nxi = D + 1
        s, ws = roots_legendre(ns); s = (s+1)/2; ws = ws/2
        xi = 2*np.pi*np.arange(nxi)/nxi
        S, X1, X2 = np.meshgrid(s, xi, xi, indexing='ij')
        W = np.meshgrid(ws, np.ones(nxi), np.ones(nxi), indexing='ij')[0] / nxi**2
        S, X1, X2, W = S.ravel(), X1.ravel(), X2.ravel(), W.ravel()
        c, sn = np.sqrt(1-S), np.sqrt(S)
        self.x = np.stack([c*np.cos(X1), c*np.sin(X1), sn*np.cos(X2), sn*np.sin(X2)], -1)
        self.w = W
        from scipy.special import eval_jacobi
        basis, deg = [], []
        for k in range(1, K+1):
            for m1 in range(-k, k+1):
                for m2 in range(-k, k+1):
                    r = k - abs(m1) - abs(m2)
                    if r < 0 or r % 2: continue
                    n = r//2
                    rad = c**abs(m1) * sn**abs(m2) * eval_jacobi(n, abs(m2), abs(m1), 1-2*S)
                    for f in ([np.cos(m1*X1+m2*X2)] + ([np.sin(m1*X1+m2*X2)] if (m1, m2) != (0, 0) else [])):
                        basis.append(rad*f); deg.append(k)
        B = np.array(basis); deg = np.array(deg)
        # orthonormalise within each degree (exact quadrature)
        out = []; od = []
        for k in range(1, K+1):
            Bk = B[deg == k]; G = (Bk*W) @ Bk.T
            e, V = np.linalg.eigh(G); keep = e > 1e-10*e.max()
            Bk = (V[:, keep]/np.sqrt(e[keep])).T @ Bk
            out.append(Bk); od += [k]*len(Bk)
        self.B = np.vstack(out); self.deg = np.array(od)
        self.BW = self.B*W
    def Q(self, F):
        """F: (..., N) values; returns sum_k |F_k|^2/k."""
        c = F @ self.BW.T
        return (c**2/self.deg).sum(-1)

def qser(c, h, m, nmax=40):
    """Q of (<x,h>-m) e^{<x,c>} on S^3, closed form by Bessel series; c (...,4)."""
    k = np.linalg.norm(c, axis=-1); k = np.maximum(k, 1e-12)
    cosphi = (c @ h)/k; sin2 = 1-cosphi**2
    n = np.arange(1, nmax+1)[:, None]
    I = iv(n+1, k[None]); Ip = 0.5*(iv(n, k[None])+iv(n+2, k[None]))
    e = 2*(n+1)*I/k; ep = 2*(n+1)*(Ip/k - I/k**2)
    t = ((cosphi*ep - m*e)**2 + sin2*e**2*n*(n+2)/(3*k**2))/n
    return t.sum(0)

def bound(b1, b2, b3, h, beta, g):
    x = g.x; w8 = g.w
    w = x                                   # w on the grid
    c1 = b1 + 4*beta*w                      # field on x1 given w
    k1 = np.linalg.norm(c1, axis=-1)
    v2 = b2 + qmul(qconj(w), qconj(b3)[None])   # field on x2 given w
    A = Z0(np.linalg.norm(v2, axis=-1))
    Z = np.sum(w8*A*Z0(k1))
    gh = Z0p(k1)*(c1 @ h)/np.maximum(k1, 1e-14)  # <h, grad Z0(c1)>
    m = np.sum(w8*A*gh)/Z
    t1 = 0.5*np.sum(w8*A*np.sqrt(qser(c1, h, m)))/Z
    gfun = gh - m*Z0(k1)                    # Pi_1 u = e^{<x2,b2>+<w,conj(x2 b3)>} g(c1)/Z
    # T_w: for each x2 on the grid, F(w) = e^{<w, conj(x2 b3)>} g(c1(w)); weight e^{<x2,b2>}/Z
    V = qconj(qmul(x, b3[None]))            # conj(x2 b3) for x2 on the grid
    E = np.exp(V @ w.T)                     # (x2, w)
    F = E*gfun[None]
    tw = 0.5*np.sum(w8*np.exp(x @ b2)*np.sqrt(np.maximum(g.Q(F), 0)))/Z
    # T_2(x2) = e^{<x2,b2>} G(x2)/Z, G(x2) = int dsigma(w) F(w)
    G = F @ w8
    t2 = 2*0.5*np.sqrt(max(g.Q(np.exp(x @ b2)*G/Z), 0))
    return t1+tw+t2, (t1, tw, t2, m)

if __name__ == "__main__":
    g = Grid(20, 8)
    print("points", len(g.w), "basis", len(g.deg))
    # sanity: Q(<x,h>) = 1/4
    h = np.array([1., 0, 0, 0]); print("Q(x1) =", g.Q(g.x @ h))
    z = np.zeros(4)
    for beta in (0.0, 0.05, 0.0833, 0.09):
        print(beta, bound(z, z, z, h, beta, g))

def bound_direct(b1, b2, b3, h, beta, g, swap=False):
    """Order (x1, x2, x3) in the original coordinates (all Lipschitz factors one); swap=True uses
    (x1, x3, x2)."""
    x = g.x; w8 = g.w
    c1 = b1 + 4*beta*x; k1 = np.linalg.norm(c1, axis=-1)
    v2 = b2 + qmul(qconj(x), qconj(b3)[None]); A = Z0(np.linalg.norm(v2, axis=-1))
    Z = np.sum(w8*A*Z0(k1))
    gh = Z0p(k1)*(c1 @ h)/np.maximum(k1, 1e-14); m = np.sum(w8*A*gh)/Z
    t1 = 0.5*np.sum(w8*A*np.sqrt(qser(c1, h, m)))/Z
    def gf(c):
        k = np.linalg.norm(c, axis=-1)
        return Z0p(k)*(c @ h)/np.maximum(k, 1e-14) - m*Z0(k)
    # w = conj(x2 x3).  Outer variable y (x3, or x2 if swap), inner variable z.
    if not swap:   # w = conj(z y) = conj(y) conj(z)
        W = qmul(qconj(x)[:, None], qconj(x)[None, :])          # [y, z]
        bin_, bout = b2, b3
    else:          # inner x3 = z, outer x2 = y: w = conj(y z) = conj(z) conj(y)
        W = qmul(qconj(x)[None, :], qconj(x)[:, None])          # [y, z]
        bin_, bout = b3, b2
    F = np.exp(x @ bin_)[None, :]*gf(b1 + 4*beta*W)               # [y, z]
    t2 = 0.5*np.sum(w8*np.exp(x @ bout)*np.sqrt(np.maximum(g.Q(F), 0)))/Z
    G = F @ w8
    t3 = 0.5*np.sqrt(max(g.Q(np.exp(x @ bout)*G/Z), 0))
    return t1+t2+t3, (t1, t2, t3, m)

def bound_anova(b1, b2, b3, h, beta, g, thetas=None):
    """ANOVA split of Pi_1 u = R + P2 + P3 (interaction plus the two main effects), with the
    interaction shared: T2 = theta R + P3(x2), T3 = (1-theta) R + P2(x3).  Returns min over theta."""
    x = g.x; w8 = g.w
    c1 = b1 + 4*beta*x; k1 = np.linalg.norm(c1, axis=-1)
    v2 = b2 + qmul(qconj(x), qconj(b3)[None]); A = Z0(np.linalg.norm(v2, axis=-1))
    Z = np.sum(w8*A*Z0(k1))
    gh = Z0p(k1)*(c1 @ h)/np.maximum(k1, 1e-14); m = np.sum(w8*A*gh)/Z
    t1 = 0.5*np.sum(w8*A*np.sqrt(qser(c1, h, m)))/Z
    def gf(c):
        k = np.linalg.norm(c, axis=-1)
        return Z0p(k)*(c @ h)/np.maximum(k, 1e-14) - m*Z0(k)
    W = qmul(qconj(x)[:, None], qconj(x)[None, :])        # [y=x3, z=x2]: w = conj(x2 x3)
    V = np.exp(x @ b3)[:, None]*np.exp(x @ b2)[None, :]*gf(b1 + 4*beta*W)/Z
    Py = V @ w8                  # function of x3 (inner x2 averaged)
    Pz = w8 @ V                  # function of x2
    R = V - Py[:, None] - Pz[None, :]
    cR2, cP2 = R @ g.BW.T, Pz @ g.BW.T            # x2-harmonic coefficients, rows y
    cR3, cP3 = R.T @ g.BW.T, Py @ g.BW.T          # x3-harmonic coefficients, rows z
    def val(t):
        a = 0.5*np.sum(w8*np.sqrt(np.maximum((((t*cR2 + cP2[None])**2)/g.deg).sum(-1), 0)))
        b = 0.5*np.sum(w8*np.sqrt(np.maximum(((((1-t)*cR3 + cP3[None])**2)/g.deg).sum(-1), 0)))
        return a + b
    ts = np.linspace(0, 1, 11) if thetas is None else thetas
    vs = [val(t) for t in ts]; i = int(np.argmin(vs))
    return t1 + vs[i], (t1, vs[i], ts[i], m)

def bound_flow(b1, b2, b3, h, beta, g, kmax=2, theta0=0.5):
    """As bound_anova, with a general doubly centred correction S in the span of products of
    harmonics of degree <= kmax in x2 and in x3 (convex minimisation)."""
    from scipy.optimize import minimize
    x = g.x; w8 = g.w
    c1 = b1 + 4*beta*x; k1 = np.linalg.norm(c1, axis=-1)
    v2 = b2 + qmul(qconj(x), qconj(b3)[None]); A = Z0(np.linalg.norm(v2, axis=-1))
    Z = np.sum(w8*A*Z0(k1))
    gh = Z0p(k1)*(c1 @ h)/np.maximum(k1, 1e-14); m = np.sum(w8*A*gh)/Z
    t1 = 0.5*np.sum(w8*A*np.sqrt(qser(c1, h, m)))/Z
    def gf(c):
        k = np.linalg.norm(c, axis=-1)
        return Z0p(k)*(c @ h)/np.maximum(k, 1e-14) - m*Z0(k)
    W = qmul(qconj(x)[:, None], qconj(x)[None, :])
    V = np.exp(x @ b3)[:, None]*np.exp(x @ b2)[None, :]*gf(b1 + 4*beta*W)/Z
    Py = V @ w8; Pz = w8 @ V
    R = V - Py[:, None] - Pz[None, :]
    C2 = (theta0*R) @ g.BW.T + (Pz @ g.BW.T)[None]           # [y, a]  x2-coeffs of T2
    C3 = ((1-theta0)*R).T @ g.BW.T + (Py @ g.BW.T)[None]     # [z, b]  x3-coeffs of T3
    sel = np.where(g.deg <= kmax)[0]; nb = len(sel)
    By = g.B[sel].T                                           # [y, b] basis values (x3)
    Bz = g.B[sel].T                                           # [z, a] basis values (x2)
    d = g.deg
    def obj(s):
        s = s.reshape(nb, nb)                                  # s[a, b]
        D2 = C2.copy(); D2[:, sel] += By @ s.T                 # S = sum s_ab B_a(x2) B_b(x3)
        D3 = C3.copy(); D3[:, sel] -= Bz @ s
        n2 = np.sqrt(np.maximum(((D2**2)/d).sum(1), 1e-30)); n3 = np.sqrt(np.maximum(((D3**2)/d).sum(1), 1e-30))
        f = 0.5*(w8 @ n2 + w8 @ n3)
        G2 = (w8/n2)[:, None]*D2[:, sel]/d[sel]; G3 = (w8/n3)[:, None]*D3[:, sel]/d[sel]
        grad = 0.5*(G2.T @ By - (G3.T @ Bz).T)   # d/ds_ab: sum_y G2[y,a] By[y,b] - sum_z G3[z,b] Bz[z,a]
        return f, grad.ravel()
    f0 = obj(np.zeros(nb*nb))[0]
    r = minimize(obj, np.zeros(nb*nb), jac=True, method='L-BFGS-B', options={'maxiter': 500})
    return t1 + r.fun, (t1, f0, r.fun)

def bound_cert(b1, b2, b3, h, beta, g, K0=2, weight='vmf'):
    """Certifiable form: K <= 1/4 + (1+beta) J, J bounded by the two martingale orders
    (x2 then x3, and x3 then x2), Cauchy-Schwarz over the outer variable with a vMF weight,
    and the degree-K0 majorant of Q."""
    x = g.x; w8 = g.w
    W = qmul(qconj(x)[:, None], qconj(x)[None, :])            # [y=x3, z=x2]
    C1 = b1 + 4*beta*W; K1 = np.linalg.norm(C1, axis=-1)
    E3 = np.exp(x @ b3); E2 = np.exp(x @ b2)
    rho = E3[:, None]*E2[None, :]*Z0(K1); Z = w8 @ rho @ w8; rho /= Z
    psi = Z0p(K1)*(C1 @ h)/np.maximum(K1, 1e-14)/Z0(K1)
    m = w8 @ (rho*psi) @ w8
    V = rho*(psi - m)                                           # density of (psi - E psi) mu23
    def Qt(F):            # rows: functions of the inner variable
        c = F @ g.BW.T
        var = (F**2) @ w8 - (F @ w8)**2
        low = c[:, g.deg <= K0]; dl = g.deg[g.deg <= K0]
        return (low**2/dl).sum(-1) + np.maximum(var - (low**2).sum(-1), 0)/(K0+1)
    out = []
    for VV, Eo in ((V, E3), (V.T, E2)):
        Po = VV @ w8                                            # function of the outer variable
        T_in = VV - Po[:, None]
        a = Eo if weight == 'vmf' else (rho if VV is V else rho.T) @ w8
        a = a/(w8 @ a)
        t_in = 0.5*np.sqrt(w8 @ (Qt(T_in)/a))
        t_out = 0.5*np.sqrt(Qt(Po[None, :])[0])
        out.append(t_in + t_out)
    J = min(out)
    return 0.25 + (1+beta)*J, (J, out, m)
