"""Checks of Theorem "effectivelocal": the effective forms S_k of the block-constrained
quadratic approximation are S_k(w) = sum_P int K_P |F_P^|^2, with
K_P = prod_{s not in P} atilde_s / Y the functions of the paper (Feynman gauge, the pole of the
fine propagator cancelled), and the constants of the proof of part (ii).
(A) K_P against the constrained minimum computed directly from the plaquette form (pure gauges
    projected out, no Feynman gauge), at random momenta, k = 1, 2, 3.
(B) K_P(0) = 1, K_P > 0 and Y >= (2/pi)^32 on a grid of the real torus.
(C) the bounds of Step 4 at random points of the strip |Im theta| <= 1/2, k = 1, 2, 3.
(D) the arithmetic of Step 4: c_1, the bound for r, M, eta = 1e-36, C = 1e26.
(E) floating point (no statement rests on it): the decay of the kernels of K_P by
    l1 distance and the range of K_P on a grid of 12^4 momenta, k = 1, 2, 3.
Run: python3 effective_locality.py (about eight minutes). Needs numpy."""
import numpy as np, itertools, math
PL = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
passed = failed = 0


def check(name, ok):
    global passed, failed
    print(('PASS ' if ok else 'FAIL ') + name)
    if ok:
        passed += 1
    else:
        failed += 1


def DL(t, L):
    t = np.asarray(t, complex)
    return sum(np.exp(1j * t * j) for j in range(L)) / L


def sharp2(f, t, L):
    """f f^sharp at t: the analytic continuation of |f|^2."""
    return f(t, L) * np.conj(f(np.conj(t), L))


def blocks(th, k, want_parts=False):
    L = 2 ** k
    th = np.asarray(th, complex)
    r = np.zeros(th.shape[:-1] + (4,), complex)
    minre = np.inf
    for n in itertools.product(range(L), repeat=4):
        q = (th + 2 * np.pi * np.array(n)) / L
        q = q - 2 * np.pi * np.round(q.real / (2 * np.pi))
        D2 = sharp2(DL, q, L)
        Psi2 = np.prod(D2, -1)
        eps = L * L * np.sum(2 - 2 * np.cos(q), -1)
        if not any(n):
            Psi0, D0, eh = Psi2, D2, eps
        else:
            r += Psi2[..., None] * D2 / eps[..., None]
            s = L * q.real
            minre = min(minre, np.min(eps.real / np.sum(s * s, -1)))
    E2 = (np.exp(1j * th) - 1) * (np.exp(-1j * th) - 1)
    A = Psi0[..., None] * D0
    at = A + eh[..., None] * r
    Y = Psi0 ** 3 * np.prod(D0, -1)
    for s in range(4):
        oth = [t for t in range(4) if t != s]
        for m in range(1, 4):
            for J in itertools.combinations(oth, m):
                term = eh ** (m - 1)
                for t in oth:
                    term = term * (r[..., t] if t in J else A[..., t])
                Y = Y + E2[..., s] * term
    if want_parts:
        return at, Y, r, A, eh, minre
    return at, Y


def K(th, k):
    at, Y = blocks(th, k)
    return np.stack([np.prod([at[..., s] for s in range(4) if s not in P], 0) / Y for P in PL], -1)


def S_from_K(th, k):
    KK = K(th, k)
    E = np.exp(1j * np.asarray(th, complex)) - 1
    out = np.zeros(E.shape[:-1] + (4, 4), complex)
    for i, (l, m) in enumerate(PL):
        c = np.zeros(E.shape, complex)
        c[..., m] = E[..., l]
        c[..., l] = -E[..., m]
        out += KK[..., i][..., None, None] * np.conj(c)[..., :, None] * c[..., None, :]
    return out


def S_direct(th, k):
    """min Q at one coarse momentum from the plaquette form itself: the fine field at the
    L^4 aliases, the constraint L^-4 sum_n c(theta_n) x_n = w, pure gauges projected out."""
    L = 2 ** k
    th = np.asarray(th, float)
    G = np.zeros((4, 4), complex)
    for n in itertools.product(range(L), repeat=4):
        q = (th + 2 * np.pi * np.array(n)) / L
        e = np.exp(1j * q) - 1
        Qm = np.sum(abs(e) ** 2) * np.eye(4) - np.outer(np.conj(e), e).T
        ee = np.outer(e, np.conj(e))
        Minv = np.linalg.inv(Qm + ee) - ee / np.sum(abs(e) ** 2) ** 2
        c = L * np.prod(DL(q, L)) * DL(q, L)
        G += (c[:, None] * Minv * np.conj(c)[None, :]) / L ** 4
    b = np.exp(1j * th) - 1
    P = np.eye(4) - np.outer(b, np.conj(b)) / np.vdot(b, b)
    # minimum over w + C b of z^* G^+ z restricted: S = (P G P)^+ on the complement of b
    return np.linalg.pinv(P @ G @ P, rcond=1e-12)


def main():
    global passed, failed
    rng = np.random.default_rng(7)
    # (A)
    for k in (1, 2, 3):
        pts = rng.uniform(-np.pi, np.pi, (6 if k < 3 else 3, 4))
        err = max(np.abs(S_from_K(p, k) - S_direct(p, k)).max() / np.abs(S_direct(p, k)).max() for p in pts)
        check('(A) k=%d: sum_P K_P c_P^* c_P equals the constrained minimum (rel. err %.1e)' % (k, err), err < 1e-8)
    # (B)
    N = 8
    xs = np.linspace(-np.pi, np.pi, N)
    X = np.stack(np.meshgrid(xs, xs, xs, xs, indexing='ij'), -1).reshape(-1, 4)
    for k in (1, 2, 3):
        at, Y = blocks(np.zeros(4), k)
        K0 = np.array([np.prod([at[s] for s in range(4) if s not in P]) / Y for P in PL])
        check('(B) k=%d: K_P(0) = 1' % k, np.abs(K0 - 1).max() < 1e-12)
        at, Y = blocks(X, k)
        KK = K(X, k)
        check('(B) k=%d: K_P real and positive, Y >= (2/pi)^32 on a grid (min Y %.3g)' % (k, Y.real.min()),
              np.abs(KK.imag).max() < 1e-9 and KK.real.min() > 0 and Y.real.min() >= (2 / np.pi) ** 32)
    # (C)
    for k in (1, 2, 3):
        pts = rng.uniform(-np.pi, np.pi, (300, 4)) + 0.5j * rng.uniform(-1, 1, (300, 4))
        at, Y, r, A, eh, minre = blocks(pts, k, True)
        ok = (np.abs(r).max() < 2.2e7 and np.abs(Y).max() < 1e28 and np.abs(at).max() < 1.6e9
              and np.abs(eh).max() < 70 and np.abs(A).max() <= math.e ** 5 and minre >= 0.28)
        check('(C) k=%d: strip bounds |r|<2.2e7, |Y|<1e28, |atilde|<1.6e9, |eps^|<70, Re L^2 eps_n >= 0.28|Re s|^2'
              ' (max|r| %.3g, min ratio %.3f)' % (k, np.abs(r).max(), minre), ok)
    # (D)
    c1 = (math.exp(0.5) + 1) * math.exp(0.25) * math.pi / 2
    rb = c1 ** 10 / (0.28 * math.pi ** 2) * (4 / 3) ** 4
    eh = 4 * (math.pi ** 2 + 0.25) * math.exp(0.5)
    EE = (math.exp(0.5) + 1) ** 2
    A = math.e ** 5
    M = math.e ** 16 + 4 * EE * (3 * rb * A ** 2 + 3 * eh * rb ** 2 * A + eh ** 2 * rb ** 3)
    low = (2 / math.pi) ** 32
    eta, C = 1e-36, 1e26
    check('(D) c1 = %.4f < 5.35, r bound %.3g < 2.2e7, M = %.3g <= 1e28' % (c1, rb, M),
          c1 < 5.35 and rb < 2.2e7 and M <= 1e28 and eh < 70 and EE < 7.1)
    check('(D) 16 eta M <= (2/pi)^32/2 and (A + 70 r)^2 / ((2/pi)^32/2) = %.3g <= C' % ((A + 70 * 2.2e7) ** 2 / (low / 2)),
          16 * eta * 1e28 <= low / 2 and A + 70 * 2.2e7 < 1.6e9 and (1.6e9) ** 2 / (low / 2) <= C)
    # (E)
    N = 12
    js = np.arange(N) * 2 * np.pi / N
    X = np.stack(np.meshgrid(js, js, js, js, indexing='ij'), -1).reshape(-1, 4)
    for k in (1, 2, 3):
        at, Y = blocks(X, k)
        KP = (np.prod(at[:, 2:], -1) / Y).real
        Kr = np.fft.ifftn(KP.reshape((N,) * 4)).real
        rr = np.indices((N,) * 4)
        rr = np.minimum(rr, N - rr).sum(0)
        mag = [np.abs(Kr)[rr == d].max() for d in range(6)]
        fac = [mag[d + 1] / mag[d] for d in range(1, 5)]
        print('     k=%d: K_P in [%.3f, %.2f]; max|kernel| by l1 distance %s; factors %s'
              % (k, KP.min(), KP.max(), ' '.join('%.3g' % m for m in mag), ' '.join('%.3f' % f for f in fac)))
    print('%d checks passed, %d failed' % (passed, failed))


if __name__ == '__main__':
    main()
