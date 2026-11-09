"""Leak-subspace angles for Lemma A4: r(x',y') = ||P_{V_x} P_{V_y}||, V_x = ran (I-P_c) P_{F,x} on C'_x.
Same abelian Gaussian model as gauss_rho.py. Vectors of the Gaussian space are div-free fine currents with
inner product <f,g> = f^T G g (G = scalar Green per component, harmonic mode dropped).
P_{F,x} u = current u' with u'_in = 0, u'_out = u_out - K_{out,in} K_in^+ u_in.
<P_c f, P_c g> = (Pth^T G f)^T M^+ (Pth^T G g), Pth = straight path currents, M = coarse_M.
Cube positions: off + n' z, z in (Z_m)^4, m = bL/n'. Reports r(0,z) for z up to axis permutation, and
R_ang = sum_z r(0,z) (translation invariance, n' multiple of b)."""
import sys, json, itertools, numpy as np, scipy.linalg as sl, scipy.sparse as sp
from gauss_rho import sites, sidx, d0T, scalar_green, coarse_M, cube_K

def nullsp(A, rcond=None):  # null space via the gesvd driver (gesdd fails to converge on some cube matrices)
    u, sv, vh = sl.svd(A, full_matrices=True, lapack_driver='gesvd')
    tol = (sv.max() if sv.size else 0.0) * (rcond if rcond is not None else max(A.shape) * np.finfo(float).eps)
    return vh[int((sv > tol).sum()):].conj().T

def full_K(N):
    X = sites(N); r, c, v = [], [], []; row = 0
    for a in range(4):
        for cc in range(a+1, 4):
            ea, ec = np.eye(4, dtype=int)[a], np.eye(4, dtype=int)[cc]
            for lk, sg in ((sidx(X, N)*4+a, 1), (sidx(X+ea, N)*4+cc, 1), (sidx(X+ec, N)*4+a, -1), (sidx(X, N)*4+cc, -1)):
                r.append(row+np.arange(len(X))); c.append(lk); v.append(np.full(len(X), sg))
            row += len(X)
    d1 = sp.csr_matrix((np.concatenate(v), (np.concatenate(r), np.concatenate(c))), shape=(row, 4*N**4))
    return (d1.T @ d1).tocsc()

class Model:
    def __init__(s, b, L, npr):
        s.b, s.L, s.npr, s.N = b, L, npr, b*L; N = s.N; n = L**4
        k = 2*np.pi*np.fft.fftfreq(N); lam = sum(np.meshgrid(*(2-2*np.cos(k),)*4, indexing='ij'))
        with np.errstate(divide='ignore'): s.ghat = np.where(lam > 1e-12, 1/lam, 0.0)
        s.K = full_K(N); s.M = coarse_M(b, L); s.Mp = np.linalg.pinv(s.M, hermitian=True); s.D = d0T(L)
        s.Z = nullsp(s.D)
        X = sites(L); r, c = [], []
        for mu in range(4):
            for i in range(b):
                f = b*X; f[:, mu] += i; r.append(sidx(f, N)*4+mu); c.append(sidx(X, L)*4+mu)
        s.Pth = sp.csr_matrix((np.ones(4*n*b), (np.concatenate(r), np.concatenate(c))), shape=(4*N**4, 4*n))

    def G(s, F):  # apply Green, F: (4N^4, m)
        N = s.N; m = F.shape[1]; A = F.reshape(N, N, N, N, 4, m)
        out = np.fft.ifftn(np.fft.fftn(A, axes=(0, 1, 2, 3)) * s.ghat[..., None, None], axes=(0, 1, 2, 3)).real
        return out.reshape(4*N**4, m)

    def harm0(s, F):  # drop the constant (harmonic) component per direction: G ignores it, K_in^+ must too
        N = s.N; A = F.reshape(N**4, 4, -1); return (A - A.mean(0, keepdims=True)).reshape(F.shape)

    def cube(s, off):
        b, L, N, n = s.b, s.L, s.N, s.L**4
        inl, Kin = cube_K(b, L, s.npr, off)
        gin = np.array([sidx(np.array(x), N)*4+mu for (x, mu) in inl])
        Pin = s.Pth[gin]; T = np.where(np.abs(Pin).sum(0).A1 > 0)[0]; O = np.setdiff1d(np.arange(4*n), T)
        zo = nullsp(s.D[:, O]); ZO = np.zeros((4*n, zo.shape[1])); ZO[O] = zo
        U = s.Z @ nullsp(ZO.T @ s.M @ s.Z)            # coarse coefficients of C'
        u = s.harm0(s.Pth @ U)                                 # fine currents, torus-harmonic part removed
        Kp = np.linalg.pinv(Kin, hermitian=True)
        corr = s.K[:, gin] @ (Kp @ u[gin])
        up = u - corr; up[gin] = 0.0                           # P_F u
        Gu = s.G(up); Y = s.Pth.T @ Gu                         # for P_c
        return dict(gin=gin, Kp=Kp, U=U, u=u, up=up, Gu=Gu, Y=Y, dimC=U.shape[1])

    def gram(s, A, B):  # <(I-P_c)P_F a_i, (I-P_c)P_F b_j>
        return A['up'].T @ B['Gu'] - A['Y'].T @ s.Mp @ B['Y']

def whiten(Gm, tol=1e-10):
    w, V = np.linalg.eigh((Gm+Gm.T)/2); keep = w > tol*w.max(); return V[:, keep] / np.sqrt(w[keep])

def run(b, L, npr, o):
    mdl = Model(b, L, npr); m = mdl.N // npr; off0 = np.array([o]*4)
    A = mdl.cube(off0); WA = whiten(mdl.gram(A, A))
    res = []
    for z in itertools.combinations_with_replacement(range(m), 4):
        if not any(z): continue
        # distance on the cube torus: min(z, m - z) per axis
        B = mdl.cube(off0 + npr*np.array(z)); WB = whiten(mdl.gram(B, B))
        r = float(np.linalg.svd(WA.T @ mdl.gram(A, B) @ WB, compute_uv=False).max())
        orbit = len(set(itertools.permutations(z)))
        res.append(dict(z=list(z), dinf=int(max(min(t, m-t) for t in z)), orbit=orbit, r=r))
        print(json.dumps(dict(b=b, L=L, nprime=npr, **res[-1])), flush=True)
    Rang = 1.0 + sum(x['orbit']*x['r'] for x in res)
    print(json.dumps(dict(b=b, L=L, nprime=npr, off=o, dimC=A['dimC'], R_ang=Rang,
                          by_dinf={d: sum(x['orbit']*x['r'] for x in res if x['dinf'] == d) for d in sorted({x['dinf'] for x in res})})), flush=True)

def selftest():
    b, L, npr = 3, 4, 3; mdl = Model(b, L, npr); off = np.array([0]*4); A = mdl.cube(off)
    # (1) Var(u) - Var(P_F u) == u_in^T K_in^+ u_in  (Pythagoras for conditional expectation)
    lhs = np.einsum('ij,ij->j', A['u'], mdl.G(A['u'])) - np.einsum('ij,ij->j', A['up'], A['Gu'])
    rhs = np.einsum('ij,ij->j', A['u'][A['gin']], A['Kp'] @ A['u'][A['gin']])
    e1 = np.abs(lhs-rhs).max()/np.abs(rhs).max(); print('pythagoras rel err', e1); assert e1 < 1e-8
    # (2) P_F u is divergence-free (gauge invariant)
    N = mdl.N; Df = sp.csr_matrix(d0T(N)) if N**4 <= 20000 else None
    if Df is not None:
        e2 = np.abs(Df @ A['up']).max()/np.abs(A['up']).max(); print('div P_F u', e2); assert e2 < 1e-8
    # (3) consistency with gauss_rho: max ||P_F u||^2/||u||^2 on C' equals rho1^2
    from gauss_rho import rho1
    Gu0 = mdl.G(A['u']); ev = sl.eigh(A['up'].T @ A['Gu'], A['u'].T @ Gu0, eigvals_only=True)
    # gauss_rho.rho1 keeps the harmonic part in u_in; report the size of that bias (not asserted)
    r1 = rho1(b, L, npr, (0, 0, 0, 0))['rho1']; print(json.dumps(dict(rho1_gauss_rho=r1, rho1_harm0=float(np.sqrt(ev.max())))), flush=True)
    assert abs(r1 - np.sqrt(ev.max())) < 5e-2
    # (4) r in [0,1] and self-angle = 1
    WA = whiten(mdl.gram(A, A)); s0 = np.linalg.svd(WA.T @ mdl.gram(A, A) @ WA, compute_uv=False).max(); print('self r', s0); assert abs(s0-1) < 1e-8
    print('SELFTEST PASS')

if __name__ == '__main__':
    if sys.argv[1] == 'selftest': selftest()
    else: run(*map(int, sys.argv[1:5]))
