"""Floating-point symbol of the effective forms S_k(w) = inf{Q(phi) : V^(k)(phi) = w}
of the block-constrained quadratic approximation, for k = 1, 2, 3, on a half-shifted grid
of 12^4 coarse momenta: the k-fold averages are the straight-path averages over blocks of
side L = 2^k (Lemma "straightpaths"), the constraint is exact, and the pure gauges are
projected out. Prints the range of S_k/Q_k on curvatures and the decay of the kernel of S_k
by coarse l-infinity distance. The paper states these values as floating-point observations
only; no statement rests on them. Run: python3 effective_symbol.py (about six minutes).
Needs numpy."""
import numpy as np, itertools
def Mmat(a):
    # a: (...,4) complex; Q = sum_{mu<nu} |a_mu phi_nu - a_nu phi_mu|^2
    sh=a.shape[:-1]; M=np.zeros(sh+(4,4),complex)
    for m in range(4):
        for n in range(m+1,4):
            r=np.zeros(sh+(4,),complex); r[...,n]=a[...,m]; r[...,m]=-a[...,n]
            M+=np.conj(r)[...,:,None]*r[...,None,:]
    return M
def D(th,L):
    th=np.asarray(th); out=np.ones_like(th,dtype=complex)*0
    for s in range(L): out=out+np.exp(1j*th*s)
    return out/L
def symbol(k,N):
    L=2**k
    js=(np.arange(N)+0.5)*2*np.pi/N
    Q=np.stack(np.meshgrid(js,js,js,js,indexing='ij'),-1).reshape(-1,4)
    G0=np.zeros((Q.shape[0],4,4),complex)
    for n in itertools.product(range(L),repeat=4):
        p=(Q+2*np.pi*np.array(n))/L
        a=np.exp(1j*p)-1
        M=Mmat(a); aa=a[...,:,None]*np.conj(a)[...,None,:]
        na=np.sum(abs(a)**2,-1)
        Minv=np.linalg.inv(M+aa)-aa/(na**2)[:,None,None]
        Dp=np.prod(D(p,L),-1)
        c=Dp[:,None]*L*D(p,L)   # c_mu
        G0+=c[:,:,None]*Minv*np.conj(c)[:,None,:]
    b=np.exp(1j*Q)-1
    Gi=np.linalg.inv(G0)
    Gb=np.einsum('qij,qj->qi',Gi,b)
    S=Gi-Gb[:,:,None]*np.conj(Gb)[:,None,:]/np.einsum('qi,qi->q',np.conj(b),Gb)[:,None,None]
    return Q,S,b
if __name__=='__main__':
    N=12
    for k in (1,2,3):
        Q,S,b=symbol(k,N); L4=16**k
        S=S*L4   # the form in coarse-lattice units
        gauge=np.max(abs(np.einsum('qij,qj->qi',S,b)))
        M1=Mmat(b); ev=[]
        for i in range(len(Q)):
            w,V=np.linalg.eigh(M1[i]); P=V[:,w>1e-9]
            ev+=list(np.linalg.eigvals(np.linalg.solve(np.diag(w[w>1e-9]),P.conj().T@S[i]@P)).real)
        K=np.fft.ifftn(S.reshape(N,N,N,N,4,4),axes=(0,1,2,3))
        r=np.indices((N,)*4); r=np.minimum(r,N-r); rinf=r.max(0)
        mag=[np.abs(K).max((-1,-2))[rinf==d].max() for d in range(N//2+1)]
        print('k = %d: S_k vanishes on gradients (%.1e); S_k/Q_k on curvatures in [%.3f, %.2f]' % (k,gauge,min(ev),max(ev)))
        print('       kernel by l-infinity distance:', ' '.join('%.3g' % m for m in mag[:6]),
              '; ratio %.3f' % (mag[4]/mag[3]))
