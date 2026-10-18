"""Exact Hessian of S0(Q) = sum_p Re Tr Q_p for SU(2) lattice gauge fields on the periodic
torus (Z/L)^d, w.r.t. right-invariant fields Q_e -> e^{tX_e} Q_e (SZZ convention), metric
<X,Y> = -Tr(XY).  Orthonormal basis of su(2): T_k = i sigma_k / sqrt 2."""
import numpy as np, itertools
s = [np.array([[0,1],[1,0]],complex), np.array([[0,-1j],[1j,0]]), np.array([[1,0],[0,-1]],complex)]
T = [1j*sk/np.sqrt(2) for sk in s]          # -Tr(T_a T_b) = delta_ab

def expsu2(x):
    """exp(sum x_k T_k)"""
    th = np.linalg.norm(x)/np.sqrt(2)        # (x.T)^2 = -(|x|^2/2) I
    M = sum(xk*Tk for xk,Tk in zip(x,T))
    if th < 1e-14: return np.eye(2)+M
    return np.cos(th)*np.eye(2) + np.sin(th)/th*M

class Lattice:
    def __init__(s_, d, L):
        s_.d, s_.L = d, L
        s_.sites = list(itertools.product(range(L), repeat=d))
        s_.idx = {}
        for x in s_.sites:
            for mu in range(d): s_.idx[(x,mu)] = len(s_.idx)
        s_.nE = len(s_.idx)
        s_.plaq = []
        for x in s_.sites:
            for mu in range(d):
                for nu in range(mu+1,d):
                    xm = tuple((x[k]+(k==mu))%L for k in range(d))
                    xn = tuple((x[k]+(k==nu))%L for k in range(d))
                    # Q_p = Q_{x,mu} Q_{x+mu,nu} Q_{x+nu,mu}^* Q_{x,nu}^*
                    s_.plaq.append((s_.idx[(x,mu)], s_.idx[(xm,nu)], s_.idx[(xn,mu)], s_.idx[(x,nu)]))

def plaquette_hessian(Q, p):
    """12x12 Hessian of Re Tr Q_p in coordinates (e1,e2,e3,e4) x (T_1,T_2,T_3)."""
    e = p
    # factor list: (kind, data): 'v',slot,sign  or 'c',matrix
    fac = [('v',0,+1),('c',Q[e[0]]),('v',1,+1),('c',Q[e[1]]),('c',Q[e[2]].conj().T),('v',2,-1),('c',Q[e[3]].conj().T),('v',3,-1)]
    def prod(ins):
        M = np.eye(2,dtype=complex)
        for f in fac:
            if f[0]=='c': M = M@f[1]
            else:
                if f[1] in ins: M = M@ins[f[1]]
        return M
    H = np.zeros((12,12))
    for i in range(4):
        si = fac[[k for k,f in enumerate(fac) if f[0]=='v' and f[1]==i][0]][2]
        for a in range(3):
            for j in range(4):
                sj = fac[[k for k,f in enumerate(fac) if f[0]=='v' and f[1]==j][0]][2]
                for b in range(3):
                    if i==j:
                        Z = 0.5*(T[a]@T[b]+T[b]@T[a])   # sign^2 = 1
                        val = np.trace(prod({i:Z})).real
                    else:
                        val = np.trace(prod({i:si*T[a], j:sj*T[b]})).real
                    H[3*i+a,3*j+b] = val
    return H

def full_hessian(lat, Q):
    n = 3*lat.nE
    H = np.zeros((n,n))
    for p in lat.plaq:
        Hp = plaquette_hessian(Q,p)
        ids = np.concatenate([np.arange(3*e,3*e+3) for e in p])
        H[np.ix_(ids,ids)] += Hp
    return H

def lam_max(lat,Q):
    return np.linalg.eigvalsh(full_hessian(lat,Q))[-1]

def frustrated(lat):
    Q=[None]*lat.nE
    for (x,mu),k in lat.idx.items():
        Q[k] = ((-1)**(sum(x[:mu])%2))*np.eye(2,dtype=complex)
    return Q

if __name__=="__main__":
    for d in (2,3,4):
        lat = Lattice(d,2)
        Qid=[np.eye(2,dtype=complex)]*lat.nE
        ev=np.linalg.eigvalsh(full_hessian(lat,Qid))
        Qf=frustrated(lat)
        evf=np.linalg.eigvalsh(full_hessian(lat,Qf))
        print(d, "identity: min,max", ev[0], ev[-1], " frustrated: min,max", evf[0], evf[-1], " SZZ 8(d-1)=",8*(d-1), " 4d=",4*d)
