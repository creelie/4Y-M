import numpy as np, itertools, sys
from scipy.linalg import expm
rng=np.random.default_rng(1)
N=int(sys.argv[1]); d=int(sys.argv[2]); L=int(sys.argv[3])
# su(N) basis orthonormal for <X,Y>=-Tr(XY)
B=[]
for a in range(N):
  for b in range(a+1,N):
    M=np.zeros((N,N),complex); M[a,b]=1; M[b,a]=-1; B.append(M/np.sqrt(2))
    M=np.zeros((N,N),complex); M[a,b]=1j; M[b,a]=1j; B.append(M/np.sqrt(2))
for k in range(1,N):
  M=np.zeros((N,N),complex)
  for j in range(k): M[j,j]=1j
  M[k,k]=-1j*k; B.append(M/np.sqrt(k*(k+1)))
B=np.array(B); g=len(B)
for X in B: assert abs(np.trace(X@X).real+1)<1e-12
sites=list(itertools.product(range(L),repeat=d))
idx={}
for s in sites:
  for mu in range(d): idx[(s,mu)]=len(idx)
def sh(s,mu): s=list(s); s[mu]=(s[mu]+1)%L; return tuple(s)
plaq=[]
for s in sites:
  for mu in range(d):
    for nu in range(mu+1,d):
      plaq.append((idx[(s,mu)],idx[(sh(s,mu),nu)],idx[(sh(s,nu),mu)],idx[(s,nu)]))
def haar():
  Z=(rng.normal(size=(N,N))+1j*rng.normal(size=(N,N)))/np.sqrt(2)
  Q,R=np.linalg.qr(Z); Q=Q@np.diag(np.diag(R)/abs(np.diag(R)))
  return Q/np.linalg.det(Q)**(1/N)
def hessian(U):
  n=len(U)*g; H=np.zeros((n,n)); h=1e-4
  for p in plaq:
    def f(t):
      M=[expm(sum(t[4*0+0:0]) if False else np.zeros((N,N))) for _ in range(0)]
      V=[expm(np.tensordot(t[i*g:(i+1)*g],B,1))@U[p[i]] for i in range(4)]
      return np.trace(V[0]@V[1]@V[2].conj().T@V[3].conj().T).real
    m=4*g; Hl=np.zeros((m,m)); f0=f(np.zeros(m))
    E=np.eye(m)*h
    for i in range(m):
      for j in range(i,m):
        v=(f(E[i]+E[j])-f(E[i]-E[j])-f(-E[i]+E[j])+f(-E[i]-E[j]))/(4*h*h)
        Hl[i,j]=Hl[j,i]=v
    ids=np.concatenate([np.arange(p[i]*g,(p[i]+1)*g) for i in range(4)])
    H[np.ix_(ids,ids)]+=Hl
  return H
nl=len(idx)
res=[]
for trial in range(3):
  U=[haar() for _ in range(nl)]
  ev=np.linalg.eigvalsh(hessian(U)); res.append((ev.min(),ev.max()))
print("random", res, "bound 4d=",4*d, "SZZ 8(d-1)=",8*(d-1))
U=[np.eye(N,dtype=complex) for _ in range(nl)]
ev=np.linalg.eigvalsh(hessian(U)); print("identity",ev.min(),ev.max())
z=np.eye(N,dtype=complex); z[0,0]=z[1,1]=-1
U=[None]*nl
for s in sites:
  for mu in range(d):
    U[idx[(s,mu)]]=np.linalg.matrix_power(z,sum(s[:mu]))
ev=np.linalg.eigvalsh(hessian(U)); print("frustrated",ev.min(),ev.max())
