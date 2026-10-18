import sys; N0,d0,L0=sys.argv[1:4]; sys.argv=['x',N0,d0,L0]
import os; exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'hess.py')).read().split('nl=len(idx)')[0])
nl=len(idx)
def DD(U):
  # D_x X = sum_mu X_{x,mu} - U_{x-mu,mu}^* X_{x-mu,mu} U_{x-mu,mu}; matrix rows: site x, components in basis
  n=nl*g; rows=[]
  Bc=B.reshape(g,-1)
  Dm=np.zeros((len(sites)*g,n))
  for si,s in enumerate(sites):
    for mu in range(d):
      e=idx[(s,mu)]
      for k in range(g): Dm[si*g+k, e*g+k]+=1
      sm=list(s); sm[mu]=(sm[mu]-1)%L; e2=idx[(tuple(sm),mu)]; Q=U[e2]
      for k in range(g):
        Y=Q.conj().T@B[k]@Q
        coeff=[-np.trace(Y@B[j]).real for j in range(g)]
        for j in range(g): Dm[si*g+j, e2*g+k]-=coeff[j]
  return Dm.T@Dm
worst=-1
for t in range(4):
  U=[haar() for _ in range(nl)]
  H=hessian(U); M=DD(U)
  worst=max(worst,np.linalg.eigvalsh(H+M).max(),np.linalg.eigvalsh(-H+M).max())
print(N0,d0,L0,"max eig of ±Hess+D*D over random:",worst,"bound",4*d)
