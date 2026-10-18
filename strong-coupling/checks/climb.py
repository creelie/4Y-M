import sys; sys.argv=['x',sys.argv[1],sys.argv[2],'2']
import os; exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'hess.py')).read().split('nl=len(idx)')[0])
nl=len(idx)
best=-1
for start in range(2):
  U=[haar() for _ in range(nl)]
  cur=np.linalg.eigvalsh(hessian(U)).max(); step=0.5
  for it in range(150):
    V=[expm(step*np.tensordot(rng.normal(size=g),B,1))@u for u in U]
    v=np.linalg.eigvalsh(hessian(V)).max()
    if v>cur: U,cur=V,v
    else: step*=0.97
  best=max(best,cur); print("start",start,"max eig",cur)
print("best",best,"4d=",4*int(sys.argv[2]))
