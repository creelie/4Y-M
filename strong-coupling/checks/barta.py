import sympy as sp
t,k,a=sp.symbols('theta kappa a',real=True)
b=2*sp.cot(t)-k*sp.sin(t)
psi=sp.sin(t)*sp.exp(a*k*sp.cos(t))
L1=sp.diff(psi,t,2)+b*sp.diff(psi,t)-2/sp.sin(t)**2*psi
M=sp.diff(psi,t,2)+b*sp.diff(psi,t)+sp.diff(b,t)*psi
e1=sp.simplify(-L1/psi-(3+(5*a+1)*k*sp.cos(t)-a*(a+1)*k**2*sp.sin(t)**2))
e2=sp.simplify(-M/psi-(3+(5*a+2)*k*sp.cos(t)-a*(a+1)*k**2*sp.sin(t)**2))
print(e1,e2)
