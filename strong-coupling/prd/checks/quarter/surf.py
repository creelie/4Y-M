import numpy as np, sys
sys.argv=['x','26']
exec(open('lip2d.py').read().split('n=int(sys.argv[1])')[0])
from scipy.special import iv
out=open('surf.dat','w')
for k in np.linspace(0,2.0,9):
    for phi in np.linspace(0,np.pi/2,7):
        r=rate(k,phi,26)
        out.write(f"{k:.4f} {phi:.4f} {r:.5f}\n"); out.flush()
    out.write("\n")
