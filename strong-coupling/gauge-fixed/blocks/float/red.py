import numpy as np
def fields(p, beta):
    """6 parameters after the SU(2)^3 symmetry: r2, r3, rho1 in [0,1] (times 20 beta), alpha, gamma, delta."""
    R = 20*beta
    r2, r3, r1, al, ga, de = p
    b2 = np.array([R*r2, 0, 0, 0.]); b3 = np.array([R*r3, 0, 0, 0.])
    b1 = R*r1*np.array([np.cos(al), np.sin(al), 0, 0.])
    h = np.array([np.cos(ga), np.sin(ga)*np.cos(de), np.sin(ga)*np.sin(de), 0.])
    return b1, b2, b3, h
LO = np.array([0, 0, 0, 0, 0, 0.]); HI = np.array([1, 1, 1, np.pi, np.pi, np.pi])
