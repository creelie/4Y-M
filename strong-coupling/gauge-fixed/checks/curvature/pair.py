import numpy as np, sys, os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
def Mp(a,b,t,th):
    c=np.cos(th); s=th/2
    z=(1+t)*np.cos(s)-1j*(1-t)*np.sin(s)
    A=a+c; B=b+c
    return np.array([[A,z,-1],[np.conj(z),B,-z],[-1,-np.conj(z),A]])
def minpair(a,b,t,n=20000):
    ths=np.linspace(0,2*np.pi,n); return min(np.linalg.eigvalsh(Mp(a,b,t,th)).min() for th in ths)
def mindiag(a,b,t):
    out=[]
    for cphi in (-1,1):
        # a(y1^2+y3^2)+b y2^2 + cphi (y1+y2-y3)^2 -t(-2y1y2+2y2y3)
        Mq=np.diag([a,b,a]).astype(float)
        v=np.array([1,1,-1.]); Mq+=cphi*np.outer(v,v)
        Mq[0,1]+=t; Mq[1,0]+=t; Mq[1,2]-=t; Mq[2,1]-=t
        out.append(np.linalg.eigvalsh(Mq).min())
    return min(out)
if __name__=="__main__":
    a,b,t=map(float,sys.argv[1:4])
    print('pair',minpair(a,b,t),'diag',mindiag(a,b,t),'h',4*a+2*b+4*t)
