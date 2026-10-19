"""ctypes wrapper of tsmul.c: product of two truncated series given the pair table (I, J, K) of a ring, with the
same product on absolute values. The library is compiled on first use with gcc -O3 -fPIC -shared."""
import ctypes, numpy as np, os, subprocess
_here=os.path.dirname(os.path.abspath(__file__)); _so=os.path.join(_here,'libtsmul.so')
if not os.path.exists(_so) or os.path.getmtime(_so)<os.path.getmtime(os.path.join(_here,'tsmul.c')):
    subprocess.check_call(['gcc','-O3','-fPIC','-shared','-o',_so,os.path.join(_here,'tsmul.c')])
_lib=ctypes.CDLL(_so)
_P=np.ctypeslib.ndpointer
_lib.tsmul.argtypes=[ctypes.c_long,_P(np.int32),_P(np.int32),_P(np.int32),ctypes.c_long,
    _P(np.complex128),_P(np.complex128),_P(np.complex128),ctypes.c_void_p,ctypes.c_void_p,ctypes.c_void_p]
def cmul(R,a,b,aa=None,ab=None):
    a=np.ascontiguousarray(a,np.complex128); b=np.ascontiguousarray(b,np.complex128); out=np.empty(R.N,np.complex128)
    if aa is None:
        _lib.tsmul(len(R.I),R.I,R.J,R.K,R.N,a,b,out,None,None,None); return out
    aa=np.ascontiguousarray(aa,np.float64); ab=np.ascontiguousarray(ab,np.float64); oa=np.empty(R.N)
    _lib.tsmul(len(R.I),R.I,R.J,R.K,R.N,a,b,out,aa.ctypes.data,ab.ctypes.data,oa.ctypes.data); return out,oa
