__all__ = ['fiducial_deconvolute', 'af_dtype']
import numpy as np
import numpy.ctypeslib as C

# define types
af_dtype = np.dtype([('key', np.float64), ('val', np.float64)], align=True)
_af_ctype =  C.ndpointer(af_dtype, ndim=1, flags='C,A')
_double_ctype = C.ndpointer(np.float64, ndim=1, flags='C,A')

# load the c library
_C_LIB = C.load_library('fiducial_deconvolute', __path__[0])
_C_LIB.convolved_fit.restype = None
_C_LIB.convolved_fit.argtypes = [_af_ctype, C.ctypes.c_int, \
        _double_ctype, _double_ctype, C.ctypes.c_int, C.ctypes.c_double,\
        C.ctypes.c_int, C.ctypes.c_double]

def fiducial_deconvolute(af, smm, mf, scatter, repeat=40, sm_step=0.01):
    if len(smm) != len(mf):
        raise ValueError('`smf` and `mf` must have the same size!')
    sm_step = np.fabs(float(sm_step))
    sm_min = min(af['key'].min(), smm.min())
    if sm_min <= 0:
        offset = sm_step-sm_min
        af['key'] += offset
        smm += offset
    _C_LIB.convolved_fit(af, len(af), smm, mf, len(mf), float(scatter), \
        int(repeat), sm_step)
    if sm_min <= 0:
        smm -= offset
        af['key'] -= offset
    return smm

