#cython: embedsignature=True
#distutils: language = c++
"""
A Python wrapper around the C DIAlign library
"""

import numpy as np
cimport numpy as np

from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from cython.operator cimport dereference as deref, preincrement as inc, address as address
from DIAlign cimport getQuantile as _getQuantile
from DIAlign cimport AffineAlignObj as _AffineAlignObj
from libcpp.memory cimport shared_ptr
from DIAlign cimport alignChromatogramsCpp as _alignChromatogramsCpp

def getQuantile(data, q):
    """
        Compute quantile q on the data
        
    """
    cdef libcpp_vector[double] c_data = data
    cdef double result = _getQuantile( c_data, q)
    return result

cdef class AffineAlignObj:
    """
    Cython implementation of _AffineAlignObj
    """
    cdef shared_ptr[_AffineAlignObj] inst

    def __dealloc__(self):
         self.inst.reset()

    def __init__(self, int row, int col):
        """Cython signature: void AffineAlignObj()"""
        self.inst = shared_ptr[_AffineAlignObj](new _AffineAlignObj(row, col))

    property score:
        def __set__(self, np.ndarray[double, ndim=1, mode="c"] score not None):
            self.inst.get().score = score
    
        def __get__(self):
            return self.inst.get().score

    property indexA_aligned:
        def __set__(self, np.ndarray[double, ndim=1, mode="c"] indexA_aligned not None):
            self.inst.get().indexA_aligned = indexA_aligned
    
        def __get__(self):
            return self.inst.get().indexA_aligned

    property indexB_aligned:
        def __set__(self, np.ndarray[double, ndim=1, mode="c"] indexB_aligned not None):
            self.inst.get().indexB_aligned = indexB_aligned
    
        def __get__(self):
            return self.inst.get().indexB_aligned

def alignChromatogramsCppSimple(AffineAlignObj obj,
                                np.ndarray[double, ndim=2, mode="c"] r1 not None,
                                np.ndarray[double, ndim=2, mode="c"] r2 not None,
                                bytes alignType,
                                np.ndarray[double, ndim=1, mode="c"] tA not None,
                                np.ndarray[double, ndim=1, mode="c"] tB not None,
                                bytes normalization,
                                bytes simType):
    
    _alignChromatogramsCpp(deref(obj.inst.get()), r1, r2, alignType, tA, tB, normalization, simType)
           

    



def alignChromatogramsCpp(      AffineAlignObj obj,
                                np.ndarray[double, ndim=2, mode="c"] r1 not None,
                                np.ndarray[double, ndim=2, mode="c"] r2 not None,
                                bytes alignType,
                                np.ndarray[double, ndim=1, mode="c"] tA not None,
                                np.ndarray[double, ndim=1, mode="c"] tB not None,
                                bytes normalization,
                                bytes simType,
                                double B1p,
                                double B2p,
                                int noBeef,
                                double goFactor,
                                double geFactor,
                                double cosAngleThresh,
                                bool OverlapAlignment,
                                double dotProdThresh,
                                double gapQuantile,
                                bool hardConstrain,
                                double samples4gradient):
    
    _alignChromatogramsCpp(deref(obj.inst.get()), r1, r2, alignType, tA, tB, normalization, simType,
                                B1p, B2p, noBeef, goFactor, geFactor, cosAngleThresh, OverlapAlignment,
                                dotProdThresh, gapQuantile, hardConstrain, samples4gradient)


