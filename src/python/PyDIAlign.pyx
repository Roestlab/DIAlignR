#cython: embedsignature=True
#distutils: language = c++
"""
A Python wrapper around the C DIAlign library
"""

# Why was I able to get OK testoutput in build_ext?
# We have to call both of them as actual functions are in regular numpy whereas cimport numpy is needed for compatibility.
import numpy as np
cimport numpy as np

# Definition of STL containers is in libcpp
from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from cython.operator cimport dereference as deref, preincrement as inc, address as address
from DIAlign cimport getQuantile as _getQuantile
from DIAlign cimport AffineAlignObj as _AffineAlignObj
# This is from .pxd file, members mentioned there only, could be accessed here. 
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
    # std::shared_ptr<_AffineAlignObj> inst = std::make_shared<_AffineAlignObj>(); This command reserves memory and define structre as well.
    # may ot pass sthe argument
    # std::shared_ptr<_AffineAlignObj> inst; This only reserves memory.
    # Is this pointer pointing towards stack or heap?
    cdef shared_ptr[_AffineAlignObj] inst

    def __dealloc__(self):
         self.inst.reset() 
         # Is it the native reset to libcpp.memory?
         # Won't shared pointer deallocate by itself? __dealloc__() works with __cinit__().

    def __init__(self, int row, int col):
        """Cython signature: void AffineAlignObj()"""
        # Why do we need to use new? I thought shared_ptr automate this task.
        self.inst = shared_ptr[_AffineAlignObj](new _AffineAlignObj(row, col))

    property score:
        def __set__(self, np.ndarray[double, ndim=1, mode="c"] score not None):
            self.inst.get().score = score
            # Is it possible to make it immutable from python side?
    
        def __get__(self):
            return self.inst.get().score
            # Is it possible to read it without property?

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


