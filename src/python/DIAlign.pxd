# distutils: language = c++
"""
The following provides a Python wrapper around the C++ DIAlign library


	DIAlign.pxd
        hannes.rost@utoronto.ca
 
	Copyright 2019 Hannes Roest

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string
from libcpp cimport bool

cdef extern from "utils.h" namespace "DIAlign":

    double getQuantile(libcpp_vector[double] vec, double quantile) nogil except +

cdef extern from "affinealignobj.h" namespace "DIAlign":

    cdef cppclass AffineAlignObj:
        AffineAlignObj(int row, int col) nogil except +
        AffineAlignObj(AffineAlignObj o) nogil except + # Using Copy constructor
        # Can we also use Copy assignment operator here?

        # These are not functions or class, therefore, nogil and except + are not needed.
        # Other members of the class are pointers and it could be tricky with Python to
        # handle them. We need to write a separate function to deal with them.
        libcpp_vector[int] indexA_aligned
        libcpp_vector[int] indexB_aligned
        libcpp_vector[double] score

cdef extern from "CppInterface.hpp" namespace "DIAlign":

    void alignChromatogramsCpp( AffineAlignObj& obj,
                                libcpp_vector[libcpp_vector[double] ] & r1,
                                libcpp_vector[libcpp_vector[double] ] & r2,
                                libcpp_string alignType,
                                libcpp_vector[double]& tA,
                                libcpp_vector[double]& tB,
                                libcpp_string & normalization,
                                libcpp_string simType) nogil except +

    void alignChromatogramsCpp( AffineAlignObj& obj,
                                libcpp_vector[libcpp_vector[double] ] & r1,
                                libcpp_vector[libcpp_vector[double] ] & r2,
                                libcpp_string alignType,
                                libcpp_vector[double]& tA,
                                libcpp_vector[double]& tB,
                                libcpp_string & normalization,
                                libcpp_string simType,
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
                                double samples4gradient) nogil except +
