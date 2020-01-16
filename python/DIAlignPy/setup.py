"""
	setup.py
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

"""
These Python bindings rely on the original C++ code for the actual function
calls. To compile the bindings, you will need Cython and the Python headers
installed on your system.

To build and test, run

$ python setup.py build_ext --inplace
$ nosetests test_pydialign.py  

"""

import os, shutil
from setuptools import setup, Extension
from Cython.Distutils import build_ext

import os
curr_dir = os.path.dirname(os.path.realpath(__file__))
curr_dir = "."
src_dir = os.path.join(curr_dir, "..", "..", "src")
src_files = [
                       "affinealignment.cpp",
                       "affinealignobj.cpp",
                       "alignment.cpp",
                       "chromSimMatrix.cpp",
                       "constrainMat.cpp",
                       "gapPenalty.cpp",
                       "utils.cpp",
             ]

src_files = [os.path.join(src_dir, f) for f in src_files]
src_files.append("PyDIAlign.pyx") 

include_dirs = [src_dir]

ext_modules = [Extension("DIAlignPy", src_files,
                     language='c++',
                     include_dirs=include_dirs,
                     extra_compile_args=["-std=c++11"], extra_link_args=["-std=c++11"]
                     )]

setup(
    ext_modules = ext_modules,
    cmdclass = {'build_ext': build_ext},

    name="PyDIAlign",

    version="0.1.0"

)

