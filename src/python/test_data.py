#!/usr/bin/python
# -*- coding: utf-8  -*-
"""
	test_data.py
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

import numpy as np

FACTOR = 16 # length 192
FACTOR = 32 # length 384
FACTOR = 128 # length 1536
FACTOR = 256 # length 3072
FACTOR = 512 # length 6144
FACTOR = 16 # length 192

data = [ 
            100.0,
            300.00005,
            400.00010,
            450.00010,
            455.00010,
            700.00010
        ]
# 6 vectors of length 196
chromatograms = np.array(data*6*FACTOR*2).reshape(6, len(data)*FACTOR*2)

data2 = [ 
            0.0,
            0.0,
            0.0,
            100.0,
            200.0,
            300.00005,
            400.00010,
            450.00010,
            100.0,
            100.0,
            455.00010,
            700.00010
        ]
# 6 vectors of length 192
chromatograms2 = np.array(data2*6*FACTOR).reshape(6, len(data2)*FACTOR)

print("length", len(data)*FACTOR*2, len(data2)*FACTOR)

rt_data = np.ascontiguousarray(data*FACTOR*2)
rt_data2 = np.ascontiguousarray(data2*FACTOR)
chromatograms = np.ascontiguousarray(chromatograms)
chromatograms2 = np.ascontiguousarray(chromatograms2)

"""
    # run as

    python -m timeit  -r 3 -s "import DIAlignPy; import test_data as d; obj = DIAlignPy.AffineAlignObj(256, 256);" 'DIAlignPy.alignChromatogramsCppSimple(obj, d.chromatograms, d.chromatograms2,b"none", d.rt_data, d.rt_data2, b"mean", b"cosineAngle")'

"""

