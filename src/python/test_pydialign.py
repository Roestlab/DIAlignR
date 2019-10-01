#!/usr/bin/python
# -*- coding: utf-8  -*-
"""
	test_pydialign.py
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

import unittest
import numpy as np
# from unittest import assertAlmostEqual        

import DIAlignPy as PyDIAlign

data = [ 
            100.0,
            200.0,
            300.00005,
            400.00010,
            450.00010,
            455.00010,
            700.00010
        ]
# 6 vectors of length 7
chromatograms = np.array(data*6).reshape(6, len(data))

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
# 5 vectors of length 12
chromatograms2 = np.array(data2*6).reshape(6, len(data2))

data = np.ascontiguousarray(data)
data2 = np.ascontiguousarray(data2)
chromatograms = np.ascontiguousarray(chromatograms)
chromatograms2 = np.ascontiguousarray(chromatograms2)

class TestMSDIAlign(unittest.TestCase):

    def setUp(self):

        self.data = data
        self.data2 = data2
        self.chromatograms = chromatograms
        self.chromatograms2 = chromatograms2

    def test_affineAlignObj(self):
        obj = PyDIAlign.AffineAlignObj(0,0)

    def test_doAlignment(self):
        obj = PyDIAlign.AffineAlignObj(256, 256); 
        PyDIAlign.alignChromatogramsCppSimple(obj, self.chromatograms, self.chromatograms,
                b"none", self.data, self.data, b"mean", b"cosineAngle")

        self.assertAlmostEqual(obj.score[0], 1.0, 4)
        self.assertAlmostEqual(obj.score[5], 6.0, 4)

        self.assertEqual(obj.indexA_aligned[0], 1)
        self.assertEqual(obj.indexA_aligned[5], 6)
        self.assertEqual(obj.indexB_aligned[0], 1)
        self.assertEqual(obj.indexB_aligned[5], 6)

    def test_doAlignment2(self):
        obj = PyDIAlign.AffineAlignObj(256, 256); 
        PyDIAlign.alignChromatogramsCppSimple(obj, self.chromatograms, self.chromatograms2,
                b"none", self.data, self.data2, b"mean", b"cosineAngle")

        self.assertEqual(len(obj.score), 12)
        self.assertEqual(obj.score[0], 0.0)
        self.assertEqual(obj.score[11], 7.0)

        self.assertEqual(obj.indexA_aligned[0], 0)
        self.assertEqual(obj.indexA_aligned[4], 0)
        self.assertEqual(obj.indexA_aligned[5], 1)
        self.assertEqual(obj.indexA_aligned[11], 7)

        # B has no gaps
        self.assertEqual(obj.indexB_aligned[0], 1)
        self.assertEqual(obj.indexB_aligned[5], 6)
        self.assertEqual(obj.indexB_aligned[11], 12)

    def test_doAlignment2(self):
        obj = PyDIAlign.AffineAlignObj(256, 256); 
        PyDIAlign.alignChromatogramsCppSimple(obj, self.chromatograms, self.chromatograms2,
                b"none", self.data, self.data2, b"mean", b"dotProductMasked")

        self.assertEqual(len(obj.score), 12)
        self.assertEqual(obj.score[0], 0.0)
        self.assertAlmostEqual(obj.score[11], 76.99906921386719, 5)

        # Now A has a gap in the middle.
        self.assertEqual(obj.indexA_aligned[0], 0)
        self.assertEqual(obj.indexA_aligned[4], 1)
        self.assertEqual(obj.indexA_aligned[5], 2)
        self.assertEqual(obj.indexA_aligned[6], 3)
        self.assertEqual(obj.indexA_aligned[7], 4)
        self.assertEqual(obj.indexA_aligned[8], 0)
        self.assertEqual(obj.indexA_aligned[11], 7)

        # B has no gaps
        self.assertEqual(obj.indexB_aligned[0], 1)
        self.assertEqual(obj.indexB_aligned[5], 6)
        self.assertEqual(obj.indexB_aligned[11], 12)

if __name__ == '__main__':
    unittest.main()

