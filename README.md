# DIAlignR
DIAlignR is an R package for retention time alignment of targeted mass spectrometric data, including DIA and SWATH-MS data. This tool works with MS2 chromatograms directly and uses dynamic programming for alignment of raw chromatographic traces. DIAlignR uses a hybrid approach of global (feature-based) and local (raw data-based) alignment to establish correspondence between peaks.

[![Travis build status](https://travis-ci.org/shubham1637/DIAlignR.svg?branch=master)](https://travis-ci.org/shubham1637/DIAlignR)

# Documentation
For documentation please see [our vignette](http://htmlpreview.github.io/?https://github.com/shubham1637/DIAlignR/blob/master/doc/DIAlignR-vignette.html).

# Developing C++ code
```
cd DIAlignR
mkdir build && cd build
cmake -B. -H.. 
make clean && make && make test
make runTest3
cd ..
```

Documenting C++ code
```
sudo apt install doxygen doxygen-gui 
sudo apt install graphviz
cd DIAlignR
cd src
doxygen doc/Doxyfile
```

# Citation
If you use the provided algorithms or the package, please cite our paper:

Gupta S, Ahadi S, Zhou W, Röst H. "DIAlignR Provides Precise Retention Time Alignment Across Distant Runs in DIA and Targeted Proteomics." Mol Cell Proteomics. 2019 Apr;18(4):806-817. doi: https://doi.org/10.1074/mcp.TIR118.001132 Epub 2019 Jan 31.

CNPN 2018 Poster doi: https://doi.org/10.6084/m9.figshare.6200837.v1     
HUPO 2018 Poster doi: https://doi.org/10.6084/m9.figshare.7121696.v2     
