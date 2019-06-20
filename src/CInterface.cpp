
#include "CInterface.h"

#include "CppInterface.hpp"

CAffineAlignObj alignChromatogramsC(
                           int nchrom,
                           double* r1, int ndata1,
                           double* r2, int ndata2,
                           const char* alignType,
                           double* tA, int ntA,
                           double* tB, int ntB,
                           const char* normalization, const char* simType,
                           double B1p, double B2p, int noBeef,
                           double goFactor, double geFactor,
                           double cosAngleThresh, bool OverlapAlignment,
                           double dotProdThresh, double gapQuantile,
                           bool hardConstrain, double samples4gradient)
{
  // Convert C data structs to C++ structs
  std::vector<std::vector<double> > d1;
  std::vector<std::vector<double> > d2;
  for (int chr = 0; chr < nchrom; chr++)
  {
    std::vector<double> vec(r1 + chr*ndata1, r1 + chr*ndata1 + ndata1);
    d1.push_back(std::move(vec));
  }
  for (int chr = 0; chr < nchrom; chr++)
  {
    std::vector<double> vec(r2 + chr*ndata2, r2 + chr*ndata2 + ndata2);
    d2.push_back(std::move(vec));
  }

  std::vector<double> timeA(tA, tA + ntA);
  std::vector<double> timeB(tB, tB + ntB);

  // Call C++ function
  auto obj = DIAlign::alignChromatogramsCpp(d1, d2, alignType, timeA, timeB, normalization, simType, 
                           B1p, B2p, noBeef, goFactor, geFactor, cosAngleThresh, OverlapAlignment,
                           dotProdThresh, gapQuantile, hardConstrain, samples4gradient);

  // Convert C++ result back to C
  CAffineAlignObj ret;

  ret.n_indexA_aligned = obj.indexA_aligned.size();
  ret.n_indexB_aligned = obj.indexB_aligned.size();
  ret.indexA_aligned = (int*) malloc(ret.n_indexA_aligned * sizeof(int));;
  ret.indexB_aligned = (int*) malloc(ret.n_indexB_aligned * sizeof(int));;

  ret.n_score = obj.score.size();
  ret.score = (double*) malloc(ret.n_score * sizeof(double));

  std::memcpy(ret.indexA_aligned, &obj.indexA_aligned[0], ret.n_indexA_aligned * sizeof(int));
  std::memcpy(ret.indexB_aligned, &obj.indexB_aligned[0], ret.n_indexB_aligned * sizeof(int));
  std::memcpy(ret.score, &obj.score[0], ret.n_score * sizeof(double));

  return ret;
}

