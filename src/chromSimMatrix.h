#ifndef CHROMSIMMATRIX_H
#define CHROMSIMMATRIX_H

#include <vector>
#include "utils.h"
#include "similarityMatrix.h"

namespace DIAlign 
{

namespace SimilarityMatrix 
{
// functor for getting sum of previous result and square of current element.
// TODO: Need to understand the implementation.
template<typename T>
struct square
{
  T operator()(const T& Left, const T& Right) const
  {
    // We use this struct as binary operation function object. It should take current accumulation value (Left) and value of current element (Right).
    return (Left + Right*Right);
  }
};

double meanVecOfVec(const std::vector<std::vector<double>>& vec);

double eucLenVecOfVec(const std::vector<std::vector<double>>& vec);

// Eucledian length at each time-point.
std::vector<double> perSampleEucLenVecOfVec(const std::vector<std::vector<double>>& vec);

std::vector<double> perSampleSqrSumVecOfVec(const std::vector<std::vector<double>>& vec);

std::vector<double> perSampleMeanVecOfVec(const std::vector<std::vector<double>>& vec);

std::vector<double> perSampleSumVecOfVec(const std::vector<std::vector<double>>& vec);

void distToSim(SimMatrix& s, double offset, double Numerator);

void clamp(std::vector<double>& vec, double minValue, double maxValue);

std::vector<std::vector<double>> meanNormalizeVecOfVec(const std::vector<std::vector<double>>& d);

std::vector<std::vector<double>> L2NormalizeVecOfVec(const std::vector<std::vector<double>>& d);

std::vector<std::vector<double>> divideVecOfVec(const std::vector<std::vector<double>>& d, double num);

void ElemWiseSumOuterProd(const std::vector<double>& d1, const std::vector<double>& d2, SimMatrix& s);

void ElemWiseSumOuterProdMeanSub(const std::vector<double>& d1, const std::vector<double>& d2, SimMatrix& s, const std::vector<double>& mean1, const std::vector<double>& mean2);

void ElemWiseSumOuterEucl(const std::vector<double>& d1, const std::vector<double>& d2, SimMatrix& s);

void ElemWiseOuterCosine(const std::vector<double>& d1, const std::vector<double>& d2, const std::vector<double>& d1_mag, const std::vector<double>& d2_mag, SimMatrix& s);

void SumOuterProd(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, SimMatrix& s);

void SumOuterCov(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, SimMatrix& s);

void SumOuterCorr(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, SimMatrix& s);

void SumOuterEucl(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, SimMatrix& s);

void SumOuterCosine(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, SimMatrix& s);

SimMatrix getSimilarityMatrix(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, const std::string SimType, double cosAngleThresh, double dotProdThresh);

} // namespace SimilarityMatrix
} // namespace DIAlign

#endif // CHROMSIMMATRIX_H
