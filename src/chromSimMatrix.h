#ifndef CHROMSIMMATRIX_H
#define CHROMSIMMATRIX_H

#include <vector>
#include "simpleFcn.h"

struct SimMatrix
{
  std::vector<double> data;
  int n_row;
  int n_col;
};

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

std::vector<double> perSampleEucLenVecOfVec(const std::vector<std::vector<double>>& vec);

void distToSim(SimMatrix& s, double offset, double Numerator);

std::vector<std::vector<double>> meanNormalizeVecOfVec(const std::vector<std::vector<double>>& d);

std::vector<std::vector<double>> L2NormalizeVecOfVec(const std::vector<std::vector<double>>& d);

std::vector<std::vector<double>> divideVecOfVec(const std::vector<std::vector<double>>& d, double num);

void ElemWiseSumOuterProd(const std::vector<double>& d1, const std::vector<double>& d2, SimMatrix& s);

void ElemWiseSumOuterEucl(const std::vector<double>& d1, const std::vector<double>& d2, SimMatrix& s);

void ElemWiseOuterCosine(const std::vector<double>& d1, const std::vector<double>& d2, const std::vector<double>& d1_mag, const std::vector<double>& d2_mag, SimMatrix& s);

void SumOuterProd(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, SimMatrix& s);

void SumOuterEucl(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, SimMatrix& s);

void SumOuterCosine(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, SimMatrix& s);

SimMatrix getSimilarityMatrix(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2, const std::string Normalization, const std::string SimType);


#endif // CHROMSIMMATRIX_H
