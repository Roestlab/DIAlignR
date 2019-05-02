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

SimMatrix SumOuterProdMeanNormFrag(const std::vector<std::vector<double>>& d1, const std::vector<std::vector<double>>& d2);

double meanVecOfVec(std::vector<std::vector<double>> d1);

std::vector<std::vector<double>> divideVecOfVec(const std::vector<std::vector<double>>& d, double num);

void ElemWiseSumOuterProd(const std::vector<double>& d1, const std::vector<double>& d2, SimMatrix& s);

#endif // CHROMSIMMATRIX_H
