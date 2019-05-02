#ifndef CHROMSIMMATRIX_H
#define CHROMSIMMATRIX_H

#include <vector>

struct SimMatrix
{
  std::vector<double> data;
  int n_row;
  int n_col;
};

SimMatrix OuterProdMeanNormAllFunc(std::vector<std::vector<double>> d1, std::vector<std::vector<double>> d2);

double meanVecOfVec(std::vector<std::vector<double>> d1);

#endif // CHROMSIMMATRIX_H
